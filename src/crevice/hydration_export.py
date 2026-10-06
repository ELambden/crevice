"""Hydration tables, all-frame histories and figures with explicit evidence limits."""
from pathlib import Path
import json,csv

from .presentation import FigureText, annotate_option

LIMITATIONS = [
 'Water-oxygen proximity is not a hydrogen-bond or hydration-free-energy calculation.',
 'Geometric probe-accessible area is not evidence that water can reach a sealed cavity.',
 'Static deposited waters are incomplete observations; absent waters do not establish dryness.',
 'Water exchange between saved trajectory frames is unresolved; sampled persistence is not a continuous-time residence lifetime.',
 'Surface quadrature, radii, contact cutoff, force field and anatomical assignment are not included in temporal sampling intervals.',
 'Hydration, geometric exposure and motion associations do not establish functional or causal control.',
]


@annotate_option
def write_hydration_bundle(analysis,output_dir,*,prefix='crevice',confidence=.95,block_length=None,replicates=2000,seed=20260914,dpi=240,roles=None,annotate=None):
    """Write the statistics, tables, figures and scenes of ``crevice hydration``.

    For every focused residue and frame, ``analysis`` holds the unique explicit
    water oxygens within the contact cutoff of any residue heavy atom and the
    Shrake-Rupley solvent-accessible area (probe plus CREVICE van der Waals
    radii). Per residue, the mean, frame quantiles and, for trajectories,
    autocorrelation-aware block-bootstrap mean intervals are reported; the
    water-contact occupancy is the fraction of frames with at least one contact
    water. Structures without explicit water report hydration as unobserved
    (null), never zero. Optional water density, typed interactions and the
    coordinated viewer scenes and overview figures are written when present in
    ``analysis``. No HTML is written.

    Parameters
    ----------
    analysis : dict
        Result of :func:`crevice.hydration.static_hydration` or
        :func:`crevice.hydration_trajectory.analyze_hydration_trajectory`.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    confidence : float, default 0.95
        Nominal level of mean intervals, in [0, 1); 0 gives descriptive statistics only.
    block_length : int, optional
        Bootstrap block length in frames.
    replicates : int, default 2000
        Bootstrap replicates (at least 200 when intervals are requested).
    seed : int, default 20260914
        Bootstrap random seed.
    dpi : int, default 240
        Figure resolution.
    roles : dict of str to str, optional
        Residue role labels (for example ``boundary_lining``) for static input.
    annotate : bool, optional
        Draw figure and panel titles here and in the density, overview and
        colour-scale figures, and the legend, colour-scale ticks, caption and
        region banner in the PyMOL/VMD/ChimeraX scenes. ``None`` (default)
        inherits the surrounding setting, which is off unless enabled with
        :func:`crevice.presentation.figure_annotations` or ``--annotate``.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Raises
    ------
    ValueError
        For an invalid confidence level, too few replicates or a bad prefix.

    Outputs
    -------
    PREFIX_hydration.json : JSON
        Settings, definitions, per-residue statistics, shared-water proximity,
        limitations.
    PREFIX_hydration_summary.csv : table
        Every interval statistic of the JSON flattened to one row (per residue
        and metric, per shared-water pair, unique focus waters, fixed-region
        waters) plus the correlations; columns as in
        :func:`crevice.io.summary_statistics_rows`.
    PREFIX_hydration_residues.csv, PREFIX_hydration_shared_water.csv : tables
        Per-residue means, SDs and intervals; residue pairs sharing a water
        (header only when no pair shares a water).
    PREFIX_hydration_frames.csv : table
        One row per frame and focused residue: ``frame_index``, ``time_ps``,
        ``residue``, water counts, SASA (``*_A2``) and, with a cavity
        trajectory, boundary area, partner flag and aligned centroid (see
        :func:`hydration_frame_tables`).
    PREFIX_hydration_region_frames.csv : table
        One row per frame: ``frame_index``, ``time_ps``,
        ``unique_focus_water_count`` and, with a cavity trajectory,
        ``region_water_count`` and ``cavity_volume_A3``.
    PREFIX_hydration_frames.npz : arrays
        The same per-frame arrays plus sparse residue-water contact identities.
    PREFIX_hydration_residues.png : figure
        Up to 24 residues: "Mean water oxygens within contact cutoff" (or
        "Observed water oxygens" for one structure) and "Geometric
        probe-accessible area (Å²)" for all heavy atoms and N/O atoms, with
        interval bars for trajectories.
    PREFIX_hydration_heatmap.png, PREFIX_hydration_timeseries.png : figures
        Trajectories only: per-residue water contacts and SASA against time
        (ns), and unique nearby waters / waters in the fixed region against time.
    PREFIX_hydration_observed_episodes.csv, PREFIX_hydration_sampled_survival.csv, PREFIX_hydration_mobility.csv : trajectory tables
        Contact episodes, sampled survival and centroid mobility (definition in
        the manifest's ``hydration_mobility`` entry).
    PREFIX_hydration_alignment.json : JSON
        Trajectories only: per-frame alignment fits (rotation, translation,
        RMSD) and periodic boxes, kept as provenance.
    PREFIX_water_density.* : density
        See :func:`crevice.water_density.write_density_bundle`.
    PREFIX_interactions.json, PREFIX_interaction_edges.csv, PREFIX_interaction_edge_frames.csv, PREFIX_residue_conformations.csv, PREFIX_residue_conformation_frames.csv, PREFIX_interaction_frames.npz : interactions
        Typed contacts and conformations, when requested (see
        :func:`crevice.interaction_analysis.write_interaction_bundle`).
    PREFIX_analysis_* : scenes, evidence table and figures
        See :func:`crevice.analysis_viewers.write_analysis_assets` and
        :func:`crevice.analysis_viewers.write_analysis_views`.
    PREFIX_hydration_manifest.json : JSON
        File index, frame and residue counts, water status, limitations and
        the mobility definition.

    Colours and representation
    --------------------------
    * Residue figure: teal (``#178c91``) water-count bars; orange
      (``#d98c33``) all-atom SASA and blue-grey (``#547c9b``) N/O SASA bars,
      with a two-entry legend; dark interval bars.
    * Heatmap: ``viridis`` (water contacts) and ``cividis`` (SASA, Å²) with
      labelled colour bars; light grey = missing frames.
    * Timeseries: teal unique-water line, blue-grey fixed-region line.
    """
    import numpy as np
    from .hydration import METRICS
    from .hydration_trajectory import sampled_survival,contact_episodes
    from .cavity_trajectory import describe_series
    from .figures import _pyplot,_save_figure
    import math
    if not math.isfinite(confidence) or not 0<=confidence<1:raise ValueError('Hydration confidence must be in [0,1)')
    if confidence and (isinstance(replicates,bool) or int(replicates)!=replicates or replicates<200):raise ValueError('Hydration bootstrap requires at least 200 replicates')
    root=Path(output_dir);root.mkdir(parents=True,exist_ok=True)
    if not prefix or Path(prefix).name!=prefix:raise ValueError('Hydration prefix must be a simple filename stem')
    static='metrics' in analysis
    if static:
        ids=list(analysis['residue_ids']);arrays={k:np.asarray(v)[None,:] for k,v in analysis['metrics'].items()}
        arrays.update(time_ps=np.array([0.]),frame_indices=np.array([0]),residue_ids=np.asarray(ids))
        histories=[analysis['memberships']];bridges=[analysis['bridges']]
    else:arrays=dict(analysis['arrays']);ids=list(arrays['residue_ids']);histories=analysis['histories'];bridges=analysis['bridge_histories']
    n=len(arrays['time_ps']);r=len(ids);times=arrays['time_ps'];regular=n<3 or np.allclose(np.diff(times),np.diff(times)[0],rtol=1e-7,atol=1e-7)
    options=dict(confidence=confidence or .95,block_length=block_length,replicates=replicates,seed=seed,regular=regular and confidence!=0)
    def describe(x):
        result=describe_series(x,**options)
        if confidence==0:result['confidence_interval']={'status':'not_requested'}
        for field in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
            interval=result['confidence_interval']
            if field in interval:interval[field]=[max(0.,v) if v is not None else None for v in interval[field]]
        return result
    def correlation(x,y):
        good=np.isfinite(x)&np.isfinite(y)
        if good.sum()<3 or min(np.std(x[good]),np.std(y[good]))<1e-10:return None
        return float(np.corrcoef(x[good],y[good])[0,1])
    records=[];roles=roles or {}
    for j,key in enumerate(ids):
        water=arrays['water_count'][:,j];occupancy=np.where(np.isfinite(water),(water>0).astype(float),np.nan)
        role=roles.get(key,'selected_residue')
        if 'boundary_area_A2' in arrays:role='ever_boundary' if np.any(arrays['boundary_area_A2'][:,j]>0) else 'ever_nonlining_partner'
        row={'residue':str(key),'role':role,'statistics':{metric:describe(arrays[metric][:,j]) for metric in METRICS},'water_contact_occupancy':describe(occupancy),
             'unique_water_ids_observed':len(set().union(*(h[j] for h in histories))),
             'hydration_sasa_correlation':correlation(water,arrays['sasa_A2'][:,j]) if not static else None}
        if static and 'boundary_area_by_residue_A2' in analysis:row['mean_boundary_area_A2']=analysis['boundary_area_by_residue_A2'].get(str(key))
        ci=row['water_contact_occupancy']['confidence_interval']
        for field in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
            if field in ci:ci[field]=[float(np.clip(v,0,1)) if v is not None else None for v in ci[field]]
        if 'cavity_volume_A3' in arrays:
            row['hydration_volume_correlation']=correlation(water,arrays['cavity_volume_A3'])
            eligible=(arrays['boundary_area_A2'][:,j]>0)&np.isfinite(water)
            row['boundary_eligible_frames']=int(eligible.sum());row['mean_water_count_when_boundary']=float(water[eligible].mean()) if eligible.any() else None
            row['mean_boundary_area_A2']=float(np.nanmean(arrays['boundary_area_A2'][:,j])) if np.isfinite(arrays['boundary_area_A2'][:,j]).any() else None
        if 'aligned_residue_centroids_A' in arrays:
            xyz=arrays['aligned_residue_centroids_A'][:,j];row['aligned_centroid_rmsf_A']=float(np.sqrt(np.mean(np.sum((xyz-xyz.mean(0))**2,axis=1))))
            row['hydration_centroid_distance_correlation']=correlation(water,arrays['centroid_distance_A'][:,j])
        records.append(row)
    pair_frames={}
    for i,frame in enumerate(bridges):
        for pair in frame:pair_frames.setdefault(pair,[]).append(i)
    bridge_records=[]
    for (i,j),frames in sorted(pair_frames.items()):
        occupancy=np.zeros(n);occupancy[frames]=1
        stat=describe(occupancy)
        for field in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
            ci=stat['confidence_interval']
            if field in ci:ci[field]=[float(np.clip(v,0,1)) if v is not None else None for v in ci[field]]
        bridge_records.append({'source':str(ids[i]),'target':str(ids[j]),'observed_bridge_frames':len(frames),
           'eligible_frames':n,'occupancy':len(frames)/n,'occupancy_statistics':stat})
    unique=np.asarray([len(set().union(*h)) for h in histories],dtype=float) if analysis['water_population'] else np.full(n,np.nan)
    arrays['unique_focus_water_count']=unique
    report={'mode':'static_observed_structure' if static else 'trajectory','frame_count':n,'residue_count':r,
      'water_population':analysis['water_population'],'water_status':'observed_explicit_water' if analysis['water_population'] else 'unavailable_no_explicit_water',
      'settings':analysis['settings'],'water_definition':'Unique water oxygen identities within cutoff of any selected residue heavy atom; backbone/sidechain/NO sets may overlap',
      'sasa_definition':'Shrake–Rupley equal-area sphere quadrature using CREVICE VDW radii plus probe; all supplied nonwater heavy atoms are obstacles',
      'hydrogen_bonds':{'status':'not_assigned','reason':'Proximity lacks validated donor/acceptor, protonation and directional hydrogen-bond geometry'},
      'confidence_definition':'Approximate autocorrelation-aware mean intervals; no interval for single/constant/inadequately sampled series; no simultaneous inference across residues',
      'unique_focus_water_count':describe(unique),'residues':records,'shared_water_proximity':bridge_records,
      'shared_water_definition':'The same water oxygen is within the cutoff of both residues in the same frame; not a directional water-mediated hydrogen bond',
      'limitations':LIMITATIONS}
    if static:report.update(source=analysis['source'],partial_occupancy_waters=analysis['partial_occupancy_waters'],static_water_policy='Positive-occupancy deposited waters counted once; partial occupancies are reported, not converted into equilibrium populations')
    else:
        report['reader']=analysis['reader'];report['sample_spacing_ps']=float(np.diff(times)[0]) if n>1 and regular else None
        if 'region_water_count' in arrays:report['reference_region_water_count']=describe(arrays['region_water_count'])
        report['sampled_survival']={str(key):sampled_survival([h[j] for h in histories],times) for j,key in enumerate(ids)}
        if 'cavity_volume_A3' in arrays:report['regional_water_volume_correlation']=correlation(arrays['region_water_count'],arrays['cavity_volume_A3'])
    files={}
    def path(key,suffix):
        p=root/(prefix+suffix);files[key]=str(p);return p
    def table(p,rows,fields=None):
        if not rows and not fields:p.write_text('');return
        with p.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(fields or rows[0]));w.writeheader();w.writerows(rows)
    path('hydration_json','_hydration.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    from .io import write_summary_statistics_csv
    write_summary_statistics_csv(report,path('hydration_summary_csv','_hydration_summary.csv'))
    table_rows=[]
    for row in records:
        entry={k:v for k,v in row.items() if not isinstance(v,(dict,list))}
        for key,stat in row['statistics'].items():
            ci=stat['confidence_interval'];entry.update({key+'_mean':stat['mean'],key+'_sd':stat['sample_sd'],key+'_observed_frames':stat['observed_frames'],
              key+'_mean_lower':ci.get('pointwise_lower',[None])[0],key+'_mean_upper':ci.get('pointwise_upper',[None])[0],key+'_CI_status':ci['status']})
        stat=row['water_contact_occupancy'];ci=stat['confidence_interval'];entry.update(water_contact_occupancy=stat['mean'],water_contact_occupancy_lower=ci.get('pointwise_lower',[None])[0],water_contact_occupancy_upper=ci.get('pointwise_upper',[None])[0],water_contact_CI_status=ci['status']);table_rows.append(entry)
    table(path('hydration_residues_csv','_hydration_residues.csv'),table_rows)
    bridge_rows=[]
    for row in bridge_records:
        ci=row['occupancy_statistics']['confidence_interval'];bridge_rows.append({k:v for k,v in row.items() if k!='occupancy_statistics'}|{'CI_status':ci['status'],
            'occupancy_lower':float(np.clip(ci['pointwise_lower'][0],0,1)) if ci.get('pointwise_lower',[None])[0] is not None else None,
            'occupancy_upper':float(np.clip(ci['pointwise_upper'][0],0,1)) if ci.get('pointwise_upper',[None])[0] is not None else None})
    table(path('hydration_shared_water_csv','_hydration_shared_water.csv'),bridge_rows,
          None if bridge_rows else ['source','target','observed_bridge_frames','eligible_frames','occupancy','CI_status','occupancy_lower','occupancy_upper'])
    water_ids=sorted(set().union(*(set().union(*h) for h in histories)));lookup={key:i for i,key in enumerate(water_ids)}
    sparse=np.asarray([(i,j,lookup[w]) for i,h in enumerate(histories) for j,waters in enumerate(h) for w in sorted(waters)],dtype=int).reshape((-1,3))
    arrays.update(contact_water_ids=np.asarray(water_ids,dtype=str),contact_frame_residue_water_indices=sparse)
    np.savez_compressed(path('hydration_frames_npz','_hydration_frames.npz'),**arrays)
    residue_frames,region_frames=hydration_frame_tables(arrays)
    table(path('hydration_frames_csv','_hydration_frames.csv'),residue_frames)
    if region_frames and len(region_frames[0])>2:
        table(path('hydration_region_frames_csv','_hydration_region_frames.csv'),region_frames)
    if not static:
        table(path('hydration_episodes_csv','_hydration_observed_episodes.csv'),contact_episodes(histories,times,ids))
        table(path('hydration_survival_csv','_hydration_sampled_survival.csv'),[{'residue':key,**row} for key,s in report['sampled_survival'].items() for row in s['rows']])
        path('hydration_alignment_json','_hydration_alignment.json').write_text(json.dumps({'fits':analysis['alignment'],'boxes':analysis['boxes']},indent=2)+'\n')
    plt=_pyplot();ranked=sorted(range(r),key=lambda j:(records[j]['role'] not in {'ever_boundary','boundary_lining'},-(records[j].get('mean_boundary_area_A2') or 0),-(records[j]['statistics']['water_count']['mean'] or 0),ids[j]))[:24]
    labels=[str(ids[j]) for j in ranked]
    text=FigureText(f"Water contacts and geometric solvent exposure of {len(ranked)} focused residues "
                    +('in one structure' if static else f'averaged over {n} frames')
                    +('' if analysis['water_population'] else '; no explicit waters were supplied, so hydration is unobserved'))
    fig,axes=plt.subplots(1,2,figsize=(12,8),sharey=True,layout='constrained')
    counts=[records[j]['statistics']['water_count']['mean'] for j in ranked]
    axes[0].barh(labels,[v or 0 for v in counts],color='#178c91');axes[0].invert_yaxis();axes[0].set_xlabel('Observed water oxygens' if static else 'Mean water oxygens within contact cutoff')
    if not analysis['water_population']:text.notice(axes[0],.5,.5,'No explicit waters supplied\nHydration is unobserved',transform=axes[0].transAxes,ha='center',bbox={'facecolor':'white','edgecolor':'none'})
    axes[1].barh(labels,[records[j]['statistics']['sasa_A2']['mean'] or 0 for j in ranked],color='#d98c33',label='All selected heavy atoms')
    axes[1].barh(labels,[records[j]['statistics']['NO_sasa_A2']['mean'] or 0 for j in ranked],color='#547c9b',label='N/O atoms')
    axes[1].set_xlabel('Geometric probe-accessible area (Å²)');axes[1].legend(frameon=False,fontsize=8)
    if not analysis['water_population']:
        axes[0].set_xticks([]);axes[0].set_xlabel('Explicit hydration is unobserved')
    if not static:
        for axis,metric in [(axes[0],'water_count'),(axes[1],'sasa_A2')]:
            for position,j in enumerate(ranked):
                stat=records[j]['statistics'][metric];ci=stat['confidence_interval'];lower=ci.get('pointwise_lower',[None])[0];upper=ci.get('pointwise_upper',[None])[0]
                if lower is not None and upper is not None:
                    axis.errorbar(stat['mean'],position,xerr=[[max(0.,stat['mean']-lower)],[max(0.,upper-stat['mean'])]],fmt='none',ecolor='#263840',capsize=2,lw=.8)
        text.title(axes[0],f'{confidence:.0%} pointwise mean intervals where estimable' if confidence else 'Descriptive means; no confidence intervals requested',fontsize=9)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    text.suptitle(fig,'Residue water contacts and solvent exposure'+(' · static snapshot' if static else f' · all {n} frames'))
    _save_figure(fig,path('hydration_residues_png','_hydration_residues.png'),dpi=dpi,text=text)
    if not static:
        text=FigureText(f"Water-oxygen contacts (top) and geometric probe-accessible area (bottom) of {len(ranked)} residues in every frame")
        fig,axes=plt.subplots(2,1,figsize=(11,10),sharex=True,layout='constrained')
        for ax,key,title,cmap in [(axes[0],'water_count','Water-oxygen contacts','viridis'),(axes[1],'sasa_A2','Geometric probe-accessible area (Å²)','cividis')]:
            # The colour bar carries the quantity; the panel title is annotation only.
            image=ax.pcolormesh(times/1000,np.arange(len(ranked)),arrays[key][:,ranked].T,shading='nearest',cmap=cmap);ax.set_yticks(np.arange(len(ranked)),labels,fontsize=8);ax.invert_yaxis();ax.set_facecolor('#dddddd');text.title(ax,title);fig.colorbar(image,ax=ax,label=title)
        axes[1].set_xlabel('Time (ns)');_save_figure(fig,path('hydration_heatmap_png','_hydration_heatmap.png'),dpi=dpi,text=text)
        text=FigureText("Unique water oxygens near the selected residues (top) and waters in the fixed reference region (bottom) against time")
        fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True,layout='constrained')
        axes[0].plot(times/1000,unique,color='#178c91',lw=.8);axes[0].set_ylabel('Unique waters\nnear selected residues')
        if np.isfinite(arrays.get('region_water_count',[])).any():axes[1].plot(times/1000,arrays['region_water_count'],color='#547c9b',lw=.8)
        else:text.notice(axes[1],.5,.5,'No fixed-region water measurement',transform=axes[1].transAxes,ha='center')
        axes[1].set_ylabel('Waters in fixed region');axes[1].set_xlabel('Time (ns)')
        for ax in axes:ax.spines[['top','right']].set_visible(False)
        text.suptitle(fig,'Explicit solvent observations · contact persistence is limited by saved-frame spacing')
        _save_figure(fig,path('hydration_timeseries_png','_hydration_timeseries.png'),dpi=dpi,text=text)
    if not static:files.update(write_hydration_mobility(arrays,root,prefix))
    if 'water_density' in analysis:
        from .water_density import write_density_bundle
        files.update(write_density_bundle(analysis['water_density'],root,prefix,level=analysis['analysis_options']['density_level'],dpi=dpi))
    if 'interactions' in analysis:
        from .interaction_analysis import write_interaction_bundle
        files.update(write_interaction_bundle(analysis['interactions'],root,prefix,times=None if static else times,
            hydration_arrays=arrays,confidence=confidence,block_length=block_length,replicates=replicates,seed=seed))
        report['hydrogen_bonds']={'status':analysis['interactions']['metadata']['hydrogen_bond_status'],
                                 'details':Path(files['interaction_json']).name}
        Path(files['hydration_json']).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if 'reference_frame' in analysis:
        from .analysis_viewers import write_analysis_assets,write_analysis_views
        files.update(write_analysis_assets(analysis,report,root,prefix))
        if not analysis.get('analysis_options',{}).get('skip_views',False):
            files.update(write_analysis_views(root,prefix,files,dpi=dpi))

    from .radii import radii_fields
    # The set the hydration analysis measured with (it may be the set recorded by --cavity-results).
    used={'radii':analysis['settings']['radii']} if 'radii' in analysis.get('settings',{}) else radii_fields()
    manifest={'files':files,'frame_count':n,'residue_count':r,'water_status':report['water_status'],'limitations':LIMITATIONS,**used}
    if 'radii_check' in analysis.get('settings',{}):manifest['radii_check']=analysis['settings']['radii_check']
    if 'hydration_mobility_csv' in files:manifest['hydration_mobility']=MOBILITY_NOTES
    path('hydration_manifest_json','_hydration_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return files


def hydration_mobility_statistics(arrays):
    """Hydration at an interval start versus the following centroid displacement.

    For each residue and each pair of adjacent saved frames, the step is the
    displacement (Å) of the aligned residue centroid; the interval is *wet at
    both ends* when the residue has a contact water in both frames, *dry at both
    ends* when it has none in either, and otherwise *changed*. Descriptive only:
    exchanges between saved frames are not resolved.

    Parameters
    ----------
    arrays : dict of str to numpy.ndarray
        Hydration trajectory arrays with ``time_ps``, ``residue_ids``,
        ``water_count`` (frames x residues) and ``aligned_residue_centroids_A``
        (frames x residues x 3).

    Returns
    -------
    list of dict
        One row per residue: ``observed_intervals``, ``mean_centroid_step_A``,
        ``wet_at_both_ends_intervals``, ``dry_at_both_ends_intervals``,
        ``changed_hydration_state_intervals``, ``mean_step_wet_at_both_ends_A``,
        ``mean_step_dry_at_both_ends_A`` and
        ``starting_hydration_following_step_correlation`` (Pearson, ``None``
        with fewer than three observed intervals or no variation). Empty when
        there are no centroids or fewer than two frames.
    """
    import numpy as np
    if 'aligned_residue_centroids_A' not in arrays or len(arrays['time_ps'])<2:return []
    steps=np.linalg.norm(np.diff(arrays['aligned_residue_centroids_A'],axis=0),axis=2)
    rows=[]
    for j,key in enumerate(arrays['residue_ids']):
        water=arrays['water_count'][:-1,j];following=arrays['water_count'][1:,j];step=steps[:,j]
        observed=np.isfinite(water)&np.isfinite(following)
        wet=observed&(water>0)&(following>0);dry=observed&(water==0)&(following==0)
        corr=None
        if observed.sum()>=3 and min(np.std(water[observed]),np.std(step[observed]))>1e-10:
            corr=float(np.corrcoef(water[observed],step[observed])[0,1])
        rows.append({'residue':str(key),'observed_intervals':int(observed.sum()),'mean_centroid_step_A':float(step.mean()),
             'wet_at_both_ends_intervals':int(wet.sum()),'dry_at_both_ends_intervals':int(dry.sum()),
             'changed_hydration_state_intervals':int(np.sum(observed&~wet&~dry)),
             'mean_step_wet_at_both_ends_A':float(step[wet].mean()) if wet.any() else None,
             'mean_step_dry_at_both_ends_A':float(step[dry].mean()) if dry.any() else None,
             'starting_hydration_following_step_correlation':corr})
    return rows


def write_hydration_mobility(arrays,output_dir,prefix='crevice'):
    """Write ``PREFIX_hydration_mobility.csv`` (one row per residue; see :func:`hydration_mobility_statistics`).

    The definition and interpretation are recorded in the hydration manifest
    (``hydration_mobility``); frame intervals follow from ``time_ps`` in
    ``PREFIX_hydration_region_frames.csv``/``PREFIX_hydration_frames.csv``.

    Parameters
    ----------
    arrays : dict of str to numpy.ndarray
        Hydration trajectory arrays (see :func:`hydration_mobility_statistics`).
    output_dir : str or Path
        Existing output directory.
    prefix : str, default "crevice"
        File-name prefix.

    Returns
    -------
    dict of str to str
        ``{"hydration_mobility_csv": path}``, or empty when no row was computed
        (nothing is written then).
    """
    rows=hydration_mobility_statistics(arrays)
    if not rows:return {}
    root=Path(output_dir);path=root/(prefix+'_hydration_mobility.csv')
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return {'hydration_mobility_csv':str(path)}


MOBILITY_NOTES={
    'definition':'Displacement of the aligned selected-heavy-atom residue centroid between adjacent saved frames; hydration at both interval endpoints is observed',
    'interpretation':'Descriptive association without p-values or causal claims; within-interval hydration exchanges and internal sidechain motion are unresolved; conditional wet/dry means can have few observations'}


def hydration_frame_tables(arrays):
    """Per-frame hydration arrays as two long tables (lists of row dicts).

    Returns ``(residue_rows, region_rows)``. ``residue_rows`` has one row per
    frame and focused residue: ``frame_index``, ``time_ps``, ``residue`` and
    every per-residue series present (``water_count``,
    ``backbone_water_count``, ``sidechain_water_count``, ``NO_water_count``,
    ``sasa_A2``, ``backbone_sasa_A2``, ``sidechain_sasa_A2``, ``NO_sasa_A2``
    and, when a cavity trajectory was supplied, ``boundary_area_A2``,
    ``nonlining_partner``, ``centroid_distance_A`` and
    ``aligned_centroid_{x,y,z}_A``). ``region_rows`` has one row per frame:
    ``frame_index``, ``time_ps`` and the per-frame series present
    (``unique_focus_water_count``, ``region_water_count``,
    ``cavity_volume_A3``). Empty cells are unobserved values (NaN in the NPZ).
    """
    import numpy as np
    n=len(arrays['time_ps']);ids=[str(v) for v in arrays['residue_ids']];r=len(ids)
    frames=np.asarray(arrays.get('frame_indices',np.arange(n)))
    def value(v):
        v=v.item() if hasattr(v,'item') else v
        return '' if isinstance(v,float) and not np.isfinite(v) else v
    skip={'time_ps','frame_indices','residue_ids','contact_water_ids','contact_frame_residue_water_indices'}
    per_residue=[k for k,a in arrays.items() if k not in skip and np.ndim(a)==2 and np.shape(a)==(n,r)]
    xyz=[k for k,a in arrays.items() if k not in skip and np.ndim(a)==3 and np.shape(a)[:2]==(n,r) and np.shape(a)[2]==3]
    per_frame=[k for k,a in arrays.items() if k not in skip and np.ndim(a)==1 and len(a)==n]
    residue_rows=[]
    for i in range(n):
        for j,key in enumerate(ids):
            row={'frame_index':int(frames[i]),'time_ps':value(arrays['time_ps'][i]),'residue':key}
            row.update({k:value(arrays[k][i,j]) for k in per_residue})
            for k in xyz:
                stem={'aligned_residue_centroids_A':'aligned_centroid'}.get(k,k.removesuffix('_A'))
                row.update({f'{stem}_{axis}_A':value(arrays[k][i,j,m]) for m,axis in enumerate('xyz')})
            residue_rows.append(row)
    region_rows=[{'frame_index':int(frames[i]),'time_ps':value(arrays['time_ps'][i]),**{k:value(arrays[k][i]) for k in per_frame}} for i in range(n)]
    return residue_rows,region_rows
