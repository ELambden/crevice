"""Paired descriptive comparisons of different named regions in one trajectory."""
from pathlib import Path
import json

from .presentation import FigureText, annotate_option


def _load_result(path, definition_path=None, input_provenance=None):
    import numpy as np
    from .region_definition import load_definition,sha256
    path=Path(path).resolve(strict=True)
    ending='_cavity_statistics.json'
    if not path.name.endswith(ending):raise ValueError('Expected a cavity_statistics.json bundle')
    prefix=path.name[:-len(ending)];root=path.parent
    report=json.loads(path.read_text());reader=json.loads((root/(prefix+'_reader.json')).read_text())
    metadata=report['settings'].get('region_definition')
    if definition_path:
        definition,_=load_definition(definition_path)
        if 'prepared' not in definition:raise ValueError('Comparison requires prepared definitions')
        if metadata and metadata['definition_sha256']!=sha256(definition_path):raise ValueError('Comparison definition differs from the recorded region')
        definition_origin={'path':str(Path(definition_path).resolve()),'sha256':sha256(definition_path),
                           'binding':'explicitly checked definition for existing observations' if not metadata else 'recorded run definition'}
    elif metadata:
        definition=metadata['definition'];definition_origin={'sha256':metadata['definition_sha256'],'binding':'recorded run definition'}
    else:raise ValueError('Legacy observations require an explicit prepared region definition')
    if not metadata:
        if input_provenance is None:raise ValueError('Legacy observations require historical input provenance with recorded SHA256 digests')
        historical=json.loads(Path(input_provenance).read_text())
        if isinstance(historical.get('files'),list):historical={r['path']:r for r in historical['files']}
        for key in ['topology','trajectory']:
            source=str(Path(reader[key]).resolve())
            if historical.get(source,{}).get('sha256')!=definition['system'][key]['sha256']:
                raise ValueError('Historical input provenance differs for '+key)
        definition_origin['historical_input_provenance']={'path':str(Path(input_provenance).resolve()),'sha256':sha256(input_provenance)}
    system=definition['system'];geometry=definition['geometry'];settings=report['settings']
    for key in ['topology','trajectory']:
        if system[key] is None or sha256(reader[key])!=system[key]['sha256']:raise ValueError('Comparison source hash differs for '+key)
    from .cavity_trajectory import CavityReference
    reference_path=report['reference']['source_path']
    if sha256(reference_path)!=definition['reference']['asset']['sha256']:raise ValueError('Result reference map differs from the region definition')
    expected=CavityReference.from_dx(reference_path,spacing=geometry['spacing_A'],axis=geometry['axis'],grid_phase=geometry['grid_phase'],profile_padding=geometry['profile_padding_A'])
    conditions=[settings['geometry_mode']=='reference-region',settings['region_margin']==geometry['region_margin_A'],
      settings['spacing_A']==geometry['spacing_A'],settings['geometry']['probe_radius']==geometry['probe_radius_A'],
      settings['alignment_residues']==definition['alignment']['residues'],reader['selection']==system['protein_selection'],
      reader['periodic_boundary_handling']==system['pbc'],reader['frame_indices'][0]==system['reference_frame'],
      (settings.get('obstacles') or {}).get('selection')==system['obstacle_selection'],
      np.allclose(report['reference']['grid_anchor_A'],expected.anchor,rtol=0,atol=1e-7),
      np.allclose(report['reference']['grid_basis'],expected.basis,rtol=0,atol=1e-7)]
    if not all(conditions):raise ValueError('Existing result settings or coordinate frame differ from the region definition')
    with np.load(root/(prefix+'_cavity_frames.npz'),allow_pickle=False) as a:arrays={k:a[k].copy() for k in a.files}
    indices=np.asarray(list(range(*reader['frame_indices'])))
    if len(indices)!=report['frame_count'] or not np.array_equal(arrays['time_ps'],reader['frame_times_ps']):raise ValueError('Result frame indices or timestamps differ from its reader record')
    hydration=None;hp=root/(prefix+'_hydration.json')
    if hp.exists():
        hydration=json.loads(hp.read_text())
        if hydration['frame_count']!=report['frame_count']:raise ValueError('Hydration frame count differs from geometry')
        for key in ['topology','trajectory']:
            if Path(hydration['reader'][key]).resolve()!=Path(reader[key]).resolve():raise ValueError('Hydration source differs from geometry')
        with np.load(root/(prefix+'_hydration_frames.npz'),allow_pickle=False) as a:
            if not np.array_equal(a['time_ps'],arrays['time_ps']) or not np.array_equal(a['frame_indices'],indices):raise ValueError('Hydration observations have different frames/times')
            arrays['region_water_count']=a['region_water_count'].copy()
    profile_rows=report['profile']['rows'];profile_confidence=report['profile']['confidence_interval'].get('confidence',.95)
    pp=root/(prefix+'_cavity_profile_statistics.json')
    profile_extra=None
    if pp.exists():
        profile_extra=json.loads(pp.read_text());profile_rows=profile_extra['profiles']['section']['rows']
        profile_confidence=profile_extra['statistical_options']['confidence']
        positions=np.asarray([r['position_A'] for r in profile_rows])
        expected_positions=np.asarray([r['position_A'] for r in report['profile']['rows']])
        means=np.asarray([np.nan if r['mean_diameter_A'] is None else r['mean_diameter_A'] for r in profile_rows])
        if not np.array_equal(positions,expected_positions) or not np.allclose(means,np.nanmean(arrays['largest_component_diameter_A'],axis=0),rtol=0,atol=1e-8,equal_nan=True):
            raise ValueError('Pointwise profile observations differ from the geometry arrays')
    return {'path':path,'prefix':prefix,'root':root,'report':report,'reader':reader,'definition':definition,
            'definition_origin':definition_origin,'arrays':arrays,'indices':indices,'hydration':hydration,
            'profile_rows':profile_rows,'profile_confidence':profile_confidence,'profile_extra':profile_extra}


@annotate_option
def _recorded_radii(report):
    """The ``radii`` record of a cavity-statistics report, or ``None`` if it has none."""
    settings=report.get('settings') or {}
    return settings.get('radii') or (settings.get('region_definition') or {}).get('radii')


def compare_regions(left_path,right_path,output_dir,*,left_definition=None,right_definition=None,left_input_provenance=None,right_input_provenance=None,annotate=None,allow_radii_mismatch=False):
    """Compare two named regions measured on the same trajectory frames.

    Both results must come from ``region-trajectory`` (or be bound to prepared
    definitions) on the same topology and trajectory, with identical frames,
    timestamps, alignment, selections, periodic handling and geometry settings;
    otherwise a ``ValueError`` is raised. Paired per-frame differences of
    volume and fixed-neighbourhood water counts, and per-residue boundary
    statistics, are summarised descriptively. Different regions are different
    observables; no functional comparison is implied.

    Parameters
    ----------
    left_path, right_path : str or Path
        ``*_cavity_statistics.json`` of each region run.
    output_dir : str or Path
        New or empty output directory.
    left_definition, right_definition : str or Path, optional
        Prepared definitions for legacy results lacking region provenance.
    left_input_provenance, right_input_provenance : str or Path, optional
        Historical input SHA256 records required when binding legacy results.
    annotate : bool, optional
        Draw the region IDs as titles of the profile panels and the figure
        title. ``None`` (default) inherits the surrounding setting.
    allow_radii_mismatch : bool, default False
        Compare even though the two runs recorded different atomic radius sets
        (different table hashes). Otherwise that is a ``ValueError``: volumes
        measured with different radii are not comparable. Runs that recorded
        no radius set (default runs before 2 Oct 2026) cannot be checked; the
        comparison proceeds and says so in ``comparison.json["radii"]``.

    Returns
    -------
    dict
        The comparison summary written to ``comparison.json``.

    Outputs
    -------
    comparison.json : JSON
        Definitions, provenance hashes, ``radii`` (status and each run's radius set), paired statistics and limitations.
    comparison_summary.csv : table
        The interval statistics of ``comparison.json`` flattened, one row per
        quantity (each region's volume and fixed-neighbourhood water count and
        the paired right-minus-left differences); columns as in
        :func:`crevice.io.summary_statistics_rows`.
    residue_comparison.csv : table
        Per residue boundary statistics in each region.
    paired_observations.csv : table
        One row per paired frame: ``source_frame_indices``, ``time_ps``,
        ``left_volume_A3``, ``right_volume_A3`` and, when both runs measured
        hydration, ``left_region_water_count``/``right_region_water_count``.
    paired_observations.npz : arrays
        The same paired columns as arrays.
    comparison.png : figure (200 dpi)
        Four panels: "Regional volume (Å³)" and "Waters in fixed neighbourhood"
        against "Time (ns)" for both regions; for each region, "Radius (Å)"
        against "Section coordinate relative to this reference (Å)".

    Colours and representation
    --------------------------
    Left region brown (``#a45c39``), right region teal (``#167e8a``) in every
    panel; the volume panel's legend names the regions (region IDs). Profile
    panels: region-coloured mean section-equivalent radius line, the same
    colour at 18% for the observed 95% frame range, and orange (``#d98c33``,
    65%) for the pointwise mean confidence interval, with a legend.
    """
    import csv
    import numpy as np
    from .cavity_trajectory import describe_series
    from .region_definition import write_json,sha256
    from .figures import _pyplot,_save_figure
    left=_load_result(left_path,left_definition,left_input_provenance);right=_load_result(right_path,right_definition,right_input_provenance)
    root=Path(output_dir).resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Comparison requires an empty output directory')
    a,b=left['arrays'],right['arrays'];ld,rd=left['definition'],right['definition']
    if ld['region_id']==rd['region_id']:raise ValueError('Different regional observables require distinct region IDs')
    for key in ['topology','trajectory']:
        if ld['system'][key]['sha256']!=rd['system'][key]['sha256']:raise ValueError('Paired regions require the same source trajectory and topology')
    if not np.array_equal(left['indices'],right['indices']) or not np.array_equal(a['time_ps'],b['time_ps']):raise ValueError('Paired region comparisons require identical frames and timestamps')
    if ld['alignment']!=rd['alignment']:raise ValueError('Paired regions require the same declared alignment')
    for key in ['protein_selection','obstacle_selection','pbc']:
        if ld['system'][key]!=rd['system'][key]:raise ValueError('Region comparison selection/PBC definitions differ')
    for key in ['spacing_A','probe_radius_A','region_margin_A','axis','grid_phase']:
        if ld['geometry'][key]!=rd['geometry'][key]:raise ValueError('Region comparison geometry settings differ: '+key)
    if not np.array_equal(a['residue_ids'],b['residue_ids']):raise ValueError('Geometry residue identities differ')
    from .radii import radius_set_label,same_radius_set
    lr,rr=_recorded_radii(left['report']),_recorded_radii(right['report'])
    if lr is None or rr is None:radii_status='unverified_not_recorded'
    elif same_radius_set(lr,rr):radii_status='same_radius_set'
    elif allow_radii_mismatch:radii_status='mismatch_overridden'
    else:raise ValueError(f'Radius set mismatch: left run used {radius_set_label(lr)}, right run used {radius_set_label(rr)}; '
                          'volumes measured with different atomic radii are not comparable. Pass --allow-radii-mismatch to compare anyway (recorded).')
    short=lambda r:None if r is None else {'name':r.get('name'),'table_sha256':r.get('table_sha256')}
    root.mkdir(parents=True,exist_ok=True)
    regular=left['report']['regular_timestamps'] and right['report']['regular_timestamps']
    summary={'schema_version':1,'comparison':'right minus left, paired saved frames from one trajectory',
      'left':{'region_id':ld['region_id'],'label':ld['label'],'review':ld['review'],'definition':left['definition_origin'],'statistics_sha256':sha256(left['path'])},
      'right':{'region_id':rd['region_id'],'label':rd['label'],'review':rd['review'],'definition':right['definition_origin'],'statistics_sha256':sha256(right['path'])},
      'radii':{'status':radii_status,'left':short(lr),'right':short(rr)},
      'frame_count':len(a['time_ps']),'time_range_ps':[float(a['time_ps'][0]),float(a['time_ps'][-1])],
      'left_volume_A3':left['report']['volume_A3'],'right_volume_A3':right['report']['volume_A3'],
      'paired_volume_difference_A3':describe_series(b['volume_A3']-a['volume_A3'],regular=regular),
      'limitations':['The regions have different spatial definitions; a volume difference is not a gating, transport-state or functional effect.',
        'All-frame water counts belong to each fixed reference neighbourhood, not the instantaneous cavity interior.',
        'These paired observations are from one trajectory; they are not independent replicas.',
        'Missing residue hydration remains unavailable rather than zero. Candidate region status is retained.']}
    if left['hydration'] and right['hydration']:
        for key in ['contact_cutoff_A','sasa_probe_A','sasa_sphere_points','water_selection','periodic_distances']:
            if left['hydration']['settings'][key]!=right['hydration']['settings'][key]:raise ValueError('Hydration settings differ: '+key)
        summary['left_region_water_count']=left['hydration']['reference_region_water_count']
        summary['right_region_water_count']=right['hydration']['reference_region_water_count']
        summary['paired_region_water_difference']=describe_series(b['region_water_count']-a['region_water_count'],regular=regular)
    else:summary['hydration_comparison_status']='unavailable_missing_hydration_bundle'
    lining_a=np.any(a['boundary_area_A2']>0,axis=0);lining_b=np.any(b['boundary_area_A2']>0,axis=0)
    summary['boundary_residue_overlap']={'left_ever':int(lining_a.sum()),'right_ever':int(lining_b.sum()),
      'intersection':int((lining_a&lining_b).sum()),'union':int((lining_a|lining_b).sum())}
    rows=[]
    hydration=[{r['residue']:r for r in item['hydration']['residues']} if item['hydration'] else {} for item in [left,right]]
    for j,key in enumerate(a['residue_ids']):
        if not (lining_a[j] or lining_b[j] or key in hydration[0] or key in hydration[1]):continue
        row={'residue':str(key),'left_mean_boundary_A2':float(np.nanmean(a['boundary_area_A2'][:,j])),
             'right_mean_boundary_A2':float(np.nanmean(b['boundary_area_A2'][:,j]))}
        for name,arr,h in [('left',a,hydration[0]),('right',b,hydration[1])]:
            observed=np.isfinite(arr['boundary_area_A2'][:,j]);row[name+'_boundary_observed_frames']=int(observed.sum())
            row[name+'_boundary_frequency']=float(np.mean(arr['boundary_area_A2'][observed,j]>0)) if observed.any() else None
            row[name+'_hydration_occupancy']=h[key]['water_contact_occupancy']['mean'] if key in h else None
            row[name+'_water_count']=h[key]['statistics']['water_count']['mean'] if key in h else None
        rows.append(row)
    summary['residues']=rows
    def top(arr):
        order=np.argsort(-np.nanmean(arr['boundary_area_A2'],axis=0))
        return [str(arr['residue_ids'][i]) for i in order if np.nanmean(arr['boundary_area_A2'][:,i])>0][:10]
    summary['boundary_area_top10']={'left':top(a),'right':top(b),'interpretation':'Separate geometric boundary contributions; not a functional importance ranking'}
    interactions=[]
    for item in [left,right]:
        p=item['root']/(item['prefix']+'_interactions.json')
        if p.exists():
            d=json.loads(p.read_text());interactions.append({(e['kind'],*sorted([e['source'],e['target']])) for e in d['edges']})
        else:interactions.append(None)
    if all(i is not None for i in interactions):
        summary['typed_contact_scope']={'left_ever_edges':len(interactions[0]),'right_ever_edges':len(interactions[1]),'shared_ever_edges':len(interactions[0]&interactions[1]),
          'meaning':'Edges observed within each region-specific residue focus; absent from a focus does not mean zero physical contact occupancy'}
    summary['profile_uncertainty']={}
    for name,item in [('left',left),('right',right)]:
        summary['profile_uncertainty'][name]={'pointwise_confidence':item['profile_confidence'],
          'section_positions_with_mean_intervals':sum(r.get('pointwise_mean_lower_A',r.get('mean_lower_A')) is not None for r in item['profile_rows']),
          'whole_profile_status':item['report']['profile']['confidence_interval']['status'],
          'interpretation':'Separate pointwise estimates where supported; no simultaneous whole-profile coverage claim'}
    write_json(root/'comparison.json',summary)
    from .io import write_summary_statistics_csv
    write_summary_statistics_csv(summary,root/'comparison_summary.csv')
    with (root/'residue_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['residue']);w.writeheader();w.writerows(rows)
    arrays={'source_frame_indices':left['indices'],'time_ps':a['time_ps'],'left_volume_A3':a['volume_A3'],'right_volume_A3':b['volume_A3']}
    if 'region_water_count' in a and 'region_water_count' in b:arrays.update(left_region_water_count=a['region_water_count'],right_region_water_count=b['region_water_count'])
    np.savez_compressed(root/'paired_observations.npz',**arrays)
    with (root/'paired_observations.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(arrays));w.writeheader()
        w.writerows({k:('' if isinstance(v,float) and not np.isfinite(v) else v) for k,v in zip(arrays,values)}
                    for values in zip(*(np.asarray(arrays[k]).tolist() for k in arrays)))
    text=FigureText(f"Paired observations of regions {ld['region_id']} (brown) and {rd['region_id']} (teal) on the same frames: "
                    "volume and fixed-neighbourhood waters against time, and each region's section-equivalent radius profile")
    plt=_pyplot();fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');colors=['#a45c39','#167e8a'];time=a['time_ps']/1000
    for item,arr,color in zip([left,right],[a,b],colors):
        label=item['definition']['region_id'];axes[0,0].plot(time,arr['volume_A3'],lw=.8,alpha=.75,color=color,label=label)
        if 'region_water_count' in arr:axes[0,1].plot(time,arr['region_water_count'],lw=.7,alpha=.7,color=color,label=label)
    axes[0,0].set(xlabel='Time (ns)',ylabel='Regional volume (Å³)');axes[0,0].legend(fontsize=8)
    axes[0,1].set(xlabel='Time (ns)',ylabel='Waters in fixed neighbourhood')
    for ax,item,color in zip(axes[1],[left,right],colors):
        rows=item['profile_rows'];x=np.asarray([r['position_A'] for r in rows])
        ax.plot(x,[r['mean_diameter_A']/2 for r in rows],color=color,label='Mean section-equivalent radius')
        ax.fill_between(x,[r['fluctuation_lower_A']/2 for r in rows],[r['fluctuation_upper_A']/2 for r in rows],color=color,alpha=.18,label='Observed 95% frame range')
        low=np.array([np.nan if r.get('pointwise_mean_lower_A',r.get('mean_lower_A')) is None else r.get('pointwise_mean_lower_A',r.get('mean_lower_A'))/2 for r in rows]);high=np.array([np.nan if r.get('pointwise_mean_upper_A',r.get('mean_upper_A')) is None else r.get('pointwise_mean_upper_A',r.get('mean_upper_A'))/2 for r in rows])
        if np.isfinite(low).any():ax.fill_between(x,low,high,color='#d98c33',alpha=.65,label=f"{item['profile_confidence']:.0%} pointwise mean CI where estimable")
        ax.set(xlabel='Section coordinate relative to this reference (Å)',ylabel='Radius (Å)');text.title(ax,item['definition']['region_id']);ax.legend(fontsize=8)
    text.suptitle(fig,'Named regions · paired trajectory observations')
    _save_figure(fig,root/'comparison.png',dpi=200,text=text)
    return summary
