"""Reusable reports from saved all-frame cavity arrays; no geometry reanalysis."""
from pathlib import Path

from .presentation import NON_CHANNEL_CAST_HEX, FigureText, annotate_option


def sectional_statistics(values, positions, **options):
    """Separate pointwise temporal estimates; no simultaneous-coverage claim.

    Each axial position is described independently with
    :func:`crevice.cavity_trajectory.describe_series`.

    Parameters
    ----------
    values : array_like, shape (frames, positions)
        Per-frame width (diameter, Å) at each axial position; NaN = missing.
    positions : sequence of float
        Axial positions, Å, one per column.
    **options
        Passed to :func:`~crevice.cavity_trajectory.describe_series`
        (``confidence``, ``block_length``, ``replicates``, ``seed``,
        ``regular``).

    Returns
    -------
    list of dict
        One row per position: ``position_A``, ``observed_frames``,
        ``coverage_fraction``, ``mean_diameter_A``, ``fluctuation_lower_A`` /
        ``fluctuation_upper_A`` (2.5/97.5% frame quantiles),
        ``pointwise_mean_lower_A`` / ``pointwise_mean_upper_A`` (clipped at 0),
        ``mean_CI_status``, ``statistical_inefficiency``, ``effective_frames`` and
        ``occupied_fraction_observed``.

    Raises
    ------
    ValueError
        If ``values`` is not 2D with one column per position.
    """
    import numpy as np
    from .cavity_trajectory import describe_series
    values=np.asarray(values,dtype=float)
    if values.ndim!=2 or values.shape[1]!=len(positions):
        raise ValueError('One profile column is required per coordinate')
    rows=[]
    for j,x in enumerate(positions):
        stat=describe_series(values[:,j],**options);ci=stat['confidence_interval']
        lo=ci.get('pointwise_lower',[None])[0];hi=ci.get('pointwise_upper',[None])[0]
        rows.append({'position_A':float(x),'observed_frames':stat['observed_frames'],
          'coverage_fraction':stat['coverage_fraction'],'mean_diameter_A':stat['mean'],
          'fluctuation_lower_A':stat['quantile_025'],'fluctuation_upper_A':stat['quantile_975'],
          'pointwise_mean_lower_A':max(0.,lo) if lo is not None else None,
          'pointwise_mean_upper_A':max(0.,hi) if hi is not None else None,
          'mean_CI_status':ci['status'],'statistical_inefficiency':stat['statistical_inefficiency'],
          'effective_frames':stat['effective_frames'],
          'occupied_fraction_observed':float(np.mean(values[np.isfinite(values[:,j]),j]>0)) if stat['observed_frames'] else None})
    return rows


@annotate_option
def write_additional_cavity_reports(report,arrays,output_dir,*,prefix='crevice',confidence=.95,
                                   block_length=None,replicates=2000,seed=20260913,dpi=240,annotate=None):
    """Write sectional-width, free-sphere, dynamic-partner and residue-heatmap reports.

    Everything is computed from the saved all-frame arrays of a cavity
    trajectory; no geometry is reanalysed. Two per-plane width measures are
    summarised along the fixed reference axis with :func:`sectional_statistics`:
    the area-equivalent radius of the largest face-connected section
    (``sqrt(A/pi)``, not an inscribed passage radius) and the largest local
    atom-clear sphere centred on a measured sample in that plane (the sphere may
    extend beyond the measured region; no transport path is implied).
    Nonlocal boundary-partner contacts are recalculated per frame and paired
    with aligned centroid-motion correlations.

    Parameters
    ----------
    report : dict
        Statistics from :func:`crevice.cavity_trajectory.summarize_cavity_trajectory`.
    arrays : dict of ndarray
        The matching all-frame arrays.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    confidence : float, default 0.95
        Nominal level of the pointwise mean intervals; 0 disables them.
    block_length : int, optional
        Bootstrap block length in frames.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260913
        Bootstrap random seed.
    dpi : int, default 240
        Figure resolution.
    annotate : bool, optional
        Draw panel titles, figure titles, the "mean intervals unavailable" note
        and the second line of the motion-correlation axis label. ``None``
        (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_cavity_section_statistics.csv, PREFIX_cavity_free_sphere_statistics.csv : tables
        Per axial position: observed frames, coverage, mean diameter, 2.5/97.5%
        frame quantiles, pointwise mean interval and its status (Å).
    PREFIX_cavity_section_radius_statistics.png, PREFIX_cavity_free_sphere_radius_statistics.png : figures
        Radius (Å) against "Position along reference axis (Å)" with a coverage
        panel ("Fraction").
    PREFIX_cavity_profile_statistics.json : JSON
        Definitions of both measures and of the bands.
    PREFIX_cavity_partner_statistics.csv, PREFIX_cavity_contact_statistics.csv : tables
        Per-residue partner frequency; per-contact occupancy with intervals.
    PREFIX_cavity_dynamic_partners.png : figure
        Top 20 nonlocal boundary-partner contacts: frequency given measured
        geometry, and aligned centroid motion correlation (-1 to 1).
    PREFIX_cavity_residue_heatmap.png : figure
        Assigned boundary area (Å²) of the top 20 residues in every frame.

    Colours and representation
    --------------------------
    * Width statistics (a non-channel cast): violet (``#a855f7``, the
      non-channel cast colour of :data:`crevice.presentation.NON_CHANNEL_CAST_HEX`)
      band (22%) = observed 2.5-97.5% frame range, violet line = observed mean, orange (``#d98c33``) band (50%) = pointwise
      mean interval; coverage panel blue-grey (``#547c9b``) = measured coverage,
      magenta (``#a34477``) = section present given measured. Legends name them.
    * Dynamic partners: magenta bars; blue-grey dots = whole-trajectory
      correlation, translucent blue-grey bars = first-half to second-half values;
      grey zero line.
    * Residue heatmap: ``cividis`` scale with colour bar "Assigned boundary area
      (Å²)"; light grey = unresolved frames.
    """
    import csv,json
    import numpy as np
    from .figures import _pyplot,_save_figure
    from .cavity_trajectory import describe_series
    root=Path(output_dir);root.mkdir(parents=True,exist_ok=True)
    if not prefix or Path(prefix).name!=prefix:raise ValueError('prefix must be a simple filename stem')
    files={};plt=_pyplot()
    def path(key,suffix):
        p=root/(prefix+suffix);files[key]=str(p);return p
    def table(p,rows):
        if not rows:return
        with p.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    options=dict(confidence=confidence,block_length=block_length,replicates=replicates,seed=seed,regular=report['regular_timestamps'])
    positions=np.array([r['position_A'] for r in report['profile']['rows']]);time=np.asarray(arrays['time_ps'])/1000
    profiles={}
    definitions={
      'section':('largest_component_diameter_A','Area-equivalent radius (Å)','Largest connected planar section',
                 'sqrt(A/pi) for the largest face-connected section; not an inscribed passage radius'),
      'free_sphere':('maximum_atom_clear_sphere_diameter_A','Local atom-clear sphere radius (Å)','Largest local sphere centre in each measured plane',
                 'Maximum atom-surface clearance at measured grid centres in a plane; sphere may extend beyond the measurement region; no connected transport path implied')}
    for name,(key,label,title,definition) in definitions.items():
        rows=sectional_statistics(arrays[key],positions,**options);profiles[name]={'definition':definition,'rows':rows}
        table(path(name+'_statistics_csv','_cavity_'+name+'_statistics.csv'),rows)
        def column(key):return np.array([r[key] if r[key] is not None else np.nan for r in rows])
        text=FigureText(title+' ('+definition+'), along the aligned reference axis: observed frame range, mean and pointwise mean interval, with coverage below')
        fig,axes=plt.subplots(2,1,figsize=(10,6.5),sharex=True,height_ratios=[4,1],layout='constrained');ax=axes[0]
        ax.fill_between(positions,.5*column('fluctuation_lower_A'),.5*column('fluctuation_upper_A'),color=NON_CHANNEL_CAST_HEX,alpha=.22,label='Observed 2.5–97.5% frame range')
        ax.plot(positions,.5*column('mean_diameter_A'),color=NON_CHANNEL_CAST_HEX,lw=1.8,label='Observed mean')
        lower=.5*column('pointwise_mean_lower_A');upper=.5*column('pointwise_mean_upper_A')
        if np.isfinite(lower).any():ax.fill_between(positions,lower,upper,color='#d98c33',alpha=.5,label=f'{confidence:.0%} pointwise mean interval, where estimable')
        else:text.text(ax,.01,.96,'Mean intervals unavailable under sampling safeguards',transform=ax.transAxes,va='top',fontsize=9)
        ax.set_ylabel(label);ax.set_ylim(bottom=0);text.title(ax,title+' · aligned reference axis',fontsize=12);ax.legend(frameon=False,fontsize=9,loc='lower left')
        axes[1].plot(positions,column('coverage_fraction'),color='#547c9b',label='Measured coverage')
        axes[1].plot(positions,column('occupied_fraction_observed'),color='#a34477',label='Section present | measured')
        axes[1].set_ylim(-.02,1.02);axes[1].set_ylabel('Fraction');axes[1].set_xlabel('Position along reference axis (Å)');axes[1].legend(frameon=False,fontsize=8,ncol=2)
        for a in axes:a.spines[['top','right']].set_visible(False)
        _save_figure(fig,path(name+'_radius_statistics_png','_cavity_'+name+'_radius_statistics.png'),dpi=dpi,text=text)
    metadata={'profiles':profiles,'empirical_band':'Fixed 2.5–97.5% quantiles of measured frames; not confidence bounds',
              'mean_band':f'{confidence:.0%} separate pointwise batch-bootstrap intervals; serial-dependence and coverage safeguards retained; not simultaneous confidence across coordinates',
              'statistical_options':options}
    path('profile_statistics_json','_cavity_profile_statistics.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    # Dynamic partner roles are recalculated in every frame, separate from the
    # reference-lining labels retained by the older motion report.
    partner_rows=[];ids=list(arrays['residue_ids']);motion={r['residue']:r for r in report['residue_dynamics']['residues']}
    for j,key in enumerate(ids):
        series=describe_series(arrays['nonlining_partner'][:,j],**options);ci=series['confidence_interval'];m=motion[key]
        partner_rows.append({'residue':str(key),'partner_frequency_observed':series['mean'],
          'observed_frames':series['observed_frames'],'partner_mean_lower':ci.get('pointwise_lower',[None])[0],
          'partner_mean_upper':ci.get('pointwise_upper',[None])[0],'mean_CI_status':ci['status'],
          'centroid_rmsf_A':m['centroid_rmsf_A']})
    for row in partner_rows:
        for key in ['partner_mean_lower','partner_mean_upper']:
            if row[key] is not None:row[key]=float(np.clip(row[key],0,1))
    table(path('partner_statistics_csv','_cavity_partner_statistics.csv'),partner_rows)
    contacts=[]
    for row in report['residue_dynamics']['contacts']:
        stat=row['occupancy_statistics'];ci=stat['confidence_interval']
        contacts.append({**{k:v for k,v in row.items() if not isinstance(v,(dict,list))},
          'occupancy_mean_lower':ci.get('pointwise_lower',[None])[0],'occupancy_mean_upper':ci.get('pointwise_upper',[None])[0],
          'occupancy_CI_status':ci['status'],'occupancy_effective_frames':stat['effective_frames']})
    table(path('contact_statistics_csv','_cavity_contact_statistics.csv'),contacts)
    selected=sorted((r for r in contacts if r['dynamic_nonlining_boundary_partner_frames']>0),
                    key=lambda r:(-r['dynamic_partner_frequency_observed'],r['source'],r['target']))[:20]
    text=FigureText(f"{len(selected)} most frequent nonlocal boundary-partner contacts: frequency given measured geometry "
                    "and aligned centroid motion correlation (dot whole trajectory, bar first to second half)")
    fig,axes=plt.subplots(1,2,figsize=(12,7),sharey=True,layout='constrained')
    labels=[r['source']+' ↔ '+r['target'] for r in selected]
    axes[0].barh(labels,[r['dynamic_partner_frequency_observed'] for r in selected],color='#a34478');axes[0].invert_yaxis();axes[0].set_xlim(0,1.02)
    axes[0].set_xlabel('Nonlocal boundary–partner frequency | measured')
    for i,row in enumerate(selected):
        if row['dccm'] is not None:axes[1].scatter(row['dccm'],i,color='#547c9b',s=30)
        if row['dccm_first_half'] is not None and row['dccm_second_half'] is not None:
            axes[1].plot([row['dccm_first_half'],row['dccm_second_half']],[i,i],color='#547c9b',alpha=.5,lw=3)
    axes[1].axvline(0,color='#aaaaaa',lw=.8);axes[1].set_xlim(-1.03,1.03);text.axis_label(axes[1].set_xlabel,'Aligned centroid motion correlation','Dot: whole trajectory · line: half-trajectory values')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    if not selected:text.notice(axes[0],.5,.5,'No observed nonlocal boundary–partner contacts',transform=axes[0].transAxes,ha='center',wrap=True)
    text.suptitle(fig,'Changing structural partners · proximity and motion association')
    _save_figure(fig,path('dynamic_partners_png','_cavity_dynamic_partners.png'),dpi=dpi,text=text)
    ranked=sorted(report['residues'],key=lambda r:(-(r['boundary_area_A2']['mean'] or 0),r['residue']))[:20]
    idx=[ids.index(r['residue']) for r in ranked];labels=[r['residue'] for r in ranked]
    text=FigureText(f"Assigned boundary area of the {len(idx)} residues with the largest mean boundary area, in every frame")
    fig,ax=plt.subplots(figsize=(11,7),layout='constrained')
    image=ax.pcolormesh(time,np.arange(len(idx)),np.asarray(arrays['boundary_area_A2'])[:,idx].T,shading='nearest',cmap='cividis')
    ax.set_yticks(np.arange(len(idx)),labels);ax.invert_yaxis();ax.set_xlabel('Time (ns)');ax.set_facecolor('#dddddd');text.title(ax,'Boundary area by residue across every frame · grey = unresolved')
    fig.colorbar(image,ax=ax,label='Assigned boundary area (Å²)')
    _save_figure(fig,path('residue_heatmap_png','_cavity_residue_heatmap.png'),dpi=dpi,text=text)
    return files
