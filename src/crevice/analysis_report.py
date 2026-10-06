"""Coordinated residue, geometry and hydration evidence table and overview figures."""
from pathlib import Path
import csv,json,shutil

from .presentation import FigureText, annotate_option

EDGE_COLORS={'hydrogen_bond':'#d49c1f','water_bridge_hbond':'#0da896','salt_bridge_candidate':'#a64da6',
             'aromatic_parallel_candidate':'#4d73b3','aromatic_edge_face_candidate':'#4d73b3'}


@annotate_option
def write_analysis_report(output_dir,prefix,files,*,dpi=240,annotate=None):
    """Write the coordinated evidence table and overview figures.

    Reads the scene definition, hydration report, per-frame hydration arrays
    and, when present, typed interactions and water density written earlier to
    ``output_dir``. Residue evidence axes (hydration, typed contacts, boundary
    footprint, exposure, conformation) are kept separate; no composite score is
    formed. No HTML report is written: the values are in CSV tables (this
    evidence table,
    ``PREFIX_hydration_frames.csv``, ``PREFIX_hydration_region_frames.csv``,
    ``PREFIX_residue_conformation_frames.csv``,
    ``PREFIX_interaction_edge_frames.csv``) or in the hydration, interaction,
    density and scene JSON records written earlier.

    Parameters
    ----------
    output_dir : str or Path
        Directory holding ``PREFIX_analysis_scene.json``, ``PREFIX_hydration.json``
        and ``PREFIX_hydration_frames.npz``.
    prefix : str
        File-name prefix.
    files : dict of str to str
        Files written so far; used to find the optional inputs.
    dpi : int, default 240
        Figure resolution.
    annotate : bool, optional
        Draw panel titles, the figure title and the colour-scale title in the
        PNGs. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_analysis_evidence.csv : table
        One row per focused residue: roles, hydration occupancy and interval
        status, ``mean_SASA_A2``, ``boundary_area_A2``, typed-edge counts and
        occupancy sums, conformation statistics.
    PREFIX_analysis_overview.png : figure (14 x 11 inch)
        Four panels: water-contact occupancy of up to 18 residues ("Frames with
        ≥1 contact water", 0-100%); persistent typed contacts (occupancy ≥ 50%)
        as a circular network; cavity volume (Å³) and fixed-neighbourhood waters
        against time (ns); mean boundary footprint (Å²) against mean geometric
        solvent exposure (Å²).
    PREFIX_hydration_legend.png : figure
        The absolute hydration colour scale, "Water-contact occupancy" 0-100%.
    PREFIX_analysis_cavity_{radius,diameter}_profile.png, PREFIX_analysis_cavity_volume_timeseries.png : figures
        Copies of the cavity-trajectory figures, when the hydration run used one.

    Colours and representation
    --------------------------
    * Hydration colour (overview bars and nodes, legend): absolute linear scale
      from red (0% of frames with a contact water) through purple to blue
      (100%); grey = no water observations. Dark lines = pointwise mean intervals.
    * Typed-contact edges: gold (``#d49c1f``) hydrogen bond, teal (``#0da896``)
      water-bridged hydrogen bond, purple (``#a64da6``) salt-bridge candidate,
      blue (``#4d73b3``) aromatic candidates; width scales with occupancy; a
      legend names the kinds shown. Larger nodes are direct boundary residues;
      node labels are residue numbers.
    * Geometry panel: grey-blue volume line, teal fixed-neighbourhood water line
      on the right axis.
    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    from .analysis_viewers import hydration_rgb
    from .figures import _pyplot,_save_figure
    root=Path(output_dir);extra={}
    def path(key,suffix):
        p=root/(prefix+suffix);extra[key]=str(p);return p
    scene=json.loads((root/(prefix+'_analysis_scene.json')).read_text())
    hydration=json.loads((root/(prefix+'_hydration.json')).read_text())
    interaction_path=root/(prefix+'_interactions.json')
    interactions=json.loads(interaction_path.read_text()) if 'interaction_json' in files and interaction_path.exists() else {'edges':[],'residues':[],'metadata':{'status':'not_requested'}}
    density_path=root/(prefix+'_water_density.json')
    density=json.loads(density_path.read_text()) if 'water_density_json' in files and density_path.exists() else {'status':'not_requested'}
    with np.load(root/(prefix+'_hydration_frames.npz')) as loaded:arrays={k:loaded[k] for k in loaded.files}
    ids=list(arrays['residue_ids']);times=arrays['time_ps']/1000
    chemicals={r['residue']:r for r in interactions['residues']}
    hydration_rows={r['residue']:r for r in hydration['residues']}
    # Keep evidence axes separate; no unvalidated composite importance score.
    rows=[]
    for r in scene['residues']:
        edges=[e for e in interactions['edges'] if r['residue'] in [e['source'],e['target']]]
        row={k:v for k,v in r.items() if k not in {'atom_serials','color_rgb','centroid_A'}}
        row['observed_typed_edge_count']=len(edges)
        row['persistent_typed_edge_count']=sum(e['observed_occupancy']>=.5 for e in edges)
        row['hydration_occupancy_CI_status']=hydration_rows[r['residue']]['water_contact_occupancy']['confidence_interval']['status']
        for kind in EDGE_COLORS:
            row[kind+'_occupancy_sum']=sum(e['observed_occupancy'] for e in edges if e['kind']==kind)
        row.update({k:v for k,v in chemicals.get(r['residue'],{}).items() if k!='residue'})
        rows.append(row)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path('analysis_evidence_csv','_analysis_evidence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    context=hydration['settings'].get('context_statistics_json')
    context_coverage=None
    if context and Path(context).is_file():
        source=Path(context);context_report=json.loads(source.read_text())
        context_coverage=context_report.get('frame_count',context_report.get('summary',{}).get('frame_count'))
        stem=source.name.removesuffix('_cavity_statistics.json')
        for suffix,key in [('_cavity_radius_profile.png','geometry_radius_png'),('_cavity_diameter_profile.png','geometry_diameter_png'),('_cavity_volume_timeseries.png','geometry_volume_png')]:
            p=source.parent/(stem+suffix)
            if p.is_file():
                dest=path(key,'_analysis'+suffix)
                if p.resolve()!=dest.resolve():shutil.copy2(p,dest)
    plt=_pyplot();cmap=LinearSegmentedColormap.from_list('hydration',[(1,0,0),(0,0,1)])
    ranked=sorted(rows,key=lambda r:(r['residue'] not in scene['display_residues'],-(r['boundary_area_A2'] or 0),r['residue']))
    shown=ranked[:min(18,len(ranked))]
    text=FigureText(f"{prefix}: water-contact occupancy of up to 18 focused residues, persistent typed contacts (occupancy at least 50%), "
                    "cavity volume and fixed-neighbourhood waters against time, and boundary footprint against solvent exposure")
    fig,axes=plt.subplots(2,2,figsize=(14,11),layout='constrained')
    ax=axes[0,0];labels=[r['residue'] for r in shown];values=[r['hydration_occupancy'] for r in shown]
    ax.barh(labels,[0 if v is None else v for v in values],color=[hydration_rgb(v) for v in values]);ax.invert_yaxis();ax.set_xlim(0,1.02)
    ax.set_xticks([0,.5,1],['0%','50%','100%']);ax.set_xlabel('Frames with ≥1 contact water');text.title(ax,'Boundary contributors · hydration')
    if hydration['frame_count']==1:text.axis_label(ax.set_xlabel,'Observed snapshot contact','(not equilibrium probability)')
    if all(v is None for v in values):text.notice(ax,.1,.5,'No explicit waters supplied\nHydration is unobserved',transform=ax.transAxes,bbox={'facecolor':'white','edgecolor':'none'})
    for i,r in enumerate(shown):
        ci=hydration_rows[r['residue']]['water_contact_occupancy']['confidence_interval'];lo=ci.get('pointwise_lower',[None])[0];hi=ci.get('pointwise_upper',[None])[0]
        if lo is not None and hi is not None:ax.plot([lo,hi],[i,i],color='#222222',lw=1.4)
    ax=axes[0,1];edges=sorted(interactions['edges'],key=lambda e:-e['observed_occupancy'])
    direct=set(scene['display_residues']);selected=[e for e in edges if e['observed_occupancy']>=.5 and (e['source'] in direct or e['target'] in direct)][:24]
    nodes=sorted({k for e in selected for k in [e['source'],e['target']]});lookup={r['residue']:r for r in rows}
    positions={k:np.asarray([np.cos(2*np.pi*i/max(1,len(nodes))),np.sin(2*np.pi*i/max(1,len(nodes)))]) for i,k in enumerate(nodes)}
    for e in selected:
        a,b=positions[e['source']],positions[e['target']];ax.plot([a[0],b[0]],[a[1],b[1]],color=EDGE_COLORS[e['kind']],lw=.5+2*e['observed_occupancy'],alpha=.65)
    for key,point in positions.items():
        ax.scatter(*point,s=70 if key in direct else 40,color=[hydration_rgb(lookup[key]['hydration_occupancy'])],edgecolors='#444444',zorder=3)
        # Node names are this panel's category labels.
        ax.text(*(point*1.16),key.split(':')[-1],ha='center',va='center',fontsize=8)
    if not nodes:text.notice(ax,0,0,'No qualifying persistent typed edges\nSee chemistry coverage and complete edge table',ha='center')
    ax.set_xlim(-1.45,1.45);ax.set_ylim(-1.45,1.45);ax.set_aspect('equal');ax.axis('off');text.title(ax,'Persistent typed contacts · occupancy ≥50%')
    from matplotlib.lines import Line2D
    if selected:ax.legend(handles=[Line2D([0],[0],color=color,label=kind.replace('_',' ')) for kind,color in EDGE_COLORS.items() if any(e['kind']==kind for e in selected)],loc='lower center',bbox_to_anchor=(.5,-.08),fontsize=7,ncol=2)
    ax=axes[1,0]
    if 'cavity_volume_A3' in arrays:
        ax.plot(times,arrays['cavity_volume_A3'],color='#718790',lw=.8);ax.set_ylabel('Cavity volume (Å³)')
    else:text.notice(ax,.5,.7,'Cavity trajectory not supplied',ha='center',transform=ax.transAxes)
    ax.set_xlabel('Time (ns)' if len(times)>1 else 'Static snapshot');text.title(ax,'Geometry and regional hydration')
    if 'region_water_count' in arrays:
        other=ax.twinx();other.plot(times,arrays['region_water_count'],color='#0da896',lw=.7,alpha=.65);other.set_ylabel('Waters in fixed reference neighbourhood',color='#087c6f')
    ax=axes[1,1]
    points=[r for r in rows if r['boundary_area_A2'] is not None and r['mean_SASA_A2'] is not None]
    if points:
        ax.scatter([r['boundary_area_A2'] for r in points],[r['mean_SASA_A2'] for r in points],c=[hydration_rgb(r['hydration_occupancy']) for r in points],edgecolors='#444444',s=40,linewidths=.4)
        ax.set_xlabel('Mean cavity boundary footprint (Å²)');ax.set_ylabel('Mean geometric solvent exposure (Å²)')
    else:text.notice(ax,.5,.5,'Boundary-footprint history not supplied',ha='center',transform=ax.transAxes)
    text.title(ax,'Separate residue evidence axes')
    text.suptitle(fig,prefix+' · coordinated cavity, residue and hydration evidence',fontsize=16)
    _save_figure(fig,path('analysis_overview_png','_analysis_overview.png'),dpi=dpi,text=text)
    text=FigureText('Absolute hydration colour scale used in every CREVICE hydration view: red 0% to blue 100% of frames with a contact water; '
                    'residues without water observations are grey')
    # A colour bar keeps its quantity label; the title is annotation only.
    fig,ax=plt.subplots(figsize=(8,1.35),layout='constrained');ax.imshow(np.linspace(0,1,256)[None,:],aspect='auto',cmap=cmap,extent=[0,1,0,1]);ax.set_yticks([]);ax.set_xticks([0,.25,.5,.75,1],['0%','25%','50%','75%','100%']);ax.set_xlabel('Water-contact occupancy');text.title(ax,'Water-contact occupancy · unobserved residues are grey',fontsize=11)
    _save_figure(fig,path('hydration_legend_png','_hydration_legend.png'),dpi=dpi,text=text)
    return extra
