"""Stream solvent from original MD coordinates, preserving periodic water IDs.

For every analysed frame, each focus residue's water contacts and SASA
(see :mod:`crevice.hydration`) are measured from the original coordinates,
with periodic minimum-image distances and water identities taken from the
topology's residue index. Protein frames are superimposed on a reference
frame, so centroid positions, water densities and region counts are in one
aligned frame. Optional extensions add a water-oxygen density map and
chemically typed interaction occupancies.

Time-series summaries
    :func:`sampled_survival` gives the fraction of waters still present
    after a lag, and :func:`contact_episodes` lists continuous contact spans
    with explicit censoring. Both use saved frames only: exchanges between
    frames are invisible, so these are not continuous-time residence
    lifetimes.

Main entry point: :func:`analyze_hydration_trajectory` (command
``crevice hydration`` with a trajectory, or ``crevice region-trajectory``).
"""
from pathlib import Path
import json
from .radii import RadiusSet, radii_option


def sampled_survival(history,times,max_lag=20):
    """Mean fraction of starting waters observed at every intervening saved frame.

    For each lag ``L`` (0 to ``min(max_lag, n - 1)``) and each starting frame
    with at least one water, the fraction of its waters present in every frame
    from the start to ``start + L`` is computed; the rows report the mean over
    starting frames.

    Parameters
    ----------
    history : sequence of set
        Water IDs in contact at each frame (for one residue or region).
    times : sequence of float
        Frame times (ps); must be regularly spaced.
    max_lag : int, default 20
        Largest lag in frames.

    Returns
    -------
    dict
        ``status`` (``"complete"``, ``"unavailable_single_frame"`` or
        ``"unavailable_irregular_sampling"``) and ``rows`` with ``lag_frames``,
        ``lag_ps``, ``eligible_origins`` and ``survival_fraction`` (``None``
        when no start is eligible), plus the definition and limitation.

    Examples
    --------
    >>> from crevice.hydration_trajectory import sampled_survival
    >>> result = sampled_survival([{"w1", "w2"}, {"w1"}, {"w1"}], [0.0, 10.0, 20.0], max_lag=1)
    >>> [(row["lag_ps"], row["survival_fraction"]) for row in result["rows"]]
    [(0.0, 1.0), (10.0, 0.75)]
    """
    import numpy as np
    times=np.asarray(times,dtype=float);n=len(history)
    if n<2:return {'status':'unavailable_single_frame','rows':[]}
    if not np.allclose(np.diff(times),np.diff(times)[0],rtol=1e-7,atol=1e-7):return {'status':'unavailable_irregular_sampling','rows':[]}
    rows=[]
    for lag in range(min(int(max_lag),n-1)+1):
        values=[]
        for start in range(n-lag):
            initial=history[start]
            if not initial:continue
            keep=set(initial)
            for step in range(1,lag+1):
                keep.intersection_update(history[start+step])
                if not keep:break
            values.append(len(keep)/len(initial))
        rows.append({'lag_frames':lag,'lag_ps':float(lag*np.diff(times)[0]),'eligible_origins':len(values),
                     'survival_fraction':float(np.mean(values)) if values else None})
    return {'status':'complete','rows':rows,'intermittency_frames':0,
            'definition':'Mean over nonempty starting frames of the fraction of starting water IDs present in every intervening saved frame',
            'limitation':'Exchanges between saved frames are unresolved; this is not a continuous-time residence lifetime; lag-specific origin populations differ'}


def contact_episodes(histories,times,residue_ids):
    """Observed spans, with explicit left/right boundary censoring; no fitted lifetime.

    An episode is a maximal run of consecutive frames in which one water is in
    contact with one residue.

    Parameters
    ----------
    histories : sequence of sequence of set
        For each frame, one set of water IDs per residue.
    times : sequence of float
        Frame times (ps).
    residue_ids : sequence of str
        Residue labels, aligned with the inner sequences.

    Returns
    -------
    list of dict
        ``residue``, ``water_id``, ``start_row``/``end_row``,
        ``observed_frames``, ``start_ps``/``end_ps``, ``observed_span_ps``
        (end minus start time, so a single-frame episode spans 0 ps), and
        ``left_censored``/``right_censored`` (the episode touches the first or
        last analysed frame, so its true length is unknown). Sorted by residue,
        start row and water ID.

    Examples
    --------
    >>> from crevice.hydration_trajectory import contact_episodes
    >>> rows = contact_episodes([[{"w1"}], [{"w1"}], [set()]], [0.0, 10.0, 20.0], ["A:SER1"])
    >>> [(r["water_id"], r["observed_frames"], r["left_censored"], r["right_censored"]) for r in rows]
    [('w1', 2, True, False)]
    """
    active={};episodes=[]
    def close(key,start,end,right):
        i,water=key
        episodes.append({'residue':residue_ids[i],'water_id':water,'start_row':start,'end_row':end,
          'observed_frames':end-start+1,'start_ps':float(times[start]),'end_ps':float(times[end]),
          'observed_span_ps':float(times[end]-times[start]),'left_censored':start==0,'right_censored':right})
    for row,groups in enumerate(histories):
        current={(i,water) for i,ids in enumerate(groups) for water in ids}
        for key in list(active):
            if key not in current:close(key,active.pop(key),row-1,False)
        for key in current:active.setdefault(key,row)
    for key,start in active.items():close(key,start,len(histories)-1,True)
    return sorted(episodes,key=lambda r:(r['residue'],r['start_row'],r['water_id']))


def _environment_with_radii(atoms,elements,indices,radii=None,source=''):
    """Environment atom indices that have a radius; warn about and drop the rest.

    Non-water heavy atoms of the whole system (lipids, ions, ligands) occlude
    SASA. Those without a radius in the set in effect
    (:meth:`crevice.radii.RadiusSet.has_radius`) are left out with an
    :class:`~crevice.radii.UnrecognisedElementWarning`.
    """
    import warnings
    import numpy as np
    from collections import Counter
    from .radii import UnrecognisedElementWarning,effective_radii,unrecognised_element_message
    active=effective_radii(radii)
    names=np.asarray(atoms.names)[indices];resnames=np.asarray(atoms.resnames)[indices]
    known=np.asarray([active.has_radius(str(n),str(r),str(e)) for n,r,e in zip(names,resnames,elements[indices])],dtype=bool)
    if not known.all():
        by_element=dict(sorted(Counter(str(e) for e in elements[indices][~known]).items()))
        by_name=dict(sorted(Counter(f'{r}:{n}' for n,r in zip(names[~known],resnames[~known])).items()))
        warnings.warn(f'{source} (hydration environment): '+unrecognised_element_message(by_element,by_name,radii=active),
                      UnrecognisedElementWarning,stacklevel=3)
    return np.asarray(indices)[known]


def _environment_radii(atoms,elements,indices,radii=None):
    """SASA obstacle radii of the full-system environment atoms, in Å.

    Each atom gets :meth:`crevice.radii.RadiusSet.radius_for_names` of the set
    in effect (:func:`crevice.radii.effective_radii`; by default
    :data:`crevice.radii.DEFAULT_RADII`, so recognised ion residues such as
    CHARMM ``SOD``/``CLA`` get their CHARMM36 radii and every other atom its
    standard-table element radius, default 1.70 Å). ``radii`` (a
    :class:`~crevice.radii.RadiusSet`) overrides the set in effect.
    """
    import numpy as np
    from .radii import effective_radii
    active=effective_radii(radii)
    names=np.asarray(atoms.names)[indices];resnames=np.asarray(atoms.resnames)[indices]
    return np.asarray([active.radius_for_names(str(n),str(r),str(e)) for n,r,e in zip(names,resnames,elements[indices])],dtype=float)


@radii_option
def analyze_hydration_trajectory(topology,trajectory,*,context_json=None,residue_ids=None,
                                  selection='protein and not name H*',water_selection=None,
                                  start=0,stop=None,stride=1,max_frames=500,pbc='check',
                                  cutoff=3.5,sasa_probe=1.4,sasa_samples=256,alignment_residues=None,
                                  progress_callback=None,align=True,analysis_options=None,
                                  radii: RadiusSet | str | None = None,allow_radii_mismatch=False):
    """Measure per-residue hydration over an MD trajectory.

    Parameters
    ----------
    topology, trajectory : str or pathlib.Path
        MD files readable by MDAnalysis.
    context_json : str or pathlib.Path, optional
        ``*_cavity_statistics.json`` from ``crevice cavity-trajectory``. When
        given, its sibling reader, frame and provenance files are checked
        against these inputs (paths, frame indices and times, and SHA-256
        hashes when recorded), and its frames, atom selection, alignment
        residues and residue set (every residue that ever lined the cavity or
        was a partner) are adopted. With reference-region geometry, waters in
        the fixed reference neighbourhood are also counted. Requires
        ``align=True``.
    residue_ids : sequence of str, optional
        Residue labels to measure; all selected protein residues by default.
    selection : str, default "protein and not name H*"
        MDAnalysis selection of protein heavy atoms.
    water_selection : str, optional
        MDAnalysis selection of water; defaults to
        :data:`crevice.hydration.DEFAULT_WATER_SELECTION`. Exactly one oxygen
        per water residue is required.
    start, stop, stride, max_frames, pbc
        Frame range and periodic handling, as for
        :func:`crevice.trajectory.load_trajectory_report` (``max_frames``
        default 500, ``pbc`` default ``"check"``; ``"none"`` disables periodic
        distances).
    cutoff : float, default 3.5
        Water-oxygen contact distance in Å.
    sasa_probe : float, default 1.4
        SASA probe radius in Å.
    sasa_samples : int, default 256
        Points per atom sphere.
    alignment_residues : sequence of str, optional
        Residues used for the rigid fit.
    progress_callback : callable, optional
        Called as ``progress_callback(row, n_frames)`` after each frame.
    align : bool, default True
        Superimpose each frame on the reference (the first analysed frame, or
        the context's first frame).
    analysis_options : dict, optional
        Extension settings (see :func:`crevice.hydration_extensions.defaults`):
        water-density grid spacing, display smoothing and isovalue, and flags to
        skip the density map, the interaction analysis or the viewer outputs.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (SASA environment radii) comes from that set; it is
        recorded as ``settings["radii"]``. With ``context_json`` and no set
        chosen explicitly, the set recorded in the cavity results is used
        instead; see :func:`crevice.radii.reconcile_recorded_radii` and
        :doc:`/methods/atomic-radii`.
    allow_radii_mismatch : bool, default False
        With ``context_json``: measure with an explicitly chosen set even if it
        differs from the one the cavity results recorded (otherwise
        ``ValueError``). The outcome is recorded as ``settings["radii_check"]``.

    Returns
    -------
    dict
        ``arrays``: per-frame, per-residue arrays for every metric in
        :data:`crevice.hydration.METRICS`, ``time_ps``, ``residue_ids``,
        ``frame_indices``, ``region_water_count`` (NaN unless measured),
        ``centroid_distance_A`` (residue centroid to the region anchor or
        reference centroid) and, when aligned, ``aligned_residue_centroids_A``;
        context cavity arrays when a context is given. Also ``histories`` and
        ``bridge_histories`` (per-frame memberships), ``alignment`` reports,
        ``boxes``, the ``reader`` report, ``water_population``, ``settings``,
        optional ``water_density`` and ``interactions`` results, and reference
        frames/coordinates used by the viewer writers.

    Raises
    ------
    ValueError
        For inconsistent context files or hashes, no matching frames, times
        that are not finite and increasing, unknown focus residues, a water
        selection without exactly one oxygen per water, non-protein atoms in the
        protein selection, protein coordinates that differ from the solvent
        source, or a missing box with periodic handling.

    Notes
    -----
    Water IDs are topology residue indices (``"water-resindex:<i>"``), which
    stay fixed however molecules are wrapped. Protein coordinates must already
    be whole across the periodic boundary.
    """
    import numpy as np
    from .trajectory import load_trajectory_report,_open_universe,_element_names,align_frame_report
    from .hydration import HydrationTopology,DEFAULT_WATER_SELECTION,minimum_vectors,METRICS
    context=None;context_arrays=None;context_reader=None;context_verification='not_applicable'
    if context_json is not None:
        if not align:raise ValueError('Cached cavity hydration requires alignment to its reference')
        path=Path(context_json).resolve(strict=True);context=json.loads(path.read_text())
        stem=path.name.removesuffix('_cavity_statistics.json')
        context_reader=json.loads((path.parent/(stem+'_reader.json')).read_text())
        for key,current in [('topology',topology),('trajectory',trajectory)]:
            if Path(context_reader[key]).resolve()!=Path(current).resolve():raise ValueError('Hydration context source differs from the requested '+key)
        provenance=path.parent.parent/'input-provenance.json'
        context_verification='source paths, selected identities, frame indices and timestamps; content hashes not supplied'
        if provenance.is_file():
            import hashlib
            hashes=json.loads(provenance.read_text());verified=0
            for current in [topology,trajectory]:
                key=str(Path(current).resolve())
                if key in hashes:
                    h=hashlib.sha256()
                    with Path(current).open('rb') as stream:
                        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
                    if h.hexdigest()!=hashes[key]['sha256']:raise ValueError('Trajectory hydration source content differs from the files recorded in the cavity-trajectory results')
                    verified+=1
            if verified==2:context_verification='topology and trajectory SHA256 verified against geometry checkpoint'
        region=context['settings'].get('region_definition')
        if region:
            from .region_definition import sha256
            for key,current in [('topology',topology),('trajectory',trajectory)]:
                if sha256(current)!=region['definition']['system'][key]['sha256']:
                    raise ValueError('Hydration source SHA256 differs from the named region geometry')
            context_verification='topology and trajectory SHA256 verified against named region definition; identities, frames and timestamps checked'
        context_indices=list(range(*context_reader['frame_indices']))
        chosen=[i for i,f in enumerate(context_indices) if f>=start and (stop is None or f<stop) and (f-start)%stride==0]
        if not chosen:raise ValueError('No context frames match the requested range')
        chosen_frames=[context_indices[i] for i in chosen]
        start,stop,stride=chosen_frames[0],chosen_frames[-1]+1,chosen_frames[1]-chosen_frames[0] if len(chosen_frames)>1 else 1
        selection=context_reader['selection']
        alignment_residues=context['settings'].get('alignment_residues')
        with np.load(path.parent/(stem+'_cavity_frames.npz')) as a:context_arrays={k:a[k].copy() for k in ['time_ps','volume_A3','residue_ids','boundary_area_A2','nonlining_partner']}
        for key in ['time_ps','volume_A3','boundary_area_A2','nonlining_partner']:context_arrays[key]=context_arrays[key][chosen]
        if residue_ids is None:residue_ids=[r['residue'] for r in context['residues'] if r['lining_frames'] or r['nonlining_partner_frames']]
    from .radii import active_radii,effective_radii,reconcile_recorded_radii
    if context is not None:
        recorded=context['settings'].get('radii') or (context['settings'].get('region_definition') or {}).get('radii')
        radius_set,radii_check=reconcile_recorded_radii(recorded,active_radii(),allow_mismatch=allow_radii_mismatch,
                                                         what='the cavity results (--cavity-results)')
    else:
        radius_set,radii_check=effective_radii(),{'status':'no_context'}
    protein,reader=load_trajectory_report(topology,trajectory,selection=selection,start=start,stop=stop,stride=stride,max_frames=max_frames,pbc=pbc)
    times=np.asarray(protein.times(),dtype=float)
    if not np.isfinite(times).all() or (len(times)>1 and np.any(np.diff(times)<=0)):raise ValueError('Hydration trajectories require actual finite increasing times')
    if context_arrays is not None and not np.array_equal(times,context_arrays['time_ps']):raise ValueError('Hydration and cavity frames/times differ')
    first=protein.frames[0]
    reference_frame=first
    if context_reader is not None and first.frame_index!=context_reader['frame_indices'][0]:
        reference_index=context_reader['frame_indices'][0]
        reference_frame=load_trajectory_report(topology,trajectory,selection=selection,start=reference_index,stop=reference_index+1,max_frames=1,pbc=pbc)[0].frames[0]
    known={a.residue_key.label for a in first.atoms};focus=known if residue_ids is None else set(residue_ids)
    if not focus or focus-known:raise ValueError('Hydration focus must contain valid selected protein residue identities')
    universe=_open_universe(topology,trajectory);atoms=universe.atoms
    selected=universe.select_atoms(selection);elements=np.asarray(_element_names(atoms))
    water_selection=water_selection or DEFAULT_WATER_SELECTION
    water_group=universe.select_atoms(water_selection)
    oxygen_indices=water_group.indices[elements[water_group.indices]=='O']
    water_resindices=np.asarray(atoms.resindices)[oxygen_indices]
    if len(set(water_resindices))!=len(water_resindices):raise ValueError('Water selection must identify exactly one oxygen per water residue')
    water_residue_set=set(np.asarray(atoms.resindices)[water_group.indices])
    env_mask=~np.isin(elements,['H','D']) & ~np.isin(atoms.resindices,list(water_residue_set))
    env_indices=_environment_with_radii(atoms,elements,np.flatnonzero(env_mask),radius_set,str(topology))
    target_rows=np.asarray([i for i,a in enumerate(first.atoms) if a.residue_key.label in focus],dtype=int)
    if any(first.atoms[i].hetero or first.atoms[i].element.upper() in {'H','D'} for i in target_rows):raise ValueError('Hydration protein selection must contain only protein heavy atoms')
    env_lookup={idx:i for i,idx in enumerate(env_indices)}
    model=HydrationTopology(tuple(first.atoms[i] for i in target_rows),_environment_radii(atoms,elements,env_indices,radius_set),
              np.asarray([env_lookup[int(selected.indices[i])] for i in target_rows]),tuple('water-resindex:'+str(i) for i in water_resindices))
    reference_points=None;anchor=np.mean([a.coord for a in reference_frame.atoms],axis=0);margin=None
    if context is not None and context['settings'].get('geometry_mode')=='reference-region':
        refpath=context.get('reference',{}).get('source_path')
        if refpath and Path(refpath).is_file():
            from .cavity_trajectory import CavityReference
            from scipy.spatial import cKDTree
            cavity_reference=CavityReference.from_dx(refpath)
            expected=context.get('reference',{}).get('source_sha256')
            if expected and cavity_reference.source_sha256!=expected:raise ValueError('Reference cavity map hash differs from the geometry analysis')
            reference_points=cavity_reference.points;tree=cKDTree(reference_points)
            anchor=np.asarray(context['profile']['axis_origin']);margin=context['settings']['region_margin']
            bound=np.linalg.norm(reference_points-anchor,axis=1).max()+margin
    n=len(protein.frames);r=len(model.residue_ids)
    arrays={key:np.full((n,r),np.nan) for key in METRICS}
    arrays.update(time_ps=times,residue_ids=np.asarray(model.residue_ids),frame_indices=np.asarray([f.frame_index for f in protein.frames]),
                  region_water_count=np.full(n,np.nan),centroid_distance_A=np.zeros((n,r)))
    positions=np.zeros((n,r,3));histories=[];bridge_histories=[];fits=[];boxes=[]
    from .hydration_extensions import defaults,reference_points as density_reference_points
    from .water_density import WaterDensity
    from .interaction_analysis import chemical_topology_from_universe,InteractionAnalysis
    extension_options=defaults(analysis_options)
    volume_path=context.get('reference',{}).get('source_path') if context is not None else None
    density=None;interactions=None
    if align and not extension_options['skip_density']:
        points,definition=density_reference_points(reference_frame,model.residue_ids,volume_path)
        density=WaterDensity(points,margin=margin if margin is not None else (2. if volume_path else cutoff),
            spacing=extension_options['density_spacing'],smoothing=extension_options['density_smoothing'],
            expected_frames=n,max_points=extension_options['density_max_points'],definition=definition)
    universe.trajectory[first.frame_index]
    if not extension_options['skip_interactions']:
        chemistry=chemical_topology_from_universe(universe,selected,first,oxygen_indices)
        interactions=InteractionAnalysis(chemistry,model.residue_ids,n)
    snapshots=[];snapshot_times=[];view_reference=None;view_coordinates=None;view_rotation=None;view_box=None

    for row,frame in enumerate(protein.frames):
        universe.trajectory[frame.frame_index]
        if not np.isclose(float(universe.trajectory.ts.time),times[row],rtol=1e-9,atol=1e-7):raise ValueError('Solvent/protein frame times differ')
        if not np.allclose(selected.positions,np.asarray([a.coord for a in frame.atoms]),atol=1e-5,rtol=0):raise ValueError('Solvent source and selected protein coordinates differ; unwrap needs a consistently reimaged input')
        if align:aligned_frame,fit=align_frame_report(frame,reference_frame,residues=alignment_residues)
        else:
            aligned_frame=frame
            fit={'status':'not_aligned','rotation_rows':np.eye(3).tolist(),'translation_A':[0.,0.,0.]}
        rotation=np.asarray(fit['rotation_rows']);translation=np.asarray(fit['translation_A'])
        box=None if pbc=='none' else universe.dimensions
        if pbc!='none' and box is None:raise ValueError('Periodic hydration needs cell dimensions; explicitly use pbc=none for nonperiodic data')
        coords=selected.positions[target_rows].astype(float);waters=atoms.positions[oxygen_indices].astype(float)
        measured=model.measure(coords,atoms.positions[env_indices],waters,box=box,orientation=rotation,cutoff=cutoff,sasa_probe=sasa_probe,sasa_samples=sasa_samples)
        for key in METRICS:arrays[key][row]=measured['metrics'][key]
        histories.append(measured['memberships']);bridge_histories.append(measured['bridges']);fits.append(fit);boxes.append(box.tolist() if box is not None else None)
        aligned=coords@rotation+translation
        for j in range(r):positions[row,j]=aligned[model.atom_residue==j].mean(axis=0)
        arrays['centroid_distance_A'][row]=np.linalg.norm(positions[row]-anchor,axis=1)
        if reference_points is not None and len(waters):
            from .water_membership import image_waters_around_anchor,fixed_reference_water_mask
            # Shared with the instantaneous-membership pass; the definition is unchanged.
            relative=image_waters_around_anchor(waters,box,rotation,translation,anchor)
            arrays['region_water_count'][row]=np.sum(fixed_reference_water_mask(relative,anchor,tree,margin,bound))
        if density is not None:density.add(waters,rotation=rotation,translation=translation,box=box)
        if interactions is not None:interactions.add(atoms.positions,box=box)
        if row==0:
            view_reference=aligned_frame
            view_coordinates=atoms.positions.astype(float)@rotation+translation
            view_rotation=rotation.copy();view_box=box.copy() if box is not None else None
        if row in {0,n//2,n-1}:
            snapshots.append(aligned_frame);snapshot_times.append(float(times[row]))
        if progress_callback:progress_callback(row,n)
    if align:arrays['aligned_residue_centroids_A']=positions
    if context_arrays is not None:
        indices=[list(context_arrays['residue_ids']).index(key) for key in model.residue_ids]
        arrays.update(cavity_volume_A3=context_arrays['volume_A3'],boundary_area_A2=context_arrays['boundary_area_A2'][:,indices],
                      nonlining_partner=context_arrays['nonlining_partner'][:,indices])
    result={'arrays':arrays,'histories':histories,'bridge_histories':bridge_histories,'alignment':fits,'boxes':boxes,
            'reader':reader,'water_population':len(oxygen_indices),'settings':{'water_selection':water_selection,'water_identity':'stable MD topology residue index; one oxygen per residue',
              'contact_cutoff_A':cutoff,'sasa_probe_A':sasa_probe,'sasa_sphere_points':sasa_samples,'environment_heavy_atoms':len(env_indices),
              'environment_policy':'all original nonwater heavy atoms, including lipids, ligands and ions',
              'periodic_distances':pbc!='none','alignment_applied':bool(align),'context_statistics_json':str(context_json) if context_json else None,
              'context_verification':context_verification,'alignment_reference_source_index':reference_frame.frame_index,
              'region_definition':context['settings'].get('region_definition') if context is not None else None,
              'focus':'all selected protein' if residue_ids is None else 'fixed requested/ever-boundary/ever-partner residue set',
              'region_water_definition':'oxygen positions inside the fixed reference neighbourhood; not membership in a recomputed instantaneous probe cast' if reference_points is not None else 'not measured',
              'region_margin_A':margin,'region_anchor_A':anchor.tolist(),
              'radii':radius_set.provenance(),'radii_check':radii_check}}

    if density is not None:result['water_density']=density.result(len(oxygen_indices))
    if interactions is not None:result['interactions']=interactions.result()
    result.update(reference_frame=view_reference,reference_atom_indices=selected.indices.copy(),
        reference_full_coordinates=view_coordinates,reference_rotation=view_rotation,reference_box=view_box,
        snapshot_frames=snapshots,snapshot_times_ps=snapshot_times,analysis_options=extension_options,
        reference_volume_path=volume_path)
    return result
