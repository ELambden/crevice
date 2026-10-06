"""Synthetic trajectory geometry/statistics contracts; no biological assignment."""
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from crevice.models import ChannelPoint,VoidComponent,VoidCast
from crevice.cavity_trajectory import CavityReference,cavity_section_profile,match_reference_component,describe_series,analyze_cavity_frame
from crevice.rolling import rolling_probe_cast
from crevice.volume_export import write_void_cast_dx
from test_channel_coordinate import rectangular_channel,transform


def component(xyz):
    points=tuple(ChannelPoint(i,tuple(map(float,p)),2.,2.,float(p[2])) for i,p in enumerate(xyz))
    return VoidComponent(1,'test',(0,0,0),len(points),points)


def test_planar_area_width_separates_disconnected_lobes_and_empty_sections():
    xyz=np.array([(x,y,z) for x in [0,1,5,6] for y in [0,1] for z in [0,1]])
    ref=CavityReference(xyz,np.eye(3),np.zeros(3),1.,np.arange(-1,3),None,None,None)
    result=cavity_section_profile(component(xyz),ref)
    assert result['total_area_A2'].tolist()==[0,8,8,0]
    assert result['largest_area_A2'].tolist()==[0,4,4,0]
    assert result['section_component_count'].tolist()==[0,2,2,0]
    assert result['largest_component_diameter_A'][1]==pytest.approx(4/np.sqrt(np.pi))
    assert result['maximum_atom_clear_sphere_diameter_A'][1]==4


def test_reference_matching_rejects_far_and_ambiguous_components():
    xyz=np.array([(x,0,0) for x in range(5)]);ref=CavityReference(xyz,np.eye(3),np.zeros(3),1.,np.arange(5),None,None,None)
    c=component(xyz);far=component(xyz+100)
    cast=VoidCast((far,),far.points,1.,0.,'rolling')
    assert match_reference_component(cast,ref)[1]['status']=='unresolved_no_reference_match'
    cast=replace(cast,components=(c,replace(c,id=2)))
    assert match_reference_component(cast,ref)[1]['status']=='unresolved_ambiguous_reference_match'
    assert match_reference_component(replace(cast,components=(c,)),ref)[0]==c


def test_missing_geometry_is_not_zero_and_constant_occupancy_has_no_false_certainty():
    stats=describe_series([10.,np.nan,20.],replicates=200)
    assert stats['mean']==15 and stats['coverage_fraction']==pytest.approx(2/3)
    assert stats['confidence_interval']['status']=='unavailable_missing_frames'
    assert describe_series(np.zeros(100),replicates=200)['confidence_interval']['status']=='unavailable_observed_constant'
    assert describe_series(np.arange(100),regular=False,replicates=200)['confidence_interval']['status']=='unavailable_irregular_timestamps'


def test_correlated_volume_series_reports_effective_sampling():
    rng=np.random.default_rng(92);x=np.zeros(2000)
    for i in range(1,len(x)):x[i]=.8*x[i-1]+rng.normal()
    stat=describe_series(x+100,replicates=200)
    assert stat['statistical_inefficiency']>3
    assert stat['effective_frames']<len(x)/3
    assert stat['sample_sd']>1


def test_fixed_lattice_frame_analysis_is_rigid_motion_invariant(tmp_path):
    frame=rectangular_channel(half_x=4,half_y=5,half_length=8,capped=True)
    config=dict(probe_radius=.8,outer_radius=6.,selection_mode='all',enclosure_fraction=0,min_component_volume=1)
    cast=rolling_probe_cast(frame,spacing=.5,grid_basis=np.eye(3),grid_anchor=(0,0,0),**config)
    main=max(cast.components,key=lambda c:c.volume);cast=replace(cast,components=(main,),points=main.points)
    dx=write_void_cast_dx(cast,tmp_path/'reference.dx')
    ref=CavityReference.from_dx(dx,axis=(0,0,1));a=analyze_cavity_frame(frame,frame,ref,geometry=config)
    moved,_,_=transform(frame);b=analyze_cavity_frame(moved,frame,ref,geometry=config)
    assert a['status']==b['status']=='matched'
    assert a['volume_A3']==b['volume_A3']
    assert b['alignment']['rmsd_after_A']<1e-10
    assert a['profile']['largest_component_diameter_A']==pytest.approx(b['profile']['largest_component_diameter_A'])
    assert [(r['residue'],r['boundary_area_A2']) for r in a['boundary']['residues']]==[(r['residue'],r['boundary_area_A2']) for r in b['boundary']['residues']]


def test_reference_region_clips_volume_and_handles_empty_probe_space(tmp_path):
    from crevice.regional_volume import regional_probe_cast
    from crevice import Atom
    from crevice.spatial import SpatialIndex
    from scipy.spatial import cKDTree
    frame=rectangular_channel(half_x=4,half_y=5,half_length=8,capped=True)
    cast=rolling_probe_cast(frame,spacing=.5,selection_mode='all',enclosure_fraction=0,probe_radius=.8,outer_radius=6,min_component_volume=1)
    main=max(cast.components,key=lambda c:c.volume);cast=replace(cast,components=(main,),points=main.points)
    ref=CavityReference.from_dx(write_void_cast_dx(cast,tmp_path/'ref.dx'),axis=(0,0,1))
    volume,region,report=regional_probe_cast(frame,ref,probe_radius=1.4,margin=1.)
    assert region.volume>0
    points=np.array([p.position for p in region.points])
    assert cKDTree(ref.points).query(points)[0].max()<=1.+1e-8
    spatial=SpatialIndex(frame.atoms)
    assert min(spatial.nearest_surface(p)[0] for p in points[::20])>=-1e-9
    fill=tuple(Atom(10000+i,'C','ALA','B',i,*map(float,p),'C') for i,p in enumerate(ref.points[::2]))
    blocked=replace(frame,atoms=frame.atoms+fill)
    _,empty,diagnostic=regional_probe_cast(blocked,ref,probe_radius=2.5,margin=0)
    assert empty.volume==0 and not empty.points
    assert diagnostic['sampled_volume_components']==0


def test_crop_faces_are_excluded_from_residue_attribution():
    from crevice import Atom,StructureFrame
    from crevice.residue_evidence import boundary_residue_evidence
    grid=np.zeros((5,5,5));grid[1:4,1:4,1:4]=1
    frame=StructureFrame((Atom(1,'CA','ALA','A',1,0,2,2,'C'),Atom(2,'CA','ALA','A',2,4,2,2,'C')))
    ordinary=boundary_residue_evidence(frame,grid,np.zeros(3),np.eye(3))
    excluded=boundary_residue_evidence(frame,grid,np.zeros(3),np.eye(3),boundary_exclusion=lambda centers:np.ones(len(centers),bool))
    assert ordinary['assigned_boundary_area_A2']>0
    assert excluded['assigned_boundary_area_A2']==0
    assert excluded['excluded_boundary_area_A2']==excluded['boundary_area_A2']


def test_full_trajectory_summary_and_export_are_repeatable(tmp_path):
    from crevice.trajectory import Trajectory
    from crevice.cavity_trajectory import analyze_cavity_trajectory,summarize_cavity_trajectory,write_cavity_trajectory_bundle
    import json
    frame=rectangular_channel(half_x=4,half_y=5,half_length=8,capped=True)
    cast=rolling_probe_cast(frame,spacing=.75,selection_mode='all',enclosure_fraction=0,min_component_volume=1)
    main=max(cast.components,key=lambda c:c.volume);cast=replace(cast,components=(main,),points=main.points)
    ref=CavityReference.from_dx(write_void_cast_dx(cast,tmp_path/'ref.dx'),axis=(0,0,1))
    moved,_,_=transform(frame)
    frames=[replace(frame,frame_index=0),replace(moved,frame_index=1),replace(frame,frame_index=2)]
    traj=Trajectory(frames,time_step=1.)
    a=analyze_cavity_trajectory(traj,ref,geometry_mode='reference-region',geometry={'probe_radius':1.4},workers=1)
    b=analyze_cavity_trajectory(traj,ref,geometry_mode='reference-region',geometry={'probe_radius':1.4},workers=2)
    assert [r['volume_A3'] for r in a['frames']]==[r['volume_A3'] for r in b['frames']]
    before=json.dumps(a['residue_dynamics'],sort_keys=True)
    report,arrays=summarize_cavity_trajectory(a,replicates=200)
    paths=write_cavity_trajectory_bundle(a,tmp_path/'bundle',replicates=200,dpi=60)
    assert json.dumps(a['residue_dynamics'],sort_keys=True)==before
    assert report['frame_count']==report['matched_frames']==3
    assert len(np.load(paths['all_frame_arrays_npz'])['time_ps'])==3
    assert Path(paths['diameter_profile_png']).is_file()


def test_saved_profile_statistics_separate_missing_constant_and_estimable_series():
    from crevice.cavity_reports import sectional_statistics
    rng=np.random.default_rng(100);values=np.column_stack([np.ones(2000),10+rng.normal(size=2000),10+rng.normal(size=2000)])
    values[35,2]=np.nan
    rows=sectional_statistics(values,[0,1,2],replicates=200)
    assert rows[0]['mean_CI_status']=='unavailable_observed_constant'
    assert rows[0]['pointwise_mean_lower_A'] is None
    assert rows[1]['pointwise_mean_lower_A']<rows[1]['mean_diameter_A']<rows[1]['pointwise_mean_upper_A']
    assert rows[1]['fluctuation_lower_A']<rows[1]['pointwise_mean_lower_A']
    assert rows[2]['mean_CI_status']=='unavailable_missing_frames'
    assert rows[2]['observed_frames']==1999


def test_reference_provenance_and_invalid_fixed_grid(tmp_path):
    import hashlib
    frame=rectangular_channel(half_x=4,half_y=5,half_length=8,capped=True)
    for basis,anchor in [(np.ones((3,3)),[0,0,0]),(np.eye(3),[0,np.nan,0])]:
        with pytest.raises(ValueError):rolling_probe_cast(frame,grid_basis=basis,grid_anchor=anchor)
    cast=rolling_probe_cast(frame,spacing=1.,min_component_volume=1,selection_mode='all',enclosure_fraction=0)
    path=write_void_cast_dx(cast,tmp_path/'ref.dx');ref=CavityReference.from_dx(path)
    assert ref.source_sha256==hashlib.sha256(Path(path).read_bytes()).hexdigest()
    assert ref.source_path==str(Path(path).resolve())


def test_cavity_cli_exposes_explicit_all_frame_guard_and_geometry_definitions():
    from crevice.cli import build_parser
    args=build_parser().parse_args(['cavity-trajectory','prod.gro','prod.xtc','--reference-volume-dx','pocket.dx','--out-dir','output',
      '--max-frames','1100','--geometry-mode','reference-region','--region-margin','2','--axis','0,0,1','--probe-radius','1.4'])
    assert args.max_frames==1100 and args.stop is None and args.pbc=='check'
    assert args.geometry_mode=='reference-region' and args.probe_radius==1.4
