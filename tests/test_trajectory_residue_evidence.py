"""Independent scientific contracts for uncertainty, alignment and residue evidence."""
from dataclasses import replace
import json
import math
import numpy as np
import pytest
from crevice import Atom,StructureFrame
from crevice.models import ChannelPoint,PoreProfile,FrameAnalysis,TrajectoryAnalysis
from crevice.trajectory import analyze_trajectory,align_frame_report,_check_periodic_pairs,Trajectory
from crevice.uncertainty import block_mean_confidence,statistical_inefficiency
from crevice.ensemble import profile_distribution
from crevice.networks import build_residue_network
from crevice.residue_dynamics import ResidueDynamics
from crevice.residue_evidence import boundary_residue_evidence,read_binary_dx


def atom(i,xyz,chain='A'):
    return Atom(i,'CA','ALA',chain,i,*xyz,'C')


def series(values):
    frames=[]
    for i,y in enumerate(np.asarray(values)):
        points=tuple(ChannelPoint(j,(0,0,float(j)),float(v),float(v),float(j)) for j,v in enumerate(y))
        p=PoreProfile((0,0,0),(0,0,1),points,method='axial-connected',metadata={'status':'resolved'})
        frames.append(FrameAnalysis(i,time=float(i),profile=p))
    return TrajectoryAnalysis(tuple(frames),('profile',),{'aligned_to_first_frame':True})


def test_joint_bootstrap_preserves_scaling_and_constant_columns():
    rng=np.random.default_rng(19);x=rng.normal(size=500)
    result=block_mean_confidence(np.column_stack([x,2*x,np.ones(len(x))]),block_length=8,replicates=300,seed=5)
    assert result['status']=='estimated'
    for key in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
        assert result[key][1]==pytest.approx(2*result[key][0])
        assert result[key][2]==pytest.approx(1)
    assert result==block_mean_confidence(np.column_stack([x,2*x,np.ones(len(x))]),block_length=8,replicates=300,seed=5)


def test_correlated_series_has_wider_uncertainty_than_iid_resampling():
    rng=np.random.default_rng(14);x=np.zeros(3000)
    for i in range(1,len(x)):x[i]=.93*x[i-1]+rng.normal()
    assert statistical_inefficiency(x)>10
    iid=block_mean_confidence(rng.permutation(x)[:,None],block_length=1,replicates=400)
    blocked=block_mean_confidence(x[:,None],block_length=60,replicates=400)
    width=lambda r:r['pointwise_upper'][0]-r['pointwise_lower'][0]
    assert width(blocked)>2*width(iid)


@pytest.mark.parametrize('kw',[{'confidence':1},{'confidence':float('nan')},{'block_length':0},{'block_length':1.5},{'replicates':10}])
def test_invalid_uncertainty_options(kw):
    with pytest.raises(ValueError):block_mean_confidence(np.ones((40,2)),**kw)


def test_short_and_missing_series_do_not_receive_confidence_claims():
    assert block_mean_confidence(np.ones((10,2)))['status']=='insufficient_independent_blocks'
    x=np.ones((100,2));x[20,1]=np.nan
    result=block_mean_confidence(x,block_length=5,replicates=200)
    assert result['pointwise_lower'][1] is None
    assert result['pointwise_lower'][0]==pytest.approx(1)
    result=block_mean_confidence(np.arange(100)[:,None],block_length=30)
    assert result['status']=='insufficient_independent_blocks'


def test_distribution_keeps_missing_frames_and_suppresses_inference():
    trajectory=series(np.ones((50,3)))
    frames=list(trajectory.frames);frames[10]=replace(frames[10],profile=None)
    distribution=profile_distribution(replace(trajectory,frames=tuple(frames)),confidence=.95,bootstrap_replicates=200)
    assert distribution['missing_profile_count']==1
    assert distribution['rows'][0]['coverage_fraction']==pytest.approx(49/50)
    assert distribution['confidence_interval']['status']=='no_complete_coverage'
    assert distribution['rows'][0]['mean_band_lower'] is None


def test_confidence_requires_alignment_and_regular_time_sampling():
    trajectory=series(np.ones((50,3)))
    with pytest.raises(ValueError,match='alignment'):
        profile_distribution(replace(trajectory,metadata={}),confidence=.95)
    frames=list(trajectory.frames);frames[20]=replace(frames[20],time=19.5)
    with pytest.raises(ValueError,match='regularly'):
        profile_distribution(replace(trajectory,frames=tuple(frames)),confidence=.95)


def test_unvalidated_axis_scans_never_receive_pore_confidence_bands():
    trajectory=series(np.ones((50,3)))
    frames=tuple(replace(f,profile=replace(f.profile,metadata={'status':'unvalidated_axis_scan'})) for f in trajectory.frames)
    d=profile_distribution(replace(trajectory,frames=frames),confidence=.95)
    assert d['confidence_interval']['status']=='unavailable_unvalidated_axis_scan'


def test_alignment_reports_proper_transform_and_fit_error():
    reference=StructureFrame(tuple(atom(i,p) for i,p in enumerate([(0,0,0),(3,0,0),(0,4,0),(0,0,5)],1)))
    moved=replace(reference,atoms=tuple(replace(a,x=-a.y+9,y=a.x-6,z=a.z+2) for a in reference.atoms))
    aligned,report=align_frame_report(moved,reference)
    assert report['rmsd_after_A']<1e-12<report['rmsd_before_A']
    assert np.linalg.det(report['rotation_rows'])==pytest.approx(1)
    np.testing.assert_allclose([a.coord for a in aligned.atoms],[a.coord for a in reference.atoms],atol=1e-12)
    with pytest.raises(ValueError,match='unknown'):
        align_frame_report(moved,reference,residues=['A:ALA999'])


def test_periodic_screen_rejects_wrapped_bond_and_accepts_whole():
    pair=np.array([[0,1]]);box=np.array([10,10,10,90,90,90],dtype=float)
    with pytest.raises(ValueError,match='periodic'):
        _check_periodic_pairs(np.array([[.5,0,0],[9.5,0,0]]),pair,box,12)
    _check_periodic_pairs(np.array([[.5,0,0],[-.5,0,0]]),pair,box,12)


def test_center_distance_is_distinct_from_legacy_surface_cutoff():
    frame=StructureFrame((atom(1,(0,0,0)),atom(2,(6,0,0))))
    assert len(build_residue_network(frame,cutoff=4.5).edges)==1
    assert not build_residue_network(frame,cutoff=4.5,distance_metric='center').edges
    graph=build_residue_network(frame,cutoff=6,distance_metric='center')
    assert graph.edges[0].distance==pytest.approx(6)
    assert graph.edges[0].metadata['closest_pair']['surface_distance']==pytest.approx(2.6)


def test_contact_occupancy_is_not_split_when_interaction_type_changes():
    frames=[StructureFrame((atom(1,(0,0,0)),atom(2,(3,0,0),'B'))) for _ in range(4)]
    accumulator=ResidueDynamics(lining_residues=['A:ALA1'],min_occupancy=.75)
    for i,frame in enumerate(frames):
        network=build_residue_network(frame,distance_metric='center')
        network=replace(network,edges=(replace(network.edges[0],interaction='typeA' if i%2 else 'typeB'),))
        accumulator.add(frame,network,frame_index=i,time=i)
    result=accumulator.finish(aligned=True);edge=result['contacts'][0]
    assert edge['occupancy']==1
    assert edge['contact_type_fractions']=={'typeA':.5,'typeB':.5}
    assert edge['nonlocal_contact'] and edge['reference_lining_contact']
    assert edge['dccm'] is None


def test_rigid_motion_does_not_become_residue_coupling():
    reference=StructureFrame(tuple(atom(i,p) for i,p in enumerate([(0,0,0),(3,0,0),(0,3,0),(0,0,3)],1)))
    frames=[]
    for i in range(20):
        c,s=math.cos(i*.13),math.sin(i*.13)
        frames.append(replace(reference,frame_index=i,atoms=tuple(replace(a,x=c*a.x-s*a.y+i,y=s*a.x+c*a.y,z=a.z-i*.2) for a in reference.atoms)))
    result=analyze_trajectory(frames,analyses=['network'],retain_networks=False,network_kwargs={'distance_metric':'center'})
    report=result.metadata['residue_dynamics']
    assert all(r['centroid_rmsf_A']<1e-10 for r in report['residues'])
    assert all(e['dccm'] is None for e in report['contacts'])


def test_missing_reference_profile_retains_other_analyses(monkeypatch):
    from crevice.sections import ChannelResolutionError
    import crevice.trajectory as module
    calls=[]
    def fail(*a,**k):calls.append(1);raise ChannelResolutionError('reference unresolved')
    monkeypatch.setattr(module,'pore_profile',fail)
    frame=StructureFrame(tuple(atom(i,p) for i,p in enumerate([(0,0,0),(3,0,0),(0,3,0)],1)))
    result=analyze_trajectory([frame,frame],analyses=['profile','network'],network_kwargs={'distance_metric':'center'})
    assert len(calls)==1 and len(result.frames)==2
    assert all(f.profile is None and f.network is not None for f in result.frames)
    assert result.metadata['reference_profile_status']=='unresolved'
    with pytest.raises(ChannelResolutionError):analyze_trajectory([frame],on_profile_error='raise')


def test_boundary_area_and_nonlocal_partners_survive_rigid_transform():
    grid=np.zeros((5,5,5));grid[1:4,1:4,1:4]=1
    origin=np.array([-1.,-1.,-1.]);deltas=np.eye(3)*.5
    frame=StructureFrame((atom(1,(-2.5,0,0)),atom(2,(2.5,0,0)),atom(3,(-5,0,0),'B')))
    before=boundary_residue_evidence(frame,grid,origin,deltas)
    rot=np.array([[0,1,0],[-1,0,0],[0,0,1]]);shift=np.array([11,-7,3])
    xyz=np.asarray([a.coord for a in frame.atoms])@rot+shift
    moved=replace(frame,atoms=tuple(replace(a,x=p[0],y=p[1],z=p[2]) for a,p in zip(frame.atoms,xyz)))
    after=boundary_residue_evidence(moved,grid,origin@rot+shift,deltas@rot)
    assert before['boundary_area_A2']==pytest.approx(after['boundary_area_A2'])
    rows={r['residue']:r for r in before['residues']};other={r['residue']:r for r in after['residues']}
    for key in rows:assert rows[key]['boundary_area_A2']==pytest.approx(other[key]['boundary_area_A2'])
    assert rows['B:ALA3']['role']=='nonlocal_lining_partner'
    assert rows['A:ALA1']['boundary_area_A2']==pytest.approx(rows['A:ALA2']['boundary_area_A2'])
    assert sum(r['boundary_area_A2'] for r in rows.values())==pytest.approx(before['assigned_boundary_area_A2'])


def test_smoothed_dx_is_rejected_for_measured_boundary_analysis(tmp_path):
    p=tmp_path/'fake.dx';p.write_text('object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 class array type double rank 0 items 8 data follows\n0 0 0 0 0 0 0 .7\nattribute "dep" string "positions"\n')
    with pytest.raises(ValueError,match='binary'):read_binary_dx(p)


def test_constriction_finds_between_sample_limiting_residues():
    from crevice.residue_evidence import profile_constriction_evidence
    frame=StructureFrame((atom(1,(2,0,2)),atom(2,(-2,0,2))))
    points=(ChannelPoint(0,(0,0,0),1.1,1.1,0),ChannelPoint(1,(0,0,4),1.1,1.1,4))
    profile=PoreProfile((0,0,0),(0,0,1),points,method='axial-connected')
    report=profile_constriction_evidence(frame,profile,tolerance=0)
    assert report['bottleneck_radius_A']==pytest.approx(.3)
    assert {r['residue'] for r in report['residues']}=={'A:ALA1','A:ALA2'}
    assert all(r['path_position_A']==pytest.approx([0,0,2]) for r in report['residues'])


def test_distribution_rejects_parallel_distinct_channels():
    trajectory=series(np.ones((50,3)))
    frames=list(trajectory.frames)
    frames[10]=replace(frames[10],profile=replace(frames[10].profile,axis_origin=(4,0,0)))
    with pytest.raises(ValueError,match='origin'):
        profile_distribution(replace(trajectory,frames=tuple(frames)),confidence=.95)


def test_dynamic_correlation_recovers_known_coupled_and_independent_motion():
    reference=StructureFrame(tuple(atom(i,p,chain=str(i)) for i,p in enumerate([(0,0,0),(3,0,0),(0,3,0)],1)))
    acc=ResidueDynamics(min_occupancy=1)
    for i in range(200):
        phase=2*np.pi*i/200
        displacement=[.1*np.cos(phase),.1*np.cos(phase),.1*np.sin(phase)]
        frame=replace(reference,atoms=tuple(replace(a,z=z) for a,z in zip(reference.atoms,displacement)))
        network=build_residue_network(frame,cutoff=4.5,distance_metric='center')
        acc.add(frame,network,frame_index=i,time=i)
    pairs={(r['source'],r['target']):r['dccm'] for r in acc.finish(aligned=True)['contacts']}
    assert pairs[('1:ALA1','2:ALA2')]==pytest.approx(1)
    assert abs(pairs[('1:ALA1','3:ALA3')])<1e-12


def test_viewer_residue_identities_restore_original_numbers_and_chains(tmp_path):
    from crevice.presentation import write_viewer_structure
    from crevice import load_structure
    from crevice.residue_evidence import restore_viewer_identities
    frame=StructureFrame(tuple(replace(atom(i,p,chain='SYSTEM'),serial=i*10,resid=i+380,icode='A' if i==1 else '') for i,p in enumerate([(0,0,0),(2,0,0),(0,3,0)],1)))
    path=write_viewer_structure(frame,tmp_path/'viewer.pdb')
    loaded=load_structure(path)
    assert loaded.atoms[0].resid==1
    restored,provenance=restore_viewer_identities(loaded,path)
    assert [a.residue_key for a in restored.atoms]==[a.residue_key for a in frame.atoms]
    assert [a.serial for a in restored.atoms]==[a.serial for a in frame.atoms]
    assert provenance['sidecar']
