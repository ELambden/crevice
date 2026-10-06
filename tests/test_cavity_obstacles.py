"""Physical-obstacle, periodic-image and lattice controls; no biological claims."""
from dataclasses import replace
from itertools import product
import numpy as np
import pytest
from crevice.models import Atom
from crevice.cavity_trajectory import CavityReference,analyze_cavity_trajectory
from crevice.cavity_obstacles import aligned_obstacle_images,ObstacleSource
from crevice.regional_volume import regional_probe_cast
from crevice.rolling import rolling_probe_cast
from crevice.volume_export import write_void_cast_dx,cast_grid
from crevice.residue_evidence import boundary_residue_evidence
from crevice.trajectory import load_trajectory_report
from test_channel_coordinate import rectangular_channel


@pytest.mark.parametrize('box',[[18,20,22,90,90,90],[18,20,22,75,80,65]])
def test_periodic_obstacles_agree_with_independent_lattice_enumeration(box):
    from MDAnalysis.lib.mdamath import triclinic_vectors
    from scipy.spatial.transform import Rotation
    cell=triclinic_vectors(box).astype(float)
    xyz=np.array([[1.,2.,3.],[16,2,4],[-1,19,7]])
    atoms=[Atom(i,'O','LIG','L',i,*p,'O',hetero=True) for i,p in enumerate(xyz)]
    R=Rotation.from_rotvec([.2,-.4,.3]).as_matrix();t=np.array([2,3,-4])
    lower=np.array([-10,-9,-8]);upper=np.array([13,14,15]);padding=2.
    got,_=aligned_obstacle_images(atoms,xyz+2*cell[0]-cell[1],lower,upper,rotation=R,translation=t,padding=padding,box=box)
    expected=[]
    for shift in product(range(-4,5),repeat=3):
        p=(xyz+np.asarray(shift)@cell)@R+t
        expected.extend((i,*np.round(q,6)) for i,q in enumerate(p) if np.all(q>=lower-padding)&np.all(q<=upper+padding))
    actual=[(a.serial,*np.round(a.coord,6)) for a in got]
    assert sorted(actual)==sorted(expected)


def cavity(tmp_path):
    frame=rectangular_channel(half_x=4,half_y=4,half_length=4,capped=True)
    cast=rolling_probe_cast(frame,spacing=.5,probe_radius=.8,selection_mode='all',enclosure_fraction=0,min_component_volume=1,grid_basis=np.eye(3),grid_anchor=(0,0,0))
    main=max(cast.components,key=lambda c:c.volume)
    dx=write_void_cast_dx(replace(cast,components=(main,),points=main.points),tmp_path/'reference.dx')
    return frame,CavityReference.from_dx(dx,axis=(0,0,1))


def test_occupied_ligand_reduces_volume_and_owns_its_boundary(tmp_path):
    frame,ref=cavity(tmp_path)
    ligand=Atom(9000,'C','LIG','L',1,0,0,0,'C',hetero=True)
    base,_,_=regional_probe_cast(frame,ref,probe_radius=.8)
    filled,_,_=regional_probe_cast(frame,ref,probe_radius=.8,obstacle_atoms=frame.atoms+(ligand,))
    assert len(filled.points)<len(base.points)
    assert min(np.linalg.norm(p.position) for p in filled.points)>=1.7-1e-8
    grid,origin=cast_grid(filled)
    report=boundary_residue_evidence(frame,grid,origin,ref.basis*ref.spacing,obstacle_atoms=frame.atoms+(ligand,))
    assert report['other_obstacle_boundary_area_A2']>0
    assert report['assigned_boundary_area_A2']==pytest.approx(report['protein_assigned_boundary_area_A2']+report['other_obstacle_boundary_area_A2'])
    assert 'L:LIG1' not in report['lining_residues']
    assert all(r['residue'].startswith('A:ALA') for r in report['residues'])
    # Additional distant physical atoms must leave the cast exactly unchanged.
    distant=replace(ligand,x=100)
    other,_,_=regional_probe_cast(frame,ref,probe_radius=.8,obstacle_atoms=frame.atoms+(distant,))
    assert [p.position for p in other.points]==[p.position for p in base.points]
    assert [p.raw_clearance for p in other.points]==[p.raw_clearance for p in base.points]


def md_pair(tmp_path):
    import MDAnalysis as mda
    frame,ref=cavity(tmp_path);n=len(frame.atoms)+1
    u=mda.Universe.empty(n,n_residues=n,atom_resindex=np.arange(n),trajectory=True)
    u.add_TopologyAttr('names',['C']*n);u.add_TopologyAttr('resnames',['ALA']*(n-1)+['LIG'])
    u.add_TopologyAttr('resids',np.arange(1,n+1));u.add_TopologyAttr('segids',['A'])
    xyz=np.vstack(([a.coord for a in frame.atoms],[0,0,0]))+20
    u.atoms.positions=xyz;u.dimensions=[50,50,50,90,90,90]
    gro=tmp_path/'system.gro';xtc=tmp_path/'run.xtc';u.atoms.write(str(gro))
    with mda.Writer(str(xtc),n_atoms=n,dt=10) as w:
        for i in range(2):
            u.atoms.positions=xyz+np.array([i,0,0]);u.trajectory.ts.frame=i;u.trajectory.ts.time=i*10;w.write(u.atoms)
    ref.points=ref.points+20;ref.anchor=ref.anchor+20;ref.origin=ref.origin+20
    traj,_=load_trajectory_report(gro,xtc,selection='protein',pbc='none')
    return gro,xtc,traj,ref


def test_streamed_source_keeps_display_separate_and_parallel_matches(tmp_path):
    gro,xtc,traj,ref=md_pair(tmp_path)
    source=ObstacleSource(gro,xtc,protein_selection='protein',obstacle_selection='all',reference_frame=traj.frames[0],frame_times=dict(zip([0,1],traj.times())))
    config=dict(geometry_mode='reference-region',geometry={'probe_radius':.8},obstacles=source)
    serial=analyze_cavity_trajectory(traj,ref,**config)
    parallel=analyze_cavity_trajectory(traj,ref,workers=2,**config)
    assert source.metadata['other_heavy_atom_count']==1
    for a,b in zip(serial['frames'],parallel['frames']):
        assert a['volume_A3']==b['volume_A3']
        assert a['boundary']['residues']==b['boundary']['residues']
        assert a['boundary']['other_obstacle_boundary_area_A2']>0
    assert serial['settings']['obstacles']['selection']=='all'
    assert serial['frames'][0]['volume_A3']==serial['frames'][1]['volume_A3']
    with pytest.raises(ValueError,match='retain every'):
        ObstacleSource(gro,xtc,protein_selection='protein',obstacle_selection='resname LIG',reference_frame=traj.frames[0],frame_times={0:0.})
    source.frame_times[0]=999
    with pytest.raises(ValueError,match='timestamps'):source.snapshot(0)
    with pytest.raises(ValueError,match='reference-region'):
        analyze_cavity_trajectory(traj,ref,obstacles=source)


def test_grid_phase_changes_lattice_without_moving_reference_region(tmp_path):
    frame,ref=cavity(tmp_path)
    phase=CavityReference.from_dx(ref.source_path,grid_phase=(.5,.25,.75))
    assert np.array_equal(ref.points,phase.points)
    assert (phase.anchor-ref.anchor)@phase.basis.T==pytest.approx(np.array([.5,.25,.75])*phase.spacing)
    for value in [(1,0,0),(-.1,0,0),(np.nan,0,0),(0,0)]:
        with pytest.raises(ValueError,match='grid_phase'):CavityReference.from_dx(ref.source_path,grid_phase=value)


def test_d_amino_acid_templates_retain_identity_and_signed_torsions():
    from crevice.models import StructureFrame
    from crevice.interaction_analysis import chemical_topology_from_frame
    names=['N','CA','CB','CG','C','O']
    xyz=np.array([[0,0,0],[1.4,0,0],[1.4,1.4,0],[1.4,1.4,1.4],[2.4,-.8,0],[3.4,-.8,0]])
    atoms=tuple(Atom(i,n,'DLE','A',1,*p,'N' if n=='N' else 'O' if n=='O' else 'C') for i,(n,p) in enumerate(zip(names,xyz)))
    model=chemical_topology_from_frame(StructureFrame(atoms),[])
    a=model.measure(xyz,['A:DLE1']);b=model.measure(xyz*np.array([1,1,-1]),['A:DLE1'])
    assert not model.unknown_residues
    assert model.acceptors.tolist()==[5]
    assert a['torsions']['A:DLE1']['chi1']==pytest.approx(-b['torsions']['A:DLE1']['chi1'])
    assert abs(a['torsions']['A:DLE1']['chi1'])==pytest.approx(90)
    assert a['secondary']['A:DLE1']=='NA'
    assert model.metadata()['modified_residue_templates'][0]['chirality'].startswith('D;')
