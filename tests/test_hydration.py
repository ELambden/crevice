"""Hydration contracts: geometry, solvent identity, sampling and integration."""
from dataclasses import replace
import json,math
from pathlib import Path
import numpy as np
import pytest
from crevice import Atom,StructureFrame
from crevice.hydration import static_hydration,sampled_sasa,HydrationTopology
from crevice.hydration_trajectory import sampled_survival,contact_episodes
from crevice.hydration_export import write_hydration_bundle


def atom(serial,name,xyz,resid=1,resname='ALA',element='C',hetero=False,occupancy=None):
    return Atom(serial,name,resname,'A',resid,*xyz,element=element,hetero=hetero,occupancy=occupancy)


def test_sasa_matches_analytic_isolated_and_overlapping_spheres():
    radius=3.1
    single=sampled_sasa([[0,0,0]],[[0,0,0]],[1.7],[0],samples=10000)
    assert single[0]==pytest.approx(4*math.pi*radius**2,rel=1e-12)
    overlapping=sampled_sasa([[0,0,0],[radius,0,0]],[[0,0,0],[radius,0,0]],[1.7,1.7],[0,1],samples=10000)
    assert overlapping==pytest.approx([3*math.pi*radius**2]*2,rel=.002)


def test_contacts_deduplicate_water_and_keep_backbone_sidechain_overlap():
    frame=StructureFrame((atom(1,'CA',(0,0,0)),atom(2,'CB',(0,1,0)),atom(3,'OG',(1,0,0),resid=2,resname='SER',element='O'),
                          atom(4,'O',(0,0,3),resid=4,resname='HOH',element='O',hetero=True),atom(5,'H1',(0,0,2),resid=4,resname='HOH',element='H',hetero=True)))
    result=static_hydration(frame,sasa_samples=64)
    assert result['metrics']['water_count'].tolist()==[1,1]
    assert result['metrics']['backbone_water_count'].tolist()==[1,0]
    assert result['metrics']['sidechain_water_count'].tolist()==[1,1]
    assert result['metrics']['NO_water_count'].tolist()==[0,1]
    assert result['bridges']=={(0,1):['A:HOH4']}
    assert result['settings']['environment_heavy_atoms']==3


def test_missing_deposited_waters_is_not_zero_hydration(tmp_path):
    result=static_hydration(StructureFrame((atom(1,'CA',(0,0,0)),)),sasa_samples=32)
    assert np.isnan(result['metrics']['water_count']).all()
    assert result['metrics']['sasa_A2'][0]>100
    paths=write_hydration_bundle(result,tmp_path,dpi=60,replicates=200)
    data=json.loads(Path(paths['hydration_json']).read_text())
    assert data['water_status']=='unavailable_no_explicit_water'
    assert data['residues'][0]['statistics']['water_count']['mean'] is None
    assert data['hydrogen_bonds']['status']=='not_assigned'
    assert all(Path(p).is_file() for p in paths.values())


@pytest.mark.parametrize('box',[[30,30,30,90,90,90],[30,32,35,80,75,60]])
def test_periodic_contacts_and_sasa_match_explicit_neighbour_image(box):
    from MDAnalysis.lib.mdamath import triclinic_vectors
    vectors=triclinic_vectors(box).astype(float)
    base=np.array([3.,3.,3.]);near=base+np.array([0,0,3.]);shifted=near+vectors[1]
    frame=StructureFrame((atom(1,'CA',base),atom(2,'O',shifted,resid=5,resname='HOH',element='O',hetero=True)))
    periodic=static_hydration(frame,box=box,sasa_samples=64)
    assert periodic['metrics']['water_count'][0]==1
    assert static_hydration(frame,sasa_samples=64)['metrics']['water_count'][0]==0
    expected=sampled_sasa([base],[base,base+[3.,0,0]],[1.7,1.7],[0],samples=256)
    actual=sampled_sasa([base],[base,base+[3.,0,0]+vectors[1]],[1.7,1.7],[0],samples=256,box=box)
    assert actual==pytest.approx(expected,rel=1e-8)


def test_sasa_uses_ligand_obstacles_and_is_rigid_invariant_with_reference_orientation():
    xyz=np.array([[0.,0.,0.],[2,1,0.]])
    rotation=np.array([[0.,-1,0],[1,0,0],[0,0,1.]])
    a=sampled_sasa(xyz[:1],xyz,[1.7,1.7],[0],samples=256)
    b=sampled_sasa((xyz@rotation+10)[:1],xyz@rotation+10,[1.7,1.7],[0],samples=256,orientation=rotation.T)
    assert a==pytest.approx(b,abs=1e-10)
    frame=StructureFrame((atom(1,'CA',xyz[0]),atom(2,'C1',xyz[1],resid=2,resname='LIG',hetero=True)))
    assert static_hydration(frame,sasa_samples=256)['metrics']['sasa_A2']==pytest.approx(a)


def test_zero_occupancy_and_ambiguous_water_identity():
    protein=atom(1,'CA',(0,0,0));water=atom(2,'O',(0,0,3),resid=2,resname='HOH',element='O',hetero=True,occupancy=0)
    assert static_hydration(StructureFrame((protein,water)),sasa_samples=32)['water_population']==0
    water=replace(water,occupancy=.5)
    assert static_hydration(StructureFrame((protein,water)),sasa_samples=32)['partial_occupancy_waters']==1
    with pytest.raises(ValueError,match='exactly one'):
        static_hydration(StructureFrame((protein,water,replace(water,serial=3))),sasa_samples=32)


def test_sampled_survival_and_episode_censoring_are_not_lifetimes():
    histories=[{'a'},{'a'},{'b'},{'b'}];times=np.arange(4)*100
    result=sampled_survival(histories,times,max_lag=3)
    assert result['rows'][1]['survival_fraction']==pytest.approx(2/3)
    assert result['rows'][2]['survival_fraction']==0
    rows=contact_episodes([[h] for h in histories],times,['A:ALA1'])
    assert [r['observed_span_ps'] for r in rows]==[100,100]
    assert rows[0]['left_censored'] and not rows[0]['right_censored']
    assert not rows[1]['left_censored'] and rows[1]['right_censored']
    assert sampled_survival(histories,[0,100,250,300])['status']=='unavailable_irregular_sampling'


def test_invalid_controls_and_unknown_residue_focus():
    frame=StructureFrame((atom(1,'CA',(0,0,0)),))
    for kw in [{'cutoff':-1},{'sasa_samples':10},{'sasa_probe':float('nan')},{'residue_ids':['A:LYS4']},{'box':[2,2,2,90,90,90]}]:
        with pytest.raises(ValueError):static_hydration(frame,**kw)


def test_hetatm_focus_residues_are_reported_untyped_not_fatal():
    """1GRM as PDB: D-amino acids are HETATM, so they are not measured as protein.

    They are kept as occluding atoms and listed explicitly; only identities that
    are absent from the structure are an error.
    """
    frame=StructureFrame((atom(1,'CA',(0,0,0),resid=1,resname='VAL'),
                          atom(2,'CA',(3.8,0,0),resid=2,resname='DLE',hetero=True),
                          atom(3,'O',(0,0,3),resid=9,resname='HOH',element='O',hetero=True)))
    result=static_hydration(frame,residue_ids=['A:VAL1','A:DLE2','A:HOH9'],sasa_samples=64)
    assert result['residue_ids']==('A:VAL1',)
    assert result['untyped_focus_residues']==[{'residue':'A:DLE2','reason':'hetero_record'},
                                              {'residue':'A:HOH9','reason':'water'}]
    assert result['settings']['untyped_focus_residues']==result['untyped_focus_residues']
    assert result['settings']['environment_heavy_atoms']==2   # DLE still occludes
    with pytest.raises(ValueError,match='not present'):
        static_hydration(frame,residue_ids=['A:VAL1','A:LYS7'],sasa_samples=64)
    with pytest.raises(ValueError,match='no protein residue'):
        static_hydration(frame,residue_ids=['A:DLE2'],sasa_samples=64)


def test_untyped_focus_residues_reach_the_hydration_json(tmp_path):
    frame=StructureFrame((atom(1,'CA',(0,0,0),resid=1,resname='VAL'),
                          atom(2,'CA',(3.8,0,0),resid=2,resname='DLE',hetero=True)))
    result=static_hydration(frame,residue_ids=['A:VAL1','A:DLE2'],sasa_samples=64)
    paths=write_hydration_bundle(result,tmp_path,dpi=60,replicates=200)
    data=json.loads(Path(paths['hydration_json']).read_text())
    assert data['settings']['untyped_focus_residues']==[{'residue':'A:DLE2','reason':'hetero_record'}]


def md_fixture(tmp_path):
    import MDAnalysis as mda
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    from crevice.presentation import write_viewer_structure
    from crevice.cli import main
    frame=StructureFrame((atom(1,'CA',(10,10,10)),atom(2,'CA',(12,10,10),resid=2),atom(3,'CA',(10,12,10),resid=3),
                          atom(4,'O',(10,10,13),resid=4,resname='HOH',element='O',hetero=True)))
    source=tmp_path/'input.pdb'
    lines=[]
    for a in frame.atoms:
        record='HETATM' if a.hetero else 'ATOM  '
        lines.append(f'{record}{a.serial:5d} {a.name:^4} {a.resname:3} A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00          {a.element:>2}')
    source.write_text('\n'.join(lines)+'\nEND\n')
    u=mda.Universe(source);trajectory=tmp_path/'input.xtc';coords=u.atoms.positions.copy()
    with mda.Writer(str(trajectory),len(u.atoms),dt=100) as writer:
        for i in range(3):
            u.atoms.positions=coords.copy()
            if i==1:u.atoms[-1].position=[20,20,20]
            if i==2:u.atoms.positions=coords+[29,0,0]
            u.dimensions=[50,50,50,90,90,90];u.trajectory.ts.frame=i;u.trajectory.ts.time=i*100;writer.write(u.atoms)
    return source,trajectory


def test_streaming_hydration_preserves_water_ids_and_periodic_alignment(tmp_path):
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    from crevice.cli import main
    source,trajectory=md_fixture(tmp_path)
    result=analyze_hydration_trajectory(source,trajectory,sasa_samples=32,max_frames=3)
    assert result['arrays']['time_ps'].tolist()==[0,100,200]
    assert result['arrays']['water_count'][:,0].tolist()==[1,0,1]
    assert result['histories'][0][0]==result['histories'][2][0]
    assert result['arrays']['sasa_A2'][0]==pytest.approx(result['arrays']['sasa_A2'][2],abs=1e-10)
    files=write_hydration_bundle(result,tmp_path/'bundle',dpi=60,replicates=200)
    assert all(Path(p).is_file() for p in files.values())
    assert main(['hydration',str(source),'--trajectory',str(trajectory),'--out-dir',str(tmp_path/'cli'),'--max-frames','3','--hydration-sasa-points','32','--dpi','60'])==0


def reference_dx(path):
    grid=np.zeros((7,7,7));grid[2:5,2:5,2:5]=1
    lines=['object 1 class gridpositions counts 7 7 7','origin 8.5 8.5 11.5','delta 0.5 0 0','delta 0 0.5 0','delta 0 0 0.5',
           'object 2 class gridconnections counts 7 7 7','object 3 class array type double rank 0 items 343 data follows']
    lines+=[' '.join(str(int(x)) for x in grid.ravel())];lines+=['attribute "dep" string "positions"']
    path.write_text('\n'.join(lines)+'\n');return path


def test_standard_static_hydration_and_explicit_skip_preserve_geometry(tmp_path):
    from crevice.cli import main
    source,_=md_fixture(tmp_path);dx=reference_dx(tmp_path/'ref.dx');before=dx.read_bytes()
    base=['residue-evidence',str(source),'--volume-dx',str(dx),'--dpi','60','--hydration-sasa-points','32']
    assert main(base+['--out-dir',str(tmp_path/'standard')])==0
    hydration=json.loads((tmp_path/'standard/crevice_hydration.json').read_text())
    assert hydration['water_population']==1 and hydration['residue_count']>0
    assert main(base+['--out-dir',str(tmp_path/'skipped'),'--skip-hydration'])==0
    assert not (tmp_path/'skipped/crevice_hydration.json').exists()
    assert dx.read_bytes()==before


def test_standard_cavity_trajectory_and_subset_keep_original_alignment_reference(tmp_path):
    from crevice.cli import main
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    source,trajectory=md_fixture(tmp_path);dx=reference_dx(tmp_path/'ref.dx');output=tmp_path/'standard'
    assert main(['cavity-trajectory',str(source),str(trajectory),'--reference-volume-dx',str(dx),'--geometry-mode','reference-region',
      '--max-frames','3','--out-dir',str(output),'--dpi','60','--bootstrap-replicates','200','--hydration-sasa-points','32'])==0
    report=json.loads((output/'crevice_hydration.json').read_text());assert report['frame_count']==3
    cached=output/'crevice_cavity_statistics.json'
    subset=analyze_hydration_trajectory(source,trajectory,context_json=cached,start=2,stop=3,max_frames=1,sasa_samples=32)
    assert subset['arrays']['time_ps'].tolist()==[200]
    assert subset['settings']['alignment_reference_source_index']==0
    assert subset['alignment'][0]['translation_A']==pytest.approx([-29,0,0],abs=1e-4)
    dx.write_text(dx.read_text()+'# changed reference provenance\n')
    with pytest.raises(ValueError,match='map hash differs'):
        analyze_hydration_trajectory(source,trajectory,context_json=cached,start=2,stop=3,max_frames=1,sasa_samples=32)


def test_hydration_mobility_uses_adjacent_aligned_frames_and_tracks_transitions():
    from crevice.hydration_export import hydration_mobility_statistics
    arrays={'time_ps':np.arange(4)*100.,'residue_ids':np.array(['A:ALA1']),
      'aligned_residue_centroids_A':np.array([[[0.,0,0]],[[1.,0,0]],[[3.,0,0]],[[3.,0,0]]]),
      'water_count':np.array([[0.],[0.],[1.],[1.]])}
    row=hydration_mobility_statistics(arrays)[0]
    assert row['mean_centroid_step_A']==1.
    assert row['mean_step_dry_at_both_ends_A']==1.
    assert row['mean_step_wet_at_both_ends_A']==0.
    assert row['changed_hydration_state_intervals']==1


def test_confidence_opt_out_and_invalid_parameters_are_explicit(tmp_path):
    from crevice.cli import main
    result=static_hydration(StructureFrame((atom(1,'CA',(0,0,0)),)),sasa_samples=32)
    files=write_hydration_bundle(result,tmp_path/'disabled',confidence=0,dpi=60)
    report=json.loads(Path(files['hydration_json']).read_text())
    assert report['residues'][0]['statistics']['sasa_A2']['confidence_interval']['status']=='not_requested'
    with pytest.raises(ValueError,match='confidence'):
        write_hydration_bundle(result,tmp_path/'bad',confidence=2,dpi=60)
    assert main(['hydration','missing.pdb','--out-dir',str(tmp_path/'early'),'--hydration-cutoff','-1'])==2
    assert not (tmp_path/'early').exists()
