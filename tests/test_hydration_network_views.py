"""Geometry, chemistry, normalization and colour contracts for integrated views."""
import ast
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import pytest
from crevice import Atom,StructureFrame
from crevice.water_density import WaterDensity,write_scalar_dx
from crevice.interaction_analysis import chemical_topology_from_frame,InteractionAnalysis,write_interaction_bundle
from crevice.interaction_chemistry import torsion_degrees
from crevice.analysis_viewers import hydration_rgb


def atom(i,name,xyz,*,resid=1,resname='LYS',element='N',hetero=False):
    return Atom(i,name,resname,'A',resid,*xyz,element=element,hetero=hetero)


def chemical(frame,waters=(),box=None):
    return chemical_topology_from_frame(frame,list(waters),box=box)


def test_density_count_occupancy_and_integral_are_distinct():
    density=WaterDensity([[0,0,0]],margin=2,spacing=1,smoothing=.5,expected_frames=2)
    density.add([[0,0,0],[.1,.1,.1]])
    density.add([[1,0,0]])
    r=density.result(2)
    assert r['counts'].sum()==3
    assert r['density_A3'].sum()==pytest.approx(1.5)
    assert r['occupied_frames'].sum()==2
    assert r['voxel_occupancy'].max()==.5
    assert r['metadata']['density_integral_water_count']==pytest.approx(1.5)
    assert r['metadata']['display_integral_water_count']==pytest.approx(1.5,rel=1e-8)
    assert r['first_half_counts'].sum()==2
    with pytest.raises(ValueError,match='More density frames'):density.add([])


@pytest.mark.parametrize('box',[[30,30,30,90,90,90],[30,32,34,80,75,60]])
def test_density_periodic_images_and_rigid_alignment(box):
    from MDAnalysis.lib.mdamath import triclinic_vectors
    cell=triclinic_vectors(box)
    points=np.asarray([[2,3,4],[3,3,4]])
    base=WaterDensity(points,margin=2,expected_frames=1)
    base.add([[2,3,4]])
    shifted=WaterDensity(points,margin=2,expected_frames=1)
    shifted.add(np.asarray([[2,3,4]])+cell[1],box=box)
    assert np.array_equal(base.result(1)['counts'],shifted.result(1)['counts'])
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    translated=WaterDensity(points,margin=2,expected_frames=1)
    shift=np.array([10,20,30])
    translated.add(np.asarray([[2,3,4]])@rotation+shift,rotation=rotation.T,translation=-shift@rotation.T)
    assert np.array_equal(base.result(1)['counts'],translated.result(1)['counts'])


def test_density_missing_water_and_zero_observed_occupancy_differ():
    d=WaterDensity([[0,0,0]],expected_frames=1);d.add([])
    missing=d.result(0);observed=d.result(10)
    assert np.isnan(missing['density_A3']).all()
    assert missing['metadata']['status']=='unavailable_no_explicit_water'
    assert np.all(observed['density_A3']==0)
    assert observed['metadata']['status']=='observed'
    with pytest.raises(ValueError,match='exceeds'):WaterDensity([[0,0,0]],max_points=1)
    with pytest.raises(ValueError,match='incomplete'):WaterDensity([[0,0,0]],expected_frames=2).result(1)


def test_dx_roundtrip_retains_density_values_centres_and_axis_order(tmp_path):
    from gridData import Grid
    values=np.arange(24).reshape(2,3,4)/17
    p=tmp_path/'density.dx';write_scalar_dx(values,[2,3,4],.5,p)
    loaded=Grid(str(p))
    assert loaded.grid==pytest.approx(values,abs=1e-10)
    assert loaded.origin==pytest.approx([2,3,4])
    assert loaded.delta==pytest.approx([.5,.5,.5])


def test_hydrogen_bond_requires_typed_acceptor_and_directional_explicit_H():
    atoms=[atom(1,'NZ',(0,0,0)),atom(2,'HZ1',(1,0,0),element='H'),
           atom(3,'OD1',(3,0,0),resid=2,resname='ASP',element='O')]
    frame=StructureFrame(tuple(atoms));model=chemical(frame)
    r=model.measure([a.coord for a in atoms],['A:LYS1','A:ASP2'])
    assert any(e['kind']=='hydrogen_bond' and e['angle_deg']==pytest.approx(180) for e in r['observations'])
    bent=np.asarray([a.coord for a in atoms],dtype=float);bent[1]=[0,1,0]
    assert not any(e['kind']=='hydrogen_bond' for e in model.measure(bent,['A:LYS1','A:ASP2'])['observations'])
    no_H=StructureFrame((atoms[0],atoms[2]));m=chemical(no_H)
    assert m.metadata()['hydrogen_bond_status']=='unavailable_no_explicit_donor_hydrogens'
    assert any(e['kind']=='salt_bridge_candidate' for e in m.measure([a.coord for a in no_H.atoms],['A:LYS1','A:ASP2'])['observations'])
    amide=StructureFrame(tuple(atoms[:2]+[replace(atoms[2],name='ND2',element='N',resname='ASN')]))
    assert not any(e['kind']=='hydrogen_bond' for e in chemical(amide).measure([a.coord for a in amide.atoms],['A:LYS1','A:ASN2'])['observations'])


def test_named_histidine_tautomer_and_protonated_acid_are_respected():
    frame=StructureFrame((atom(1,'ND1',(0,0,0),resname='HSD'),atom(2,'HD1',(1,0,0),resname='HSD',element='H'),
                          atom(3,'NE2',(0,2,0),resname='HSD')))
    model=chemical(frame)
    assert model.acceptors.tolist()==[2]
    assert model.donor_pairs.tolist()==[[0,1]]
    generic=StructureFrame((replace(frame.atoms[0],resname='HIS'),replace(frame.atoms[2],resname='HIS')))
    assert not len(chemical(generic).acceptors)
    acid=StructureFrame((atom(1,'OD1',(0,0,0),resname='ASH',element='O'),atom(2,'HD1',(1,0,0),resname='ASH',element='H'),
                         atom(3,'OD2',(0,2,0),resname='ASH',element='O')))
    model=chemical(acid)
    assert model.acceptors.tolist()==[2] and not len(model.negative)


def test_water_bridge_requires_two_directional_bonds():
    atoms=(atom(1,'OG',(-2.8,0,0),resname='SER',element='O'),atom(2,'HG',(-1.8,0,0),resname='SER',element='H'),
           atom(3,'O',(0,0,0),resid=3,resname='HOH',element='O',hetero=True),
           atom(4,'H1',(1,0,0),resid=3,resname='HOH',element='H',hetero=True),
           atom(5,'H2',(-.3,.95,0),resid=3,resname='HOH',element='H',hetero=True),
           atom(6,'OD1',(2.8,0,0),resid=2,resname='ASN',element='O'))
    frame=StructureFrame(atoms);model=chemical(frame,[2]);xyz=np.asarray([a.coord for a in atoms])
    edges=model.measure(xyz,['A:SER1','A:ASN2'])['observations']
    bridge=next(e for e in edges if e['kind']=='water_bridge_hbond')
    assert bridge['atom_indices']==[5,2,0] or bridge['atom_indices']==[0,2,5]
    assert bridge['distance_A']==pytest.approx(2.8)
    xyz[3]=[0,1,0]
    assert not any(e['kind']=='water_bridge_hbond' for e in model.measure(xyz,['A:SER1','A:ASN2'])['observations'])


def ring_frame(shift):
    names=['CG','CD1','CE1','CZ','CE2','CD2'];points=np.c_[1.4*np.cos(np.arange(6)*np.pi/3),1.4*np.sin(np.arange(6)*np.pi/3),np.zeros(6)]
    return StructureFrame(tuple(atom(r*6+i+1,name,point+np.asarray(shift)*r,resid=r+1,resname='PHE',element='C') for r in range(2) for i,(name,point) in enumerate(zip(names,points))))


def test_aromatic_geometry_rejects_large_offset_and_nonplanar_rings():
    frame=ring_frame([0,0,3.6]);model=chemical(frame);xyz=np.asarray([a.coord for a in frame.atoms]);focus=['A:PHE1','A:PHE2']
    assert any(e['kind']=='aromatic_parallel_candidate' for e in model.measure(xyz,focus)['observations'])
    xyz[6:,0]+=4
    assert not model.measure(xyz,focus)['observations']
    xyz=np.asarray([a.coord for a in frame.atoms]);xyz[6,2]+=1.5
    assert not model.measure(xyz,focus)['observations']


def test_torsion_matches_independent_MDA_geometry_and_breaks_are_missing():
    from MDAnalysis.lib.distances import calc_dihedrals
    points=np.asarray([[1,0,0],[0,0,0],[0,1,0],[0,1,1]],dtype=float)
    assert torsion_degrees(points)==pytest.approx(np.degrees(calc_dihedrals(*points)))
    points[-1]+=[10,0,0]
    assert np.isnan(torsion_degrees(points))


def test_interaction_occupancy_uses_all_frames_and_circular_statistics(tmp_path):
    frame=ring_frame([0,0,3.6]);model=chemical(frame);analysis=InteractionAnalysis(model,['A:PHE1','A:PHE2'],3)
    xyz=np.asarray([a.coord for a in frame.atoms]);analysis.add(xyz);far=xyz.copy();far[6:,0]+=10;analysis.add(far);analysis.add(xyz)
    r=analysis.result();files=write_interaction_bundle(r,tmp_path,'test',times=np.array([0,100,200]),replicates=200)
    report=json.loads(Path(files['interaction_json']).read_text())
    assert report['edges'][0]['observed_occupancy']==pytest.approx(2/3)
    assert report['edges'][0]['denominator_frames']==3
    assert all(Path(p).is_file() for p in files.values())


def test_hydration_palette_is_absolute_and_missing_is_grey():
    assert hydration_rgb(0)==[1,0,0]
    assert hydration_rgb(1)==[0,0,1]
    assert hydration_rgb(.5)==[.5,0,.5]
    assert hydration_rgb(None)==[.55,.55,.55]
    for value in [float('nan'),-.1,1.1]:
        with pytest.raises(ValueError):hydration_rgb(value)


from pathlib import Path


def test_ambiguous_residue_chemistry_is_reported_without_breaking_other_analyses():
    frame=StructureFrame((atom(1,'NZ',(0,0,0)),atom(2,'NZ',(2,0,0))))
    model=chemical(frame)
    assert model.metadata()['ambiguous_residues_excluded']==['A:LYS1']
    assert not len(model.donor_pairs)
    analysis=InteractionAnalysis(model,['A:LYS1'],1);analysis.add([a.coord for a in frame.atoms])
    assert np.isnan(analysis.result()['hydrogen_bonded_water_count']).all()


def test_protonated_acid_without_its_hydrogen_does_not_guess_acceptor_oxygen():
    frame=StructureFrame((atom(1,'OD1',(0,0,0),resname='ASH',element='O'),atom(2,'OD2',(0,2,0),resname='ASH',element='O')))
    model=chemical(frame)
    assert not len(model.acceptors)


def test_complete_static_views_and_explicit_opt_out_do_not_reuse_stale_chemistry(tmp_path):
    from crevice.hydration_workflow import static_hydration_bundle
    atoms=[atom(1,'N',(0,0,0),resname='ALA'),atom(2,'CA',(1.45,0,0),resname='ALA',element='C'),
           atom(3,'C',(2,1.4,0),resname='ALA',element='C'),atom(4,'O',(1.8,2.5,0),resname='ALA',element='O'),
           atom(5,'O',(0,0,3),resid=2,resname='HOH',element='O',hetero=True)]
    frame=StructureFrame(tuple(atoms));prefix='toy sample'
    files=static_hydration_bundle(frame,tmp_path,prefix=prefix,sasa_samples=32,dpi=60)
    for name in ['hydration_pml','hydration_vmd','hydration_tcl','hydration_cxc','network_pml','network_vmd','network_cxc','water_density_dx','analysis_evidence_csv',
                 'hydration_frames_csv','residue_conformation_frames_csv','interaction_edge_frames_csv']:
        assert Path(files[name]).is_file()
    assert not list(tmp_path.glob('*.html')) and 'analysis_html' not in files and 'analysis_report_data_json' not in files
    assert Path(files['analysis_evidence_csv']).name.startswith('toy sample')
    report=json.loads((tmp_path/(prefix+'_interactions.json')).read_text())
    assert report['residues'][0]['mean_hydrogen_bonded_waters'] is None
    files=static_hydration_bundle(frame,tmp_path,prefix=prefix,sasa_samples=32,dpi=60,analysis_options={'skip_interactions':True,'skip_density':True})
    # The stale interaction and density files of the first run are not reused.
    assert 'interaction_json' not in files and 'water_density_json' not in files
    scene=json.loads(Path(files['analysis_scene_json']).read_text())
    assert scene['contacts']==[] and not scene['surfaces']
    import csv
    with open(files['analysis_evidence_csv'],newline='',encoding='utf-8') as handle:
        rows=list(csv.DictReader(handle))
    assert rows and all(r['observed_typed_edge_count']=='0' for r in rows)
