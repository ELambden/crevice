"""Named-region provenance, independent geometry checks and CLI contracts."""
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pytest
from crevice.region_definition import (file_record,load_definition,validate_definition,
    frame_contract,resolve_landmarks,prepare_region,trajectory_arguments,write_json)
from crevice.cli import build_parser,main
from crevice.water_density import write_scalar_dx
from test_cavity_obstacles import md_pair


def definition(tmp_path):
    gro,xtc,traj,ref=md_pair(tmp_path)
    dx=tmp_path/'reference.dx';write_scalar_dx(ref.grid,ref.origin,ref.spacing,dx,basis=ref.basis)
    frame=traj.frames[0]
    # Opposite protein corners have a midpoint inside this symmetric synthetic cavity.
    xyz=np.array([a.coord for a in frame.atoms]);center=np.array([20.,20.,20.])
    i=int(np.argmin(xyz.sum(1)));j=int(np.argmin(np.linalg.norm(xyz-(2*center-xyz[i]),axis=1)))
    rows=[{'residue':frame.atoms[k].residue_key.label,'atoms':['C'],'evidence':'Synthetic opposite-wall geometry control'} for k in [i,j]]
    data={'schema_version':1,'region_id':'synthetic','label':'Synthetic candidate',
        'review':{'status':'candidate','rationale':'Software geometry control only','sources':[]},
        'system':{'topology':file_record(gro,tmp_path),'trajectory':file_record(xtc,tmp_path),'reference_frame':0,
                  'protein_selection':'protein','obstacle_selection':'protein','pbc':'none','assembly':'synthetic','preparation':'synthetic','membrane':'not applicable'},
        'alignment':{'residues':None,'rationale':'Synthetic matching heavy atoms'},
        'geometry':{'spacing_A':.5,'probe_radius_A':.8,'region_margin_A':1.,'axis':[0,0,1],
                    'grid_phase':[0,0,0],'profile_padding_A':2.,'max_grid_points':1000000},
        'reference':{'kind':'existing_dx','asset':file_record(dx,tmp_path)},'landmarks':rows}
    path=tmp_path/'definition.json';write_json(path,data)
    return path,data,traj


@pytest.mark.parametrize('change,match',[
    (lambda d:d.update(schema_version=True),'schema_version'),
    (lambda d:d.update(unknown=1),'unknown fields'),
    (lambda d:d.update(region_id='../escape'),'region_id'),
    (lambda d:d['geometry'].update(grid_phase=[1,0,0]),'grid_phase'),
    (lambda d:d['geometry'].update(probe_radius_A=float('nan')),'finite'),
    (lambda d:d['review'].update(status='reviewed'),'reviewer'),
    (lambda d:d['alignment'].update(residues=['A:ALA1','A:ALA1']),'unique'),
])
def test_invalid_definition_fails_before_geometry(tmp_path,change,match):
    _,data,_=definition(tmp_path);change(data)
    with pytest.raises(ValueError,match=match):validate_definition(data)


def test_input_hash_and_duplicate_json_keys_rejected(tmp_path):
    path,data,_=definition(tmp_path)
    data['system']['topology']['sha256']='0'*64;write_json(path,data)
    with pytest.raises(ValueError,match='SHA256'):load_definition(path)
    path.write_text('{"schema_version":1,"schema_version":1}')
    with pytest.raises(ValueError,match='Duplicate JSON'):load_definition(path)


def test_exact_landmark_identity_no_residue_number_fallback(tmp_path):
    _,data,traj=definition(tmp_path)
    rows=resolve_landmarks(traj.frames[0],data['landmarks'])
    assert rows[0]['residue_atom_serials']
    wrong=deepcopy(data['landmarks']);wrong[0]['residue']=wrong[0]['residue'].replace('ALA','GLY')
    with pytest.raises(ValueError,match='resolves to 0'):resolve_landmarks(traj.frames[0],wrong)
    from dataclasses import replace
    frame=traj.frames[0];ambiguous=replace(frame,atoms=frame.atoms+(frame.atoms[rows[0]['anchor_atom_serials'][0]-1],))
    with pytest.raises(ValueError,match='resolves to 2'):resolve_landmarks(ambiguous,data['landmarks'])


def test_prepared_region_roundtrip_and_reference_frame_binding(tmp_path):
    path,data,traj=definition(tmp_path)
    out=tmp_path/'review with spaces';files=prepare_region(path,out)
    prepared=Path(files['definition_json']);resolved,assets=load_definition(prepared)
    assert resolved['review']['status']=='candidate'
    assert resolved['reference']['asset']['sha256']==data['reference']['asset']['sha256']
    argv,metadata=trajectory_arguments(prepared,tmp_path/'analysis')
    args=build_parser().parse_args(argv)
    assert args.spacing==.5 and args.region_margin==1. and args.geometry_mode=='reference-region'
    assert metadata['definition']['prepared']['frame_contract']==frame_contract(traj.frames[0])
    scene=json.loads(Path(files['review_scene_json']).read_text())
    assert set(i for r in scene['landmarks'] for i in r['residue_atom_serials'])==set(i for r in resolve_landmarks(traj.frames[0],data['landmarks']) for i in r['residue_atom_serials'])
    for ext in ['pml','tcl','vmd','cxc']:assert Path(files['review_'+ext]).is_file()
    with pytest.raises(ValueError,match='empty output'):prepare_region(path,out)
    resolved['prepared']['frame_contract']['coordinates_float32_sha256']='0'*64;write_json(prepared,resolved)
    with pytest.raises(ValueError,match='coordinate frame'):trajectory_arguments(prepared,tmp_path/'analysis')


def test_landmark_proposal_is_connected_and_excludes_atom_cores(tmp_path):
    path,data,traj=definition(tmp_path)
    data['reference']={'kind':'landmark_sphere','radius_A':3.,'max_seed_distance_A':2.,'ambiguity_distance_A':.25}
    write_json(path,data);files=prepare_region(path,tmp_path/'proposal')
    from crevice.residue_evidence import read_binary_dx
    from scipy.ndimage import label
    from crevice.radii import atom_vdw_radius
    from crevice.rolling import _SphereQueries
    grid,origin,deltas=read_binary_dx(files['reference_dx']);xyz=np.argwhere(grid)@deltas+origin
    atoms=traj.frames[0].atoms
    q=_SphereQueries(np.array([a.coord for a in atoms]),np.array([atom_vdw_radius(a) for a in atoms]))
    assert q.points(xyz)[0].min()>=-1e-8
    report=json.loads(Path(files['review_json']).read_text());anchor=np.asarray(report['proposal']['anchor_A'])
    assert np.linalg.norm(xyz-anchor,axis=1).max()<=3.+1e-8
    assert label(grid)[1]==1
    assert report['review']['status']=='candidate'
    assert report['proposal']['selected']['nearest_distance_A']<=2


def test_region_cli_runs_each_frame_and_retains_region_provenance(tmp_path):
    path,_,_=definition(tmp_path);files=prepare_region(path,tmp_path/'review')
    out=tmp_path/'run'
    code=main(['region-trajectory',files['definition_json'],'--out-dir',str(out),'--skip-hydration','--max-frames','2','--confidence','.95','--dpi','60'])
    assert code==0
    report=json.loads((out/'synthetic_cavity_statistics.json').read_text())
    assert report['frame_count']==2
    assert report['settings']['region_definition']['region_id']=='synthetic'
    assert report['settings']['region_definition']['review']['status']=='candidate'
    a=np.load(out/'synthetic_cavity_frames.npz')
    assert a['volume_A3'][0]==a['volume_A3'][1]
    assert main(['region-trajectory',files['definition_json'],'--out-dir',str(out),'--skip-hydration'])==2


def test_static_mmcif_preserves_canonical_identity_and_nonprotein_obstacles(tmp_path):
    import gemmi
    from crevice.presentation import write_viewer_structure
    from crevice.parser import load_structure
    path,data,traj=definition(tmp_path)
    complete=load_structure(tmp_path/'system.gro',md_selection='all')
    pdb=write_viewer_structure(complete,tmp_path/'model.pdb')
    cif=tmp_path/'model.cif';model=gemmi.read_structure(str(pdb));model.setup_entities();model.assign_label_seq_id();model.make_mmcif_document().write_file(str(cif))
    static=load_structure(cif)
    data['system']['topology']=file_record(cif,tmp_path);data['system']['trajectory']=None
    data['system']['obstacle_selection']='all'
    data['landmarks']=[{'residue':static.atoms[0].residue_key.label,'atoms':['C'],'evidence':'Synthetic mmCIF identity control'}]
    write_json(path,data);files=prepare_region(path,tmp_path/'mmcif-review')
    prepared,_=load_definition(files['definition_json'])
    report=json.loads(Path(files['review_json']).read_text())
    assert prepared['landmarks'][0]['residue']==static.atoms[0].residue_key.label
    assert report['reference_frame']['atom_count']==len(static.atoms)-1
    assert report['reference_samples_inside_obstacle_spheres']>0
    assert 'context_pdb' in files
    assert report['obstacles']['periodic'] is False
    with pytest.raises(ValueError,match='requires a declared trajectory'):trajectory_arguments(files['definition_json'],tmp_path/'run')


def test_proposal_rejects_a_seed_neighbourhood_inside_an_atom(tmp_path):
    path,data,_=definition(tmp_path)
    data['landmarks']=data['landmarks'][:1]
    data['reference']={'kind':'landmark_sphere','radius_A':2.,'max_seed_distance_A':.2,'ambiguity_distance_A':0.}
    write_json(path,data)
    with pytest.raises(ValueError,match='No admissible probe center'):prepare_region(path,tmp_path/'blocked-seed')
    assert not list((tmp_path/'blocked-seed').glob('*_reference.dx'))


def test_comparison_rejects_unpaired_frames_and_preserves_region_identity(tmp_path):
    from crevice.region_comparison import compare_regions
    path,data,_=definition(tmp_path)
    files=[]
    for name,stop in [('left',2),('right',2),('short',1)]:
        data['region_id']=name;write_json(path,data)
        prepared=prepare_region(path,tmp_path/(name+'-review'))
        out=tmp_path/name
        assert main(['region-trajectory',prepared['definition_json'],'--out-dir',str(out),'--skip-hydration','--stop',str(stop),'--dpi','60'])==0
        files.append(out/(name+'_cavity_statistics.json'))
    result=compare_regions(files[0],files[1],tmp_path/'comparison')
    assert result['frame_count']==2
    assert result['left']['region_id']=='left' and result['right']['region_id']=='right'
    assert result['paired_volume_difference_A3']['mean']==0
    assert result['paired_volume_difference_A3']['confidence_interval']['status']=='unavailable_observed_constant'
    assert result['hydration_comparison_status']=='unavailable_missing_hydration_bundle'
    assert result['profile_uncertainty']['left']['section_positions_with_mean_intervals']==0
    from crevice.region_comparison import _load_result
    profile_file=files[1].with_name('right_cavity_profile_statistics.json')
    profile=json.loads(profile_file.read_text());original=deepcopy(profile)
    profile['profiles']['section']['rows'][0]['position_A']+=1;write_json(profile_file,profile)
    with pytest.raises(ValueError,match='Pointwise profile observations'):_load_result(files[1])
    write_json(profile_file,original)

    with pytest.raises(ValueError,match='identical frames'):compare_regions(files[0],files[2],tmp_path/'mismatch')
    from crevice.region_comparison import _load_result
    report=json.loads(files[0].read_text());report['settings'].pop('region_definition');write_json(files[0],report)
    prepared=tmp_path/'left-review/left_region.json'
    with pytest.raises(ValueError,match='historical input provenance'):_load_result(files[0],prepared)
    d,assets=load_definition(prepared)
    provenance={str(assets[k]):d['system'][k] for k in ['topology','trajectory']}
    record=tmp_path/'historical-inputs.json';write_json(record,provenance)
    assert _load_result(files[0],prepared,record)['report']['frame_count']==2
    provenance[str(assets['trajectory'])]['sha256']='0'*64;write_json(record,provenance)
    with pytest.raises(ValueError,match='Historical input provenance'):_load_result(files[0],prepared,record)



def test_equally_near_disconnected_probe_components_remain_ambiguous():
    from dataclasses import replace
    from types import SimpleNamespace
    from crevice.models import Atom
    from crevice.region_definition import propose_landmark_region
    from crevice.radii import atom_vdw_radius
    from test_channel_coordinate import rectangular_channel
    frame=rectangular_channel(half_x=8,half_y=8,half_length=8,capped=True)
    plane=tuple(Atom(10000+i,'C','ALA','B',i,0.,float(y),float(z),'C')
                for i,(y,z) in enumerate((y,z) for y in range(-4,5,2) for z in range(-4,5,2)))
    frame=replace(frame,atoms=frame.atoms+plane)
    coords=np.asarray([a.coord for a in frame.atoms]);radii=np.array([atom_vdw_radius(a) for a in frame.atoms])
    source=SimpleNamespace(atoms=frame.atoms,radii=radii,snapshot=lambda i:(coords,None))
    landmarks=[{'centroid_A':[-8,0,0]},{'centroid_A':[8,0,0]}]
    g={'spacing_A':.25,'probe_radius_A':.8,'axis':[0,0,1],'grid_phase':[0,0,0],'max_grid_points':1000000}
    settings={'radius_A':3.5,'max_seed_distance_A':3.5,'ambiguity_distance_A':.25}
    with pytest.raises(ValueError,match='Ambiguous landmark seed'):
        propose_landmark_region(frame,landmarks,settings,g,source)
