"""Synthetic access/identity contracts; no biological substrate validation."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pytest
from scipy.spatial import cKDTree
from crevice import Atom,StructureFrame
from crevice.rolling import rolling_probe_cast
from crevice.residue_evidence import boundary_residue_evidence,read_binary_dx,write_residue_evidence_bundle
from crevice.volume_export import write_volume_viewer_bundle
from test_channel_coordinate import transform


@pytest.fixture(scope='module')
def two_rooms():
    xyz=set()
    for x in range(-8,9):
        for y in range(-4,5):
            for z in range(-5,6):
                if abs(x)==8 or abs(y)==4 or abs(z)==5 or (x==0 and y*y+z*z>=10):xyz.add((x,y,z))
    return StructureFrame(tuple(Atom(i+1,'C','ALA','A',i+1,*p,'C') for i,p in enumerate(sorted(xyz))))


def seeded(frame,probe=.8,seed=(-4,0,0)):
    return rolling_probe_cast(frame,spacing=.5,probe_radius=probe,outer_radius=6,
                              seed=seed,min_component_volume=1,selection_mode='all',enclosure_fraction=0)


def xyz(cast):return np.array([p.position for c in cast.components for p in c.points])


def test_seed_retains_closed_room_and_respects_narrow_neck(two_rooms):
    small=seeded(two_rooms);large=seeded(two_rooms,1.8)
    assert xyz(small)[:,0].max()>4
    assert xyz(large)[:,0].max()<1
    assert cKDTree(xyz(large)).query([[-4,0,0]])[0][0]<.8
    assert len(large.components)==1
    assert large.metadata['seed_selection']['snap_distance_A']<=.75


def test_seeded_cast_rigid_transform(two_rooms):
    moved,rotation,shift=transform(two_rooms)
    a=seeded(two_rooms,1.8);b=seeded(moved,1.8,np.array([-4,0,0])@rotation.T+shift)
    assert a.total_volume==b.total_volume
    assert cKDTree(xyz(a)@rotation.T+shift).query(xyz(b))[0].max()<1e-7


@pytest.mark.parametrize('seed,match',[((0,4,0),'accommodate'),((100,100,100),'close enough'),((0,float('nan'),0),'finite')])
def test_invalid_or_unreachable_seed_fails(two_rooms,seed,match):
    with pytest.raises(ValueError,match=match):seeded(two_rooms,seed=seed)


def test_measured_map_without_scene_creates_all_contexts_and_preserves_samples(two_rooms,tmp_path):
    from crevice.presentation import write_viewer_structure
    cast=seeded(two_rooms,1.8)
    structure=write_viewer_structure(two_rooms,tmp_path/'source.pdb')
    files=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=tmp_path/'base',prefix='test',frame=two_rooms)
    grid,origin,delta=read_binary_dx(files['volume_dx'])
    evidence=boundary_residue_evidence(two_rooms,grid,origin,delta)
    paths=write_residue_evidence_bundle(evidence,two_rooms,tmp_path/'context',prefix='test',volume_path=files['volume_dx'],dpi=60)
    actual,aorigin,adelta=read_binary_dx(paths['volume_dx'])
    before=np.argwhere(grid)@delta+origin;after=np.argwhere(actual)@adelta+aorigin
    assert len(before)==len(after)
    assert cKDTree(before).query(after)[0].max()<1e-7
    for extension in ['pml','vmd','tcl','cxc']:
        assert Path(paths['residue_context_'+extension]).is_file()
    report=json.loads(Path(paths['residue_evidence_json']).read_text())
    assert report['overlay_representation']=='cartoon_with_evidence_sticks'
    assert report['display_atom_serials']['lining']


def test_cast_cli_passes_explicit_seed_and_probe(two_rooms,tmp_path):
    from crevice.cli import main
    from crevice.presentation import write_viewer_structure
    structure=write_viewer_structure(two_rooms,tmp_path/'rooms.pdb')
    out=tmp_path/'cli'
    assert main(['cast',str(structure),'--out-dir',str(out),'--prefix','rooms',
                 '--seed=-4,0,0','--probe-radius','1.8','--cavity-selection','all',
                 '--enclosure-fraction','0','--min-component-volume','1','--dpi','60'])==0
    cast=json.loads((out/'rooms_void_cast.json').read_text())
    assert cast['metadata']['seed_selection']['coordinate_A']==[-4.,0.,0.]
    assert cast['total_volume']==pytest.approx(seeded(two_rooms,1.8).total_volume)
    for extension in ['pml','vmd','tcl','cxc']:assert (out/f'rooms_pore_cast.{extension}').is_file()
