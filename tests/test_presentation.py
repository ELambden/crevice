"""Display extensions, geometric mouths, GRO filtering and native scene contracts."""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
from dataclasses import replace
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from crevice import Atom, StructureFrame, pore_profile, void_cast
from crevice.parser import load_structure_report, load_structure
from crevice.presentation import extend_channel_cast, orthogonal_map_cast, write_viewer_structure
from crevice.spatial import SpatialIndex
from crevice.volume_export import cast_grid, write_volume_viewer_bundle
from test_channel_coordinate import rectangular_channel, transform


@pytest.fixture(scope='module')
def channel_display():
    frame = rectangular_channel()
    profile = pore_profile(frame, enclosure_radius=0.8, search_radius=8, samples=41)
    cast = void_cast(frame, profile=profile, spacing=.5, min_radius=0, max_grid_points=200000)
    return frame, profile, cast


def test_mouth_brackets_capture_end_sphere_enclosure(channel_display):
    frame, profile, _ = channel_display
    mouths = profile.metadata['channel_mouths']
    # End-wall spheres (r + q = 2.5), 1 A atom separation: last touching disks.
    expected = 10+math.sqrt(2.5**2-.5**2)
    for label, sign in [('lower',-1),('upper',1)]:
        mouth = mouths[label]
        assert mouth['axial_bracket_width_A'] <= .025
        assert mouth['bracket_A'][0] <= mouth['coordinate_A'] <= mouth['bracket_A'][1]
        assert mouth['coordinate_A'] == pytest.approx(sign*expected, abs=.04)


def test_extension_does_not_change_measurement_and_is_atom_clear(channel_display):
    frame, profile, cast = channel_display
    before = cast.to_dict(include_points=True)
    shown = extend_channel_cast(cast, frame, 2)
    assert shown.total_volume > cast.total_volume
    assert cast.to_dict(include_points=True) == before
    assert shown.metadata['measured_volume_A3'] == cast.total_volume
    assert shown.metadata['extension_actual_A'] == [2,2]
    original = {p.position for c in cast.components for p in c.points}
    assert original <= {p.position for c in shown.components for p in c.points}
    spatial = SpatialIndex(frame.atoms)
    for c in shown.components:
        for p in c.points:
            assert spatial.nearest_surface(p.position)[0] >= -1e-10
            if p.position not in original:
                sign = 1 if p.t > profile.points[-1].t else -1
                predecessor = tuple(x-sign*cast.spacing*a for x,a in zip(p.position,profile.axis_direction))
                assert spatial.segment_clearance(predecessor,p.position) >= -1e-10


def test_extension_respects_zero_crop_and_invalid_length(channel_display):
    frame, _, cast = channel_display
    assert extend_channel_cast(cast,frame,0) is cast
    crop = replace(cast,metadata={**cast.metadata,'focus_points':[(0,0,0)]})
    assert extend_channel_cast(crop,frame,2) is crop
    for length in (-1,float('nan'),float('inf')):
        with pytest.raises(ValueError,match='cast_extension'):
            extend_channel_cast(cast,frame,length)


def test_positive_pymol_grid_preserves_every_occupied_sample(channel_display):
    frame, profile, cast = channel_display
    moved, rotation, shift = transform(frame)
    p = pore_profile(moved, enclosure_radius=0.8, search_radius=8, samples=41)
    other = void_cast(moved,profile=p,spacing=.5,min_radius=0,max_grid_points=200000)
    local = orthogonal_map_cast(other)
    assert np.array_equal(cast_grid(other)[0],cast_grid(local)[0])
    basis = np.asarray(other.metadata['grid_basis'])
    assert np.allclose(np.asarray([v.position for c in local.components for v in c.points]) @ basis,
                       [v.position for c in other.components for v in c.points])


def test_viewer_structure_preserves_selected_model_and_chain_identity(tmp_path):
    frame = StructureFrame((Atom(17,'CA','ALA','long_chain',12,1.23456,2,3,'C'),
                            Atom(98,'CA','GLY','A',12,5,6,7,'C')))
    path = write_viewer_structure(frame,tmp_path/'viewer.pdb')
    read = load_structure(path)
    assert len(read.atoms) == 2
    assert read.atoms[0].chain_id != read.atoms[1].chain_id
    assert np.allclose([a.coord for a in read.atoms],[a.coord for a in frame.atoms],atol=.0005)
    mapping = json.loads(path.with_suffix('.identities.json').read_text())
    assert mapping['original_serials_in_viewer_order'] == [17,98]
    assert mapping['chain_mapping']['A'] == 'A'


def mixed_gro(path):
    rows=[(1,'ALA','CA',1,(.1,.2,.3)),(1,'ALA','HA',2,(.2,.2,.3)),
          (2,'TIP3','OH2',3,(.3,.2,.3)),(3,'SOD','SOD',4,(.4,.2,.3)),
          (4,'DOPC','C1',5,(.5,.2,.3))]
    path.write_text('mixed selection\n5\n'+''.join(
        f'{r:5d}{res:<5}{name:>5}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n'
        for r,res,name,i,(x,y,z) in rows)+'   3.0 3.0 3.0\n')
    return path


def test_gro_filters_solvent_lipid_ion_and_converts_units(tmp_path):
    frame, report = load_structure_report(mixed_gro(tmp_path/'mixed.gro'))
    assert [a.resname for a in frame.atoms] == ['ALA','ALA']
    assert len(frame.selected_atoms()) == 1
    assert frame.atoms[0].coord == pytest.approx((1,2,3))
    assert report['excluded_residue_atom_counts'] == {'DOPC':1,'SOD':1,'TIP3':1}
    all_frame, _ = load_structure_report(tmp_path/'mixed.gro',md_selection='all')
    assert len(all_frame.atoms) == 5
    assert len(all_frame.selected_atoms(include_hetero=False)) == 1
    with pytest.raises(ValueError,match='matched no atoms'):
        load_structure_report(tmp_path/'mixed.gro',md_selection='resname XXX')


@pytest.mark.parametrize('split_regions',[False,True])
def test_native_pymol_relocated_rotated_cast(channel_display,tmp_path,split_regions):
    pytest.importorskip('pymol',reason='Native PyMOL not installed in this interpreter')
    frame, profile, cast = channel_display
    moved, _, _ = transform(frame)
    p = pore_profile(moved, enclosure_radius=0.8,search_radius=8,samples=21)
    c = void_cast(moved,profile=p,spacing=.5,min_radius=0,max_grid_points=200000)
    initial=tmp_path/'initial'
    initial.mkdir()
    structure=write_viewer_structure(moved,initial/'viewer.pdb')
    shown=extend_channel_cast(c,moved,2)
    points=tuple(point for component in shown.components for point in component.points)
    cut=(min(point.t for point in points)+max(point.t for point in points))/2
    regions={'lower':tuple(point for point in points if point.t<=cut),
             'upper':tuple(point for point in points if point.t>cut)} if split_regions else None
    files=write_volume_viewer_bundle(c,structure_path=structure,output_dir=initial,
                                    frame=moved,profile=p,regions=regions)
    relocated=tmp_path/'relocated space'
    shutil.move(str(initial),relocated)
    runner=Path(__file__).resolve().parents[1]/'scripts/check_pymol.py'
    done=subprocess.run([sys.executable,str(runner),str(relocated/'crevice_volume.pml'),
                         '--output',str(tmp_path/'render')],capture_output=True,text=True,timeout=200)
    assert done.returncode == 0, done.stdout+done.stderr
    report=json.loads((tmp_path/'render/pymol_verification.json').read_text())
    assert report['source_grid_occupied_samples'] > c.point_count
    assert report['helper_grid_objects']==[]
    assert report['protein_representation_atoms']['sticks']==0
    assert report['protein_representation_atoms']['lines']==0
    assert report['native_expected_surface_max_error_A']<.002
    assert report['alignment_max_error_A'] < .001


def test_native_geometry_rejects_matching_wrong_transform_metadata(channel_display,tmp_path):
    pytest.importorskip('pymol',reason='Native PyMOL not installed in this interpreter')
    frame,profile,cast=channel_display
    moved,_,_=transform(frame)
    p=pore_profile(moved, enclosure_radius=0.8,search_radius=8,samples=21)
    c=void_cast(moved,profile=p,spacing=.5,min_radius=0,max_grid_points=200000)
    root=tmp_path/'wrong orientation'
    root.mkdir()
    structure=write_viewer_structure(moved,root/'viewer.pdb')
    files=write_volume_viewer_bundle(c,structure_path=structure,output_dir=root,frame=moved,profile=p)
    metadata_path=Path(files['volume_metadata_json'])
    metadata=json.loads(metadata_path.read_text())
    old=metadata['pymol_local_to_input_matrix']
    wrong=np.eye(4).ravel().tolist()
    scene=Path(files['volume_pml'])
    scene.write_text(scene.read_text().replace(repr(old),repr(wrong)))
    metadata['pymol_local_to_input_matrix']=wrong
    metadata_path.write_text(json.dumps(metadata))
    runner=Path(__file__).resolve().parents[1]/'scripts/check_pymol.py'
    done=subprocess.run([sys.executable,str(runner),str(scene),'--output',str(tmp_path/'render')],capture_output=True,text=True,timeout=200)
    report=json.loads((tmp_path/'render/pymol_verification.json').read_text())
    assert done.returncode==1
    assert report['alignment_max_error_A']<.001  # the former check is fooled
    assert report['extent_max_error_A']<.001
    assert report['native_surface_sample_error_A']>c.spacing


def test_display_smoothing_preserves_samples_necks_and_empty_pockets():
    from crevice.volume_export import presentation_grid
    from scipy.ndimage import label
    grid=np.zeros((19,19,19),dtype=np.float32)
    grid[2:8,5:14,5:14]=1; grid[11:17,5:14,5:14]=1
    grid[7:12,9,9]=1  # one-voxel neck
    grid[4,8,8]=0      # one-voxel empty pocket
    for width in (.3,.6,1.2):
        shown=presentation_grid(grid,.5,width)
        assert np.array_equal(shown>.5,grid>.5)
        assert label(shown>.5)[1]==1
        assert shown[4,8,8]<.5 and shown[9,9,9]>.5
        assert np.any((shown>0)&(shown<1))
    assert np.array_equal(presentation_grid(grid,.5,0),grid)
    for bad in (-1,float('nan'),float('inf')):
        with pytest.raises(ValueError,match='surface_smoothing'):
            presentation_grid(grid,.5,bad)


def test_display_smoothing_leaves_measured_dx_identical(channel_display,tmp_path):
    frame,profile,cast=channel_display
    structure=write_viewer_structure(frame,tmp_path/'viewer.pdb')
    before=cast.to_dict(include_points=True)
    raw=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=tmp_path/'raw',
                                  frame=frame,profile=profile,surface_smoothing=0)
    smooth=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=tmp_path/'smooth',
                                     frame=frame,profile=profile,surface_smoothing=.6)
    assert Path(raw['volume_dx']).read_bytes()==Path(smooth['volume_dx']).read_bytes()
    assert Path(raw['display_dx']).read_bytes()!=Path(smooth['display_dx']).read_bytes()
    assert cast.to_dict(include_points=True)==before


def test_native_residue_overlay_survives_atom_sorting_and_relocation(channel_display,tmp_path):
    pytest.importorskip('pymol',reason='Native PyMOL not installed in this interpreter')
    from dataclasses import replace
    from crevice.residue_evidence import read_binary_dx,boundary_residue_evidence,write_residue_evidence_bundle
    frame,profile,cast=channel_display
    frame=replace(frame,atoms=tuple(replace(a,serial=a.serial*10) for a in frame.atoms))
    root=tmp_path/'initial';root.mkdir()
    structure=write_viewer_structure(frame,root/'viewer.pdb')
    files=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=root,frame=frame,profile=profile)
    # Deliberately reverse analysis atom order; the PDB serial identities match.
    frame=replace(frame,atoms=tuple(reversed(frame.atoms)))
    g,o,d=read_binary_dx(files['volume_dx']);report=boundary_residue_evidence(frame,g,o,d)
    bundle=write_residue_evidence_bundle(report,frame,root,prefix='crevice',scene_path=files['volume_pml'],dpi=70)
    relocated=tmp_path/'relocated overlay';shutil.move(str(root),relocated)
    runner=Path(__file__).resolve().parents[1]/'scripts/check_pymol.py'
    done=subprocess.run([sys.executable,str(runner),str(relocated/'crevice_residue_context.pml'),'--output',str(tmp_path/'render')],capture_output=True,text=True,timeout=200)
    assert done.returncode==0,done.stdout+done.stderr
    result=json.loads((tmp_path/'render/pymol_verification.json').read_text())
    assert result['residue_highlights_passed']
    assert result['residue_sticks_passed']
    assert result['residue_highlights']['lining']['expected_atoms']>0


def test_native_vmd_residue_sticks_relocate_and_reject_coordinate_mismatch(channel_display,tmp_path):
    if shutil.which('vmd') is None:pytest.skip('Native VMD not installed')
    from crevice.residue_evidence import read_binary_dx,boundary_residue_evidence,write_residue_evidence_bundle
    frame,profile,cast=channel_display
    root=tmp_path/'initial';root.mkdir()
    structure=write_viewer_structure(frame,root/'viewer.pdb')
    files=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=root,frame=frame,profile=profile)
    g,o,d=read_binary_dx(files['volume_dx']);report=boundary_residue_evidence(frame,g,o,d)
    write_residue_evidence_bundle(report,frame,root,prefix='crevice',scene_path=files['volume_pml'],dpi=60)
    relocated=tmp_path/'relocated context';shutil.move(str(root),relocated)
    runner=Path(__file__).resolve().parents[1]/'scripts/check_vmd.py'
    command=[sys.executable,str(runner),str(relocated/'crevice_residue_context.vmd'),'--output',str(tmp_path/'render')]
    done=subprocess.run(command,capture_output=True,text=True,timeout=180)
    assert done.returncode==0,done.stdout+done.stderr
    native=json.loads((tmp_path/'render/vmd_verification.json').read_text())
    assert native['residue_sticks_passed']
    reference=relocated/'crevice_residue_reference.tsv';lines=reference.read_text().splitlines()
    fields=lines[0].split();fields[1]=str(float(fields[1])+2);lines[0]='\t'.join(fields);reference.write_text('\n'.join(lines)+'\n')
    command[-1]=str(tmp_path/'rejected')
    done=subprocess.run(command,capture_output=True,text=True,timeout=60)
    assert done.returncode!=0
    assert 'Residue reference coordinates differ' in (tmp_path/'rejected/vmd.log').read_text()


def test_mouth_guides_are_opt_in_and_non_channel_casts_are_violet(channel_display, tmp_path):
    from crevice.presentation import (CHANNEL_CAST_RGB, NON_CHANNEL_CAST_RGB, NON_CHANNEL_CAST_HEX,
                                      cast_rgb, display_guides)
    frame, profile, cast = channel_display
    assert cast_rgb(cast) == CHANNEL_CAST_RGB
    plain = write_volume_viewer_bundle(cast, structure_path=tmp_path/'p.pdb', output_dir=tmp_path/'plain',
                                       frame=frame, profile=profile)
    scene = json.loads(Path(plain['scene_json']).read_text())
    assert scene['mouth_rings'] and not scene['mouth_rings_drawn']
    assert 'radius 0.065' not in Path(plain['volume_tcl']).read_text()
    assert '.cylinder' not in Path(plain['mouth_guides_bild']).read_text()
    assert 'mouths.bild' not in Path(plain['volume_cxc']).read_text()
    meta = json.loads(Path(plain['volume_metadata_json']).read_text())
    assert meta['channel_mouths'] and meta['mouth_guides_drawn'] is False
    with display_guides(mouth_guides=True):
        guided = write_volume_viewer_bundle(cast, structure_path=tmp_path/'p.pdb', output_dir=tmp_path/'rings',
                                            frame=frame, profile=profile)
    assert json.loads(Path(guided['scene_json']).read_text())['mouth_rings_drawn']
    assert 'radius 0.065' in Path(guided['volume_tcl']).read_text()
    assert 'mouths.bild' in Path(guided['volume_cxc']).read_text()
    # Every non-channel cast takes the violet default, in all three viewers.
    other = replace(cast, mode='rolling')
    assert cast_rgb(other) == NON_CHANNEL_CAST_RGB
    files = write_volume_viewer_bundle(other, structure_path=tmp_path/'p.pdb', output_dir=tmp_path/'other',
                                      frame=frame)
    assert f'color #2 {NON_CHANNEL_CAST_HEX}' in Path(files['volume_cxc']).read_text()
    assert repr(list(NON_CHANNEL_CAST_RGB)) in Path(files['volume_pml']).read_text()
    assert 'color change rgb 10 0.659 0.333 0.969' in Path(files['volume_tcl']).read_text()


def test_profile_mouth_brackets_are_opt_in(channel_display, tmp_path):
    from PIL import Image
    from crevice.figures import plot_profile_radius
    from crevice.presentation import display_guides, figure_description
    _, profile, _ = channel_display

    def orange(path):
        image = np.asarray(Image.open(path).convert('RGB')).astype(int)
        return int(((abs(image[..., 0]-0xd9) < 12) & (abs(image[..., 1]-0x82) < 12) & (abs(image[..., 2]-0x2b) < 12)).sum())
    plot_profile_radius(profile, tmp_path/'plain.png', dpi=60)
    assert 'mouth guides not drawn' in figure_description(tmp_path/'plain.png')['Description']
    with display_guides(mouth_guides=True):
        plot_profile_radius(profile, tmp_path/'guides.png', dpi=60)
    assert orange(tmp_path/'plain.png') == 0 < orange(tmp_path/'guides.png')
