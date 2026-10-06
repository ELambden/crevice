"""Render native PyMOL scenes and verify surfaces, coordinates and presentation.

Runs a cast or residue-evidence scene (``PREFIX_volume.pml`` or
``PREFIX_residue_context.pml``) in PyMOL (``python -m pymol -cq``, so the
``pymol`` module must be importable), ray-traces an image and checks the
displayed triangles against the expected surface (``*_geometry_reference.npz``
or the PyMOL mesh), records the protein representations, checks the stick residues
and colours of an evidence overlay, the measured grid, and that the image
shows cast pixels (teal channel or violet non-channel cast).

Usage::

    python scripts/check_pymol.py OUT/PREFIX_volume.pml --output check/

Writes ``pymol_native.png``, ``native_geometry.npz`` and
``pymol_verification.json``; exits with status 1 unless the checks pass.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    """Parse the command line, run PyMOL and write ``pymol_verification.json``.

    Options: ``scene``, ``--output`` (directory), ``--timeout`` (seconds, default
    300) and ``--width`` (pixels, default 1200).

    Raises
    ------
    SystemExit
        With status 1 when any check fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scene', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--timeout', type=float, default=300)
    parser.add_argument('--width', type=int, default=1200)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    scene = args.scene.resolve(strict=True)
    overlay_path=scene.parent/(scene.stem.removesuffix('_residue_context')+'_residue_evidence.json')
    overlay=json.loads(overlay_path.read_text()) if scene.stem.endswith('_residue_context') and overlay_path.exists() else None
    with tempfile.TemporaryDirectory(prefix='crevice-pymol-') as tmp:
        tmp = Path(tmp)
        image = tmp/'native.png'
        report = tmp/'scene.json'
        native_geometry = tmp/'native_geometry.npz'
        program = f'''
import json
import numpy as np
from pymol import cmd
cmd.set("max_threads", 4)
cmd.rebuild()
objects = cmd.get_names('objects')
result = {{'version': cmd.get_version()[0], 'atom_count': cmd.count_atoms('crevice_protein'),
          'objects': objects, 'object_types': {{name:cmd.get_type(name) for name in objects}},
          'helper_grid_objects': [name for name in objects if name == 'crevice_grid' or name.startswith('crevice_region_grid')],
          'protein_representation_atoms': {{rep:cmd.count_atoms('crevice_protein and rep '+rep) for rep in ['cartoon','sticks','lines','spheres','surface'] }},
          'cartoon_transparency': cmd.get_setting_float('cartoon_transparency','crevice_protein'),
          'camera_rotation': np.asarray(cmd.get_view()[:9]).reshape(3,3,order='F').tolist(),
          'volume_extent': cmd.get_extent('crevice_volume'),
          'volume_matrix': cmd.get_object_ttt('crevice_volume'),
          'protein_extent': cmd.get_extent('crevice_protein')}}
if {bool(overlay)!r}:
    colors = []
    cmd.iterate('crevice_protein', 'colors.append((ID,color))', space={{'colors':colors}})
    color_map = dict(colors)
    membership = {overlay.get('display_atom_serials',{}) if overlay else {}!r}
    palette = {{'lining':cmd.get_color_index('crevice_boundary_evidence'), 'partners':cmd.get_color_index('crevice_partner_evidence')}}
    result['residue_highlights'] = {{role:{{'expected_atoms':len(indices),'matched_atoms':sum(color_map.get(i)==palette[role] for i in indices),'color_index':palette[role]}}
                                    for role,indices in membership.items()}}
    sticks = []
    cmd.iterate('crevice_protein and rep sticks', 'sticks.append(ID)', space={{'sticks':sticks}})
    expected_sticks = sorted(set(i for indices in membership.values() for i in indices))
    result['residue_sticks_passed'] = sorted(sticks)==expected_sticks if {bool(overlay and overlay.get('overlay_representation')=='cartoon_with_evidence_sticks')!r} else len(sticks)==0
    result['residue_highlights_passed'] = result['residue_sticks_passed'] and bool(membership) and all(palette[role]>=0 and all(color_map.get(i)==palette[role] for i in indices) for role,indices in membership.items())
if 'crevice_grid' in objects:
    field = np.asarray(cmd.get_volume_field('crevice_grid', copy=1))
    result.update(grid_shape=list(field.shape), grid_occupied_samples=int((field>.5).sum()),
                  grid_values=np.unique(field)[:16].tolist())
# OBJ exports the actual displayed triangles. An identity camera makes this
# check independent of the presentation camera and map/object metadata.
view = cmd.get_view()
enabled = cmd.get_names('objects', enabled_only=1)
cmd.disable('all')
for name in enabled:
    if name.startswith('crevice_volume'):
        cmd.enable(name)
cmd.set_view([1,0,0,0,1,0,0,0,1,0,0,-100,0,0,0,1,200,1])
_, obj = cmd.get_mtl_obj()
vertices = np.asarray([[float(v) for v in line.split()[1:4]] for line in obj.splitlines() if line.startswith('v ')])
np.savez_compressed({str(native_geometry)!r}, vertices=vertices,
                    protein_coords=cmd.get_coords('crevice_protein', state=1))
result['triangle_vertices'] = len(vertices)
cmd.set_view(view)
cmd.disable('all')
for name in enabled:
    cmd.enable(name)
open({str(report)!r}, 'w').write(json.dumps(result, indent=2))
cmd.png({str(image)!r}, width={args.width}, height={round(args.width*.75)}, dpi=200, ray=1)
# ray=1 writes synchronously; sync() here can wait on this PML block's own queue.
assert __import__('os').path.getsize({str(image)!r}) > 0
cmd.quit()
'''
        check_script = tmp/'inspect_scene.py'
        check_script.write_text(program)
        # Startup .py arguments can run while a PML Python block is still using
        # the API. Execute inspection at the end of the same block instead.
        source=scene.read_text()
        if source.count('python end')!=1:
            raise ValueError('Native checker requires one CREVICE Python scene block')
        source=source.replace('from pymol import cmd',
                              'from pymol import cmd\ncmd._pymol.__script__ = '+repr(str(scene)),1)
        combined=tmp/'scene_and_inspection.pml'
        combined.write_text(source.replace('python end',program+'\npython end',1))
        environment = dict(os.environ)
        if environment.get('PYTHONPATH'):
            environment['PYTHONPATH'] = os.pathsep.join(str(Path(p or '.').resolve()) for p in environment['PYTHONPATH'].split(os.pathsep))
        try:
            done = subprocess.run([sys.executable, '-m', 'pymol', '-cq', str(combined)], cwd=tmp,
                                  capture_output=True, text=True, timeout=args.timeout, env=environment)
        except subprocess.TimeoutExpired as error:
            def decoded(value):
                return value.decode(errors='replace') if isinstance(value,bytes) else (value or '')
            (root/'pymol.log').write_text(decoded(error.stdout)+decoded(error.stderr))
            result=json.loads(report.read_text()) if report.is_file() else {}
            result.update(status='timeout',timeout_seconds=args.timeout,scene=str(scene))
            (root/'pymol_verification.json').write_text(json.dumps(result,indent=2)+'\n')
            raise SystemExit('Native PyMOL render exceeded its timeout; log and inspection retained')
        log = done.stdout+done.stderr
        (root/'pymol.log').write_text(log)
        result = json.loads(report.read_text()) if report.exists() else {}
        ok = done.returncode == 0 and image.exists() and not any(
            marker in log for marker in ('Traceback', ' Error:', 'SyntaxError:', 'FileNotFoundError:', 'cmd.sync() timed out', 'ObjectSurfaceUpdate-Error'))
        if overlay is not None:
            ok=ok and result.get('residue_highlights_passed',False)
        if ok:
            import itertools
            import numpy as np
            metadata_files = list(scene.parent.glob('*_volume_metadata.json'))
            if len(metadata_files) == 1:
                metadata = json.loads(metadata_files[0].read_text())
                scene_geometry = json.loads((scene.parent/Path(metadata['files']['scene_json']).name).read_text())
                camera_expected = np.asarray(scene_geometry['camera_rotation'])
                if scene.stem.endswith('_end_view'):
                    angle = np.deg2rad(scene_geometry.get('end_view_turn_degrees',74))
                    camera_expected = np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],[0,np.sin(angle),np.cos(angle)]])@camera_expected
                result['camera_rotation_max_error'] = float(np.abs(camera_expected-np.asarray(result['camera_rotation'])).max())
                expected = np.asarray(metadata['pymol_local_to_input_matrix']).reshape(4,4)
                actual = np.asarray(result['volume_matrix']).reshape(4,4)
                result['alignment_max_error_A'] = float(np.max(np.abs(actual-expected)))
                local_path = scene.parent / Path(metadata['files']['pymol_local_dx']).name
                lines = local_path.read_text().splitlines()
                shape = tuple(int(v) for v in lines[1].split()[-3:])
                origin = np.array([float(v) for v in lines[2].split()[1:]])
                deltas = np.array([[float(v) for v in line.split()[1:]] for line in lines[3:6]])
                field = np.fromstring(local_path.read_text().split('data follows\n',1)[1].split('attribute',1)[0], sep=' ')
                # A supersampled display field contains every measured sample at
                # indices divisible by its factor; count only those samples.
                display_smoothing = metadata.get('display_smoothing') or {}
                factor = int(display_smoothing.get('supersample_factor') or 1)
                measured_nodes = field.reshape(shape)[::factor, ::factor, ::factor] if len(field) == int(np.prod(shape)) else field
                source_count = int((measured_nodes>.5).sum())
                result['display_supersample_factor'] = factor
                result['source_grid_occupied_samples'] = source_count
                result['source_grid_shape'] = list(shape)
                standalone = metadata.get('pymol_representation') == 'standalone_triangle_mesh'
                if standalone:
                    mesh_path = scene.parent/Path(metadata['files']['pymol_mesh_npz']).name
                    with np.load(mesh_path,allow_pickle=False) as mesh:
                        corners = mesh['volume_vertices']
                    result['grid_data_origin'] = 'source DX; no live PyMOL map'
                    presentation_ok = (not result['helper_grid_objects']
                        and result['protein_representation_atoms']['cartoon'] > 0
                        and all(result['protein_representation_atoms'][rep] == 0 for rep in (['lines','spheres','surface'] if overlay and overlay.get('overlay_representation')=='cartoon_with_evidence_sticks' else ['sticks','lines','spheres','surface'])))
                else:
                    corners = np.array([origin+np.array(idx)@deltas for idx in itertools.product(*[(0,n-1) for n in shape])])
                    result['grid_data_origin'] = 'source DX and live PyMOL map'
                    presentation_ok = result.get('grid_occupied_samples') == source_count
                # PyMOL reports CGO extents in object-local coordinates even
                # when a TTT is active. Native OBJ triangles below independently
                # verify the transformed surface against input-frame geometry.
                extent_points = corners if standalone else corners@expected[:3,:3].T+expected[:3,3]
                result['extent_coordinate_frame'] = 'object local' if standalone else 'input world'
                bounds = np.array([extent_points.min(axis=0),extent_points.max(axis=0)])
                result['extent_max_error_A'] = float(np.max(np.abs(bounds-np.array(result['volume_extent']))))
                result['presentation_passed'] = presentation_ok
                ok = (result['camera_rotation_max_error'] < 1e-5 and result['alignment_max_error_A'] < .001 and result['extent_max_error_A'] < .001
                      and source_count == metadata['display_points'] and len(field)==int(np.prod(shape)) and presentation_ok)
        if ok:
            from scipy.spatial import cKDTree
            reference_files = list(scene.parent.glob('*_geometry_reference.npz'))
            ok = len(reference_files) == 1 and native_geometry.exists()
            if ok:
                reference = np.load(reference_files[0])
                native = np.load(native_geometry)
                vertices = native['vertices']
                spacing = metadata['spacing_A']
                surface_error = float(cKDTree(reference['points']).query(vertices)[0].max()) if len(vertices) else float('inf')
                coverage_error = float(cKDTree(vertices).query(reference['boundary_points'])[0].max()) if len(vertices) else float('inf')
                protein_error = 0.0
                if len(reference['protein_coords']):
                    protein_error = (max(float(cKDTree(reference['protein_coords']).query(native['protein_coords'])[0].max()),
                                         float(cKDTree(native['protein_coords']).query(reference['protein_coords'])[0].max()))
                                     if reference['protein_coords'].shape == native['protein_coords'].shape else float('inf'))
                result.update(native_surface_sample_error_A=surface_error,
                              native_surface_coverage_error_A=coverage_error,
                              native_protein_coordinate_error_A=protein_error)
                # Occupancy-preserving interpolation can move a crossing along
                # its full occupied-to-empty lattice edge (length h).
                # The constrained supersampled surface stays inside measured
                # cells containing the unsmoothed boundary: at most one cell
                # diagonal (sqrt(3) h) from an occupied sample.
                if display_smoothing.get('method') == 'constrained_gaussian_twicing':
                    bound = spacing*3**.5+.001
                else:
                    bound = spacing*(1.0 if metadata.get('smoothing') else .9)+.001
                ok = surface_error <= bound and coverage_error <= spacing*1.8+.001 and protein_error <= .001
                result['native_surface_sample_bound_A'] = bound
                if 'surface_vertices' in reference:
                    expected_vertices = reference['surface_vertices']
                    mesh_error = max(float(cKDTree(expected_vertices).query(vertices)[0].max()),
                                     float(cKDTree(vertices).query(expected_vertices)[0].max()))
                    result['native_expected_surface_max_error_A'] = mesh_error
                    ok = ok and mesh_error <= .002
                import shutil
                shutil.copy2(native_geometry, root/'native_geometry.npz')
            else:
                result['native_geometry_error'] = 'Missing independent geometry reference'
        if ok:
            from PIL import Image
            import numpy as np
            with Image.open(image) as im:
                rgb = np.asarray(im.convert('RGB'))
                teal = (rgb[:,:,1] > rgb[:,:,0]*1.3) & (rgb[:,:,2] > rgb[:,:,0]*1.3) & (rgb[:,:,0] < 180)
                # Non-channel casts are violet (crevice.presentation.NON_CHANNEL_CAST_RGB).
                violet = (rgb[:,:,2] > rgb[:,:,1]*1.3) & (rgb[:,:,0] > rgb[:,:,1]*1.2) & (rgb[:,:,1] < 180)
                teal = teal | violet
                result.update(width=im.width, height=im.height, teal_fraction=float(teal.mean()))
                ok = int(teal.sum()) >= 25 and result.get('source_grid_occupied_samples',0) > 0
                im.save(root/'pymol_native.png')
        result.update(status='passed' if ok else 'failed', scene=str(scene), returncode=done.returncode)
        (root/'pymol_verification.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2))
        if not ok:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
