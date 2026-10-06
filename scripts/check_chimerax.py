"""Run ChimeraX scenes and verify native coordinates, residue sticks and surfaces.

Opens a cast or residue-evidence scene (``PREFIX_volume.cxc`` or
``PREFIX_residue_context.cxc``) in ChimeraX (offscreen, or under ``xvfb-run``
with ``--xvfb``), saves a supersampled PNG and checks: protein coordinates
against ``PREFIX_residue_reference.tsv`` (0.001 Å), the stick residues and
their gold/magenta colours against ``PREFIX_residue_evidence.json``, the
displayed surface against ``*_geometry_reference.npz``, the camera against
the scene JSON, and that the image has content and a clean border.

Usage::

    python scripts/check_chimerax.py OUT/PREFIX_volume.cxc --output check/

Writes ``chimerax_native.png``, ``native_inspection.json``,
``native_geometry.npz``, ``chimerax.log`` and ``chimerax_verification.json``;
exits with status 1 unless the checks pass.
"""
from pathlib import Path
import argparse,json,subprocess,tempfile,shutil
from workspace_paths import local_chimerax


def main():
    """Parse the command line, run ChimeraX and write ``chimerax_verification.json``.

    Options: ``scene``, ``--output`` (directory), ``--executable`` (default
    ``chimerax`` on ``PATH``, else ``workspace_paths.local_chimerax()``),
    ``--xvfb`` (virtual X display instead of offscreen rendering), ``--width``
    (pixels, default 900) and ``--timeout`` (seconds, default 180).

    Raises
    ------
    SystemExit
        With status 1 when any check fails.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scene',type=Path);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',default=shutil.which('chimerax'))
    parser.add_argument('--xvfb',action='store_true',help='Use a virtual X display instead of offscreen rendering')
    parser.add_argument('--width',type=int,default=900);parser.add_argument('--timeout',type=float,default=180)
    args=parser.parse_args();scene=args.scene.resolve(strict=True);out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    executable=args.executable or str(local_chimerax())
    prefix=scene.stem.removesuffix('_residue_context');evidence=scene.parent/(prefix+'_residue_evidence.json')
    overlay=json.loads(evidence.read_text()) if evidence.is_file() else None
    inspection=f'''
from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure,Atom
from chimerax.map import Volume
from pathlib import Path
import numpy as np,json,traceback
out=Path({str(out)!r})
try:
    run(session,'open '+{json.dumps(str(scene))!r})
    run(session,'save '+{json.dumps(str(out/'chimerax_native.png'))!r}+' width {args.width} height {round(args.width*.75)} supersample 2')
    atoms=session.models.list(type=AtomicStructure)[0].atoms
    result={{'atom_count':len(atoms),'displayed_atom_count':int(atoms.displays.sum()),'model_ids':[m.id_string for m in session.models.list()], 'surface_vertices':0}}
    serials=np.asarray(atoms.serial_numbers);coords=np.asarray(atoms.coords)
    result['displayed_stick_serials']=sorted(int(i) for i in serials[atoms.displays & (atoms.draw_modes==Atom.STICK_STYLE)])
    result['atom_serials']=serials.tolist();result['atom_colors']=atoms.colors.tolist()
    vertices=[]
    for volume in session.models.list(type=Volume):
        if volume.display:
            for surface in volume.surfaces:
                if surface.display and surface.vertices is not None:vertices.append(surface.scene_position*surface.vertices)
    joined=np.concatenate(vertices) if vertices else np.empty((0,3))
    result['surface_vertices']=len(joined)
    result['camera_matrix']=session.main_view.camera.position.matrix.tolist()
    np.savez_compressed(out/'native_geometry.npz',protein_coords=coords,atom_serials=serials,surface_vertices=joined)
    result['native_execution_passed']=True
    (out/'native_inspection.json').write_text(json.dumps(result,indent=2))
except Exception as error:
    (out/'native_inspection.json').write_text(json.dumps({{'native_execution_passed':False,'error':str(error),'traceback':traceback.format_exc()}},indent=2))
    raise
'''
    with tempfile.TemporaryDirectory(prefix='crevice-chimerax-') as directory:
        script=Path(directory)/'inspect.py';script.write_text(inspection)
        command=[executable,'--silent','--exit',str(script)]
        if args.xvfb:command=['xvfb-run','-a','-s','-screen 0 1280x1024x24',*command]
        else:command.insert(1,'--offscreen')
        try:
            done=subprocess.run(command,capture_output=True,text=True,timeout=args.timeout)
            log=done.stdout+done.stderr;code=done.returncode
        except subprocess.TimeoutExpired as error:
            def text(x):return x.decode(errors='replace') if isinstance(x,bytes) else (x or '')
            log=text(error.stdout)+text(error.stderr);code=-1
    (out/'chimerax.log').write_text(log)
    path=out/'native_inspection.json';result=json.loads(path.read_text()) if path.is_file() else {'native_execution_passed':False}
    ok=code==0 and result['native_execution_passed'] and (out/'chimerax_native.png').is_file()
    if ok:
        import numpy as np
        from scipy.spatial import cKDTree
        native=np.load(out/'native_geometry.npz')
        if overlay is not None:
            expected=sorted(i for group in overlay['display_atom_serials'].values() for i in group)
            result['residue_sticks_passed']=expected==result['displayed_stick_serials'] and result['displayed_atom_count']==len(expected)
            colors=dict(zip(result['atom_serials'],result['atom_colors']))
            result['residue_colors_passed']=all(colors[i]==palette for role,palette in [('lining',[217,140,51,255]),('partners',[163,68,120,255])] for i in overlay['display_atom_serials'][role])
            reference=np.loadtxt(scene.parent/(prefix+'_residue_reference.tsv'))
            actual=dict(zip(native['atom_serials'],native['protein_coords']))
            result['protein_coordinate_error_A']=float(np.max(np.abs(np.array([actual[int(i)] for i in reference[:,0]])-reference[:,1:])))
            ok=ok and result['residue_sticks_passed'] and result['residue_colors_passed'] and result['protein_coordinate_error_A']<=.001
        references=list(scene.parent.glob('*_geometry_reference.npz'))
        if len(references)==1 and len(native['surface_vertices']):
            ref=np.load(references[0]);points=ref['surface_vertices'] if 'surface_vertices' in ref else ref['boundary_points']
            error=max(cKDTree(points).query(native['surface_vertices'])[0].max(),cKDTree(native['surface_vertices']).query(points)[0].max())
            result['surface_reference_error_A']=float(error)
            meta=json.loads(next(scene.parent.glob('*_volume_metadata.json')).read_text())
            # Different contour triangulations need not have identical vertices.
            spacing=float(meta.get('spacing',meta.get('spacing_A',.5)))
            result['surface_reference_bound_A']=spacing*1.8+.001
            ok=ok and error<=result['surface_reference_bound_A']
        else:ok=False;result['surface_reference_error']='Missing surface or unique independent reference'
        # Full analysis bundles also contain an analysis_scene.json; the volume
        # manifest identifies the camera belonging to this cast unambiguously.
        metadata_files=list(scene.parent.glob('*_volume_metadata.json'))
        camera_files=[]
        if len(metadata_files)==1:
            metadata=json.loads(metadata_files[0].read_text())
            camera_name=metadata.get('files',{}).get('scene_json')
            if camera_name:
                camera_path=scene.parent/Path(camera_name).name
                if camera_path.is_file():camera_files=[camera_path]
        if len(camera_files)==1:
            geometry=json.loads(camera_files[0].read_text())
            expected_rotation=np.asarray(geometry['camera_rotation']).T
            result['camera_rotation_error']=float(np.max(np.abs(np.asarray(result['camera_matrix'])[:,:3]-expected_rotation)))
            ok=ok and result['camera_rotation_error']<1e-8
        else:ok=False;result['camera_rotation_error']='Missing unique scene camera reference'
        from PIL import Image
        rgb=np.asarray(Image.open(out/'chimerax_native.png').convert('RGB'));border=np.concatenate([rgb[0],rgb[-1],rgb[:,0],rgb[:,-1]])
        result['foreground_fraction']=float(np.any(rgb<245,axis=2).mean());result['border_foreground_fraction']=float(np.any(border<245,axis=1).mean())
        ok=ok and result['foreground_fraction']>.001 and result['border_foreground_fraction']<.01
    result.update(status='passed' if ok else 'failed',returncode=code,scene=str(scene),executable=executable)
    (out/'chimerax_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in {'atom_serials','atom_colors','displayed_stick_serials'}},indent=2))
    if not ok:raise SystemExit(1)


if __name__=='__main__':main()
