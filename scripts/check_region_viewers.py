"""Run native region review scenes and independently inspect their visible data.

Loads a region-review scene from ``crevice region-prepare``
(``PREFIX_review.pml``, ``.vmd`` or ``.cxc``) twice in the real viewer and
checks: protein coordinates (0.001 Å), landmark stick residues and colour,
transferred coordinate landmarks, the world-space region surface against
``PREFIX_geometry_geometry_reference.npz``, the camera, that reloading does
not duplicate the protein, the absence of helper grids, and the rendered
image. ChimeraX runs under ``xvfb-run``.

Usage::

    python scripts/check_region_viewers.py REVIEW/PREFIX_review.pml --viewer pymol --output check/

The output directory must not exist. Writes ``native.png``,
``inspection.json``, ``geometry.npz`` and ``verification.json``; exits with
status 1 unless every check passes.
"""
from pathlib import Path
import argparse,json,os,subprocess,sys
import numpy as np
from crevice.volume_export import _tcl_word
from workspace_paths import local_chimerax


def main():
    """Parse the command line, run the viewer and write ``verification.json``.

    Options: ``scene``, ``--viewer {pymol,vmd,chimerax}``, ``--output`` (new
    directory) and ``--width`` (pixels, height = 0.75 x width; default 1200).

    Raises
    ------
    SystemExit
        With status 1 when any check fails.
    """
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scene',type=Path);p.add_argument('--viewer',required=True,choices=['pymol','vmd','chimerax'])
    p.add_argument('--output',required=True,type=Path);p.add_argument('--width',type=int,default=1200)
    args=p.parse_args();scene=args.scene.resolve(strict=True);root=scene.parent;prefix=scene.stem.removesuffix('_review')
    spec=json.loads((root/(prefix+'_review_scene.json')).read_text());out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    width=args.width;height=round(width*.75)
    if args.viewer=='pymol':
        code=r'''import json,numpy as np
from pathlib import Path
from pymol import cmd
out=Path(__OUT__)
cmd.set('max_threads',4)
model=cmd.get_model('crevice_protein');serials=np.asarray([a.id for a in model.atom])
colors=[];cmd.iterate('crevice_protein','colors.append((ID,color))',space={'colors':colors})
result={'sticks':sorted(a.id for a in cmd.get_model('crevice_protein and rep sticks').atom),
 'colors':{str(i):list(cmd.get_color_tuple(c)) for i,c in colors},
 'camera_rotation':np.asarray(cmd.get_view()[:9]).reshape(3,3,order='F').tolist(),
 'helper_grids':[n for n in cmd.get_names() if n=='crevice_grid' or n.startswith('crevice_region_grid')],
 'protein_objects':sum(n=='crevice_protein' for n in cmd.get_names('objects'))}
arrays={'serials':serials,'protein_coords':np.array([a.coord for a in model.atom]),
 'markers':cmd.get_coords('crevice_review_points') if 'crevice_review_points' in cmd.get_names('objects') else np.empty((0,3))}
cmd.png(str(out/'native.png'),width=__WIDTH__,height=__HEIGHT__,dpi=240,ray=1)
view=cmd.get_view();cmd.disable('all');cmd.enable('crevice_volume');cmd.set_view([1,0,0,0,1,0,0,0,1,0,0,-100,0,0,0,1,200,1])
_,obj=cmd.get_mtl_obj();arrays['vertices']=np.asarray([[float(v) for v in l.split()[1:4]] for l in obj.splitlines() if l.startswith('v ')])
np.savez_compressed(out/'geometry.npz',**arrays);(out/'inspection.json').write_text(json.dumps(result,indent=2,allow_nan=False))
cmd.quit()
'''
        code=code.replace('__OUT__',repr(str(out))).replace('__WIDTH__',str(width)).replace('__HEIGHT__',str(height))
        script=out/'inspect.pml';script.write_text('@'+str(scene)+'\n@'+str(scene)+'\npython\n'+code+'\npython end\n')
        command=[sys.executable,'-m','pymol','-cq',str(script)]
    elif args.viewer=='chimerax':
        code=r'''import json,numpy as np
from pathlib import Path
from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure,Atom
out=Path(__OUT__)
run(session,'open '+__SCENE__);run(session,'open '+__SCENE__)
proteins=[m for m in session.models.list(type=AtomicStructure) if m.name=='crevice_protein'];protein=proteins[0]
serials=protein.atoms.serial_numbers
result={'sticks':sorted(map(int,serials[protein.atoms.displays & (protein.atoms.draw_modes==Atom.STICK_STYLE)])),
 'colors':{str(i):list(map(float,c[:3]/255)) for i,c in zip(serials,protein.atoms.colors)},
 'camera_rotation':session.main_view.camera.position.matrix[:,:3].T.tolist(),
 'protein_objects':len(proteins),'helper_grids':[]}
surface=next(m for m in session.models.list() if m.name=='crevice_volume')
markers=next(m for m in session.models.list() if m.name=='crevice_review_points')
np.savez_compressed(out/'geometry.npz',serials=serials,protein_coords=protein.atoms.coords,
 vertices=surface.scene_position*surface.vertices,markers=markers.atoms.coords)
run(session,'save '+json.dumps(str(out/'native.png'))+' width __WIDTH__ height __HEIGHT__ supersample 2')
(out/'inspection.json').write_text(json.dumps(result,indent=2,allow_nan=False))
'''.replace('__OUT__',repr(str(out))).replace('__SCENE__',repr(json.dumps(str(scene)))).replace('__WIDTH__',str(width)).replace('__HEIGHT__',str(height))
        script=out/'inspect.py';script.write_text(code)
        command=['xvfb-run','-a','-s','-screen 0 1280x1024x24',str(local_chimerax().resolve()),'--silent','--exit',str(script)]
    else:
        q=_tcl_word
        code=['if {[catch {','source '+q(str(scene)),'source '+q(str(scene)),'set out '+q(str(out)),r'''
set f [open [file join $out inspection.tsv] w]
set a [atomselect $crevice_protein all]
foreach serial [$a get serial] xyz [$a get {x y z}] {puts $f "atom|$serial|$xyz"}
$a delete
set sticks {}
for {set i 0} {$i<[molinfo $crevice_protein get numreps]} {incr i} {
    set style [lindex [molinfo $crevice_protein get [list [list rep $i]]] 0]
    if {[string match "Licorice*" $style]} {
        set selection [lindex [molinfo $crevice_protein get [list [list selection $i]]] 0]
        set a [atomselect $crevice_protein $selection];set sticks [concat $sticks [$a get serial]];$a delete
    }
}
puts $f "sticks|[lsort -integer $sticks]"
puts $f "color|[colorinfo rgb 17]"
puts $f "camera|[lindex [molinfo $crevice_protein get rotate_matrix] 0]"
puts $f "molecules|[llength [molinfo list]]|[llength $crevice_review_owned]"
foreach id [graphics $crevice_review_graphics list] {
    set item [graphics $crevice_review_graphics info $id]
    if {[lindex $item 0] eq "sphere"} {puts $f "marker|[lindex $item 1]"}
}
close $f
''',f'display resize {width} {height}','catch {display redraw}','display update','render TachyonInternal [file join $out native.tga]',
'foreach m [molinfo list] {mol off $m}','mol on $crevice_volume',
'foreach matrix {center_matrix rotate_matrix scale_matrix global_matrix} {molinfo $crevice_volume set $matrix [list [transidentity]]}',
'render Wavefront [file join $out surface.obj]','} err]} {puts stderr "CREVICE_CHECK_ERROR: $err";puts stderr $::errorInfo;exit 1}','quit']
        script=out/'inspect.tcl';script.write_text('\n'.join(code));command=['vmd','-dispdev','text','-e',str(script)]
    env=dict(os.environ);env['VMDFORCECPUCOUNT']='2'
    if env.get('PYTHONPATH'):env['PYTHONPATH']=os.pathsep.join(str(Path(x or '.').resolve()) for x in env['PYTHONPATH'].split(os.pathsep))
    with (out/'native.log').open('w') as log:
        try:
            with subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT,env=env) as process:
                try:returncode=process.wait(timeout=600)
                except subprocess.TimeoutExpired:process.kill();process.wait();returncode=-1
        except OSError as error:log.write(str(error));returncode=-1
    if args.viewer=='vmd' and (out/'inspection.tsv').exists():
        serials=[];coords=[];markers=[];r={};sticks=[]
        for line in (out/'inspection.tsv').read_text().splitlines():
            v=line.split('|')
            if v[0]=='atom':serials.append(int(v[1]));coords.append(list(map(float,v[2].split())))
            elif v[0]=='sticks':sticks=list(map(int,v[1].split()));r['sticks']=sticks
            elif v[0]=='color':r['colors']={str(i):list(map(float,v[1].split())) for i in sticks}
            elif v[0]=='camera':r['camera_rotation']=np.array(list(map(float,v[1].replace('{','').replace('}','').split()))).reshape(4,4)[:3,:3].tolist()
            elif v[0]=='molecules':r['protein_objects']=1 if v[1]==v[2] else 0
            elif v[0]=='marker':markers.append(list(map(float,v[1].strip('{}').split())))
        obj=out/'surface.obj';vertices=np.asarray([[float(x) for x in l.split()[1:4]] for l in obj.read_text().splitlines() if l.startswith('v ')]) if obj.exists() else np.empty((0,3))
        np.savez_compressed(out/'geometry.npz',serials=serials,protein_coords=coords,vertices=vertices,markers=np.asarray(markers).reshape(-1,3))
        (out/'inspection.json').write_text(json.dumps(r,indent=2))
        if (out/'native.tga').exists():
            from PIL import Image
            Image.open(out/'native.tga').save(out/'native.png')
    checks={};r={};log=(out/'native.log').read_text()
    checks['native_execution']=returncode==0 and (out/'inspection.json').exists() and not any(s in log for s in ['Traceback','CREVICE_CHECK_ERROR','SyntaxError',' Error:'])
    if checks['native_execution']:
        from scipy.spatial import cKDTree
        from PIL import Image
        r=json.loads((out/'inspection.json').read_text());a=np.load(out/'geometry.npz')
        expected=np.asarray([[float(l[k:k+8]) for k in [30,38,46]] for l in (root/spec['protein_pdb']).read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))])
        checks['protein_coordinates']=bool(np.allclose(a['protein_coords'][np.argsort(a['serials'])],expected,rtol=0,atol=.001))
        ids=sorted(i for row in spec['landmarks'] for i in row['residue_atom_serials'])
        checks['landmark_sticks']=r['sticks']==ids
        checks['landmark_colors']=all(np.max(np.abs(np.asarray(r['colors'][str(i)])-spec['landmark_rgb']))<.004 for i in ids)
        coords=np.array([p for row in spec['coordinate_landmarks'] for p in row['points_A']]).reshape(-1,3)
        checks['transferred_coordinates']=bool(a['markers'].shape==coords.shape and np.allclose(a['markers'],coords,atol=.001,rtol=0))
        ref=np.load(root/(prefix+'_geometry_geometry_reference.npz'))
        verts=a['vertices'];wanted=ref['surface_vertices']
        error=max(cKDTree(verts).query(wanted)[0].max(),cKDTree(wanted).query(verts)[0].max()) if len(verts) else float('inf')
        meta=json.loads((root/(prefix+'_geometry_volume_metadata.json')).read_text())
        r['surface_error_A']=float(error) if np.isfinite(error) else None;checks['world_surface']=error<(.005 if args.viewer!='vmd' else 1.8*meta['spacing_A']+.005)
        checks['camera']=bool(np.allclose(r['camera_rotation'],spec['camera']['camera_rotation'],atol=1e-6))
        checks['reload_no_duplicate_protein']=r['protein_objects']==1
        checks['no_helper_grids']=not r.get('helper_grids',[])
        rgb=np.asarray(Image.open(out/'native.png').convert('RGB'));border=np.concatenate([rgb[0],rgb[-1],rgb[:,0],rgb[:,-1]])
        r['foreground_fraction']=float(np.any(rgb<245,axis=2).mean());r['border_foreground_fraction']=float(np.any(border<245,axis=1).mean())
        checks['rendered_image']=r['foreground_fraction']>.001 and r['border_foreground_fraction']<.01
    checks={k:bool(v) for k,v in checks.items()}
    result={k:v for k,v in r.items() if k not in {'colors','sticks'}}
    result.update(status='passed' if all(checks.values()) else 'failed',checks=checks,command=command,returncode=returncode)
    (out/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2,allow_nan=False))
    if result['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
