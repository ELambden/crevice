"""Native validation of coordinated hydration scenes, identities and controls.

Opens a unified analysis scene (``PREFIX_analysis.pml``, ``.vmd`` or
``.cxc`` written by ``crevice hydration``, ``publish``, ``cavity-trajectory``
or ``region-trajectory``) in the real viewer, renders it and checks what the
viewer actually shows against the scene JSON: residue colours, the stick
residues of each view mode (``crevice_view``/``creviceView``), snapshot
coordinates (0.001 Å), surface vertices against the reference meshes, the
camera, the absence of helper grids, and that drawn text matches the
annotation setting. Requires the viewer to be installed (ChimeraX is located
with ``workspace_paths.local_chimerax``).

Usage::

    python scripts/check_analysis_viewers.py OUT/PREFIX_analysis.pml --viewer pymol --output check/

Writes ``native.png`` (and ``network.png`` with ``--render-network``),
``native_inspection.json``, ``native_geometry.npz``, ``native.log`` and
``verification.json`` to ``--output``; exits with status 1 unless every check
passes.
"""
from pathlib import Path
import argparse,json,os,subprocess,sys
import numpy as np
from crevice.volume_export import _tcl_word
from workspace_paths import local_chimerax


def main():
    """Parse the command line, run the viewer and write ``verification.json``.

    Options: ``scene`` (scene script path), ``--viewer {pymol,vmd,chimerax}``,
    ``--output`` (directory, created), ``--width`` (image width in pixels,
    height = 0.75 x width; default 1200), ``--timeout`` (seconds, default 900),
    ``--render-network`` (also render the network view) and ``--reload`` (load
    the scene twice to check that reloading does not duplicate objects).

    Raises
    ------
    SystemExit
        With status 1 when any check fails.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scene',type=Path);parser.add_argument('--viewer',choices=['pymol','vmd','chimerax'],required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--width',type=int,default=1200)
    parser.add_argument('--timeout',type=int,default=900)
    parser.add_argument('--render-network',action='store_true')
    parser.add_argument('--reload',action='store_true')
    args=parser.parse_args();scene_path=args.scene.resolve(strict=True);out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    prefix=scene_path.stem.removesuffix('_network').removesuffix('_analysis')
    scene=json.loads((scene_path.parent/(prefix+'_analysis_scene.json')).read_text())
    w,h=args.width,round(args.width*.75);log=out/'native.log'
    if args.viewer=='pymol':
        program='''
import json,numpy as np
from pathlib import Path
out=Path(__OUT__)
cmd.set('max_threads',4)
result={'viewer':'pymol','version':cmd.get_version()[0],'modes':{},'colors':{},'snapshots':{}}
model=cmd.get_model('crevice_protein');actual={a.id:a for a in model.atom}
atom_colors=[]
cmd.iterate('crevice_protein','atom_colors.append((ID,color))',space={'atom_colors':atom_colors})
for i,color in atom_colors:result['colors'][str(i)]=list(cmd.get_color_tuple(color))
arrays={'serials':np.array(list(actual)),'protein_coords':np.array([a.coord for a in actual.values()])}
for mode in ['hydration','network','partners','geometry']:
    crevice_view(mode);result['modes'][mode]=sorted(a.id for a in cmd.get_model('crevice_protein and rep sticks').atom)
crevice_view('hydration')
cmd.png(str(out/'native.png'),width=__WIDTH__,height=__HEIGHT__,dpi=240,ray=1)
crevice_view('network')
if __NETWORK__:cmd.png(str(out/'network.png'),width=__WIDTH__,height=__HEIGHT__,dpi=240,ray=1)
for i,snapshot in enumerate(crevice_scene['snapshots'],1):
    crevice_frame(i);m=cmd.get_model('crevice_protein');arrays['snapshot_'+str(i)]=np.asarray([a.coord for a in sorted(m.atom,key=lambda a:a.id)])
    result['snapshots'][str(i)]={'reference_edges_visible':'crevice_interactions' in cmd.get_names('objects',enabled_only=1)}
crevice_frame(1);crevice_view('hydration');view=cmd.get_view();enabled=cmd.get_names('objects',enabled_only=1)
result['text_objects']={'label_atoms':cmd.count_atoms('rep labels'),
    'legend_objects':[n for n in cmd.get_names('all') if n in ('crevice_hydration_legend','crevice_hydration_labels','crevice_region_status')]}
for surface in crevice_scene['surfaces']:
    cmd.disable('all');cmd.enable(surface['name']);cmd.set_view([1,0,0,0,1,0,0,0,1,0,0,-100,0,0,0,1,200,1])
    _,obj=cmd.get_mtl_obj();arrays['surface_'+surface['key']]=np.array([[float(v) for v in line.split()[1:4]] for line in obj.splitlines() if line.startswith('v ')])
cmd.set_view(view)
result['camera_rotation']=np.asarray(view[:9]).reshape(3,3,order='F').tolist()
result['helper_grids']=[x for x in cmd.get_names() if x=='crevice_grid' or x.startswith('crevice_region_grid')]
np.savez_compressed(out/'native_geometry.npz',**arrays)
(out/'native_inspection.json').write_text(json.dumps(result,indent=2))
cmd.quit()
'''.replace('__RELOAD__',str(args.reload)).replace('__OUT__',repr(str(out))).replace('__NETWORK__',str(args.render_network)).replace('__WIDTH__',str(w)).replace('__HEIGHT__',str(h))
        source=scene_path.read_text().replace('from pymol import cmd','from pymol import cmd\ncmd._pymol.__script__ = '+repr(str(scene_path)),1)
        script=out/'inspect.pml'
        if args.reload:
            script.write_text('@'+str(scene_path)+'\n@'+str(scene_path)+'\npython\n'+program+'\npython end\n')
        else:script.write_text(source.replace('python end',program+'\npython end',1))
        env=dict(os.environ)
        if env.get('PYTHONPATH'):env['PYTHONPATH']=os.pathsep.join(str(Path(p or '.').resolve()) for p in env['PYTHONPATH'].split(os.pathsep))
        command=[sys.executable,'-m','pymol','-cq',str(script)]
    elif args.viewer=='chimerax':
        program='''
from pathlib import Path
import json,numpy as np,traceback
from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure,Atom
from chimerax.map import Volume
out=Path(__OUT__)
try:
    run(session,'open '+__SCENE__)
    if __RELOAD__:run(session,'open '+__SCENE__)
    protein=next(m for m in session.models.list(type=AtomicStructure) if m.name=='crevice_protein')
    serials=np.asarray(protein.atoms.serial_numbers);arrays={'serials':serials,'protein_coords':protein.atoms.coords.copy()}
    result={'viewer':'chimerax','colors':{str(i):list(map(float,c[:3]/255)) for i,c in zip(serials,protein.atoms.colors)},'modes':{},'snapshots':{}}
    for mode in ['hydration','network','partners','geometry']:
        run(session,'creviceView '+mode);result['modes'][mode]=sorted(map(int,serials[protein.atoms.displays & (protein.atoms.draw_modes==Atom.STICK_STYLE)]))
    run(session,'creviceView hydration')
    run(session,'save '+json.dumps(str(out/'native.png'))+' width __WIDTH__ height __HEIGHT__ supersample 2')
    run(session,'creviceView network')
    if __NETWORK__:run(session,'save '+json.dumps(str(out/'network.png'))+' width __WIDTH__ height __HEIGHT__ supersample 2')
    for i in range(1,__SNAPSHOTS__+1):
        run(session,'creviceFrame '+str(i));arrays['snapshot_'+str(i)]=protein.atoms.coords[np.argsort(serials)].copy()
        result['snapshots'][str(i)]={'reference_edges_visible':any(m.name=='crevice_interactions' and m.display for m in session.models.list())}
    run(session,'creviceFrame 1; creviceView hydration')
    from chimerax.label.label2d import session_labels
    from chimerax.label.label3d import ObjectLabels
    from chimerax.color_key.model import get_model as key_model
    labels=session_labels(session,create=False)
    result['text_objects']={'labels_2d':[l.name for l in labels.all_labels] if labels else [],
        'label_3d_models':sum(isinstance(m,ObjectLabels) for m in session.models.list()),
        'colour_key':key_model(session,create=False) is not None}
    for volume in session.models.list(type=Volume):
        key='water' if volume.name=='crevice_water_density' else 'cavity'
        arrays['surface_'+key]=np.concatenate([s.scene_position*s.vertices for s in volume.surfaces if s.vertices is not None])
        arrays['field_'+key]=volume.data.full_matrix()
    result['camera_rotation']=session.main_view.camera.position.matrix[:,:3].T.tolist()
    result['helper_grids']=[m.name for m in session.models.list() if m.name=='crevice_grid' or m.name.startswith('crevice_region_grid')]
    np.savez_compressed(out/'native_geometry.npz',**arrays)
    (out/'native_inspection.json').write_text(json.dumps(result,indent=2))
except Exception as error:
    (out/'native_failure.json').write_text(json.dumps({'error':str(error),'traceback':traceback.format_exc()},indent=2))
    raise
'''.replace('__RELOAD__',str(args.reload)).replace('__OUT__',repr(str(out))).replace('__SCENE__',repr(json.dumps(str(scene_path)))).replace('__NETWORK__',str(args.render_network)).replace('__WIDTH__',str(w)).replace('__HEIGHT__',str(h)).replace('__SNAPSHOTS__',str(len(scene['snapshots'])))
        script=out/'inspect.py';script.write_text(program);env=None
        command=['xvfb-run','-a','-s','-screen 0 1280x1024x24',str(local_chimerax().resolve()),'--silent','--exit',str(script)]
    else:
        q=_tcl_word
        lines=['if {[catch {','source '+q(str(scene_path)),('source '+q(str(scene_path))) if args.reload else '# single load', 'set out '+q(str(out)), r'''
set report [open [file join $out native.tsv] w]
set a [atomselect $crevice_protein all]
foreach serial [$a get serial] xyz [$a get {x y z}] value [$a get user] {
    if {$value<0} {set color [colorinfo rgb 3]} else {set color [colorinfo rgb [expr {$crevice_palette_start+int($value*($crevice_palette_size-1))}]]}
    puts $report "atom|$serial|$xyz|$color"
}
$a delete
foreach mode {hydration network partners geometry} {
    crevice_view $mode
    set ids {}
    foreach entry $crevice_stick_reps {
        lassign $entry rep condition
        set text [lindex [molinfo $crevice_protein get [list [list selection $rep]]] 0]
        set a [atomselect $crevice_protein $text];set ids [concat $ids [$a get serial]];$a delete
    }
    puts $report "mode|$mode|[lsort -integer $ids]"
}
for {set i 1} {$i<=[llength $crevice_snapshot_files]} {incr i} {
    crevice_view network;crevice_frame $i
    set a [atomselect $crevice_protein all]
    foreach serial [$a get serial] xyz [$a get {x y z}] {puts $report "snapshot|$i|$serial|$xyz"}
    $a delete
    puts $report "edges|$i|[molinfo $crevice_contacts get drawn]"
}
crevice_frame 1;crevice_view hydration
set crevice_text_items 0;set crevice_legend_items 0
foreach m [molinfo list] {foreach g [graphics $m list] {
    if {[string match "text*" [graphics $m info $g]]} {incr crevice_text_items}
    if {$m==$crevice_legend} {incr crevice_legend_items}
}}
puts $report "text|$crevice_text_items|$crevice_legend_items"
puts $report "camera|[lindex [molinfo $crevice_protein get rotate_matrix] 0]"
close $report
''',f'display resize {w} {h}','catch {display redraw}','display update','render TachyonInternal [file join $out native.tga]',
               'crevice_view network','catch {display redraw}','display update','render TachyonInternal [file join $out network.tga]' if args.render_network else '# network state checked without extra image',
               'foreach molecule [molinfo list] {mol off $molecule}',
               'foreach key [array names crevice_surface] {',
               'set m $crevice_surface($key);mol on $m',
               'foreach matrix {center_matrix rotate_matrix scale_matrix global_matrix} {molinfo $m set $matrix [list [transidentity]]}',
               'render Wavefront [file join $out surface_$key.obj]','mol off $m','}',
               '} err]} {puts stderr "CREVICE_CHECK_ERROR: $err";puts stderr $::errorInfo;exit 1}','quit']
        script=out/'inspect.tcl';script.write_text('\n'.join(lines));env=None;command=['vmd','-dispdev','text','-e',str(script)]
    try:
        with log.open('w') as output:
            with subprocess.Popen(command,stdin=subprocess.PIPE,stdout=output,stderr=subprocess.STDOUT,env=env) as process:
                try:code=process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:process.kill();code=process.wait()
    except OSError as error:log.write_text(str(error));code=-1
    result={};checks={};message=log.read_text()
    if args.viewer=='vmd' and (out/'native.tsv').exists():
        result={'colors':{},'modes':{},'snapshots':{}};serials=[];coords=[];snapshots={}
        for line in (out/'native.tsv').read_text().splitlines():
            row=line.split('|')
            if row[0]=='atom':serials.append(int(row[1]));coords.append(list(map(float,row[2].split())));result['colors'][row[1]]=list(map(float,row[3].split()))
            elif row[0]=='mode':result['modes'][row[1]]=list(map(int,row[2].split()))
            elif row[0]=='snapshot':snapshots.setdefault(row[1],{})[int(row[2])]=list(map(float,row[3].split()))
            elif row[0]=='edges':result['snapshots'][row[1]]={'reference_edges_visible':bool(int(row[2])) and bool(scene['contacts'])}
            elif row[0]=='text':result['text_objects']={'text_graphics':int(row[1]),'legend_graphics':int(row[2])}
            elif row[0]=='camera':result['camera_rotation']=np.asarray(list(map(float,row[1].replace('{','').replace('}','').split()))).reshape(4,4)[:3,:3].tolist()
        arrays={'serials':np.asarray(serials),'protein_coords':np.asarray(coords)}
        for key,value in snapshots.items():arrays['snapshot_'+key]=np.asarray([value[i] for i in sorted(value)])
        for surface in scene['surfaces']:
            p=out/('surface_'+surface['key']+'.obj')
            if p.is_file():arrays['surface_'+surface['key']]=np.asarray([[float(v) for v in line.split()[1:4]] for line in p.read_text().splitlines() if line.startswith('v ')])
        np.savez_compressed(out/'native_geometry.npz',**arrays)
        from PIL import Image
        for name in ['native','network']:
            p=out/(name+'.tga')
            if p.exists():Image.open(p).save(out/(name+'.png'))
        (out/'native_inspection.json').write_text(json.dumps(result,indent=2))
    if (out/'native_inspection.json').exists():result=json.loads((out/'native_inspection.json').read_text())
    checks['execution']=code==0 and bool(result) and not any(x in message for x in ['Traceback','CREVICE_CHECK_ERROR','SyntaxError:',' Error:'])
    if checks['execution']:
        from scipy.spatial import cKDTree
        from PIL import Image
        ref=np.load(scene_path.parent/scene['geometry_npz']);actual=np.load(out/'native_geometry.npz')
        coords=actual['protein_coords'][np.argsort(actual['serials'])]
        checks['protein_coordinates']=bool(coords.shape==ref['protein_coords_A'].shape and np.allclose(coords,ref['protein_coords_A'],atol=.001,rtol=0))
        checks['hydration_colors']=all(np.max(np.abs(np.asarray(result['colors'].get(str(i),[-9]*3))-r['color_rgb']))<.004 for r in scene['residues'] for i in r['atom_serials'])
        for mode in ['hydration','network','partners','geometry']:
            shown=set(scene['display_residues'])
            if mode=='partners':shown.update(scene['partner_residues'])
            if mode=='network':shown.update(k for e in scene['contacts'] for k in [e['source'],e['target']])
            wanted=sorted(i for r in scene['residues'] if r['residue'] in shown for i in r['atom_serials']) if mode!='geometry' else []
            checks['sticks_'+mode]=result['modes'][mode]==wanted
        for i,snapshot in enumerate(scene['snapshots'],1):
            expected=np.asarray([[float(line[k:k+8]) for k in [30,38,46]] for line in (scene_path.parent/snapshot['pdb']).read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))])
            checks['snapshot_'+str(i)]=bool(np.allclose(actual['snapshot_'+str(i)],expected,atol=.001,rtol=0))
            checks['reference_edges_'+str(i)]=result['snapshots'][str(i)]['reference_edges_visible']==bool(i==1 and scene['contacts'])
        for s in scene['surfaces']:
            key=s['key'];a=actual['surface_'+key] if 'surface_'+key in actual else np.empty((0,3));b=ref[key+'_vertices']
            error=float(max(cKDTree(a).query(b)[0].max(),cKDTree(b).query(a)[0].max())) if len(a) else None
            result['surface_'+key+'_error_A']=error;checks['surface_'+key]=error is not None and error < (0.005 if args.viewer=='pymol' else .9)
        checks['camera_rotation']=bool(np.allclose(result['camera_rotation'],scene['camera']['camera_rotation'],atol=1e-6))
        checks['no_helper_grids']=not result.get('helper_grids',[])
        drawn=bool((scene.get('annotations') or {'drawn':True})['drawn'])
        text=result.get('text_objects',{})
        count=sum(len(v) if isinstance(v,list) else int(v) for v in text.values())
        result['text_object_count']=count;result['annotations_drawn']=drawn
        checks['text_matches_annotations']=bool(text) and (count>0)==drawn
        for name in (['native','network'] if args.render_network else ['native']):
            p=out/(name+'.png');checks['image_'+name]=p.exists()
            if p.exists():
                a=np.asarray(Image.open(p).convert('RGB'));border=np.concatenate([a[0],a[-1],a[:,0],a[:,-1]])
                result[name+'_foreground_fraction']=float(np.any(a<245,axis=2).mean());result[name+'_border_foreground_fraction']=float(np.any(border<245,axis=1).mean())
                checks['image_'+name]=result[name+'_foreground_fraction']>.001 and result[name+'_border_foreground_fraction']<.01
    summary={k:v for k,v in result.items() if k not in {'colors','modes'}}
    summary.update(status='passed' if all(checks.values()) else 'failed',checks=checks,returncode=code,command=command,scene=str(scene_path))
    (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if summary['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
