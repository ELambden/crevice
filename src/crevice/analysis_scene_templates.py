"""Native viewer scripts using one scene definition and absolute hydration palette."""
from pathlib import Path
import json
import math

from .presentation import annotate_option, annotations_enabled, scene_annotation_record

PML = r'''# CREVICE coordinated hydration and residue evidence. Keep companion assets together.
python
from pathlib import Path
import inspect,json
import numpy as np
from pymol import cmd
from pymol.cgo import BEGIN,END,TRIANGLES,NORMAL,VERTEX,CYLINDER
_crevice_root=Path(getattr(cmd._pymol,'__script__',inspect.currentframe().f_code.co_filename)).resolve().parent
if not (_crevice_root/__SCENE__).is_file():_crevice_root=Path(__FALLBACK__)
crevice_scene=json.loads((_crevice_root/__SCENE__).read_text())
crevice_geometry=np.load(_crevice_root/crevice_scene['geometry_npz'],allow_pickle=False)
for name in ['crevice_protein','crevice_volume','crevice_cavity','crevice_water_density','crevice_interactions','crevice_hydration_legend','crevice_hydration_labels','crevice_region_status','crevice_grid','crevice_member_waters']:
    cmd.delete(name)
cmd.delete('crevice_region_grid*')
cmd.load(str(_crevice_root/crevice_scene['protein_pdb']),'crevice_protein',state=1,zoom=0)
crevice_model=cmd.get_model('crevice_protein',state=1)
crevice_by_id={a.id:a.coord for a in crevice_model.atom}
crevice_expected=crevice_geometry['protein_coords_A']
if set(crevice_by_id)!=set(range(1,len(crevice_expected)+1)) or not np.allclose([crevice_by_id[i] for i in range(1,len(crevice_expected)+1)],crevice_expected,atol=.001,rtol=0):
    raise ValueError('CREVICE reference atom identities or coordinates differ')
cmd.set_color('crevice_context',[.66,.70,.74])
cmd.set('cartoon_color','crevice_context','crevice_protein')
def crevice_style_protein():
    cmd.hide('everything','crevice_protein');cmd.dss('crevice_protein');cmd.show('cartoon','crevice_protein')
    cmd.set('cartoon_transparency',0,'crevice_protein')
    for setting,value in [('cartoon_sampling',16),('cartoon_loop_quality',16),('cartoon_oval_quality',20),('cartoon_tube_quality',16),('cartoon_fancy_helices',1),('cartoon_fancy_sheets',1),('cartoon_smooth_loops',1),('cartoon_flat_sheets',1),('cartoon_loop_radius',.18),('cartoon_oval_length',1.25),('cartoon_oval_width',.24),('stick_radius',.20),('stick_quality',20),('stick_transparency',0)]:
        cmd.set(setting,value,'crevice_protein')
    for i,row in enumerate(crevice_scene['residues']):
        name='crevice_hydration_'+str(i);cmd.set_color(name,row['color_rgb'])
        selection='crevice_protein and id '+'+'.join(map(str,row['atom_serials']))
        cmd.color(name,selection)
crevice_style_protein()
for surface in crevice_scene['surfaces']:
    key=surface['key']
    if not surface['mesh_available']:continue
    faces=crevice_geometry[key+'_faces'];body=np.empty((len(faces),3,8))
    body[:,:,0]=NORMAL;body[:,:,1:4]=crevice_geometry[key+'_normals'][faces]
    body[:,:,4]=VERTEX;body[:,:,5:8]=crevice_geometry[key+'_vertices'][faces]
    cmd.load_cgo([BEGIN,TRIANGLES,*body.ravel().tolist(),END],surface['name'],state=1,zoom=0)
    cmd.set_color('crevice_'+key+'_color',surface['color_rgb']);cmd.color('crevice_'+key+'_color',surface['name'])
    cmd.set('cgo_transparency',1-surface['opacity'],surface['name'])
crevice_edges=[]
for edge in crevice_scene['contacts']:
    for a,b in zip(edge['points_A'],edge['points_A'][1:]):
        a,b=np.asarray(a),np.asarray(b);segments=max(2,int(np.linalg.norm(b-a)/.25))
        for i in range(0,segments,2):
            p=a+(b-a)*i/segments;q=a+(b-a)*min(i+1,segments)/segments
            crevice_edges.extend([CYLINDER,*p,*q,.065,*edge['color_rgb'],*edge['color_rgb']])
if crevice_edges:cmd.load_cgo(crevice_edges,'crevice_interactions',state=1,zoom=0)
def crevice_member_waters_load(index=1):
    # Own object, shown by default and never toggled by crevice_view; hide with: disable crevice_member_waters
    info=crevice_scene.get('member_waters')
    if not info:return
    shown='crevice_member_waters' not in cmd.get_names('objects') or 'crevice_member_waters' in cmd.get_names('objects',enabled_only=1)
    cmd.delete('crevice_member_waters');entry=info['snapshots'][int(index)-1]
    if entry['pdb']:
        cmd.load(str(_crevice_root/entry['pdb']),'crevice_member_waters',state=1,zoom=0)
        cmd.hide('everything','crevice_member_waters');cmd.alter('crevice_member_waters','vdw='+repr(info['sphere_radius_A']));cmd.rebuild('crevice_member_waters')
        cmd.set_color('crevice_member_water_color',info['color_rgb']);cmd.color('crevice_member_water_color','crevice_member_waters');cmd.show('spheres','crevice_member_waters')
    else:cmd.create('crevice_member_waters','none')
    (cmd.enable if shown else cmd.disable)('crevice_member_waters')
crevice_member_waters_load(1)
crevice_camera=crevice_scene['camera'];crevice_R=np.asarray(crevice_camera['camera_rotation']);crevice_center=np.asarray(crevice_camera['center']);crevice_span=crevice_camera['vertical_span']
__LEGEND__crevice_current_snapshot=1
crevice_current_view='hydration'
def crevice_view(mode='hydration'):
    global crevice_current_view
    mode=str(mode).strip()
    if mode not in {'hydration','network','partners','geometry'}:raise ValueError('Choose hydration, network, partners or geometry')
    crevice_current_view=mode
    shown=set(crevice_scene['display_residues'])
    if mode=='partners':shown.update(crevice_scene['partner_residues'])
    if mode=='network':shown.update(k for e in crevice_scene['contacts'] for k in [e['source'],e['target']])
    cmd.hide('sticks','crevice_protein')
    if mode!='geometry':
        ids=[i for r in crevice_scene['residues'] if r['residue'] in shown for i in r['atom_serials']]
        if ids:cmd.show('sticks','crevice_protein and id '+'+'.join(map(str,ids)))
    for name in ['crevice_water_density','crevice_hydration_legend','crevice_hydration_labels']:
        (cmd.disable if mode=='geometry' else cmd.enable)(name)
    (cmd.enable if mode=='network' and crevice_current_snapshot==1 else cmd.disable)('crevice_interactions')
    cmd.set('cgo_transparency',.45 if mode=='geometry' else .82,'crevice_cavity')
def crevice_frame(index=1):
    global crevice_current_snapshot
    index=int(index)
    if not 1<=index<=len(crevice_scene['snapshots']):raise ValueError('Snapshot index outside supplied range')
    snapshot=crevice_scene['snapshots'][index-1]
    cmd.delete('crevice_protein');cmd.load(str(_crevice_root/snapshot['pdb']),'crevice_protein',state=1,zoom=0)
    cmd.set('cartoon_color','crevice_context','crevice_protein');crevice_style_protein();crevice_member_waters_load(index)
    crevice_current_snapshot=index;crevice_view(crevice_current_view)
    print('CREVICE source frame',snapshot['source_frame_index'],'time ps',snapshot['time_ps'],'; colours and density describe all analysed frames')
cmd.extend('crevice_view',crevice_view);cmd.extend('crevice_frame',crevice_frame)
for setting,value in [('two_sided_lighting',1),('ray_opaque_background',1),('orthoscopic',1),('depth_cue',0),('ray_shadows',0),('ambient',.35),('direct',.65),('light_count',3),('transparency_mode',2),('ray_transparency_oblique',.2),('ray_trace_mode',0),('ambient_occlusion_mode',2),('ambient_occlusion_scale',12),('ambient_occlusion_smooth',10),('specular',.08),('shininess',18),('antialias',2),('ray_max_passes',8)]:cmd.set(setting,value)
cmd.bg_color('white');cmd.set('field_of_view',25)
crevice_distance=crevice_span/(2*.22169466264293988)
cmd.set_view([*crevice_R.T.ravel().tolist(),0,0,-crevice_distance,*crevice_center,crevice_distance-crevice_camera['radius']*1.5,crevice_distance+crevice_camera['radius']*1.5,1])
cmd.viewport(1200,900)
crevice_view('__MODE__')
def crevice_render(path=None):
    target=path or str(_crevice_root/(crevice_scene['prefix']+'_analysis_pymol.png'))
    cmd.png(target,width=3200,height=2400,dpi=400,ray=1)
cmd.extend('crevice_render',crevice_render)
__REGION_STATUS__crevice_geometry.close()
python end
'''

# Annotation-only PyMOL blocks: colour-scale bar with tick labels and caption,
# and the region status banner. Omitted from default scenes.
PML_LEGEND = r'''def crevice_legend_position(x,y):return (crevice_center+np.array([x,y,0])*crevice_span@crevice_R).tolist()
crevice_legend=[]
for i in range(100):
    p=i/99.;a=crevice_legend_position(-.34+.35*i/100,-.425);b=crevice_legend_position(-.34+.35*(i+1)/100,-.425)
    crevice_legend.extend([CYLINDER,*a,*b,crevice_span*.008,1-p,0,p,1-p,0,p])
cmd.load_cgo(crevice_legend,'crevice_hydration_legend',zoom=0)
for x,label in [(-.34,'0%'),(-.165,'50%'),(.01,'100%')]:
    cmd.pseudoatom('crevice_hydration_labels',pos=crevice_legend_position(x,-.465),label=label)
caption=crevice_scene['annotations']['texts']['colour_scale_caption']
cmd.pseudoatom('crevice_hydration_labels',pos=crevice_legend_position(.22,-.425),label=caption)
cmd.hide('everything','crevice_hydration_labels');cmd.show('labels','crevice_hydration_labels')
cmd.set('label_color','black','crevice_hydration_labels');cmd.set('label_outline_color','white','crevice_hydration_labels');cmd.set('label_size',-crevice_span*.018,'crevice_hydration_labels');cmd.set('label_font_id',7,'crevice_hydration_labels')
'''

PML_REGION_STATUS = r'''if crevice_scene.get('region_definition'):
    cmd.pseudoatom('crevice_region_status',pos=crevice_legend_position(0,.43),label=crevice_scene['annotations']['texts']['region_status'])
    cmd.hide('everything','crevice_region_status');cmd.show('labels','crevice_region_status')
    cmd.set('label_color','black','crevice_region_status');cmd.set('label_outline_color','white','crevice_region_status');cmd.set('label_size',-crevice_span*.017,'crevice_region_status')
'''

CXC_PY = r'''from pathlib import Path
import json,numpy as np
from chimerax.core.commands import run,register,CmdDesc,StringArg,IntArg
from chimerax.atomic import AtomicStructure,Atom
root=Path(__file__).resolve().parent
scene=json.loads((root/__SCENE__).read_text())
owner='crevice_analysis_'+scene['prefix']
for old in list(session.models.list()):
    if getattr(old,'crevice_owner',None)==owner:session.models.close([old])
next_id=max([m.id[0] for m in session.models.list() if m.id]+[0])+1
protein_id=next_id;next_id+=1
run(session,'open '+json.dumps(str(root/scene['protein_pdb']))+' id #'+str(protein_id))
protein=next(m for m in session.models.list(type=AtomicStructure) if m.id==(protein_id,))
protein.name='crevice_protein';protein.crevice_owner=owner
expected=np.load(root/scene['geometry_npz'],allow_pickle=False)['protein_coords_A']
serials=np.asarray(protein.atoms.serial_numbers)
if set(serials)!=set(range(1,len(expected)+1)) or not np.allclose(protein.atoms.coords,expected[serials-1],atol=.001,rtol=0):raise ValueError('CREVICE reference atom identities or coordinates differ')
run(session,f'hide #{protein_id} atoms; cartoon #{protein_id}; color #{protein_id} #a8b3bd target c; transparency #{protein_id} 0 target c')
# Assign the array through the property to notify ChimeraX of colour changes.
colors=protein.atoms.colors.copy()
for row in scene['residues']:colors[np.isin(serials,row['atom_serials'])]=np.rint(np.asarray([*row['color_rgb'],1])*255).astype(np.uint8)
protein.atoms.colors=colors
surface_ids={}
for surface in scene['surfaces']:
    model_id=next_id;next_id+=1
    run(session,'open '+json.dumps(str(root/surface['local_dx']))+' id #'+str(model_id))
    model=next(m for m in session.models.list() if m.id==(model_id,));model.name=surface['name'];model.crevice_owner=owner
    transform=np.zeros((3,4));transform[:,:3]=surface['transform']
    run(session,'view matrix models #'+str(model_id)+','+','.join(f'{v:.12g}' for v in transform.ravel()))
    run(session,f'volume #{model_id} style surface level {surface["level"]} step 1 surfaceSmoothing false')
    color='#'+''.join(f'{round(255*v):02x}' for v in surface['color_rgb'])
    run(session,f'color #{model_id} {color}; transparency #{model_id} {100*(1-surface["opacity"])} target s')
    surface_ids[surface['key']]=model_id
contacts_id=None
if scene['contacts']:
    contacts_id=next_id;next_id+=1
    run(session,'open '+json.dumps(str(root/(scene['prefix']+'_analysis_contacts.bild')))+' id #'+str(contacts_id))
    for model in session.models.list():
        if model.id and model.id[0]==contacts_id:model.name='crevice_interactions';model.crevice_owner=owner
crevice_member_model=None
def crevice_member_waters_load(session,index):
    # Own model, shown by default and never toggled by creviceView; hide with the Models panel or: hide #<id> models
    global crevice_member_model
    info=scene.get('member_waters')
    if not info:return
    shown=True
    if crevice_member_model is not None and not crevice_member_model.deleted:
        shown=crevice_member_model.display;session.models.close([crevice_member_model])
    entry=info['snapshots'][index-1]
    if entry['pdb']:
        before=set(session.models.list())
        run(session,'open '+json.dumps(str(root/entry['pdb'])))
        model=next(m for m in session.models.list(type=AtomicStructure) if m not in before)
        model.atoms.draw_modes=Atom.SPHERE_STYLE;model.atoms.radii=np.full(len(model.atoms),info['sphere_radius_A'],dtype=np.float32)
        model.atoms.colors=np.tile(np.rint(np.asarray([*info['color_rgb'],1])*255).astype(np.uint8),(len(model.atoms),1));model.atoms.displays=True
    else:
        from chimerax.core.models import Model
        model=Model('crevice_member_waters',session);session.models.add([model])
    model.name='crevice_member_waters';model.crevice_owner=owner;model.display=shown
    crevice_member_model=model
crevice_member_waters_load(session,1)
run(session,'set bgColor white; camera ortho; lighting soft; material dull; graphics quality 3; graphics silhouettes false; windowsize 1200 900')
camera=scene['camera'];rotation=np.asarray(camera['camera_rotation']);center=np.asarray(camera['center'])
camera_matrix=np.zeros((3,4));camera_matrix[:,:3]=rotation.T;camera_matrix[:,3]=center+rotation[2]*camera['radius']*3
run(session,'view matrix camera '+','.join(f'{v:.12g}' for v in camera_matrix.ravel()))
session.main_view.camera.field_width=camera['vertical_span']*4/3
# Remove text and colour key left by an earlier annotated load of this scene.
from chimerax.label.label2d import session_labels,label_delete
from chimerax.color_key.model import get_model as crevice_key_model
labels=session_labels(session,create=False)
for name in ['creviceHydrationLabel','creviceRegionStatus']:
    if labels and labels.named_label(name):label_delete(session,[labels.named_label(name)])
if crevice_key_model(session,create=False) is not None:session.models.close([crevice_key_model(session,create=False)])
__LEGEND__crevice_current_view='hydration';crevice_current_snapshot=1
crevice_direct=set(scene['display_residues'])
def crevice_view(session,mode='hydration'):
    global crevice_current_view
    if mode not in {'hydration','network','partners','geometry'}:raise ValueError('Choose hydration, network, partners or geometry')
    crevice_current_view=mode;shown=set(crevice_direct)
    from chimerax.color_key.model import get_model
    legend=get_model(session,create=False)
    if legend is not None:legend.display=mode!='geometry'
__LABEL_VISIBILITY__    if mode=='partners':shown.update(scene['partner_residues'])
    if mode=='network':shown.update(k for e in scene['contacts'] for k in [e['source'],e['target']])
    ids=[i for r in scene['residues'] if r['residue'] in shown for i in r['atom_serials']]
    protein.atoms.displays=np.isin(serials,ids) if mode!='geometry' else np.zeros(len(serials),dtype=bool)
    protein.atoms.draw_modes=Atom.STICK_STYLE
    run(session,f'size #{protein_id} stickRadius 0.20')
    if 'water' in surface_ids:run(session,('hide' if mode=='geometry' else 'show')+' #'+str(surface_ids['water'])+' models')
    if contacts_id is not None:run(session,('show' if mode=='network' and crevice_current_snapshot==1 else 'hide')+' #'+str(contacts_id)+' models')
    if 'cavity' in surface_ids:run(session,'transparency #'+str(surface_ids['cavity'])+(' 45' if mode=='geometry' else ' 82')+' target s')
def crevice_frame(session,index):
    global crevice_current_snapshot
    if not 1<=index<=len(scene['snapshots']):raise ValueError('Snapshot index outside supplied range')
    snapshot=scene['snapshots'][index-1];coordinates={}
    for line in (root/snapshot['pdb']).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):coordinates[int(line[6:11])]=[float(line[k:k+8]) for k in [30,38,46]]
    protein.atoms.coords=np.asarray([coordinates[int(i)] for i in serials]);crevice_member_waters_load(session,index)
    crevice_current_snapshot=index;crevice_view(session,crevice_current_view)
    session.logger.info('CREVICE source frame '+str(snapshot['source_frame_index'])+'; colours and density describe all analysed frames')
register('creviceView',CmdDesc(optional=[('mode',StringArg)]),crevice_view,logger=session.logger)
register('creviceFrame',CmdDesc(required=[('index',IntArg)]),crevice_frame,logger=session.logger)
crevice_view(session,'__MODE__')
'''

# Annotation-only ChimeraX blocks: colour key with tick labels, caption and
# region status banner. Omitted from default scenes.
CXC_LEGEND = r'''run(session,'key #ff0000:0% #800080:50% #0000ff:100% pos 0.10,0.15 size 0.28,0.035 fontSize 14')
run(session,'2dlabels create creviceHydrationLabel text '+json.dumps(scene['annotations']['texts']['colour_scale_caption'])+' color black xpos 0.43 ypos 0.16 size 16')
if scene.get('region_definition'):
    run(session,'2dlabels create creviceRegionStatus text '+json.dumps(scene['annotations']['texts']['region_status'])+' color black xpos 0.15 ypos 0.94 size 17')
'''
CXC_LABEL_VISIBILITY = "    run(session,'2dlabels change creviceHydrationLabel visibility '+('false' if mode=='geometry' else 'true'))\n"


def _tcl(scene,root,mode,annotate=False):
    import numpy as np
    from .volume_export import _tcl_word as q
    lines=['# CREVICE coordinated hydration; absolute occupancy scale.',
           'set crevice_root [file dirname [file normalize [info script]]]',
           'if {[info exists crevice_owned_molecules]} {foreach m $crevice_owned_molecules {if {[lsearch -exact [molinfo list] $m]>=0} {mol delete $m}}}',
           'catch {array unset crevice_surface}',
           'set crevice_protein [mol new [file join $crevice_root '+q(scene['protein_pdb'])+'] waitfor all]',
           'mol delrep 0 $crevice_protein',
           'color change rgb 2 0.66 0.70 0.74','color change rgb 3 0.55 0.55 0.55',
           'if {[lsearch -exact [material list] CREVICEProtein] < 0} {material add CREVICEProtein}',
           'material change opacity CREVICEProtein 1.0','material change ambient CREVICEProtein 0.30','material change diffuse CREVICEProtein 0.70','material change specular CREVICEProtein 0.08',
           'mol representation NewCartoon 0.25 24 4.0 1','mol color ColorID 2','mol selection all','mol material CREVICEProtein','mol addrep $crevice_protein',
           'if {[lsearch -exact [material list] CREVICESticks] < 0} {material add CREVICESticks}',
           'material change opacity CREVICESticks 1.0','material change specular CREVICESticks 0.08',
           'set crevice_all [atomselect $crevice_protein all]','$crevice_all set user -1',
           'set crevice_palette_start [colorinfo num]','set crevice_palette_size [expr {[colorinfo max]-[colorinfo num]}]',
           'color scale method RGB','color scale min 0','color scale max 1',
           'for {set i 0} {$i<$crevice_palette_size} {incr i} {set p [expr {double($i)/($crevice_palette_size-1)}]; color change rgb [expr {$crevice_palette_start+$i}] [expr {1-$p}] 0 $p}']
    expected=root/(scene['prefix']+'_analysis_reference.tsv')
    geometry=np.load(root/scene['geometry_npz']);expected.write_text(''.join(f'{i}\t{p[0]:.9f}\t{p[1]:.9f}\t{p[2]:.9f}\n' for i,p in enumerate(geometry['protein_coords_A'],1)));geometry.close()
    lines += ['set crevice_check [open [file join $crevice_root '+q(expected.name)+'] r]',
              'set crevice_expected [dict create]',
              'foreach line [split [read $crevice_check] "\\n"] {if {[llength $line]==4} {dict set crevice_expected [lindex $line 0] [lrange $line 1 3]}}',
              'close $crevice_check',
              'if {[$crevice_all num]!=[dict size $crevice_expected]} {error "CREVICE reference atom count mismatch"}',
              'foreach serial [$crevice_all get serial] xyz [$crevice_all get {x y z}] {if {![dict exists $crevice_expected $serial]} {error "CREVICE serial mismatch"}; foreach a $xyz b [dict get $crevice_expected $serial] {if {abs($a-$b)>0.001} {error "CREVICE reference coordinate mismatch"}}}']
    groups={'hydration':set(scene['display_residues']),'partners':{r['residue'] for r in scene['residues']}}
    groups['network']=groups['hydration']|{k for e in scene['contacts'] for k in [e['source'],e['target']]}
    for name,keys in groups.items():
        values=[i for r in scene['residues'] if r['residue'] in keys for i in r['atom_serials']]
        lines.append('set crevice_membership('+name+') {'+' '.join(map(str,values))+'}')
    for row in scene['residues']:
        if row['hydration_occupancy'] is not None:
            lines+=['set sel [atomselect $crevice_protein {serial '+' '.join(map(str,row['atom_serials']))+'}]',
                    f'$sel set user {row["hydration_occupancy"]:.12g}','$sel delete']
    lines += ['$crevice_all delete','set crevice_stick_reps {}']
    for condition,color in [('user >= 0','User'),('user < 0','ColorID 3')]:
        lines+=['mol representation Licorice 0.20 24 24','mol selection {none}','mol color '+color,
                'mol material CREVICESticks','mol addrep $crevice_protein',
                'set rep [expr {[molinfo $crevice_protein get numreps]-1}]','mol scaleminmax $crevice_protein $rep 0 1',
                'lappend crevice_stick_reps [list $rep {'+condition+'}]']
    lines+=['set crevice_surface_molecules {}','array set crevice_surface {}']
    for i,surface in enumerate(scene['surfaces'],10):
        key=surface['key'];color=' '.join(map(str,surface['color_rgb']))
        lines+=['set m [mol new [file join $crevice_root '+q(surface['dx'])+'] type dx waitfor all]',
                'mol delrep 0 $m',f'color change rgb {i} {color}',
                'if {[lsearch -exact [material list] CREVICE_'+key+']<0} {material add CREVICE_'+key+'}',
                f'material change opacity CREVICE_{key} {surface["opacity"]}',f'material change specular CREVICE_{key} 0.08',
                f'mol representation Isosurface {surface["level"]} 0 0 0 1 1',f'mol color ColorID {i}',
                'mol selection all',f'mol material CREVICE_{key}','mol addrep $m',
                'set crevice_surface('+key+') $m','lappend crevice_surface_molecules $m']
    lines+=['set crevice_contacts [mol new]','mol rename $crevice_contacts crevice_interactions']
    for i,edge in enumerate(scene['contacts']):
        color_id=20+['hydrogen_bond','water_bridge_hbond','salt_bridge_candidate','aromatic_parallel_candidate','aromatic_edge_face_candidate'].index(edge['kind'])
        lines += [f'color change rgb {color_id} '+' '.join(map(str,edge['color_rgb'])),f'graphics $crevice_contacts color {color_id}']
        for left,right in zip(edge['points_A'],edge['points_A'][1:]):
            a,b=np.asarray(left),np.asarray(right);segments=max(2,int(np.linalg.norm(b-a)/.25))
            for j in range(0,segments,2):
                p=a+(b-a)*j/segments;v=a+(b-a)*min(j+1,segments)/segments
                lines += ['graphics $crevice_contacts cylinder {'+' '.join(f'{x:.9g}' for x in p)+'} {'+' '.join(f'{x:.9g}' for x in v)+'} radius 0.065 resolution 12 filled yes']
    camera=scene['camera'];R=np.asarray(camera['camera_rotation']);center=np.asarray(camera['center']);span=camera['vertical_span']
    def position(x,y):return center+np.asarray([x,y,0])*span@R
    # The legend molecule always exists (the view procedures toggle it) but is
    # empty unless annotating.
    lines+=['set crevice_legend [mol new]','mol rename $crevice_legend crevice_hydration_legend']
    for i in range(100 if annotate else 0):
        color=f'[expr {{$crevice_palette_start+int({i}/99.0*($crevice_palette_size-1))}}]'
        a,b=position(-.34+.35*i/100,-.425),position(-.34+.35*(i+1)/100,-.425)
        lines+=['graphics $crevice_legend color '+color,
                'graphics $crevice_legend cylinder {'+' '.join(map(str,a))+'} {'+' '.join(map(str,b))+'} radius '+str(span*.008)+' resolution 12 filled yes']
    if annotate:
        texts=scene['annotations']['texts']
        lines+=['graphics $crevice_legend color black']
        for x,text in [(-.34,'0%'),(-.165,'50%'),(.01,'100%')]:
            lines+=['graphics $crevice_legend text {'+' '.join(map(str,position(x,-.465)))+'} '+q(text)+' size 0.65 thickness 1.8']
        lines+=['graphics $crevice_legend text {'+' '.join(map(str,position(.06,-.425)))+'} '+q(texts['colour_scale_caption'])+' size 0.65 thickness 1.8']
        if texts.get('region_status'):
            lines+=['graphics $crevice_legend text {'+' '.join(map(str,position(-.34,.43)))+'} '+q(texts['region_status'])+' size 0.75 thickness 1.8']
    matrices={'center':np.eye(4),'rotate':np.eye(4),'scale':np.eye(4)}
    matrices['center'][:3,3]=-center;matrices['rotate'][:3,:3]=R;matrices['scale'][:3,:3]*=2/span
    lines += ['color Display Background white','display projection Orthographic','display depthcue off','display antialias on','display shadows off',
              'render aasamples TachyonInternal 8','render aosamples TachyonInternal 8',
              'display ambientocclusion on','display aoambient 0.8','display aodirect 0.4','axes location Off','display resetview',
              'foreach molecule [concat [list $crevice_protein $crevice_contacts $crevice_legend] $crevice_surface_molecules] {']
    for key,matrix in matrices.items():
        value='{'+' '.join('{'+' '.join(f'{v:.12g}' for v in row)+'}' for row in matrix)+'}'
        lines+=['molinfo $molecule set '+key+'_matrix [list '+value+']']
    lines+=['}','display height 4.0','display resize 1200 900','set crevice_current_snapshot 1','set crevice_current_view hydration',
            'set crevice_owned_molecules [concat [list $crevice_protein $crevice_contacts $crevice_legend] $crevice_surface_molecules]',
            'set crevice_snapshot_files [list '+' '.join(q(s['pdb']) for s in scene['snapshots'])+']',
            'set crevice_member_info '+('1' if scene.get('member_waters') else '0'),
            'set crevice_member_files [list '+' '.join(q(s['pdb'] or '') for s in (scene.get('member_waters') or {}).get('snapshots',[]))+']',
            'color change rgb 30 '+' '.join(map(str,(scene.get('member_waters') or {}).get('color_rgb',[1.,.5,0.]))),
            'set crevice_member_radius '+str((scene.get('member_waters') or {}).get('sphere_radius_A',.7)),r'''
proc crevice_member_waters_load {index} {
    # Own molecule, shown by default and never toggled by crevice_view; hide with: mol off $crevice_member_waters
    global crevice_root crevice_member_info crevice_member_files crevice_member_waters crevice_member_radius crevice_protein crevice_owned_molecules
    if {!$crevice_member_info} {return}
    set shown 1
    set saved [dict create]
    foreach m [molinfo list] {dict set saved $m [molinfo $m get {center_matrix rotate_matrix scale_matrix global_matrix}]}
    if {[info exists crevice_member_waters] && [lsearch -exact [molinfo list] $crevice_member_waters]>=0} {
        set shown [molinfo $crevice_member_waters get displayed];mol delete $crevice_member_waters
        set crevice_owned_molecules [lsearch -all -inline -not -exact $crevice_owned_molecules $crevice_member_waters]
    }
    set f [lindex $crevice_member_files [expr {$index-1}]]
    if {$f eq ""} {set crevice_member_waters [mol new]} else {
        set crevice_member_waters [mol new [file join $crevice_root $f] waitfor all];mol delrep 0 $crevice_member_waters
        set s [atomselect $crevice_member_waters all];$s set radius $crevice_member_radius;$s delete
        mol representation VDW 1.0 20;mol color ColorID 30;mol selection all;mol material CREVICESticks;mol addrep $crevice_member_waters
    }
    mol rename $crevice_member_waters crevice_member_waters
    dict for {m matrices} $saved {if {[lsearch -exact [molinfo list] $m]>=0} {molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $matrices}}
    molinfo $crevice_member_waters set {center_matrix rotate_matrix scale_matrix global_matrix} [molinfo $crevice_protein get {center_matrix rotate_matrix scale_matrix global_matrix}]
    if {!$shown} {mol off $crevice_member_waters}
    mol top $crevice_protein
    lappend crevice_owned_molecules $crevice_member_waters
}
proc crevice_view {{mode hydration}} {
    global crevice_protein crevice_contacts crevice_legend crevice_surface crevice_membership crevice_stick_reps crevice_current_snapshot crevice_current_view
    if {$mode ni {hydration network partners geometry}} {error "Choose hydration, network, partners or geometry"}
    set crevice_current_view $mode
    foreach entry $crevice_stick_reps {
        lassign $entry rep condition
        if {$mode eq "geometry"} {set selection none} else {set selection "($condition) and serial [join $crevice_membership($mode) { }]"}
        mol modselect $rep $crevice_protein $selection
    }
    if {[info exists crevice_surface(water)]} {if {$mode eq "geometry"} {mol off $crevice_surface(water)} else {mol on $crevice_surface(water)}}
    if {$mode eq "geometry"} {mol off $crevice_legend} else {mol on $crevice_legend}
    if {$mode eq "network" && $crevice_current_snapshot==1} {mol on $crevice_contacts} else {mol off $crevice_contacts}
    if {[info exists crevice_surface(cavity)]} {material change opacity CREVICE_cavity [expr {$mode eq "geometry" ? .55 : .18}]}
    display update
}
proc crevice_frame {index} {
    global crevice_root crevice_protein crevice_snapshot_files crevice_current_snapshot crevice_current_view
    if {$index<1 || $index>[llength $crevice_snapshot_files]} {error "Snapshot index outside supplied range"}
    set crevice_saved_matrices [dict create]
    foreach m [molinfo list] {dict set crevice_saved_matrices $m [molinfo $m get {center_matrix rotate_matrix scale_matrix global_matrix}]}
    set temporary [mol new [file join $crevice_root [lindex $crevice_snapshot_files [expr {$index-1}]]] waitfor all]
    set source [atomselect $temporary all];set target [atomselect $crevice_protein all]
    if {[$source get serial] ne [$target get serial]} {error "Snapshot atom identities differ"}
    $target set {x y z} [$source get {x y z}]
    $source delete;$target delete;mol delete $temporary
    dict for {m matrices} $crevice_saved_matrices {molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $matrices}
    crevice_member_waters_load $index
    set crevice_current_snapshot $index;crevice_view $crevice_current_view
}
proc crevice_render {} {
    global crevice_root
    render aasamples TachyonInternal 12;render aosamples TachyonInternal 12
    display resize 2400 1800;catch {display redraw};display update
    render TachyonInternal [file join $crevice_root __RENDER__]
}
'''.replace('__RENDER__',q(scene['prefix']+'_analysis_vmd.tga')), 'crevice_member_waters_load 1','crevice_view '+mode]
    return '\n'.join(lines)+'\n'


def scene_texts(scene):
    """Text elements of a coordinated hydration scene, drawn only when annotating."""
    caption='Water-contact occupancy' if scene['frame_count']>1 else 'Observed snapshot contact'
    if scene['water_status']=='unavailable_no_explicit_water':caption='Water contacts unobserved (grey)'
    texts={'colour_scale_caption':caption,'colour_scale_ticks':'0% (red), 50% (purple), 100% (blue); grey = unobserved'}
    region=scene.get('region_definition')
    if region:texts['region_status']=region['review']['status'].upper()+' REGION | '+region['region_id']
    return texts


@annotate_option
def write_scene_scripts(output_dir,prefix,*,annotate=None):
    """Write the coordinated hydration scenes for PyMOL, VMD and ChimeraX.

    Every script reads the same ``PREFIX_analysis_scene.json`` (written by
    :func:`crevice.analysis_viewers.write_analysis_assets`) and geometry NPZ,
    checks that the loaded protein atom identities and coordinates match the
    reference to 0.001 Å, and builds the same scene: cartoon protein, residues
    coloured on the absolute hydration scale, boundary-contributor sticks, the
    reference cavity, the water-density surface, typed-contact lines and, with
    ``--water-membership``, the member-water group. Each viewer gains view and
    snapshot commands (``crevice_view``/``crevice_frame`` in PyMOL and VMD,
    ``creviceView``/``creviceFrame`` in ChimeraX) and a render command. By
    default the scenes contain no text; the colour-scale caption, tick labels
    and region banner are recorded in the scene JSON ``annotations`` entry.

    Parameters
    ----------
    output_dir : str or Path
        Directory containing the scene JSON and its assets.
    prefix : str
        File-name prefix.
    annotate : bool, optional
        Also draw the colour-scale bar with 0/50/100% tick labels and its caption
        (PyMOL CGO bar and pseudoatom labels, VMD graphics, ChimeraX ``key`` and
        2D label) and, for named regions, the "STATUS REGION | id" banner.
        ``None`` (default) inherits the surrounding setting, which is off unless
        enabled with :func:`crevice.presentation.figure_annotations` or ``--annotate``.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_analysis.pml, PREFIX_analysis_network.pml : PyMOL
        Start in hydration and network view respectively.
    PREFIX_analysis.vmd, PREFIX_analysis.tcl (and ``_network`` versions) : VMD
        Identical Tcl scripts under both extensions.
    PREFIX_analysis.cxc, PREFIX_analysis_chimerax.py (and ``_network`` versions) : ChimeraX
        The ``.cxc`` opens the Python loader.
    PREFIX_analysis_contacts.bild : BILD
        Typed-contact dashes for ChimeraX.
    PREFIX_analysis_render.cxc : ChimeraX
        Saves ``PREFIX_analysis_chimerax.png`` (2400 x 1800, supersampled).
    PREFIX_analysis_reference.tsv : table
        Reference coordinates checked by the VMD script.

    Colours and representation
    --------------------------
    * Protein: opaque light-grey cartoon (``[0.66, 0.70, 0.74]``).
    * Residue hydration colour (atoms of focused residues): absolute linear
      scale red (0% of frames with a contact water) to blue (100%),
      ``[1-p, 0, p]``; grey (``[0.55, 0.55, 0.55]``) = no water observations.
    * Sticks (radius 0.20 Å, opaque): direct boundary contributors in hydration
      view; plus all focused partners in ``partners`` view; plus typed-contact
      endpoints in ``network`` view; none in ``geometry`` view.
    * Reference cavity: pale grey-green surface (``[0.69, 0.73, 0.73]``),
      18% opaque (45% transparent in geometry view).
    * Water density: teal surface (``[0.08, 0.66, 0.59]``), 42% opaque, at the
      display isovalue; hidden in geometry view.
    * Typed contacts (network view, first snapshot only): dashed 0.065 Å
      cylinders, gold hydrogen bonds, teal water-bridged hydrogen bonds, purple
      salt-bridge candidates, blue aromatic candidates.
    * Member waters (``crevice_member_waters``): orange (``[1.0, 0.5, 0.0]``)
      spheres of radius 0.7 Å at the displayed snapshot's member oxygens; own
      group, shown by default, never toggled by view commands, reloaded by
      snapshot commands keeping a hidden state.
    * When annotating: a red-purple-blue bar with 0%/50%/100% labels and caption
      below the protein, and the region banner above it.
    * White background, orthographic camera, ambient occlusion.
    """
    import numpy as np
    root=Path(output_dir).resolve();scene_path=root/(prefix+'_analysis_scene.json');scene=json.loads(scene_path.read_text())
    annotate=annotations_enabled()
    scene['annotations']=scene_annotation_record(scene_texts(scene),drawn=annotate)
    scene_path.write_text(json.dumps(scene,indent=2,allow_nan=False)+'\n')
    files={}
    for mode,ending in [('hydration','_analysis'),('network','_analysis_network')]:
        p=root/(prefix+ending+'.pml')
        pml=PML.replace('__LEGEND__',PML_LEGEND if annotate else '').replace('__REGION_STATUS__',PML_REGION_STATUS if annotate else '')
        p.write_text(pml.replace('__SCENE__',repr(scene_path.name)).replace('__FALLBACK__',repr(str(root))).replace('__MODE__',mode));files[mode+'_pml']=str(p)
        text=_tcl(scene,root,mode,annotate)
        for extension in ['vmd','tcl']:
            p=root/(prefix+ending+'.'+extension);p.write_text(text);files[mode+'_'+extension]=str(p)
        loader=root/(prefix+ending+'_chimerax.py')
        cxc=CXC_PY.replace('__LEGEND__',CXC_LEGEND if annotate else '').replace('__LABEL_VISIBILITY__',CXC_LABEL_VISIBILITY if annotate else '')
        loader.write_text(cxc.replace('__SCENE__',repr(scene_path.name)).replace('__MODE__',mode));files[mode+'_chimerax_python']=str(loader)
        p=root/(prefix+ending+'.cxc');p.write_text('open '+json.dumps(loader.name)+'\n');files[mode+'_cxc']=str(p)
    lines=[]
    for edge in scene['contacts']:
        lines.append('.color '+' '.join(map(str,edge['color_rgb'])))
        for left,right in zip(edge['points_A'],edge['points_A'][1:]):
            a,b=np.asarray(left),np.asarray(right);segments=max(2,int(np.linalg.norm(b-a)/.25))
            for i in range(0,segments,2):
                p=a+(b-a)*i/segments;q=a+(b-a)*min(i+1,segments)/segments
                lines.append('.cylinder '+' '.join(f'{v:.9g}' for v in [*p,*q])+' 0.065')
    p=root/(prefix+'_analysis_contacts.bild');p.write_text('\n'.join(lines)+'\n');files['analysis_contacts_bild']=str(p)
    files['analysis_reference_tsv']=str(root/(prefix+'_analysis_reference.tsv'))
    render=root/(prefix+'_analysis_render.cxc');render.write_text('open '+json.dumps(prefix+'_analysis.cxc')+'\nsave '+json.dumps(prefix+'_analysis_chimerax.png')+' width 2400 height 1800 supersample 3\n');files['analysis_chimerax_render']=str(render)
    return files
