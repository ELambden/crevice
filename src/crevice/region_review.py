"""Native review overlays for named regions, separate from hydration colours."""
from pathlib import Path
import json

from .presentation import annotate_option, annotations_enabled, scene_annotation_record


@annotate_option
def write_region_review(frame, reference, report, obstacles, output_dir, prefix, *, annotate=None):
    """Write PyMOL, VMD and ChimeraX review scenes and a landmark table for a named region.

    The declared reference region (a binary DX in the reference frame) is shown
    as a surface inside the protein, with the definition's residue landmarks as
    sticks, transferred coordinate landmarks as small spheres and nearby
    non-protein obstacle atoms (within 6 Å of the reference samples) as thin
    sticks, so that a person can review the region before it is used. The
    review status (candidate, reviewed, ...) is recorded in the scene JSON and,
    with ``annotate``, drawn as a caption; nothing here promotes a candidate to an anatomical or functional assignment.

    Parameters
    ----------
    frame : StructureFrame
        Protein reference frame.
    reference : CavityReference
        Declared reference region.
    report : dict
        Review diagnostics from :func:`crevice.region_definition.prepare_region`.
    obstacles : sequence of Atom
        Physical-obstacle atoms; nearby non-protein ones are shown as context.
    output_dir : str or Path
        Output directory.
    prefix : str
        Region ID used as file prefix.
    annotate : bool, optional
        Also draw residue-landmark labels, the "STATUS REGION | id" caption and
        the colour legend line in the viewer scenes. ``None`` (default) inherits
        the surrounding setting, which is off unless enabled with
        :func:`crevice.presentation.figure_annotations` or ``--annotate``.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_review.pml, PREFIX_review.vmd, PREFIX_review.tcl, PREFIX_review.cxc, PREFIX_review_chimerax.py : scenes
        Review scenes built on the region's volume scene.
    PREFIX_review_scene.json : JSON
        Landmarks, colours, camera, caption, legend, warnings and the
        ``annotations`` record of the scene text.
    PREFIX_protein.pdb, PREFIX_context.pdb : structures
        Viewer protein and nearby non-protein obstacle atoms.
    PREFIX_geometry_* : volume scene
        See :func:`crevice.volume_export.write_volume_viewer_bundle`.
    PREFIX_review_landmarks.csv : table
        One row per residue landmark: ``residue``, ``anchor_atoms``
        (semicolon-separated), ``nearest_reference_sample_A`` and
        ``centroid_{x,y,z}_A``. Status, rationale, reference volume (Å³),
        obstacle overlap and warnings are in ``PREFIX_review.json``; no HTML
        page is written.

    Colours and representation
    --------------------------
    * Protein: grey-blue cartoon (``#8fa1b3``), opaque.
    * Residue landmarks: gold/orange sticks (``[0.88, 0.53, 0.12]``, radius 0.21 Å).
    * Transferred coordinate landmarks: purple spheres (``[0.60, 0.34, 0.72]``).
    * Nearby non-protein obstacles: green thin sticks (``[0.32, 0.56, 0.40]``).
    * Reference region: teal surface (``[0.08, 0.58, 0.63]``), 35% transparent
      in PyMOL, 65% opaque in VMD and ChimeraX. Teal is kept here although the
      region is not a channel: the non-channel violet would be hard to tell
      from the purple transferred landmarks.
    * When annotating: black residue names at landmark centroids, the status
      caption above and the legend line below the protein.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    from .models import ChannelPoint,VoidComponent,VoidCast,StructureFrame
    from .presentation import write_viewer_structure
    from .volume_export import write_volume_viewer_bundle,_tcl_word
    from .region_definition import write_json
    root=Path(output_dir)
    structure=write_viewer_structure(frame,root/(prefix+'_protein.pdb'))
    points=tuple(ChannelPoint(i,tuple(map(float,p)),0.,0.,float((p-reference.anchor)@reference.basis[2])) for i,p in enumerate(reference.points))
    h=float(np.linalg.norm(reference.deltas[0]));basis=reference.deltas/h
    if not np.allclose(basis@basis.T,np.eye(3),atol=1e-7):raise ValueError('Region review requires equally spaced orthogonal reference DX samples')
    if np.linalg.det(basis)<0:basis[0]*=-1
    component=VoidComponent(1,'declared_reference',tuple(reference.points.mean(0)),len(points)*h**3,points)
    cast=VoidCast((component,),points,h,0.,'rolling',{'grid_basis':basis.tolist(),'algorithm':'declared_reference'})
    files=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=root,prefix=prefix+'_geometry',frame=frame,cast_extension=0,
                                     cast_color=(.08,.58,.63))  # teal: violet would clash with the purple site points
    files['protein_pdb']=str(structure)
    context=[a for a,d in zip(obstacles,cKDTree(reference.points).query([a.coord for a in obstacles])[0]) if a.hetero and d<=6]
    context_path=None
    if context:
        context_path=write_viewer_structure(StructureFrame(tuple(context),source=frame.source),root/(prefix+'_context.pdb'))
        files['context_pdb']=str(context_path)
    scene={'schema_version':1,'region_id':prefix,'label':report['label'],'review_status':report['review']['status'],
      'protein_pdb':structure.name,'context_pdb':context_path.name if context_path else None,
      'context_definition':'Nonprotein physical obstacle atom images within 6 A of reference samples; local context, not a complete membrane or ligand model',
      'landmarks':report['landmarks'],'coordinate_landmarks':report['coordinate_landmarks'],
      'landmark_rgb':[.88,.53,.12],'coordinate_landmark_rgb':[.60,.34,.72],
      'context_rgb':[.32,.56,.40],'cavity_rgb':[.08,.58,.63],
      'camera':json.loads(Path(files['scene_json']).read_text()),
      'caption':report['review']['status'].upper()+' REGION | '+prefix,
      'legend':'Gold: residue landmarks   Purple: transferred site points   Teal: reference region',
      'warnings':report['warnings']}
    annotate=annotations_enabled()
    scene['annotations']=scene_annotation_record({'caption':scene['caption'],'legend':scene['legend'],
        'residue_labels':', '.join(row['residue'] for row in report['landmarks'])},drawn=annotate)
    scene_path=root/(prefix+'_review_scene.json');write_json(scene_path,scene);files['review_scene_json']=str(scene_path)
    ids=sorted(i for row in report['landmarks'] for i in row['residue_atom_serials'])
    mapping=json.loads(structure.with_suffix('.identities.json').read_text())
    mapped={r['original']:(r['viewer_chain'],r['viewer_resid']) for r in mapping['residue_mapping']}
    # The underlying volume scenes supply the same world-coordinate cast and camera.
    pml=Path(files['volume_pml']).read_text()+f'''\npython
review=json.loads(Path(crevice_input({scene_path.name!r},{str(scene_path)!r})).read_text())
for name in ['crevice_review_labels','crevice_review_points','crevice_review_context']:
    cmd.delete(name)
crevice_review_ids={ids!r}
cmd.set('cartoon_transparency',0,'crevice_protein')
cmd.set('cartoon_color','crevice_context','crevice_protein')
cmd.set_color('crevice_landmark',review['landmark_rgb'])
if crevice_review_ids:
    selection='crevice_protein and id '+'+'.join(map(str,crevice_review_ids))
    cmd.show('sticks',selection);cmd.color('crevice_landmark',selection)
cmd.set('stick_radius',.21,'crevice_protein');cmd.set('stick_quality',24)
cmd.set('cgo_transparency',.35,'crevice_volume')
view=cmd.get_view()
if review['context_pdb']:
    cmd.load(str(Path(crevice_input({scene_path.name!r},{str(scene_path)!r})).parent/review['context_pdb']),'crevice_review_context',zoom=0)
    cmd.hide('everything','crevice_review_context');cmd.show('sticks','crevice_review_context')
    cmd.set_color('crevice_other_obstacle',review['context_rgb']);cmd.color('crevice_other_obstacle','crevice_review_context');cmd.set('stick_radius',.12,'crevice_review_context')
crevice_point_index=0
for row in review['coordinate_landmarks']:
    for p in row['points_A']:
        crevice_point_index+=1
        cmd.pseudoatom('crevice_review_points',pos=p,resi=str(crevice_point_index),name='PT')
if review['coordinate_landmarks']:
    cmd.set_color('crevice_transferred_site',review['coordinate_landmark_rgb']);cmd.color('crevice_transferred_site','crevice_review_points')
    cmd.hide('everything','crevice_review_points');cmd.show('spheres','crevice_review_points');cmd.set('sphere_scale',.16,'crevice_review_points')
__REVIEW_LABELS__cmd.set_view(view)
python end
'''.replace('__REVIEW_LABELS__',PML_REVIEW_LABELS if annotate else '')
    target=root/(prefix+'_review.pml');target.write_text(pml);files['review_pml']=str(target)
    q=_tcl_word
    tcl='''if {[info exists crevice_review_owned]} {
    foreach m $crevice_review_owned {if {[lsearch -exact [molinfo list] $m]>=0} {mol delete $m}}
}
'''+Path(files['volume_tcl']).read_text()
    tcl+='\nmaterial change opacity CREVICEContext 1.0\nmaterial change opacity CREVICECast 0.65\n'
    tcl+='set crevice_review_ids {'+' '.join(map(str,ids))+'}\n'
    tcl+='color change rgb 17 .88 .53 .12\ncolor change rgb 18 .60 .34 .72\ncolor change rgb 19 .32 .56 .40\n'
    tcl+='mol representation Licorice 0.21 24 24\nmol selection {serial '+' '.join(map(str,ids))+'}\nmol color ColorID 17\nmol material Opaque\nmol addrep $crevice_protein\n' if ids else ''
    tcl+='set crevice_review_owned [list $crevice_protein $crevice_volume]\n'
    tcl+='set crevice_saved [molinfo $crevice_protein get {center_matrix rotate_matrix scale_matrix global_matrix}]\n'
    if context_path:
        tcl+='set crevice_review_context [mol new [file join $crevice_root '+q(context_path.name)+'] waitfor all]\nmol delrep 0 $crevice_review_context\nmol representation Licorice 0.12 16 16\nmol selection all\nmol color ColorID 19\nmol material Opaque\nmol addrep $crevice_review_context\nlappend crevice_review_owned $crevice_review_context\n'
    tcl+='set crevice_review_graphics [mol new]\nmol rename $crevice_review_graphics crevice_review_points_and_labels\nlappend crevice_review_owned $crevice_review_graphics\ngraphics $crevice_review_graphics color 18\n'
    for row in scene['coordinate_landmarks']:
        for p in row['points_A']:tcl+='graphics $crevice_review_graphics sphere {'+' '.join(map(str,p))+'} radius 0.25 resolution 20\n'
    if annotate:
        tcl+='graphics $crevice_review_graphics color black\n'
        for row in scene['landmarks']:
            tcl+='graphics $crevice_review_graphics text {'+' '.join(map(str,row['centroid_A']))+'} '+q(row['residue'])+' size 0.65 thickness 1.8\n'
        R=np.asarray(scene['camera']['camera_rotation']);center=np.asarray(scene['camera']['center']);span=scene['camera']['vertical_span']
        for y,caption in [(.45,scene['caption']),(-.45,scene['legend'])]:
            # VMD graphics text is left-anchored.
            p=center+np.array([-.43,y,0])*span@R
            tcl+='graphics $crevice_review_graphics text {'+' '.join(map(str,p))+'} '+q(caption)+' size 0.65 thickness 1.8\n'
    tcl+='foreach m $crevice_review_owned {molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $crevice_saved}\ndisplay update\n'
    for ext in ['tcl','vmd']:
        target=root/(prefix+'_review.'+ext);target.write_text(tcl);files['review_'+ext]=str(target)
    # A Python loader avoids fixed model IDs and makes actual native reload safe.
    cxc_py=r'''from pathlib import Path
import json,numpy as np
from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure,Atom
from chimerax.geometry import Place
from chimerax.core.models import Surface
from chimerax.markers import MarkerSet
root=Path(__file__).resolve().parent
review=json.loads((root/__SCENE__).read_text())
owner='crevice_region_review'
for m in list(session.models.list()):
    if getattr(m,'crevice_region_owner',None)==owner:session.models.close([m])
protein=run(session,'open '+json.dumps(str(root/review['protein_pdb'])))[0]
protein.name='crevice_protein';protein.crevice_region_owner=owner
run(session,'hide #'+protein.id_string+' atoms; cartoon #'+protein.id_string)
protein.atoms.colors=np.array([143,161,179,255],np.uint8)
protein.residues.ribbon_colors=np.array([143,161,179,255],np.uint8)
ids=set(__IDS__)
chosen=np.asarray([int(i) in ids for i in protein.atoms.serial_numbers])
protein.atoms.displays=chosen;protein.atoms.draw_modes=Atom.STICK_STYLE
protein.atoms[chosen].colors=np.array([224,135,31,255],np.uint8)
__ATOM_LABELS__mesh=np.load(root/__MESH__,allow_pickle=False)
basis=np.asarray(__BASIS__)
surface=Surface('crevice_volume',session);surface.crevice_region_owner=owner
surface.set_geometry((mesh['volume_vertices']@basis).astype(np.float32),(mesh['volume_normals']@basis).astype(np.float32),mesh['volume_faces'])
surface.color=(20,148,161,166);session.models.add([surface]);mesh.close()
if review['context_pdb']:
    context=run(session,'open '+json.dumps(str(root/review['context_pdb'])))[0]
    context.name='crevice_review_context';context.crevice_region_owner=owner
    run(session,'hide #'+context.id_string+' cartoons; show #'+context.id_string+' atoms; style #'+context.id_string+' stick')
    context.atoms.colors=np.array([82,143,102,255],np.uint8)
points=MarkerSet(session,'crevice_review_points');points.crevice_region_owner=owner
for row in review['coordinate_landmarks']:
    for p in row['points_A']:points.create_marker(p,(153,87,184,255),.25)
session.models.add([points])
run(session,'set bgColor white; camera ortho; lighting soft; material dull; graphics quality 3; graphics silhouettes false')
R=np.asarray(review['camera']['camera_rotation']);center=np.asarray(review['camera']['center']);span=review['camera']['vertical_span']
matrix=np.column_stack((R.T,center+R[2]*review['camera']['radius']*3))
session.main_view.camera.position=Place(matrix=matrix);session.main_view.camera.field_width=span*4/3
from chimerax.label.label2d import session_labels,label_delete
labels=session_labels(session,create=False)
for name in ['crevice_region_title','crevice_region_legend']:
    if labels and labels.named_label(name):label_delete(session,[labels.named_label(name)])
__TITLE_LABELS__run(session,'windowsize 1200 900')
'''.replace('__ATOM_LABELS__',CXC_REVIEW_ATOM_LABELS if annotate else '').replace('__TITLE_LABELS__',CXC_REVIEW_TITLE_LABELS if annotate else '').replace('__SCENE__',repr(scene_path.name)).replace('__IDS__',repr(ids)).replace('__MESH__',repr(Path(files['pymol_mesh_npz']).name)).replace('__BASIS__',repr(basis.tolist()))
    loader=root/(prefix+'_review_chimerax.py');loader.write_text(cxc_py);files['review_chimerax_python']=str(loader)
    target=root/(prefix+'_review.cxc');target.write_text('open '+json.dumps(loader.name)+'\n');files['review_cxc']=str(target)
    import csv
    target=root/(prefix+'_review_landmarks.csv')
    with target.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['residue','anchor_atoms','nearest_reference_sample_A','centroid_x_A','centroid_y_A','centroid_z_A'])
        writer.writeheader()
        for r in scene['landmarks']:
            writer.writerow({'residue':r['residue'],'anchor_atoms':';'.join(r['atoms']),
                             'nearest_reference_sample_A':r['minimum_distance_to_reference_samples_A'],
                             **{f'centroid_{axis}_A':value for axis,value in zip('xyz',r['centroid_A'])}})
    files['review_landmarks_csv']=str(target)
    return files


# Annotation-only review text: residue names, caption and legend line.
PML_REVIEW_LABELS = '''for row in review['landmarks']:
    cmd.pseudoatom('crevice_review_labels',pos=row['centroid_A'],label=row['residue'])
R=np.asarray(review['camera']['camera_rotation']);center=np.asarray(review['camera']['center']);span=review['camera']['vertical_span']
for y,text in [(.45,review['caption']),(-.45,review['legend'])]:
    cmd.pseudoatom('crevice_review_labels',pos=(center+np.array([0,y,0])*span@R).tolist(),label=text)
cmd.hide('everything','crevice_review_labels');cmd.show('labels','crevice_review_labels')
cmd.set('label_color','black','crevice_review_labels');cmd.set('label_outline_color','white','crevice_review_labels');cmd.set('label_size',-span*.015,'crevice_review_labels');cmd.set('label_font_id',7,'crevice_review_labels')
'''
CXC_REVIEW_ATOM_LABELS = '''for row in review['landmarks']:
    serial=row['anchor_atom_serials'][0]
    run(session,'label #'+protein.id_string+'@@serial_number='+str(serial)+' atoms text '+json.dumps(row['residue'])+' height 0.8 color black bgColor white')
'''
CXC_REVIEW_TITLE_LABELS = '''run(session,'2dlabels create crevice_region_title text '+json.dumps(review['caption'])+' xpos 0.12 ypos 0.95 size 17 color black')
run(session,'2dlabels create crevice_region_legend text '+json.dumps(review['legend'])+' xpos 0.06 ypos 0.03 size 14 color black')
'''

