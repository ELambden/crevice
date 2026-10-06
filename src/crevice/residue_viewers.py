"""Matched residue-stick scenes for PyMOL, VMD and ChimeraX."""
from pathlib import Path
import json
from .radii import RadiusSet, radii_option


@radii_option
def scene_from_binary_volume(frame, volume_path, output_dir, prefix, *, radii: RadiusSet | str | None = None):
    """Write volume scenes from an existing binary DX without re-detecting a cavity.

    Every occupied sample of the measured grid becomes a cast sample (with its
    atom-surface clearance); face-connected groups become separate display
    regions when there is more than one. The scenes are those of
    :func:`crevice.volume_export.write_volume_viewer_bundle` with no display
    extension. No text is drawn.

    Parameters
    ----------
    frame : StructureFrame
        Structure in the coordinate frame of the DX.
    volume_path : str or Path
        Binary DX (1 inside, 0 outside).
    output_dir : str or Path
        Output directory.
    prefix : str
        File-name prefix.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict of str to str
        Output keys and paths of the volume scene bundle.

    Outputs
    -------
    PREFIX_viewer.pdb and PREFIX_volume.* scene files : scenes
        See :func:`crevice.volume_export.write_volume_viewer_bundle`.

    Colours and representation
    --------------------------
    Grey-blue 25%-transparent cartoon; opaque cast surface, or one colour per
    connected region (cast colour, purple, blue, mauve). The cast colour is the
    one recorded for the map (``cast_color_rgb`` in the sibling
    ``*_metadata.json`` of a CREVICE ``PREFIX_volume.dx``; teal when that file
    records channel mouths); otherwise the map is treated as a non-channel cast
    and drawn violet (:data:`crevice.presentation.NON_CHANNEL_CAST_RGB`).
    """
    import numpy as np
    from scipy import ndimage
    from .residue_evidence import read_binary_dx
    from .models import ChannelPoint,VoidComponent,VoidCast
    from .presentation import write_viewer_structure
    from .volume_export import write_volume_viewer_bundle
    from .rolling import _SphereQueries
    from .radii import atom_vdw_radius
    grid,origin,deltas=read_binary_dx(volume_path)
    spacing=float(np.linalg.norm(deltas[0]));basis=deltas/spacing
    if not np.allclose(basis@basis.T,np.eye(3),atol=1e-8):
        raise ValueError('Automatic residue scenes require an orthogonal, equally spaced measured map')
    labels,n=ndimage.label(grid);components=[]
    atoms=frame.selected_atoms(include_hydrogen=False,include_hetero=False)
    query=_SphereQueries(np.asarray([a.coord for a in atoms]),np.asarray([atom_vdw_radius(a) for a in atoms]))
    for i in range(1,n+1):
        xyz=np.argwhere(labels==i)@deltas+origin;clearance,nearest=query.points(xyz)
        points=tuple(ChannelPoint(j,tuple(map(float,p)),max(0.,float(r)),float(r),float(p@basis[2]),
                                  atoms[int(k)].serial,atoms[int(k)].residue_key.label)
                     for j,(p,r,k) in enumerate(zip(xyz,clearance,nearest)))
        components.append(VoidComponent(i,'imported_binary_region',tuple(map(float,xyz.mean(0))),len(points)*spacing**3,points))
    if not components:raise ValueError('Cannot generate a scene from an empty measured map')
    points=tuple(p for c in components for p in c.points)
    cast=VoidCast(tuple(components),points,spacing,0.,'rolling',{
        'algorithm':'imported_measured_binary_samples','grid_basis':basis.tolist(),'grid_origin':origin.tolist(),
        'axis_origin':np.asarray([a.coord for a in frame.atoms]).mean(0).tolist(),'axis_direction':basis[2].tolist(),
        'export_mode':'imported_binary_samples','source_volume':str(Path(volume_path).resolve()),
        'connectivity':'face-connected sampled map components; no biological channel assignment'})
    root=Path(output_dir);structure=write_viewer_structure(frame,root/f'{prefix}_viewer.pdb')
    regions={f'region_{c.id}':c.points for c in components} if len(components)>1 else None
    color=None
    sidecar=Path(volume_path).with_name(Path(volume_path).stem+'_metadata.json')
    if sidecar.is_file():
        try:
            recorded=json.loads(sidecar.read_text())
        except ValueError:
            recorded={}
        from .presentation import CHANNEL_CAST_RGB
        color=recorded.get('cast_color_rgb') or (CHANNEL_CAST_RGB if recorded.get('channel_mouths') else None)
    return write_volume_viewer_bundle(cast,structure_path=structure,output_dir=root,prefix=prefix,
                                     frame=frame,cast_extension=0,regions=regions,cast_color=color)


def write_companion_contexts(scene, root, prefix, protein_path, membership, viewer_ids, frame,
                             *, cast_opacity=.55):
    """Write VMD and ChimeraX residue-evidence overlays on an existing volume scene.

    The base scene's own ``.tcl`` and ``.cxc`` are reused unchanged; the VMD
    overlay first checks every protein atom serial and coordinate against the
    evidence reference (0.001 Å), then adds opaque stick representations for
    the two evidence roles and makes the cast translucent. No labels or text
    are added.

    Parameters
    ----------
    scene : str or Path
        Base volume scene (``*_volume.pml``); its ``.tcl`` and ``.cxc`` siblings
        must exist.
    root : str or Path
        Output directory.
    prefix : str
        File-name prefix.
    protein_path : Path
        Viewer PDB of the protein.
    membership : dict of str to list of int
        Viewer atom serials of ``lining`` and ``partners`` residues.
    viewer_ids : sequence of int
        Viewer serial of each frame atom, in frame order.
    frame : StructureFrame
        Reference coordinates for the check.
    cast_opacity : float, default 0.55
        Cast surface opacity in the overlays.

    Returns
    -------
    dict of str to str
        ``residue_context_{vmd,tcl,cxc}`` and ``residue_reference_tsv`` paths.

    Raises
    ------
    ValueError
        If the base ``.tcl`` or ``.cxc`` scene is missing.

    Outputs
    -------
    PREFIX_residue_context.vmd, PREFIX_residue_context.tcl, PREFIX_residue_context.cxc : scenes
        Overlays to open from the output directory.
    PREFIX_residue_reference.tsv : table
        Viewer serial and reference coordinates (Å) checked by VMD.

    Colours and representation
    --------------------------
    Gold boundary residues (``[0.85, 0.55, 0.20]``; ChimeraX ``#d98c33``) and
    magenta nonlocal partners (``[0.64, 0.27, 0.47]``; ``#a34478``) as opaque
    0.19 Å sticks over the base cartoon; the cast 55% opaque.
    """
    import re
    from .volume_export import _tcl_word
    scene=Path(scene);root=Path(root)
    tcl=scene.with_suffix('.tcl');cxc=scene.with_suffix('.cxc')
    if not tcl.is_file() or not cxc.is_file():
        raise ValueError('Residue contexts require companion CREVICE .tcl and .cxc volume scenes')
    coordinates=root/f'{prefix}_residue_reference.tsv'
    coordinates.write_text(''.join(f'{serial}\t{a.x:.9f}\t{a.y:.9f}\t{a.z:.9f}\n' for a,serial in zip(frame.atoms,viewer_ids)))
    definitions='\n'.join(f'set crevice_evidence_serials({role}) {{{" ".join(map(str,ids))}}}' for role,ids in membership.items())
    extra=f'''
# Matched highlighted residues: opaque sticks, transparent cavity context.
{definitions}
set crevice_reference_file [open [file join $crevice_root {_tcl_word(coordinates.name)}] r]
set crevice_expected_coordinates [dict create]
foreach crevice_line [split [read $crevice_reference_file] "\\n"] {{
    if {{[llength $crevice_line] == 4}} {{dict set crevice_expected_coordinates [lindex $crevice_line 0] [lrange $crevice_line 1 3]}}
}}
close $crevice_reference_file
set crevice_all_atoms [atomselect $crevice_protein all]
if {{[$crevice_all_atoms num] != [dict size $crevice_expected_coordinates]}} {{error "Residue reference atom count mismatch"}}
foreach crevice_serial [$crevice_all_atoms get serial] crevice_xyz [$crevice_all_atoms get {{x y z}}] {{
    if {{![dict exists $crevice_expected_coordinates $crevice_serial]}} {{error "Residue reference atom serial mismatch"}}
    foreach crevice_a $crevice_xyz crevice_b [dict get $crevice_expected_coordinates $crevice_serial] {{
        if {{abs($crevice_a-$crevice_b)>0.001}} {{error "Residue reference coordinates differ"}}
    }}
}}
$crevice_all_atoms delete
color change rgb 30 0.85 0.55 0.20
color change rgb 31 0.64 0.27 0.47
if {{[lsearch -exact [material list] CREVICEResidues] < 0}} {{material add CREVICEResidues}}
material change opacity CREVICEResidues 1.0
material change diffuse CREVICEResidues 0.7
material change specular CREVICEResidues 0.12
foreach crevice_role {{lining partners}} crevice_color {{30 31}} {{
    if {{[llength $crevice_evidence_serials($crevice_role)]}} {{
        mol representation Licorice 0.19 24 24
        mol selection "serial [join $crevice_evidence_serials($crevice_role) {{ }}]"
        mol color ColorID $crevice_color
        mol material CREVICEResidues
        mol addrep $crevice_protein
    }}
}}
material change opacity CREVICECast {cast_opacity:g}
# No wireframe, labels or helper-grid representations are introduced.
display update
'''
    text=tcl.read_text()+extra
    paths={}
    for extension in ['vmd','tcl']:
        p=root/f'{prefix}_residue_context.{extension}';p.write_text(text);paths[f'residue_context_{extension}']=str(p)
    # CXC interprets file names relative to the command script's own directory.
    commands=cxc.read_text().splitlines()
    commands[0]='open '+json.dumps(protein_path.name)+' id #1'
    volume_ids=[]
    for line in commands:
        match=re.match(r'volume (#\d+) style surface',line)
        if match:volume_ids.append(match[1])
    commands += ['# Evidence categories use original IDs in the companion report; selection uses verified viewer serials.']
    viewer_residues={}
    for line in protein_path.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            viewer_residues[int(line[6:11])]=(line[21:22].strip(),line[22:27].strip())
    for role,color in [('partners','#a34478'),('lining','#d98c33')]:
        ids=membership[role]
        if ids:
            residues=sorted({viewer_residues[i] for i in ids})
            spec='|'.join('#1'+('/'+chain if chain else '')+':'+resid for chain,resid in residues)
            commands += [f'show {spec} atoms',f'style {spec} stick',f'size {spec} stickRadius 0.19',
                         f'color {spec} {color} target a',f'transparency {spec} 0 target a']
    commands += [f'transparency {model} {100*(1-cast_opacity):g} target s' for model in volume_ids]
    commands += ['view', 'zoom 0.9']
    p=root/f'{prefix}_residue_context.cxc';p.write_text('\n'.join(commands)+'\n');paths['residue_context_cxc']=str(p)
    paths['residue_reference_tsv']=str(coordinates)
    return paths
