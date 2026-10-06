"""Software/synthetic geometry contracts; no biological channel assignment."""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
from dataclasses import replace
import json
import math
from pathlib import Path
import numpy as np
import pytest
from crevice import Atom, StructureFrame
from crevice.rolling import rolling_probe_cast, select_cast_atoms, _SphereQueries, _connected
from crevice.spatial import SpatialIndex
from test_channel_coordinate import rectangular_channel, transform


def test_analytic_edges_cannot_jump_through_atom():
    coords=np.array([[1.,0,0]])
    queries=_SphereQueries(coords,np.array([.5]))
    mask=np.ones((2,1,1),dtype=bool)
    clearance=np.full(mask.shape,.5)
    labels=_connected(mask,clearance,np.zeros(3),2.,.1,queries)
    assert labels[0,0,0] != labels[1,0,0]
    assert queries.segments(np.array([[1.,0,0],[1.,2.,0]]),0,1.,.1).tolist()==[False,True]


def test_batched_clearance_handles_unequal_radii():
    atoms=(Atom(1,'C','ALA','A',1,0,0,0,'C'),Atom(2,'O','ALA','A',2,3,1,0,'O'))
    from crevice.radii import atom_vdw_radius
    queries=_SphereQueries(np.array([a.coord for a in atoms]),np.array([atom_vdw_radius(a) for a in atoms]))
    points=np.random.default_rng(17).normal(size=(50,3))*4
    expected=SpatialIndex(atoms)
    values,_=queries.points(points)
    assert values==pytest.approx([expected.nearest_surface(p)[0] for p in points])
    for axis in range(3):
        delta=np.eye(3)[axis]*.35
        actual=queries.segments(points,axis,.35,.8)
        assert actual.tolist()==[expected.segment_clearance(p-delta,p+delta)>=.8 for p in points]


@pytest.fixture(scope='module')
def pocket():
    frame=rectangular_channel(half_x=4,half_y=5,half_length=8,capped=True)
    cast=rolling_probe_cast(frame,spacing=.5,outer_radius=6,min_component_volume=30)
    return frame,cast


def test_one_ended_pocket_and_open_channel_are_retained(pocket):
    frame,cast=pocket
    main=cast.components[0]
    assert main.volume>100
    assert any(math.dist(p.position,(0,0,0))<1 for p in main.points)
    assert min(p.position[2] for p in main.points)<-5
    assert max(p.position[2] for p in main.points)<8
    assert cast.metadata['measurement_point_count']==sum(c.point_count for c in cast.components)
    assert cast.total_volume==sum(c.point_count for c in cast.components)*.5**3
    assert all(p.raw_clearance>=0 for c in cast.components for p in c.points)
    json.dumps(cast.to_dict())
    open_cast=rolling_probe_cast(rectangular_channel(half_x=4,half_y=5,half_length=8),spacing=.5,outer_radius=6)
    assert any(p.position[2]>5 for p in open_cast.components[0].points)
    assert any(p.position[2]<-5 for p in open_cast.components[0].points)


def test_closed_cavity_is_retained_without_bulk_access():
    frame=rectangular_channel(half_x=4,half_y=5,half_length=6,capped=True)
    extra=tuple(Atom(10000+i,'C','ALA','A',10000+i,x,y,-6,'C')
                for i,(x,y) in enumerate((x,y) for x in range(-4,5) for y in range(-5,6)))
    cast=rolling_probe_cast(replace(frame,atoms=frame.atoms+extra),spacing=.5,outer_radius=5)
    assert any(math.dist(p.position,(0,0,0))<1 for c in cast.components for p in c.points)


def test_rigid_transform_preserves_entire_cast(pocket):
    frame,cast=pocket
    moved,rotation,shift=transform(frame)
    other=rolling_probe_cast(moved,spacing=.5,outer_radius=6,min_component_volume=30)
    assert other.total_volume==cast.total_volume
    from scipy.spatial import cKDTree
    points=np.array([p.position for c in cast.components for p in c.points])@rotation.T+shift
    actual=np.array([p.position for c in other.components for p in c.points])
    assert cKDTree(points).query(actual)[0].max()<1e-7


def test_three_lumens_and_interfacial_void_are_retained():
    centers=[(9.6*math.cos(a),9.6*math.sin(a)) for a in (0,2*math.pi/3,4*math.pi/3)]
    atoms=[]
    for chain,(cx,cy) in zip('ABC',centers):
        for z in range(-8,9):
            for angle in np.linspace(0,2*math.pi,32,endpoint=False):
                atoms.append(Atom(len(atoms)+1,'C','ALA',chain,len(atoms)+1,cx+4.5*math.cos(angle),cy+4.5*math.sin(angle),z,'C'))
    frame=StructureFrame(tuple(atoms))
    cast=rolling_probe_cast(frame,spacing=.6,outer_radius=6,min_component_volume=30)
    assigned=[]
    for x,y in centers+[(0,0)]:
        hits=[c.id for c in cast.components if any(math.dist(p.position,(x,y,0))<1 for p in c.points)]
        assert len(hits)==1
        assigned.extend(hits)
    assert len(set(assigned))==4
    # The original fixture has wide lateral gaps (only 18 blocked rays at
    # its interface). Narrow those gaps for the default 24-ray core policy.
    by_chain=dict(zip('ABC',centers))
    enclosed=replace(frame,atoms=tuple(replace(a,
        x=by_chain[a.chain_id][0]+(a.x-by_chain[a.chain_id][0])*4/3,
        y=by_chain[a.chain_id][1]+(a.y-by_chain[a.chain_id][1])*4/3) for a in frame.atoms))
    dominant=rolling_probe_cast(enclosed,spacing=.5,outer_radius=6,min_component_volume=30,
                                enclosure_fraction=.9,selection_mode='dominant',max_components=8)
    assert len(dominant.components)==4
    for x,y in centers+[(0,0)]:
        assert sum(any(math.dist(p.position,(x,y,0))<1 for p in c.points)
                   for c in dominant.components)==1
    from crevice.rolling import chain_enclosed_regions
    before=cast.to_dict()
    regions,summary=chain_enclosed_regions(frame,cast)
    assert {'chain_A','chain_B','chain_C','shared'} <= regions.keys()
    assert sum(len(points) for points in regions.values())==sum(c.point_count for c in cast.components)
    assert sum(row['volume_A3'] for row in summary)==pytest.approx(cast.total_volume)
    assert cast.to_dict()==before
    assert len({p.position for points in regions.values() for p in points})==sum(len(points) for points in regions.values())


def test_explicit_selection_and_resolution_limits(pocket):
    frame,_=pocket
    water=Atom(99999,'O','HOH','W',1,0,0,0,'O',hetero=True)
    selected,report=select_cast_atoms(replace(frame,atoms=frame.atoms+(water,)))
    assert selected.atoms==frame.atoms
    assert report['excluded_residue_atom_counts']=={'HOH':1}
    assert select_cast_atoms(replace(frame,atoms=frame.atoms+(water,)),'all')[0].atoms[-1]==water
    with pytest.raises(ValueError,match='explicitly choose a coarser spacing'):
        rolling_probe_cast(frame,spacing=.1,max_grid_points=1000)
    with pytest.raises(ValueError,match='outer_radius must exceed'):
        rolling_probe_cast(frame,probe_radius=2,outer_radius=1)
    with pytest.raises(ValueError,match='min_depth'):
        rolling_probe_cast(frame,min_depth=-1)


def test_directional_enclosure_allows_two_channel_exits():
    from crevice.rolling import _ray_enclosure
    blocked=np.zeros((9,9,15),dtype=bool)
    blocked[0,:,:]=blocked[-1,:,:]=blocked[:,0,:]=blocked[:,-1,:]=True
    assert _ray_enclosure(blocked)[4,4,7]==24
    blocked[:,:,0]=True
    assert _ray_enclosure(blocked)[4,4,7]==25


def test_cli_cast_does_not_require_a_through_channel(pocket,tmp_path,capsys):
    from crevice.cli import main
    from crevice.presentation import write_viewer_structure
    frame,_=pocket
    source=write_viewer_structure(frame,tmp_path/'input.pdb')
    main(['cast',str(source),'--out-dir',str(tmp_path/'bundle'),'--spacing','0.5','--min-component-volume','30'])
    manifest=json.loads((tmp_path/'bundle/input_manifest.json').read_text())
    assert manifest['void_cast']['total_volume']>100
    assert manifest['channel_profile_status']=='not_assigned_by_3D_void_detection'
    assert 'profile_json' not in manifest['files']
    assert Path(manifest['files']['pore_cast_pml']).exists()
    assert Path(manifest['files']['end_view_tcl']).exists()


def test_publication_auto_falls_back_without_inventing_a_profile(pocket,tmp_path):
    from crevice.figures import write_static_publication_bundle
    from crevice.presentation import write_viewer_structure
    frame,_=pocket
    source=write_viewer_structure(frame,tmp_path/'pocket.pdb')
    with pytest.warns(RuntimeWarning,match='No unique through-channel'):
        files=write_static_publication_bundle(frame, enclosure_radius=0.8,structure_path=source,output_dir=tmp_path/'published',
            cast_spacing=.5,cast_max_grid_points=500000,cast_min_component_volume=30,
            include_hetero=False,include_cavities=False,include_tunnels=False,include_network=False,dpi=80)
    manifest=json.loads(Path(files['manifest_json']).read_text())
    assert manifest['single_channel_unresolved_reason']
    assert manifest['void_cast']['mode']=='rolling'
    assert manifest['void_cast']['total_volume']>100
    assert 'profile_json' not in files


def test_dominant_limit_does_not_split_two_chain_single_pore(pocket):
    frame,_=pocket
    frame=replace(frame,atoms=tuple(replace(a,chain_id='A' if a.z<0 else 'B') for a in frame.atoms))
    casts=[rolling_probe_cast(frame,spacing=.5,outer_radius=6,min_component_volume=30,
            enclosure_fraction=.9,selection_mode='dominant',max_components=n) for n in (1,8)]
    assert [len(c.components) for c in casts]==[1,1]
    assert casts[0].components==casts[1].components
    from crevice.volume_export import cast_grid
    from scipy.ndimage import label
    assert label(cast_grid(casts[0])[0]>.5)[1]==1
    from crevice.spatial import SpatialIndex
    spatial=SpatialIndex(frame.atoms)
    assert min(spatial.nearest_surface(p.position)[0] for p in casts[0].components[0].points)>=-1e-9


def test_variable_sphere_union_fills_between_qualified_centers():
    from crevice.cavity_selection import _paint_balls
    scores=np.full((17,17,17),-np.inf,dtype=np.float32)
    labels=np.zeros(scores.shape,dtype=np.int32)
    centers=np.array([[6,8,8],[10,8,8]])
    _paint_balls(centers,np.array([2.,2.]),.5,scores,labels,1)
    grid=np.indices(scores.shape).transpose(1,2,3,0)
    expected=np.any(np.sum(((grid[:,:,:,None,:]-centers)*.5)**2,axis=-1)<=4,axis=-1)
    assert np.array_equal(labels>0,expected)
    assert labels[8,8,8]==1
    # A point farther than a fixed .8 A probe, yet inside an atom-clear ball.
    assert labels[6,11,8]==1


def test_new_cast_cli_defaults_and_explicit_enumeration():
    from crevice.cli import build_parser
    parser=build_parser()
    args=parser.parse_args(['cast','input.pdb','--out-dir','out'])
    assert args.max_components==1 and args.cavity_selection=='dominant'
    assert args.surface_smoothing is None and args.smooth is None  # --smooth default applies
    args=parser.parse_args(['cast','input.pdb','--out-dir','out','--max-cavities','8','--cavity-selection','all'])
    assert args.max_components==8 and args.cavity_selection=='all'
    args=parser.parse_args(['publish','input.pdb','--out-dir','out','--cast-max-cavities','4'])
    assert args.cast_max_components==4 and args.cast_selection=='dominant'


def test_dominant_rejects_invalid_selection_settings(pocket):
    frame,_=pocket
    with pytest.raises(ValueError,match='selection_mode'):
        rolling_probe_cast(frame,selection_mode='invent')
    with pytest.raises(ValueError,match='core_fraction'):
        rolling_probe_cast(frame,selection_mode='dominant',core_fraction=0)


def test_dominant_fill_is_rigid_transform_equivariant(pocket):
    frame,_=pocket
    moved,rotation,shift=transform(frame)
    settings=dict(spacing=.5,outer_radius=6,min_component_volume=30,
                  enclosure_fraction=.9,selection_mode='dominant',max_components=8)
    original=rolling_probe_cast(frame,**settings)
    rotated=rolling_probe_cast(moved,**settings)
    assert [c.volume for c in original.components]==[c.volume for c in rotated.components]
    from scipy.spatial import cKDTree
    expected=np.asarray([p.position for c in original.components for p in c.points])@rotation.T+shift
    actual=np.asarray([p.position for c in rotated.components for p in c.points])
    assert cKDTree(expected).query(actual)[0].max()<1e-8


def test_dominant_empty_input_geometry_returns_no_invented_cavity():
    frame=StructureFrame((Atom(1,'C','ALA','A',1,0,0,0,'C'),))
    cast=rolling_probe_cast(frame,spacing=.5,outer_radius=6,enclosure_fraction=.9,
                            selection_mode='dominant',max_components=8)
    assert cast.components==() and cast.total_volume==0
    assert cast.metadata['dominance']['qualified_count']==0
    assert cast.metadata['dominance']['reason']
    json.dumps(cast.to_dict())
