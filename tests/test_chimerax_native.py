"""Optional native positive and negative controls for oriented ChimeraX maps."""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
from pathlib import Path
import json,shutil,subprocess,sys
import pytest
from crevice import pore_profile,void_cast
from crevice.presentation import write_viewer_structure
from crevice.volume_export import write_volume_viewer_bundle
from crevice.residue_evidence import read_binary_dx,boundary_residue_evidence,write_residue_evidence_bundle
from test_channel_coordinate import rectangular_channel,transform
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from workspace_paths import local_chimerax  # noqa: E402


@pytest.mark.parametrize('split_regions',[False,True])
def test_native_chimerax_rotated_relocated_map_and_residue_context(tmp_path,split_regions):
    root=Path(__file__).resolve().parents[1];executable=local_chimerax(root)
    if not executable.is_file() or shutil.which('xvfb-run') is None:pytest.skip('Local ChimeraX runtime and Xvfb required')
    frame,_,_=transform(rectangular_channel());profile=pore_profile(frame, enclosure_radius=0.8,search_radius=8,samples=21)
    cast=void_cast(frame,profile=profile,spacing=.5,min_radius=0,max_grid_points=200000)
    folder=tmp_path/'initial';folder.mkdir();structure=write_viewer_structure(frame,folder/'viewer.pdb')
    points=cast.points;cut=(min(p.t for p in points)+max(p.t for p in points))/2
    regions={'lower':tuple(p for p in points if p.t<=cut),'upper':tuple(p for p in points if p.t>cut)} if split_regions else None
    files=write_volume_viewer_bundle(cast,structure_path=structure,output_dir=folder,frame=frame,profile=profile,regions=regions,cast_extension=0)
    g,o,d=read_binary_dx(files['volume_dx']);evidence=boundary_residue_evidence(frame,g,o,d)
    write_residue_evidence_bundle(evidence,frame,folder,scene_path=files['volume_pml'],dpi=60)
    dest=tmp_path/'relocated scene';shutil.move(folder,dest);scene=dest/'crevice_residue_context.cxc'
    command=[sys.executable,str(root/'scripts/check_chimerax.py'),str(scene),'--xvfb','--width','500','--output',str(tmp_path/'native')]
    done=subprocess.run(command,capture_output=True,text=True,timeout=200)
    assert done.returncode==0,done.stdout+done.stderr
    result=json.loads((tmp_path/'native/chimerax_verification.json').read_text())
    assert result['residue_sticks_passed'] and result['surface_reference_error_A']<.01
    if not split_regions:
        lines=scene.read_text().splitlines();lines=['view matrix models #2,1,0,0,0,0,1,0,0,0,0,1,0' if line.startswith('view matrix models #2,') else line for line in lines]
        scene.write_text('\n'.join(lines)+'\n');command[-1]=str(tmp_path/'rejected')
        done=subprocess.run(command,capture_output=True,text=True,timeout=200)
        assert done.returncode!=0
        result=json.loads((tmp_path/'rejected/chimerax_verification.json').read_text())
        assert result['native_execution_passed'] and result['surface_reference_error_A']>result['surface_reference_bound_A']
