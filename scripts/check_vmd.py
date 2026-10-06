"""Run an exported VMD scene headlessly and retain a native render and log.

Sources a cast or residue-evidence scene (``PREFIX_volume.vmd``/``.tcl`` or
``PREFIX_residue_context.vmd``) in VMD's text mode, renders it with
TachyonInternal and checks that the image has content, cast pixels (teal channel or violet non-channel) and a
clean border; for an evidence overlay it also compares VMD's representations
(stick residues and colours) with ``PREFIX_residue_evidence.json``. Requires
``vmd`` on ``PATH``.

Usage::

    python scripts/check_vmd.py OUT/PREFIX_volume.vmd --output check/

Writes ``vmd_native.tga``, ``vmd_native.png``, ``vmd.log`` and
``vmd_verification.json``; exits with status 1 unless the checks pass.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from crevice.volume_export import _tcl_word


def main():
    """Parse the command line, run VMD and write ``vmd_verification.json``.

    Options: ``scene``, ``--output`` (directory), ``--rotation RX RY RZ``
    (inspection rotations in degrees), ``--timeout`` (seconds, default 360) and
    ``--width`` (pixels, at least 100; default 1200).

    Raises
    ------
    SystemExit
        When VMD is missing, or with status 1 when a check fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scene", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rotation", type=float, nargs=3, default=None,
                        help="Inspection-view rotations around x, y and z in degrees")
    parser.add_argument("--timeout", type=float, default=360, help="Native-render timeout in seconds")
    parser.add_argument("--width", type=int, default=1200, help="Native inspection image width")
    args = parser.parse_args()
    if args.width<100:parser.error("width must be at least 100")
    executable = shutil.which("vmd")
    if executable is None:
        raise SystemExit("VMD is not installed")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    scene = args.scene.resolve(strict=True)
    image = root / "vmd_native.tga"
    rotation = ""
    if args.rotation is not None:
        rotation = "\nset crevice_rotation [transmult " + " ".join(
        f"[transaxis {axis} {angle}]" for axis, angle in zip("xyz", args.rotation)) + "]"
        rotation += "\nforeach crevice_molecule [molinfo list] {molinfo $crevice_molecule set rotate_matrix [list $crevice_rotation]}"
    overlay_path=scene.parent/(scene.stem.removesuffix('_residue_context')+'_residue_evidence.json')
    overlay=json.loads(overlay_path.read_text()) if scene.stem.endswith('_residue_context') and overlay_path.is_file() else None
    timed_out = False
    with tempfile.TemporaryDirectory(prefix="crevice-vmd-") as tmp:
        fresh_image = Path(tmp) / "render.tga"
        native_reps=Path(tmp)/'representations.tsv'
        inspection=''
        if overlay is not None:
            inspection = '\nset crevice_inspection [open CREVICE_NATIVE_REPS w]\nfor {set crevice_rep 0} {$crevice_rep < [molinfo $crevice_protein get numreps]} {incr crevice_rep} {\n    set crevice_style [lindex [molinfo $crevice_protein get [list [list rep $crevice_rep]]] 0]\n    set crevice_color [lindex [molinfo $crevice_protein get [list [list color $crevice_rep]]] 0]\n    set crevice_seltext [lindex [molinfo $crevice_protein get [list [list selection $crevice_rep]]] 0]\n    set crevice_sel [atomselect $crevice_protein $crevice_seltext]\n    puts $crevice_inspection "$crevice_style|$crevice_color|[lsort -integer [$crevice_sel get serial]]"\n    $crevice_sel delete\n}\nclose $crevice_inspection\n'.replace('CREVICE_NATIVE_REPS',_tcl_word(str(native_reps)))
        on_error = 'crevice_error]} {puts stderr "CREVICE_RENDER_ERROR: $crevice_error"; exit 1}\n'
        wrapper = (
            "if {[catch {source " + _tcl_word(str(scene)) + "} " + on_error +
            inspection + rotation + f"\ndisplay resize {args.width} {round(args.width*.75)}\ndisplay redraw\ndisplay update\n" +
            "if {[catch {render TachyonInternal " + _tcl_word(str(fresh_image)) + "} " + on_error + "quit\n")
        path = Path(tmp) / "render.tcl"
        path.write_text(wrapper)
        log_path = Path(tmp) / "vmd.log"
        # Keep stdin open: text-mode VMD otherwise exits on EOF between queued
        # scene commands, before the view changes and rendering are complete.
        with log_path.open("w") as output:
            with subprocess.Popen([executable, "-dispdev", "text", "-e", str(path)],
                                  stdin=subprocess.PIPE, stdout=output,
                                  stderr=subprocess.STDOUT, text=True) as process:
                try:
                    returncode = process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    process.kill()
                    returncode = process.wait()
                    timed_out = True
        log = log_path.read_text()
        representations=native_reps.read_text() if native_reps.is_file() else ''
        generated = fresh_image.is_file()
        if generated:
            shutil.copyfile(fresh_image, image)
    (root / "vmd.log").write_text(log)
    ok = not timed_out and returncode == 0 and generated and "CREVICE_RENDER_ERROR" not in log
    report = {"status": "passed" if ok else ("timeout" if timed_out else "failed"),
              "timeout_seconds": args.timeout, "executable": executable,
              "scene": str(scene), "returncode": returncode, "native_image": str(image),
              "inspection_rotation_degrees": args.rotation}
    if overlay is not None:
        rows=[row.split('|') for row in representations.splitlines()]
        actual={}
        unsupported=[]
        for style,color,serials in rows:
            kind=style.split()[0]
            if kind=='Licorice':
                actual.setdefault(color.strip(),[]).extend(int(i) for i in serials.split())
            elif kind!='NewCartoon':unsupported.append(style)
        expected={f'ColorID {number}':sorted(overlay['display_atom_serials'][role])
                  for role,number in [('lining',30),('partners',31)] if overlay['display_atom_serials'][role]}
        actual={color:sorted(ids) for color,ids in actual.items()}
        report['residue_sticks_passed']=actual==expected and not unsupported
        report['native_stick_serials_by_color']=actual
        report['coordinate_mapping_verified_during_scene_load']=ok and bool(representations)
        report['protein_representations']=[row[0] for row in rows]
        ok=ok and report['residue_sticks_passed']
        report['status']='passed' if ok else ('timeout' if timed_out else 'failed_residue_representations')
        (root/'native_representations.tsv').write_text(representations)
    if ok:
        from PIL import Image
        import numpy as np
        with Image.open(image) as rendered:
            rgb = np.asarray(rendered.convert("RGB"))
            foreground_fraction = float(np.any(rgb < 245, axis=2).mean())
            teal = (rgb[:,:,1] > rgb[:,:,0]*1.3) & (rgb[:,:,2] > rgb[:,:,0]*1.3) & (rgb[:,:,0] < 180)
            # Non-channel casts are violet (crevice.presentation.NON_CHANNEL_CAST_RGB).
            violet = (rgb[:,:,2] > rgb[:,:,1]*1.3) & (rgb[:,:,0] > rgb[:,:,1]*1.2) & (rgb[:,:,1] < 180)
            teal = teal | violet
            border = np.concatenate((rgb[0],rgb[-1],rgb[:,0],rgb[:,-1]))
            border_fraction = float(np.any(border < 245,axis=1).mean())
            ok = bool(rgb.std() > 1 and foreground_fraction > 0.001 and teal.sum() >= 25
                      and border_fraction < .01)
            report.update(teal_pixels=int(teal.sum()), border_foreground_fraction=border_fraction)
            report["status"] = "passed" if ok else "failed_blank_image"
            preview = root / "vmd_native.png"
            rendered.save(preview)
            report.update(width=rendered.width, height=rendered.height, pixel_std=float(rgb.std()),
                          foreground_fraction=foreground_fraction,
                          preview_png=str(preview))
    (root / "vmd_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
