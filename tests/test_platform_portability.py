"""Behaviour that differs between Linux, macOS and Windows, checked on any platform.

Software tests only: output-stream encodings narrower than UTF-8 (Windows
redirects output in cp1252), UTF-8 text files whatever the locale, and
relative paths that cannot exist (Windows files on another drive).
"""
import os
import subprocess
import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parents[1] / "src")


def test_help_with_symbols_survives_a_cp1252_output_stream():
    env = {**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONPATH": SRC + os.pathsep + os.environ.get("PYTHONPATH", "")}
    done = subprocess.run([sys.executable, "-m", "crevice", "trajectory", "--help"], env=env,
                          capture_output=True, timeout=120)
    assert done.returncode == 0, done.stderr.decode("cp1252", "replace")
    assert b"\\u2264 0.75" in done.stdout


def test_radius_file_is_read_as_utf8_whatever_the_locale(tmp_path):
    from crevice.radii import load_radius_file
    path = tmp_path / "radii.csv"
    path.write_bytes("# radii in Å – ≤ test\nkind,key,radius\nelement,C,1.9\n".encode("utf-8"))
    assert load_radius_file(path).element_radii["C"] == 1.9


def test_region_file_paths_fall_back_to_absolute_without_a_relative_path(tmp_path, monkeypatch):
    from crevice import region_definition
    target = tmp_path / "data" / "map.dx"
    target.parent.mkdir()
    target.write_bytes(b"x")
    assert region_definition.file_record(target, tmp_path.resolve())["path"] == "data/map.dx"

    def other_drive(path, start):
        raise ValueError("path is on mount 'D:', start on mount 'C:'")
    monkeypatch.setattr(os.path, "relpath", other_drive)
    assert region_definition.file_record(target, tmp_path.resolve())["path"] == str(target.resolve())
