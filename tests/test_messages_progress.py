"""Software tests of user-facing error messages and stderr progress feedback.

All offline: downloads are replaced by functions that raise the errors urllib
would raise. These check message wording and that progress never reaches
standard output or output files; they are not real-system validation.
"""
import io
import json
import urllib.error
import urllib.request

import pytest

from crevice import examples, pdb, progress
from crevice.cli import UNRESOLVED_CHANNEL_HINT, main
from crevice.sections import ChannelResolutionError


def _http_404(url, *args, **kwargs):
    target = url.full_url if isinstance(url, urllib.request.Request) else url
    raise urllib.error.HTTPError(target, 404, "Not Found", {}, None)


def _offline(url, *args, **kwargs):
    raise urllib.error.URLError(OSError(-3, "Temporary failure in name resolution"))


@pytest.fixture
def one_file_dataset(monkeypatch):
    dataset = examples.ExampleDataset(
        name="fake", description="test data", licence="CC0",
        files=(examples.ExampleFile("t.xtc", "0" * 64, 1, "plain"),))
    monkeypatch.setitem(examples.EXAMPLES, "fake", dataset)
    for name in (examples.BASE_URL_ENV, examples.SOURCE_DIR_ENV, examples.CACHE_DIR_ENV):
        monkeypatch.delenv(name, raising=False)
    return dataset


# --- download messages -------------------------------------------------------

def test_example_404_names_the_file_url_and_release(tmp_path, one_file_dataset, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _http_404)
    with pytest.raises(RuntimeError) as info:
        examples.fetch_example("fake", cache_dir=tmp_path)
    message = str(info.value)
    assert "t.xtc" in message and examples.DEFAULT_BASE_URL + "t.xtc" in message
    assert "HTTP 404" in message and "not attached to the GitHub release v0.1.1" in message
    assert "--source-dir" in message
    assert not list((tmp_path / "fake").glob("*.part"))


def test_example_network_failure_suggests_offline_use(tmp_path, one_file_dataset, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", _offline)
    with pytest.raises(RuntimeError, match=r"t\.xtc.*name resolution.*offline.*--source-dir"):
        examples.fetch_example("fake", cache_dir=tmp_path)


def test_example_missing_from_source_dir_is_explained(tmp_path, one_file_dataset):
    with pytest.raises(RuntimeError, match=r"t\.xtc is not in the source directory"):
        examples.fetch_example("fake", cache_dir=tmp_path / "cache", source_dir=tmp_path)


def test_cli_fetch_example_404_is_one_error_line(tmp_path, one_file_dataset, monkeypatch, capsys):
    monkeypatch.setattr(urllib.request, "urlopen", _http_404)
    assert main(["fetch-example", "fake", "--cache-dir", str(tmp_path)]) == 2
    err = capsys.readouterr().err
    lines = err.strip().splitlines()
    assert len(lines) == 1 and lines[0].startswith("crevice: error: Example file t.xtc was not found")
    assert "Traceback" not in err


def test_rcsb_404_asks_to_check_the_identifier(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb, "urlopen", _http_404)
    with pytest.raises(RuntimeError) as info:
        pdb.fetch_structure("9ZZZ", tmp_path)
    message = str(info.value)
    assert "9ZZZ" in message and "https://files.rcsb.org/" in message
    assert "HTTP 404" in message and "check the identifier" in message


def test_alphafold_404_mentions_the_uniprot_accession(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb, "urlopen", _http_404)
    with pytest.raises(RuntimeError, match=r"HTTP 404.*UniProt accession"):
        pdb.fetch_alphafold_model("AF-P99999-F1", tmp_path)


def test_structure_network_failure_suggests_offline(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb, "urlopen", _offline)
    with pytest.raises(RuntimeError, match=r"Could not fetch 1GRM .*name resolution.*--offline"):
        pdb.fetch_structure("1GRM", tmp_path)


# --- CLI error handling ------------------------------------------------------

def test_key_error_message_is_not_quoted(capsys):
    assert main(["fetch-example", "no-such-dataset"]) == 2
    assert capsys.readouterr().err.startswith("crevice: error: Unknown example dataset 'no-such-dataset'")


def test_debug_prints_the_traceback(capsys, monkeypatch):
    monkeypatch.delenv("CREVICE_DEBUG", raising=False)
    assert main(["fetch-example", "no-such-dataset"]) == 2
    assert "Traceback" not in capsys.readouterr().err
    assert main(["fetch-example", "no-such-dataset", "--debug"]) == 2
    assert "Traceback" in capsys.readouterr().err
    monkeypatch.setenv("CREVICE_DEBUG", "1")
    assert main(["fetch-example", "no-such-dataset"]) == 2
    assert "Traceback" in capsys.readouterr().err


def test_unresolved_channel_error_gets_a_hint(monkeypatch, capsys):
    import crevice.cli as cli

    def unresolved(args):
        raise ChannelResolutionError("No enclosure probe resolved a stable channel (automatic choice).")

    monkeypatch.setattr(cli, "_cmd_profile", unresolved)
    assert main(["profile", "1GRM", "-o", "unused.json"]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert err[0].startswith("crevice: error: No enclosure probe resolved")
    assert err[1] == f"crevice: hint: {UNRESOLVED_CHANNEL_HINT}"
    for option in ("--axis", "--lateral-exits", "--enclosure-radius", "crevice cast", "crevice cavities"):
        assert option in UNRESOLVED_CHANNEL_HINT


def test_quiet_and_debug_are_accepted_before_or_after_the_command():
    from crevice.cli import build_parser
    parser = build_parser()
    for argv in (["--quiet", "profile", "1GRM", "-o", "x.json"], ["profile", "1GRM", "-o", "x.json", "-q"]):
        args = parser.parse_args(argv)
        assert args.quiet is True and args.debug is False


# --- progress ----------------------------------------------------------------

class _Terminal(io.StringIO):
    def isatty(self):
        return True


def _run(stream, steps=3, **kwargs):
    with progress.ProgressBar(steps, "frames", stream=stream, **kwargs) as bar:
        for _ in range(steps):
            bar.update()


def test_progress_is_silent_unless_enabled(monkeypatch):
    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    stream = _Terminal()
    _run(stream)
    assert stream.getvalue() == ""


def test_progress_is_quiet_when_not_a_terminal(monkeypatch):
    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    stream = io.StringIO()
    with progress.progress_reporting():
        _run(stream)
    assert stream.getvalue() == ""


def test_progress_redraws_one_line_on_a_terminal(monkeypatch):
    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    stream = _Terminal()
    with progress.progress_reporting():
        _run(stream, min_interval=0)
    text = stream.getvalue()
    assert "\r" in text and "frames: [" in text and "3/3 (100%)" in text and text.endswith("\n")
    assert text.count("\n") == 1


def test_progress_env_forces_lines_or_silence(monkeypatch):
    stream = io.StringIO()
    monkeypatch.setenv(progress.PROGRESS_ENV, "1")
    with progress.progress_reporting():
        _run(stream, steps=20)
    lines = stream.getvalue().splitlines()
    assert len(lines) == 10 and lines[-1].startswith("frames: [") and "20/20" in lines[-1]
    monkeypatch.setenv(progress.PROGRESS_ENV, "0")
    terminal = _Terminal()
    with progress.progress_reporting():
        _run(terminal)
    assert terminal.getvalue() == ""


def test_nested_and_worker_bars_are_silent(monkeypatch):
    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    outer, inner = _Terminal(), _Terminal()
    with progress.progress_reporting():
        with progress.ProgressBar(2, "outer", stream=outer, min_interval=0) as bar:
            _run(inner)
            bar.update(2)
        monkeypatch.setattr(progress.os, "getpid", lambda: -1)  # as in a forked worker
        worker = _Terminal()
        _run(worker)
    assert inner.getvalue() == "" and worker.getvalue() == ""
    assert "outer" in outer.getvalue()


def test_trajectory_progress_goes_to_stderr_and_never_changes_outputs(tmp_path, monkeypatch, capsys):
    from test_trajectory_input import build_universe
    import MDAnalysis as mda
    import numpy as np

    universe = build_universe()
    topology, trajectory = tmp_path / "system.gro", tmp_path / "run.xtc"
    universe.atoms.write(str(topology))
    base = np.arange(len(universe.atoms) * 3, dtype=np.float32).reshape(-1, 3)
    with mda.Writer(str(trajectory), n_atoms=len(universe.atoms)) as writer:
        for index in range(5):
            universe.atoms.positions = base + index * 10.0
            writer.write(universe.atoms)

    def run(name, *extra):
        out = tmp_path / f"{name}.json"
        assert main(["trajectory", str(trajectory), "--topology", str(topology), "--selection", "name CA",
                     "-o", str(out), "--samples", "5", "--search-radius", "0", *extra]) == 0
        return out, capsys.readouterr()

    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    quiet_out, quiet = run("default")                  # stderr is not a terminal: no progress
    assert quiet.err == ""
    monkeypatch.setenv(progress.PROGRESS_ENV, "1")
    logged_out, logged = run("logged")
    assert "trajectory frames: [" in logged.err and "5/5 (100%)" in logged.err
    assert logged.out.replace("logged", "default") == quiet.out
    silenced_out, silenced = run("silenced", "--quiet")
    assert silenced.err == ""

    def payload(path):
        data = json.loads(path.read_text(encoding="utf-8"))
        data.pop("metadata", None)
        return data

    assert payload(logged_out) == payload(quiet_out) == payload(silenced_out)
