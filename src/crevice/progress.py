"""Progress feedback for long runs, written to standard error.

A tiny built-in progress printer (no extra dependency). It reports only while
:func:`progress_reporting` is active, which the ``crevice`` command line turns
on for every command unless ``--quiet`` is given; Python callers see nothing
unless they opt in. Progress never goes to standard output and never changes
any output file.

On a terminal, a single line is redrawn in place, at most a few times per
second. When standard error is not a terminal (a log file, a pipe, a batch
job) nothing is printed, unless the environment variable ``CREVICE_PROGRESS``
is ``1``: then a plain line is written every 10% of the work. ``CREVICE_PROGRESS=0``
silences progress everywhere.

Worker processes never print: the trajectory commands update the bar in the
parent process as each frame's result arrives, so forked pools are safe.

Examples
--------
>>> import io
>>> from crevice.progress import ProgressBar, progress_reporting
>>> stream = io.StringIO()            # not a terminal: quiet by default
>>> with progress_reporting(), ProgressBar(3, "frames", stream=stream) as bar:
...     for _ in range(3):
...         bar.update()
>>> stream.getvalue()
''
"""
from __future__ import annotations

import os
import sys
import time
from contextlib import contextmanager
from typing import Iterator, TextIO

PROGRESS_ENV = "CREVICE_PROGRESS"
"""Environment variable that forces plain progress lines (``1``) or silences progress (``0``)."""

_ACTIVE = False
_OWNER_PID: int | None = None  # forked workers inherit _ACTIVE but never print
_OPEN = 0  # bars currently open; a bar nested inside another stays silent


@contextmanager
def progress_reporting(enabled: bool = True) -> Iterator[None]:
    """Enable (or disable) progress reporting inside a ``with`` block.

    Parameters
    ----------
    enabled : bool, default True
        Whether progress bars created inside the block may print. They still
        print only to a terminal unless ``CREVICE_PROGRESS=1``.

    Yields
    ------
    None
    """
    global _ACTIVE, _OWNER_PID
    previous = _ACTIVE, _OWNER_PID
    _ACTIVE, _OWNER_PID = bool(enabled), os.getpid()
    try:
        yield
    finally:
        _ACTIVE, _OWNER_PID = previous


def progress_mode(stream: TextIO | None = None) -> str | None:
    """Return how progress would be shown on ``stream`` right now.

    Parameters
    ----------
    stream : typing.TextIO, optional
        Destination; default ``sys.stderr``.

    Returns
    -------
    str or None
        ``"tty"`` (redrawn line), ``"lines"`` (plain lines, forced by
        ``CREVICE_PROGRESS=1``) or None (silent: reporting is off, this is a
        worker process, the stream is not a terminal, or ``CREVICE_PROGRESS=0``).
    """
    if not _ACTIVE or os.getpid() != _OWNER_PID:
        return None
    setting = os.environ.get(PROGRESS_ENV, "").strip().lower()
    if setting in {"0", "off", "no", "false"}:
        return None
    stream = sys.stderr if stream is None else stream
    try:
        tty = stream.isatty()
    except (AttributeError, ValueError):
        tty = False
    if tty:
        return "tty"
    if setting in {"1", "on", "yes", "true", "always"}:
        return "lines"
    return None


def _clock(seconds: float) -> str:
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    return f"{hours}:{rest // 60:02d}:{rest % 60:02d}" if hours else f"{rest // 60}:{rest % 60:02d}"


class ProgressBar:
    """A counter that reports ``done/total`` with elapsed and remaining time.

    Parameters
    ----------
    total : int
        Number of steps expected (0 or less: an open-ended count).
    label : str
        What is being counted, such as ``"profile frames"``.
    stream : typing.TextIO, optional
        Destination; default ``sys.stderr`` (looked up when the bar is made).
    min_interval : float, default 0.25
        Minimum seconds between redraws on a terminal.

    Notes
    -----
    Create and update the bar in the parent process only. ``update`` is cheap
    when progress is silent, so it can be called once per step. A bar created
    while another is open stays silent, so nested steps (such as the probe
    ladder inside a trajectory) never garble the outer line. Always
    :meth:`close` a bar, or use it as a context manager.
    """

    def __init__(self, total: int, label: str, *, stream: TextIO | None = None, min_interval: float = 0.25):
        self.total = int(total)
        self.label = label
        self.stream = sys.stderr if stream is None else stream
        global _OPEN
        self.mode = progress_mode(self.stream) if _OPEN == 0 else None
        _OPEN += 1
        self._open = True
        self.min_interval = float(min_interval)
        self.done = 0
        self._start = time.monotonic()
        self._last = -1e9
        self._next_decile = 1
        self._drawn = False

    def __enter__(self) -> "ProgressBar":
        """Return the bar itself.

        Returns
        -------
        ProgressBar
        """
        return self

    def __exit__(self, *exc) -> None:
        """Finish the line (see :meth:`close`).

        Parameters
        ----------
        *exc
            Exception information, ignored.
        """
        self.close()

    def _text(self) -> str:
        elapsed = time.monotonic() - self._start
        if self.total > 0:
            fraction = min(1.0, self.done / self.total)
            width = 20
            filled = int(fraction * width)
            text = (f"{self.label}: [{'#' * filled}{'.' * (width - filled)}] {self.done}/{self.total} "
                    f"({fraction:.0%}) {_clock(elapsed)} elapsed")
            if 0 < self.done < self.total:
                text += f", ~{_clock(elapsed / self.done * (self.total - self.done))} left"
            return text
        return f"{self.label}: {self.done} done, {_clock(elapsed)} elapsed"

    def update(self, n: int = 1) -> None:
        """Advance by ``n`` steps and redraw if due.

        Parameters
        ----------
        n : int, default 1
            Steps completed since the last update.
        """
        self.done += n
        if self.mode is None:
            return
        now = time.monotonic()
        if self.mode == "tty":
            if now - self._last >= self.min_interval or (self.total > 0 and self.done >= self.total):
                self._last = now
                self.stream.write("\r\033[K" + self._text())
                self.stream.flush()
                self._drawn = True
        elif self.total > 0:
            if self.done * 10 >= self._next_decile * self.total:
                while self.done * 10 >= self._next_decile * self.total:
                    self._next_decile += 1
                self.stream.write(self._text() + "\n")
                self.stream.flush()

    def close(self) -> None:
        """End the redrawn terminal line with a newline (no-op otherwise)."""
        global _OPEN
        if self._open:
            _OPEN -= 1
            self._open = False
        if self.mode == "tty" and self._drawn:
            self.stream.write("\r\033[K" + self._text() + "\n")
            self.stream.flush()
            self._drawn = False


def note(message: str, *, stream: TextIO | None = None) -> None:
    """Write one progress note line, under the same rules as :class:`ProgressBar`.

    Parameters
    ----------
    message : str
        Text of the line, such as ``"automatic enclosure probe: trying 0.50-3.00 A"``.
    stream : typing.TextIO, optional
        Destination; default ``sys.stderr``.
    """
    stream = sys.stderr if stream is None else stream
    if progress_mode(stream) is not None:
        stream.write(message + "\n")
        stream.flush()
