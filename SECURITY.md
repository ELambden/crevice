# Security policy

CREVICE is research software for analysing molecular structures and trajectories.
It downloads structures from the RCSB PDB when asked, reads user-supplied files
and writes scripts for molecular viewers (PyMOL, VMD, ChimeraX) that users run
themselves.

## Supported versions

Only the latest released version (and `main`) receives fixes. Please upgrade to
the latest release before reporting.

## Reporting a vulnerability

Please do **not** open a public issue for a suspected vulnerability (for
example, a crafted input file that causes code execution, or a generated viewer
script that runs unintended commands). Instead, report it privately through GitHub's private vulnerability
reporting ("Report a vulnerability" on the repository's Security tab).

Include the CREVICE version or commit, the command you ran and a minimal input
that reproduces the problem. We aim to acknowledge reports within two weeks;
this is a small academic project without a dedicated security team.
