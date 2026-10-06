"""Run the CREVICE command line with ``python -m crevice``.

Equivalent to the ``crevice`` console script: importing this module calls
:func:`crevice.cli.main` and exits with its return code. It is therefore
excluded from the API documentation and from doctest collection.
"""

from .cli import main

raise SystemExit(main())
