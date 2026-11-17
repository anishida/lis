#!/usr/bin/env python3
"""Backward-compatible entrypoint for the generic FEA operator auditor."""
from fea_matrix_audit import *  # noqa: F401,F403 - keep the original module API

if __name__ == "__main__":
    import sys
    sys.exit(main())
