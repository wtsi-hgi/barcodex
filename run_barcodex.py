#!/usr/bin/env python3
"""
Launcher:  python run_barcodex.py [options]  is the same as  python -m barcodex [options],
but works from any directory (the package barcodex/ lives next to this file).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from barcodex.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
