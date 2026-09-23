#!/usr/bin/env python3
"""Compare CPU and GPU CIA.cmcrt files for a 1D atmosphere."""

import sys

from compare_rayleigh import main


if __name__ == "__main__":
    sys.exit(main("CIA"))
