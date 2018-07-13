#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 12/06/2018

import sys
import sortscaffolds


__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use: ./gene_over.py genes.gff")
    else:
        process(open(sys.argv[1], "r"))
