#!/usr/bin/env python3

import urllib.request
import sys
import os

proteinFtpDir = sys.argv[1]
organism = sys.argv[2]

os.makedirs(organism, exist_ok=True)

urllib.request.urlretrieve(proteinFtpDir, filename = f"{organism}/proteinCodingGenes.tsv.gz")