#!/usr/bin/env python3

import fileinput
import re

firstLine = None
for line in fileinput.input():
    firstLine=line
    break

## sample: : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 6.0.1+652ffd1
version = re.search(r'Version (\d.\d.\d)+', firstLine).group(1)
print(f"{version}")

