#!/home/admin/Desktop/Bioformatics-toolbox/env/bin/python3
# -*- coding: utf-8 -*-
import re
import sys
from biorun.methods.comm import run
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(run())
