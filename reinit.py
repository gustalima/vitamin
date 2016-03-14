#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os

files = glob.glob('*')

files.remove('vit.py')
files.remove('vit.glade')
files.remove('src')
files.remove('gui')
files.remove('reinit.py')

for file in files:
        os.system('rm -rf %s' %file)
