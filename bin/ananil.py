#!/usr/bin/env python

# test dataset: /neurospin/unicog/protocols/IRMf/Dighiero_Dehaene_7TLinguistics_2023/

import sys
# remove our executable dir from import path
# it would reload the executable as a module instead of looking for the real
# module
del sys.path[0]

from ananil import main

main()
