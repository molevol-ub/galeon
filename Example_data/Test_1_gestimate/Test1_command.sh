#!/bin/bash

which GALEON_ControlScript.py 

# Run using one g value
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 100

# Test several g values
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 150,200,300,400 -outdir g_estimation_Results_Directory_Several_gValues