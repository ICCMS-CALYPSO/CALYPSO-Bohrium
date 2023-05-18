#!/bin/bash

# nohup python -u calypso_bohrium.py > out 2>&1 &
nohup run_calypso --dft qe > out 2>&1 &
