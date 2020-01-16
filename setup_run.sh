#!/bin/bash
# Setup python paths to run codes
PYTHONPATH=${PYTHONPATH}:`pwd`
export PYTHONPATH
cd common_python
source setup_run.sh
cd ..
