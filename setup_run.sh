#!/bin/bash
# Setup python paths to run codes
PYTHONPATH=${PYTHONPATH}:`pwd`
export PYTHONPATH
cd common_python
bash setup_run.sh
cd ..
