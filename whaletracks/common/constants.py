"""Constants used in the whale project."""

import os
import sys

PROJECT_NAME = "whaletracks"

# PATHS
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
PROJECT_CODE = os.path.join(PROJECT_DIR, PROJECT_NAME)
COMMON_CODE = os.path.join(PROJECT_DIR, "common_python")

# Add search paths
sys.path.insert(0, PROJECT_CODE)
sys.path.insert(0, COMMON_CODE)
