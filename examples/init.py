import sys
import os

PROJECT_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
sys.path.append(PROJECT_DIR)

def addSubmodule(name, project_dir=PROJECT_DIR):
  new_dir = os.path.join(project_dir, name)
  sys.path.append(new_dir)

addSubmodule("common_python")
