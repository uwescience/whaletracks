import sys
import os

PROJECT_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))

def addSubmodule(name, project_dir=PROJECT_DIR):
  new_dir = os.path.join(project_dir, name)
  sys.path.append(new_dir)

addSubmodule(PROJECT_DIR)
addSubmodule("common_python")
