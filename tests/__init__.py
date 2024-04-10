import os, sys


### preamble
# Get the current script's directory
current_dir = os.path.dirname(os.path.abspath('__init__.py'))
# Get the parent directory by going one level up
parent_dir = os.path.dirname(current_dir)
top_parent_dir = os.path.dirname(parent_dir)
# Add the parent directory to sys.path
sys.path.append(top_parent_dir)

from chem461 import pchem
