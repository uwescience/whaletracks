#! /bin/bash
# Sets up the conda environment whaletracks
NAME="whalet"
# Setup for project
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info --all
conda env create --quiet --name ${NAME} --file environment.yml
conda install jupyter notebook
# Install the database browser
sudo apt-get install sqlitebrowser
echo "Use 'conda activate ${NAME}' to enter environment."
echo "Use 'conda deactivate' to exit environment."
