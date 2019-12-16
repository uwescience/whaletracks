#! /bin/bash
# arg1: name
if [ -z "$1" ]
then
      echo "Usage: 'setup.sh <env name>'"
      echo "       Defaulting to 'whaletracks'"
      NAME="whaletracks"
else
      NAME=$1
fi
echo ${NAME}
exit
# Setup for project
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info --all
conda env create --quiet --name ${NAME} --file environment.yml
conda install jupyter notebook
echo "Use 'conda activate ${NAME}' to enter environment."
echo "Use 'conda deactivate' to exit environment."
