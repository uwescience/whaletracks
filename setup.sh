# Setup for project
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info --all
conda env create --quiet --name whaletracks --file environment.yml
echo "Now do: conda activate whaletracks"
