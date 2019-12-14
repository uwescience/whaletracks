conda config --set channel_priority strict
conda config --add channels conda-forge
conda install -c conda-forge libiconv
conda install geopandas descartes
echo "Now do: conda activate base"
