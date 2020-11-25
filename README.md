![](https://travis-ci.com/uwescience/whaletracks.svg?branch=master)

# whaletracks
Analysis of sound data to detect and track whales.

The data are structured as ...

## Helpful hints

- Google will search for GPS coordinates in decimal

## Resources

- [IRIS rest interface](https://service.iris.edu/irisws/fedcatalog/1/) (includes querying stations by geographic coordinates)

- [Map with US counties](https://www.randymajors.com/p/countygmap.html)
- [Query map by lat/lon](http://ds.iris.edu/gmap/#maxlat=50&maxlon=-124&minlat=38&minlon=-132&network=*&drawingmode=box&planet=earth)
- [IRIS metadata aggragator](http://ds.iris.edu/mda/7D/FC03D/?starttime=2014-09-07T00:00:00&endtime=2015-10-02T23:59:59)


## Installation
Install [github desktop](https://desktop.github.com/) (or you can use git in the terminal if you prefer)
Install [Anaconda](https://www.anaconda.com/products/individual) for your operating system


Assumes that ``anaconda`` and ``git`` are installed.

- ``git clone --recurse-submodules https://github.com/uwescience/whaletracks.git``
- ``cd whaletracks``
- ``bash setup.sh``



## Running codes
You must setup the PYTHONPATH before running codes.
- `cd whaletracks`
- `source setup_run.sh`
