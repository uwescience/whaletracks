"""Acquires data used from obspy and builds the database."""

from common_python.database import database_util as util
import whaletracks.common. constants as cn

import matplotlib.pyplot as plt
import obspy
import os
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from shapely.geometry import Point, Polygon
import pandas as pd
import numpy as np


class DBBuilder(object):

  def __init__(self, network_code="7D",
      start_time= UTCDateTime("2001-01-01"), 
      end_time=UTCDateTime("2020-01-02"), db_pth=cn.DB_PTH):
    """
    :param str network_code:
    :param UTCDateTime start_time: starting time of data
    :param UTCDateTime end_time: starting time of data
    Usage
      builder = DBUilder()
      builder.build()
    """
    self._network_code = network_code
    self._start_time = start_time
    self._end_time = end_time
    self._db_pth = db_pth

  def build(self):
    """
    Builds all of the database tables.
    """
    util.updateDBTable(self._getStationDF(), self._db_pth,
        cn.SCM_STATION.tablename)

  ########## Table building methods. ##############
  def _getStationDF(self):
    """
    Constructs a dataframe of station metadata.
    """
    client = Client("IRIS")
    inventory = client.get_stations(network=self._network_code, station="*",
        starttime=self._start_time, endtime=self._end_time)
    network = inventory[0]
    station_count = network.total_number_of_stations
    station_dct = {k: [] for k in cn.SCM_STATION.columns}
    station_dct["network"] = list(np.repeat(self._network_code,
        station_count))
    for idx in range(station_count):
        station = network[idx]
        for key in cn.SCM_STATION.columns:
            station_dct[key].append(str(station.__getattribute__(key)))
    return pd.DataFrame(station_dct)
