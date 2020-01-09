"""Acquires data used from obspy."""

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


class DataProvider(object):

  def __init__(self, network_code="7D",
      start_time= UTCDateTime("2001-01-01"), 
      end_time=UTCDateTime("2020-01-02")):
    """
    :param str network_code:
    :param UTCDateTime start_time: starting time of data
    :param UTCDateTime end_time: starting time of data
    """
    self._network_code = network_code
    self._start_time = start_time
    self._end_time = end_time
    #
    self.df_station = self._getStationDF()

  def _getStationDF(self):
    """
    Constructs a dataframe of station metadata.
    :return pd.DataFrame:
        station_code, channel_count, start_time, end_time, access, latitude, longitude, elevation
    """
    client = Client("IRIS")
    inventory = client.get_stations(network=self._network_code, station="*",
        starttime=self._start_time, endtime=self._end_time)
    network = inventory[0]
    station_count = network.total_number_of_stations
    station_dct = {k: [] for k in cn.STATION_SCM.columns}
    station_dct["network"] = list(np.repeat(self._network_code,
        station_count))
    for idx in range(station_count):
        station = network[idx]
        for key in cn.STATION_SCM.columns:
            station_dct[key].append(str(station.__getattribute__(key)))
    return pd.DataFrame(station_dct)

  def updateDB(self):
    util.updateDBTable(self.df_station, cn.DB_PTH,
        cn.STATION_SCM.tablename)
