"""Access methods for acquiring data from the Whale Database"""

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


#Tables
TABLE_STATIONS = "stations"
TABLE_DETECTIONS = "detections"
TABLE_CHANNELS = "channels"


class DBAccessor(object):

  def __init__(self):
    """
    :param str network_code:
    :param UTCDateTime start_time: starting time of data
    :param UTCDateTime end_time: starting time of data
    """
    self._conn = None
    self._tables = {}  # Dataframes for tables

  @property
  def df_station(self):
    return self._readTable(self, cn.SCM_STATION.tablename)

  @property
  def df_channel(self):
    return self._readTable(self, cn.SCM_CHANNEL.tablename)

  @property
  def df_detection(self):
    return self._readTable(self, cn.SCM_DETECTION.tablename)

  def _connect(self):
    self._conn = sqlite3.connect(cn.DB_PTH)

  def _close(self):
    if self._conn is not None:
      self._conn.close()
      self._conn = None

  def _readTable(self, table_name):
    """
    Reads the specified table from the whales database.
    :param str table_name:
    :return pd.DataFrame:
    """
    if not table_name in self._tables.keys():
      self._connect()
      self._tables[table_name] = pd.read_table(tablename, self._conn)
      self._close()
    return self._tables[table_name]
