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

################### HELPER FUNCTIONS ########################
def _makeChannelID(network_code, station_code, channel_code):
  return "%s.%s.%s" % (network_code, station_code, channel_code)

def _makeStationID(network_code, station_code):
  return "%s.%s" % (network_code, station_code)

def _stringifyComplexs(complexs):
  return ";".join([str((x.real, x.imag)) for x in complexs])

def _getPolesZeroes(channel):
  try:
      zeroes = _stringifyComplexs(channel.response.get_paz().zeros)
  except:
      zeroes = ""
  try:
      poles = _stringifyComplexs(channel.response.get_paz().poles)
  except:
      poles = ""
  return poles, zeroes


################### CLASS ########################
class DBBuilder(object):
  """
  Builds the whaletracks database. Useage:
    builder = DBBuilder()
    builder.build()
  Database path is specified in constants.py as DB_PTH.
  """

  def __init__(self, network_code="7D", channel="*Z",
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
    #
    self.client = Client("IRIS")
    self.inventory = self.client.get_stations(
        network=self._network_code,
        station="*", channel=channel, level="response", 
        starttime=self._start_time, endtime=self._end_time)
    self.network = self.inventory[0]

  def build(self):
    """
    Builds all of the database tables.
    """
    util.updateDBTable(self._makeStationDF(), self._db_pth,
        cn.SCM_STATION.tablename)
    df = self._makeChannelDF()
    df[cn.END_DATE] = [str(d) for d in df[cn.END_DATE]]
    df[cn.START_DATE] = [str(d) for d in df[cn.START_DATE]]
    util.updateDBTable(df, self._db_pth,
        cn.SCM_CHANNEL.tablename)

  ########## Table building methods. ##############
  def _makeStationDF(self):
    """
    Constructs a dataframe of station metadata.
    """
    CODE = "code"
    attributes = list(cn.SCM_STATION.columns)
    attributes.remove(cn.NETWORK_CODE)
    attributes.remove(cn.STATION_CODE)
    attributes.append(CODE)
    #
    station_count = self.network.total_number_of_stations
    station_dct = {k: [] for k in attributes}
    for station in self.network:
      for key in attributes:
        station_dct[key].append(str(station.__getattribute__(key)))
    station_dct[cn.STATION_CODE] = station_dct[CODE]
    station_dct[cn.NETWORK_CODE] = np.repeat(self._network_code,
        len(station_dct[cn.STATION_CODE]))
    del station_dct[CODE]
    return pd.DataFrame(station_dct)
        
  def _makeChannelDF(self):
    """
    :param obspy.inventory.network.Network network:
    :return pd.DataFrame: channels, poles, zeros
      channel dataframe: cn.SCM_CHANNEL.columns
      poles, zeros
          
    """
    channel_dct = {k: [] for k in cn.SCM_CHANNEL.columns}
    for station in self.network:
      station_id = _makeStationID(self.network.code, station.code)
      for channel in station:
        poles, zeroes = _getPolesZeroes(channel)
        channel_dct[cn.AZIMUTH].append(channel.azimuth)
        channel_dct[cn.CHANNEL_ID].append(_makeChannelID(
            self.network.code, station.code, channel.code))
        channel_dct[cn.CHANNEL_TYPES].append(" ".join(channel.types))
        channel_dct[cn.DIP].append(channel.dip)
        channel_dct[cn.END_DATE].append(channel.end_date)
        channel_dct[cn.POLES].append(poles)
        channel_dct[cn.SENSITIVITY_FREQUENCY].append(channel.response.instrument_sensitivity.frequency)
        channel_dct[cn.SENSITIVITY_VALUE].append(channel.response.instrument_sensitivity.value)
        channel_dct[cn.SENSOR].append(channel.sensor.description)
        channel_dct[cn.START_DATE].append(channel.start_date)
        channel_dct[cn.STATION_ID].append(_makeStationID(
            self.network.code, station.code))
        channel_dct[cn.ZEROES].append(zeroes)
    return pd.DataFrame(channel_dct)


if __name__ == '__main__':
  # Create the database
  builder = DBBuilder()
  builder.build()
