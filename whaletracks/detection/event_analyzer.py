#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:15:18 2020

@author: wader
"""

import whaletracks.common.constants as cn

import pandas as pd
import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
from whaletracks.common import util

SECONDS_IN_MINUTE = 60
EXCLUDED_COLUMNS = [cn.THRESHOLD, cn.STATION_CODE, cn.NETWORK_CODE]

class EventAnalyzer(object):
    
    def __init__(self, times, values, start_chunk, dur=10, prominence=.5):
        """
        :param list-float times: offsets in seconds
        :param list-float values: values at times
        :param UTCdatetime: start time of chunk
        :param float dur: duration in seconds of call (default for blue whale)
        """
        self.times = [start_chunk + t for t in times]
        self.values = values
        self.start_chunk = start_chunk
        peak_indicies, peak_properties=sig.find_peaks(self.values,
            distance=dur*(1/(self.times[1]-self.times[0])),
            width=(dur/2)*(1/(self.times[1]-self.times[0])),
            prominence=prominence,
            wlen=SECONDS_IN_MINUTE*(1/(self.times[1]-self.times[0])),
            rel_height=.7)
        self.df = self._makeDetectionDF(peak_indicies, peak_properties)
        self.df[cn.THRESHOLD] = prominence
     
    def _makeDetectionDF(self, peak_indicies, peak_properties):
        """
        :param int index: index of peak
        :return pd.DataFrame: all columns, except for EXCLUDED_COLUMNS
        """
        dct = {k: [] for k in cn.SCM_DETECTION.columns 
               if not k in EXCLUDED_COLUMNS}
        for index in range(len(peak_indicies)):
            dct[cn.PEAK_TIME].append(self.times[peak_indicies[index]])
            dct[cn.PEAK_SIGNAL].append(peak_properties["prominences"][index])
            dct[cn.START_TIME].append(
                    self.times[peak_properties["left_ips"].astype(int)[index]])
            dct[cn.END_TIME].append(
                    self.times[peak_properties["right_ips"].astype(int)[index]])
            dct[cn.MIN_SIGNAL].append(peak_properties["width_heights"][index])
            dct[cn.END_EPOCH].append(util.datetimeToEpoch(
                    self.times[peak_properties["right_ips"].astype(int)[index]]))
            dct[cn.DURATION].append(self.times[peak_properties["right_ips"].astype(int)[index]]-
                                 self.times[peak_properties["left_ips"].astype(int)[index]])
            dct[cn.PEAK_EPOCH] = list(np.repeat(None, len(dct[cn.PEAK_TIME])))
            dct[cn.START_EPOCH] = list(np.repeat(None, len(dct[cn.START_TIME])))
            dct[cn.END_EPOCH] = list(np.repeat(None, len(dct[cn.END_TIME])))
        return pd.DataFrame(dct)
        
        
    def plot(self, is_plot=True):
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(self.times,self.values)
        ax.plot(self.df.peak_time,self.df.peak_signal,'x')
        plt.hlines(self.df.min_signal,self.df.start_time,
                   self.df.end_time,color="C2")
        if is_plot:
            plt.show(block=True)
