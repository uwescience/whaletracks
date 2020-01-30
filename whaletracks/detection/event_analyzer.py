#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:15:18 2020

@author: wader
"""

import pandas as pd
import scipy.signal as sig
import matplotlib.pyplot as plt

PEAK_TIME = "peak_time"
PEAK_SIGNAL = "peak_signal"
DURATION = "duration"
START_TIME = "start_time"
END_TIME = "end_time"
THRESHOLD = "threshold"
MIN_SIGNAL = "min_signal"

COLUMNS = [PEAK_TIME, PEAK_SIGNAL, DURATION, START_TIME, END_TIME, 
           THRESHOLD, MIN_SIGNAL]
       

SECONDS_IN_MINUTE = 60

class EventAnalyzer(object):
    
    def __init__(self, times, values, start_chunk, dur=10, prominence=1.5):
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
            rel_height=.9)
        self.df = self._makePeakDF(peak_indicies, peak_properties)
        self.df[THRESHOLD] = prominence
     
    def _makePeakDF(self, peak_indicies, peak_properties):
        """
        :param int index: index of peak
        :return pd.DataFrame: COLUMNS, except THRESHOLD
        """
        dct = {k: [] for k in COLUMNS if k != THRESHOLD}
        for index in range(len(peak_indicies)):
            dct[PEAK_TIME].append(self.times[peak_indicies[index]])
            dct[PEAK_SIGNAL].append(peak_properties["prominences"][index])
            dct[START_TIME].append(
                    self.times[peak_properties["left_ips"].astype(int)[index]])
            dct[END_TIME].append(
                    self.times[peak_properties["right_ips"].astype(int)[index]])
            dct[MIN_SIGNAL].append(peak_properties["width_heights"][index])
            dct[DURATION].append(self.times[peak_properties["right_ips"].astype(int)[index]]-
                                 self.times[peak_properties["left_ips"].astype(int)[index]])
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