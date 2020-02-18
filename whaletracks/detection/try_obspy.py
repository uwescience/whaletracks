#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 11:10:59 2020

@author: wader
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read, read_inventory
import matplotlib.pyplot as plt
import os
import scipy.io.wavfile as siow
import numpy as np
import scipy.signal as sig
import matplotlib.colors as color
import matplotlib.animation as animation
import datetime
import whaletracks.detection.detect_calls as detect
from whaletracks.detection.event_analyzer import EventAnalyzer
import pandas as pd
import whaletracks.common.constants as cn


f0 = 16 #average start frequency
f1 = 14.6 #average end frequency
bdwdth = 1 # average bandwidth
dur = 10 #average duration
thresh = 2;
kernel_lims=[f1-bdwdth*3, f0+bdwdth*3]
CHUNK_FILE = "analyzers.csv"

client_code = 'IRIS'
client = Client('IRIS')

starttime=("2011-12-14T12:00:00.000")
endtime=("2011-12-14T12:20:00.000")

#starttime=("2011-10-01T12:00:00.000")
#endtime=("2012-07-01T12:00:00.000")

UTCstart=UTCDateTime(starttime)
UTCend=UTCDateTime(endtime)

UTCstart_chunk=UTCstart
UTCend_chunk=UTCstart+1800

# Check for existing file
if os.path.isfile(CHUNK_FILE):
    analyzers = [pd.read_csv(CHUNK_FILE)]
else:
    analyzers = []

while UTCend > UTCstart_chunk:

    st_raw=client.get_waveforms(network="7D",station='*',location='*',
                                channel='BHZ,HHZ',starttime=UTCstart_chunk,endtime=UTCend_chunk,attach_response=True)
    st_raw.detrend(type="demean")
    st_raw.detrend(type="linear")
    st_raw.remove_response(output='VEL', pre_filt=(.2,.5,22,24))
    st_raw.remove_response()
    st_raw.remove_sensitivity()
    inventory=client.get_stations(network="7D",station='*',channel='BHZ,HHZ',
                                  level="response",starttime=UTCstart,endtime=UTCend)
    
    num_sta=len(st_raw)

    analyzers_chunk=[]
    for j in range(2,num_sta):

        tr=st_raw[j]
        

        tr_filt=tr.copy()
        #tr_filt.detrend(type="demean")
        #tr_filt.remove_response(output='VEL', pre_filt=(.002,.005,23,25))
        #tr_filt.remove_sensitivity()
        
        #if len(tr_filt.data) < tr_filt.stats.sampling_rate*10:
            #continue
    
        [f,t,Sxx]=detect.plotwav(tr_filt.stats.sampling_rate, tr_filt.data, window_size=5, overlap=.95)
    
        [tvec, fvec, BlueKernel, freq_inds]=detect.buildkernel(f0, f1, bdwdth, dur, f, t, tr_filt.stats.sampling_rate)
        
        #subset spectrogram to be in same frequency range as kernel
        Sxx_sub=Sxx[freq_inds,:][0]
        f_sub=f[freq_inds]
        
        [times, values]=detect.xcorr(t,f_sub,Sxx_sub,tvec,fvec,BlueKernel)
        
        
        analyzer_j = EventAnalyzer(times, values, UTCstart_chunk)
        analyzers_chunk.append(analyzer_j.df)
        
        analyzer_j.plot()
        
        
    UTCstart_chunk=UTCstart_chunk+1800
    UTCend_chunk=UTCend_chunk+1800
    
    
    analyzers.extend(analyzers_chunk)
    new_df = pd.concat(analyzers)
    new_df.to_csv(CHUNK_FILE)
    
    
final_analyzer=pd.concat(analyzers)
final_analyzer.to_csv(cn.DETECTION_PATH)
    
    #tr_filt.filter('bandpass',freqmin=5,freqmax=22,corners=2,zerophase=True)

    #plt.plot(tr_filt.data[100:3000])

    #tr_filt.spectrogram(log=False, title='spect' + str(tr_filt.stats.starttime))

