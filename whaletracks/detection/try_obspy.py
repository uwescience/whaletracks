#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 11:10:59 2020

@author: wader
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import os
import scipy.io.wavfile as siow
import numpy as np
import scipy.signal as sig
import matplotlib.colors as color
import matplotlib.animation as animation
import datetime
import detect_calls as detect

f0 = 16 #average start frequency
f1 = 14.6 #average end frequency
bdwdth = 1 # average bandwidth
dur = 10 #average duration


client_code = 'IRIS'
client = Client('IRIS')

starttime=("2011-12-14T12:00:00.000")
endtime=("2011-12-14T12:30:00.000")

UTCstart=UTCDateTime(starttime)
UTCend=UTCDateTime(endtime)

st_raw=client.get_waveforms(network="7D",station='*',location='*',channel='BHZ,HHZ',starttime=UTCstart,endtime=UTCend)
  
num_sta=len(st_raw)


for j in range(14,num_sta):

    tr=st_raw[j]

    tr_filt=tr.copy()
    
    [f,t,Sxx]=detect.plotwav(tr_filt.stats.sampling_rate, tr_filt.data, window_size=5, overlap=.95)
    
    [tvec, fvec, BlueKernel]=detect.buildkernel(f0, f1, bdwdth, dur, f, t, tr_filt.stats.sampling_rate)
    
    [t_scale, CorrVal_scale]=detect.xcorr(t,f,Sxx,tvec,fvec,BlueKernel)
 
    
    #tr_filt.filter('bandpass',freqmin=5,freqmax=22,corners=2,zerophase=True)

    #plt.plot(tr_filt.data[100:3000])

    #tr_filt.spectrogram(log=False, title='spect' + str(tr_filt.stats.starttime))

