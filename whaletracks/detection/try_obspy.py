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
import detect_calls as detect
from whaletracks.detection.event_analyzer import EventAnalyzer

f0 = 16 #average start frequency
f1 = 14.6 #average end frequency
bdwdth = 1 # average bandwidth
dur = 10 #average duration
thresh = 2;

client_code = 'IRIS'
client = Client('IRIS')

starttime=("2011-12-14T12:00:00.000")
endtime=("2011-12-14T12:20:00.000")

starttime=("2011-10-01T12:00:00.000")
endtime=("2012-07-01T12:00:00.000")

UTCstart=UTCDateTime(starttime)
UTCend=UTCDateTime(endtime)

UTCstart_chunk=UTCstart
UTCend_chunk=UTCstart+1800

while UTCend > UTCstart_chunk:

    st_raw=client.get_waveforms(network="7D",station='*',location='*',
                                channel='BHZ,HHZ',starttime=UTCstart_chunk,endtime=UTCend_chunk,attach_response=True)
    
    inventory=client.get_stations(network="7D",station='*',channel='BHZ,HHZ',
                                  level="response",starttime=UTCstart,endtime=UTCend)
    
    num_sta=len(st_raw)


    for j in range(1,num_sta):

        tr=st_raw[j]
        inv=inventory[j]

        tr_filt=tr.copy()
        
        if len(tr_filt.data) < tr_filt.stats.sampling_rate*10:
            continue
    
        [f,t,Sxx]=detect.plotwav(tr_filt.stats.sampling_rate, tr_filt.data, window_size=5, overlap=.95)
    
        [tvec, fvec, BlueKernel]=detect.buildkernel(f0, f1, bdwdth, dur, f, t, tr_filt.stats.sampling_rate)
    
        [times, values]=detect.xcorr(t,f,Sxx,tvec,fvec,BlueKernel)
        
        
        
        
        
        some_object=detect.pick_peaks(t_scale,CorrVal_scale,dur)
        
        corr_peaks, peak_properties=sig.find_peaks(CorrVal_scale,distance=dur*(1/(t_scale[1]-t_scale[0])),
                                                   width=(dur/2)*(1/(t_scale[1]-t_scale[0])),
                                                   prominence=1.5,wlen=60*(1/(t_scale[1]-t_scale[0])),rel_height=.9)
   
        plt.plot(t_scale,CorrVal_scale)
        plt.plot(t_scale[corr_peaks],CorrVal_scale[corr_peaks],"x")
        plt.hlines(peak_properties["width_heights"],t_scale[peak_properties["left_ips"].astype(int)],
                                   t_scale[peak_properties["right_ips"].astype(int)],color="C2")
        
        det_times
        det_strength
        det_start
        det_end
        det_length
        det_thresh
        det_width_height
        
        
    UTCstart_chunk=UTCstart_chunk+1800
    UTCend_chunk=UTCend_chunk+1800
    
    #tr_filt.filter('bandpass',freqmin=5,freqmax=22,corners=2,zerophase=True)

    #plt.plot(tr_filt.data[100:3000])

    #tr_filt.spectrogram(log=False, title='spect' + str(tr_filt.stats.starttime))

