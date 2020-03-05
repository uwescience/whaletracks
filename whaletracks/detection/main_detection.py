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
import whaletracks.detection.detect_calls as detect
from whaletracks.detection.event_analyzer import EventAnalyzer
import pandas as pd
import whaletracks.common.constants as cn
from datetime import datetime

F0 = 15.7 #average start frequency
F1 = 14.6 #average end frequency
BDWDTH = .7 # average bandwidth
DUR = 10 #average duration


CHUNK_FILE = "analyzers.csv"

CLIENT_CODE = 'IRIS'
PLOTFLAG = False

#TEST STARTTIME = ("2011-12-14T12:00:00.000")
#TEST ENDTIME = ("2011-12-14T12:20:00.000")

STARTTIME = ("2012-02-07T02:30:00.000")
ENDTIME = ("2012-07-01T00:00:00.000")

HALF_HOUR = 1800  # in seconds
CHUNK_LENGTH=HALF_HOUR   #secnods

#starttime=("2011-10-01T12:00:00.000")
#endtime=("2012-07-01T12:00:00.000")


def main(STARTTIME, ENDTIME,
         client_code=CLIENT_CODE, f0=F0,
         f1=F1,bdwdth=BDWDTH,dur=DUR,
         detection_pth=cn.SCM_DETECTION.csv_path, 
         chunk_pth=CHUNK_FILE,station_ids="*",
         is_restart=True):
    """
    :param UTCDateTime starttime:
    :param UTCDateTime endtime:
    """
    # Check for existing file
    if os.path.isfile(CHUNK_FILE) and is_restart:
        analyzers = [pd.read_csv(CHUNK_FILE)]
    else:
        analyzers = []
        
    client = Client(client_code)
    utcstart = UTCDateTime(STARTTIME)
    utcend = UTCDateTime(ENDTIME)

    utcstart_chunk = utcstart
    utcend_chunk = utcstart + CHUNK_LENGTH
    
    while utcend > utcstart_chunk:
        #import pdb; pdb.set_trace()
        print(utcstart_chunk)
        
        retry=0
        st_raw_exist=False
        
        while st_raw_exist == False and retry < 5:
            try:
                st_raw=client.get_waveforms(network="7D", station=station_ids, location='*',
                                            channel='BHZ,HHZ', starttime=utcstart_chunk,
                                            endtime=utcend_chunk, attach_response=True)
                st_raw_exist=True
                retry=5
            except:
                retry=retry+1
                st_raw_exist=False
                print("Client failed: Retry " + str(retry) + " of 5 attempts")
                
                
        if st_raw_exist == False:
            print("WARNING: no data available from input station/times")
            utcstart_chunk=utcstart_chunk+CHUNK_LENGTH
            utcend_chunk=utcend_chunk+CHUNK_LENGTH
            continue
                
            
        st_raw.detrend(type="demean")
        st_raw.detrend(type="linear")
        st_raw.remove_response(output='VEL',pre_filt=[.2,.5,22,24])
        st_raw.remove_sensitivity()
        
        num_sta=len(st_raw)
        
        
    
        analyzers_chunk=[]
        for idx in range(1, num_sta+1):
    
            j = idx - 1
            tr=st_raw[j]
            
    
            tr_filt=tr.copy()
            
            if len(tr_filt.data) < tr_filt.stats.sampling_rate*60: #skip if less than 1 min of data
                continue
            if tr_filt.data[0]==tr_filt.data[1]: #skip if data is bad (indicated by constant data)
                continue
        
            [f,t,Sxx]=detect.plotwav(tr_filt.stats.sampling_rate, tr_filt.data, window_size=5, overlap=.95, plotflag=PLOTFLAG)
        
            [tvec, fvec, BlueKernel, freq_inds]=detect.buildkernel(f0, f1, bdwdth, dur, f, t, tr_filt.stats.sampling_rate, plotflag=PLOTFLAG)
            
            #subset spectrogram to be in same frequency range as kernel
            Sxx_sub=Sxx[freq_inds,:][0]
            f_sub=f[freq_inds]
            
            [times, values]=detect.xcorr_log(t,f_sub,Sxx_sub,tvec,fvec,BlueKernel, plotflag=PLOTFLAG)
            
            #import pdb; pdb.set_trace()
            analyzer_j = EventAnalyzer(times, values, utcstart_chunk)
            
            station_codes = np.repeat(tr_filt.stats.station,analyzer_j.df.shape[0])
            network_codes = np.repeat(tr_filt.stats.network,analyzer_j.df.shape[0])
            analyzer_j.df[cn.STATION_CODE] = station_codes
            analyzer_j.df[cn.NETWORK_CODE] = network_codes
            
            analyzers_chunk.append(analyzer_j.df)
            
            #analyzer_j.plot()
            
            
        utcstart_chunk=utcstart_chunk+CHUNK_LENGTH
        utcend_chunk=utcend_chunk+CHUNK_LENGTH
        
        
        analyzers.extend(analyzers_chunk)
        new_df = pd.concat(analyzers)
        new_df.to_csv(chunk_pth, index=False)
        
        
        
    #import pdb; pdb.set_trace()  
    
    if len(analyzers) == 0:
        print('WARNING: detections dataframe empty')
        final_analyzer_df = []
        
    else:
        final_analyzer_df = pd.concat(analyzers)
        final_analyzer_df.to_csv(detection_pth,index=False)
        
    return final_analyzer_df
    
    #tr_filt.filter('bandpass',freqmin=5,freqmax=22,corners=2,zerophase=True)

    #plt.plot(tr_filt.data[100:3000])

    #tr_filt.spectrogram(log=False, title='spect' + str(tr_filt.stats.starttime))

if __name__ == "__main__":
    main(STARTTIME, ENDTIME)
