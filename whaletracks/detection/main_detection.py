#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 11:10:59 2020

This code detects fin or blue whale calls using spectrogram cross-correlation and
stores detection metrics in comma-separated variable files.

@author: wader
"""

import sys

sys.path.append('/Users/wader/Desktop/whaletracks/') 

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read, read_inventory
import matplotlib.pyplot as plt
import os
import scipy.io.wavfile as siow
import numpy as np
import scipy.signal as sig
import matplotlib.colors as color
import detect_calls as detect
from whaletracks.detection.event_analyzer import EventAnalyzer
from whaletracks.common.util import datetimeToEpoch
import pandas as pd
import whaletracks.common.constants as cn
from datetime import datetime
from whaletracks.common.util import datetimeToEpoch
from scipy.signal import hilbert

FINFLAG = True #True if detecting fins, False if detecting blues
#CHUNK_FILE = "Blues_chunk_test.csv"
CHUNK_FILE = "Fins_chunk_marianas.csv" #Name of saved call file
#FIN_DET_SERIES = "fin_series.csv"
PLOTFLAG = False #Use if troubleshooting and want to see plots.
MP_FLAG = True #Use if storing Fin multipath info
CLIENT_CODE = 'IRIS'
network="XF" #Network name "OO" for OOI, "7D" for Cascadia, "XF" for marianas
station="B19" #Specific station, or '*' for all available stations
location='*'  # '*' for all available locations
channel='BHZ,HHZ,ELZ' #Choose channels,  you'll want 'BHZ,HHZ' for Cascadia
                      #Check http://ds.iris.edu/mda/OO/ for OOI station channels

#DET_PATH=cn.SCM_DETECTION.csv_path
DET_PATH="Fins_final_marianas.csv" #Final save file

if FINFLAG == False:
    #Build blue whale B-call characteristics - wide
    F0 = 16 #average start frequency
    F1 = 14.3 #average end frequency
    BDWDTH = 1 # average bandwidth
    DUR = 15 #average duration

if FINFLAG:
    #Build fin whale call characteristics
    F0 = 25 #average start frequency
    F1 = 15 #average end frequency
    BDWDTH = 3 # average bandwidth
    DUR = 2 #average duration

#Blue whale A-call characteristics
#F0 = 14.5 #average start frequency
#F1 = 14.2 #average end frequency
#BDWDTH = .5 # average bandwidth
#DUR = 20 #average duration

#Blue whale B-call characteristics - narrow
#F0 = 15.7 #average start frequency
#F1 = 14.6 #average end frequency
#BDWDTH = .7 # average bandwidth
#DUR = 10 #average duration


STARTTIME = ("2012-02-02T00:00:00.000") #for marianas fins
ENDTIME = ("2013-02-06T00:00:00.000")

#STARTTIME = ("2012-01-09T04:10:00.000") # for blue whale freq testing FN14A
#ENDTIME = ("2012-01-09T04:20:00.000")

#STARTTIME = ("2012-03-30T21:37:00.000") # for fin max call testing marianas
#ENDTIME = ("2012-03-30T22:38:00.000")

#STARTTIME = ("2018-08-25T13:07:00.000") #for testing on FN14A fins
#ENDTIME = ("2018-08-25T15:37:00.000")

HALF_HOUR = 1800  # in seconds
CHUNK_LENGTH=HALF_HOUR  #secnods

#starttime=("2011-10-01T12:00:00.000")
#endtime=("2012-07-01T12:00:00.000")

#Fin instruments: network ="XF" and station="B19"

def main(STARTTIME, ENDTIME,
         client_code=CLIENT_CODE, f0=F0,
         f1=F1,bdwdth=BDWDTH,dur=DUR,
         detection_pth=DET_PATH, 
         chunk_pth=CHUNK_FILE,station_ids=station,
         is_restart=True):
    """
    :param UTCDateTime starttime: ex. STARTTIME = ("2012-03-30T21:37:00.000")
    :param UTCDateTime endtime: ex. ENDTIME = ("2012-03-30T22:38:00.000")
    """

    if os.path.isfile(CHUNK_FILE) and is_restart:
        analyzers = [pd.read_csv(CHUNK_FILE)]
    else:
        analyzers = []
        
    client = Client(client_code)
    utcstart = UTCDateTime(STARTTIME)
    utcend = UTCDateTime(ENDTIME)

    utcstart_chunk = utcstart
    utcend_chunk = utcstart + CHUNK_LENGTH
    
    #Loop through times between starttime and endtime
    while utcend > utcstart_chunk:
        #import pdb; pdb.set_trace()
        print(utcstart_chunk)
        
        retry=0
        st_raw_exist=False

        #Attempt to get waveforms
        while st_raw_exist == False and retry < 5:
            try:
                st_raw=client.get_waveforms(network=network, station=station_ids, location=location,
                                            channel=channel, starttime=utcstart_chunk,
                                            endtime=utcend_chunk, attach_response=True)
                st_raw_exist=True
                retry=5
            except:
                retry=retry+1
                st_raw_exist=False
                print("Client failed: Retry " + str(retry) + " of 5 attempts")
                
        #Check if waveform data exists        
        if st_raw_exist == False:
            print("WARNING: no data available from input station/times")
            utcstart_chunk=utcstart_chunk+CHUNK_LENGTH
            utcend_chunk=utcend_chunk+CHUNK_LENGTH
            continue
                

        #Remove sensitivity and response, and filter data
        st_raw.detrend(type="demean")
        st_raw.detrend(type="linear")
        st_raw.remove_response(output='VEL',pre_filt=[1,3,40,45])
        st_raw.remove_sensitivity()
        

        num_sta=len(st_raw)
        analyzers_chunk=[]
        #Run detector on each station
        for idx in range(1, num_sta+1):
    
            j = idx - 1
            tr=st_raw[j]
            
    
            tr_filt=tr.copy()
            
            if len(tr_filt.data) < tr_filt.stats.sampling_rate*60*5: #skip if less than 1 min of data
                continue
            if tr_filt.data[0]==tr_filt.data[1]: #skip if data is bad (indicated by constant data)
                continue



            #Build detection metrics for either fin or blue whale calls
            if FINFLAG:

                #Spectrogram metrics
                window_size=2
                overlap=.95
                freqlim=[6, 34]
                #SNR metrics
                snr_limits=[15, 25]
                snr_calllength=1
                snr_freqwidth=5
                #Event metrics
                prominence=.6 #min threshold
                event_dur= .5 #minimum width of detection
                distance=18 #minimum distance between detections
                rel_height=.8



            if FINFLAG == False:
                #Spectrogram metrics
                window_size=5
                overlap=.95
                freqlim=[10, 20]
                #SNR metrics
                snr_limits=[14, 16]
                snr_calllength=4
                snr_freqwidth=.6
                #Event metrics
                prominence=.5 #min threshold
                event_dur= 5 #minimum width of detection
                distance=18 #minimum distance between detections
                rel_height=.7
                
             
            #Make spectrogram
            [f,t,Sxx]=detect.plotwav(tr_filt.stats.sampling_rate, tr_filt.data, window_size=window_size, overlap=overlap, plotflag=PLOTFLAG,filt_freqlim=freqlim,ylim=freqlim)
            
            #Make detection kernel
            [tvec, fvec, BlueKernel, freq_inds]=detect.buildkernel(f0, f1, bdwdth, dur, f, t, tr_filt.stats.sampling_rate, plotflag=PLOTFLAG, kernel_lims=detect.defaultKernelLims)
            
            #subset spectrogram to be in same frequency range as kernel
            Sxx_sub=Sxx[freq_inds,:][0]
            f_sub=f[freq_inds]
            
            #Run detection using built kernel and spectrogram
            [times, values]=detect.xcorr_log(t,f_sub,Sxx_sub,tvec,fvec,BlueKernel, plotflag=PLOTFLAG,ylim=freqlim)
            
           #Pick detections using EventAnalyzer class
            analyzer_j = EventAnalyzer(times, values, utcstart_chunk, dur=event_dur, prominence=prominence, distance=distance, rel_height=rel_height)
            #analyzer_j.plot()

            #get multipath data for fins
            if MP_FLAG: #if MP_FLAG is True
                mp_df = pd.DataFrame(columns=cn.SCM_MULTIPATHS.columns)
                samples =list(range(0,len(tr_filt.data)))
                sos = sig.butter(4, np.array([15, 20]), 'bp', fs=tr_filt.stats.sampling_rate, output = 'sos') 
                filtered_data = sig.sosfiltfilt(sos, tr_filt.data)
                amplitude_envelope = abs(hilbert(filtered_data))
                sos = sig.butter(4, 3, 'lp', fs=tr_filt.stats.sampling_rate, output = 'sos')
                filtered_env = sig.sosfiltfilt(sos, amplitude_envelope)
                #seconds=[s/tr_filt.stats.sampling_rate for s in samples] 

                utctimes=[utcstart_chunk + t for t in times]
                dt_up=10
                dt_down=25
                
                for det in range(0,len(analyzer_j.df)):
                    master_det = analyzer_j.df['start_time'][det]
                    master_det_sec = (master_det - utcstart_chunk)*int(tr_filt.stats.sampling_rate)

                    start_sample=max([master_det_sec-dt_up*tr_filt.stats.sampling_rate,min(samples)])
                    end_sample=min([master_det_sec+dt_down*tr_filt.stats.sampling_rate,max(samples)])

                    #start_seq=max([master_det-dt_up,min(utctimes)])
                    #end_seq=min([master_det+dt_down,max(utctimes)])

                    #seq_utc = [utctime for utctime in utctimes if start_seq < utctime < end_seq]
                    seq_samples = [s for s in samples if start_sample < s < end_sample]
                    seq_seconds = [s/tr_filt.stats.sampling_rate for s in seq_samples]
                    #samp_seconds = [s-seq_seconds[0] for s in seq_samples]
                    seq_data = filtered_env[seq_samples]

                    #xy, x_ind, seq_inds = np.intersect1d(seq_utc, utctimes, return_indices=True)
                    #seq_vals = values[seq_inds]

                    #mp_event=analyzer_j.mp_picker(seq_utc, seq_vals, master_det)
                    
                    mp_event=analyzer_j.mp_picker(seq_seconds, seq_data, utcstart_chunk, prominence=10**-17)
                    mp_df=mp_df.append(mp_event,ignore_index=True)
                    
                    
                   
                #import pdb; pdb.set_trace() 
                analyzer_j.df= pd.concat([analyzer_j.df, mp_df], axis=1)
                if PLOTFLAG:
                    seconds=[s/100 for s in samples]
                    utctimes=[utcstart_chunk+s for s in seconds]
                    plt.plot(utctimes,filtered_data)
                    plt.plot(utctimes,filtered_env)
                    plt.scatter(mp_df['arrival_1'],mp_df['amp_1'])
                    plt.scatter(mp_df['arrival_2'],mp_df['amp_2'])
                    plt.scatter(mp_df['arrival_3'],mp_df['amp_3'])
                    plt.scatter(mp_df['arrival_4'],mp_df['amp_4'])
                    plt.scatter(mp_df['arrival_5'],mp_df['amp_5'])
                    plt.xlabel('seconds')
                    plt.ylabel('amplitude')
                    plt.show()

                #import pdb; pdb.set_trace()

            #Calculate SNR info
            [snr,ambient_snr] = detect.get_snr(analyzer_j, t, f_sub, Sxx_sub, utcstart_chunk,snr_limits=snr_limits,snr_calllength=snr_calllength,snr_freqwidth=snr_freqwidth,dur=dur)

            #Add freq info for blue whales
            if FINFLAG == False:
                #These freqency metrics only work for blue whale calls
                [peak_freqs,start_freqs,end_freqs,peak_stds,start_stds,end_stds] = detect.freq_analysis(analyzer_j,t,f_sub,Sxx_sub,utcstart_chunk)
            if FINFLAG:
                #Fill frequency metric columns in csv with Nones if Fin calls
                peak_freqs=list(np.repeat(None, len(snr)))
                start_freqs=list(np.repeat(None, len(snr)))
                end_freqs=list(np.repeat(None, len(snr)))
                peak_stds=list(np.repeat(None, len(snr)))
                start_stds=list(np.repeat(None, len(snr)))
                end_stds=list(np.repeat(None, len(snr)))

            #Make dataframe with detections from current time chunk
            station_codes = np.repeat(tr_filt.stats.station,analyzer_j.df.shape[0])
            network_codes = np.repeat(tr_filt.stats.network,analyzer_j.df.shape[0])
            peak_epoch=datetimeToEpoch(analyzer_j.df['peak_time'])
            end_epoch=datetimeToEpoch(analyzer_j.df['end_time'])
            start_epoch=datetimeToEpoch(analyzer_j.df['start_time'])
            analyzer_j.df[cn.STATION_CODE] = station_codes
            analyzer_j.df[cn.NETWORK_CODE] = network_codes
            analyzer_j.df[cn.SNR] = snr
            analyzer_j.df['ambient_snr'] = ambient_snr
            analyzer_j.df['peak_frequency'] = peak_freqs
            analyzer_j.df['start_frequency'] = start_freqs
            analyzer_j.df['end_frequency'] = end_freqs
            analyzer_j.df['peak_frequency_std'] = peak_stds
            analyzer_j.df['start_frequency_std'] = start_stds
            analyzer_j.df['end_frequency_std'] = end_stds
            analyzer_j.df['peak_epoch']=peak_epoch
            analyzer_j.df['start_epoch']=start_epoch
            analyzer_j.df['end_epoch']=end_epoch
            analyzers_chunk.append(analyzer_j.df)

                      
        utcstart_chunk=utcstart_chunk+CHUNK_LENGTH
        utcend_chunk=utcend_chunk+CHUNK_LENGTH
        
        #Extend final dataframe with detections from current time chunk
        analyzers.extend(analyzers_chunk)
        new_df = pd.concat(analyzers)
        new_df.to_csv(chunk_pth, index=False)
        

    
    if len(analyzers) == 0:
        print('WARNING: detections dataframe empty')
        final_analyzer_df = []
        
    else:
        final_analyzer_df = pd.concat(analyzers)
        final_analyzer_df.to_csv(detection_pth,index=False)
    #import pdb; pdb.set_trace()    
    return final_analyzer_df
    
 

if __name__ == "__main__":
    main(STARTTIME, ENDTIME)
