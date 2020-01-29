#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:48:49 2020

@author: wader
"""

#Needs inputs of data (just timeseries) sample rate
#The spectrogram work is done by scipy.signal.spectrogram
#Check out the man pages and figure out what kind of windows and overlap you'd like for the best looking spectrogram.

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
import math

# Constants
PLT_TIMESERIES = 1
PLT_SPECTROGRAM = 2
PLT_KERNEL = 3
PLT_SCORE = 4
FIGSIZE = [9, 3]

FILTER_OFFSET = 10

# Helper functions
# TODO: change this to median?
def defaultScaleFunction(Sxx):
    vmin=np.median(20*np.log10(Sxx))+0*np.std(20*np.log10(Sxx)) 
    vmax=np.median(20*np.log10(Sxx))+3*np.std(20*np.log10(Sxx)) 
    return vmin, vmax

def plotwav(samp, data, filt_type='bandpass', filt_freqlim=[8, 20], 
            filt_order=4, window_size=4, overlap=.95, window_type='hann',
            plotflag=True, scale_func=defaultScaleFunction,ylim=[12, 18]): 
    """
    Calculate spectogram and plot.
    :param float samp: sampling rate
    :param numpy.array data: data to process
    :param string filt_type: filter type (highpass,lowpass,bandpass etc.)
    :param tuple filt_freqlim: frequency limits of filter
    :param int filt_order: order of filter 
    :param float window_size: seconds spectrogram window size
    :param float overlap: ratio of overlap of spectrogram window
    :param string window_type: spectrogram window type
    :param bool plotflag: If True, makes plots. If False, no plot.
    :param function scal_func: Single argument Sxx, returns tuple vmin and vmax floats
    :param tuple ylim: lower and uppet bounds of frequency for spectrogram plot
    :param numpy.array data: vector of min and max spectrogram frequency 
    :return numpy.array, numpy.array, 2-d numpy.array: 
        vector of frequency, vector of seconds, matrix of power
    Key variables
      f - frequency
      t - time (seconds)
      Sxx - Spectrogram of amplitudes
    """
    #filter data to spectral bands where B-call is
    [b, a] = sig.butter(filt_order, np.array(filt_freqlim)/samp, filt_type, 'ba') 
    filtered_data = sig.filtfilt(b, a, data)

    
    datalength = data.size
    times = (np.arange(datalength)/samp) 

    #plot timeseries on upper axis
    if plotflag == True:
        plt.figure(PLT_TIMESERIES, figsize=FIGSIZE)
        plt.subplot(211)
        plt.plot(times[FILTER_OFFSET:],filtered_data[FILTER_OFFSET:])
        plt.axis([min(times), max(times), min(filtered_data[FILTER_OFFSET:]), 
                  max(filtered_data[FILTER_OFFSET:])])
        plt.xlabel('Seconds')
        plt.ylabel('Amplitude')

    #plot spectrogram on lower axis
    [f, t, Sxx] = sig.spectrogram(filtered_data, int(samp), 
    window_type,int(samp*window_size),int(samp*window_size*overlap))
    if plotflag == True:
        cmap = plt.get_cmap('magma')
        vmin, vmax = scale_func(Sxx)
        norm = color.Normalize(vmin=vmin, vmax=vmax)
        plt.subplot(212)
        plt.pcolormesh(t, f, 20*np.log10(Sxx), cmap=cmap, norm=norm)    
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.ylim(ylim)
    return [f, t, Sxx]


#Define function to build 2-D kernel linear sweep to cross-correlate with spectrograms
def buildkernel(f0, f1, bdwdth, dur, f, t, samp, plotflag=True):
    
    """
    Calculate kernel and plot
    :param float f0: starting frequency
    :param float f1: ending frequency
    :param float bdwidth: frequency width of call
    :param float dur: call length (seconds)
    :param np.array f: vector of frequencies returned from plotwav
    :param np.array t: vector of times returned from plotwav
    :param float samp: sample rate
    :param bool plotflag: If True, plots kernel. If False, no plot.
    :return numpy.array, numpy.array, 2-d numpy.array: 
        vector of kernel times, vector of kernel frequencies, matrix of kernel values
    Key variables
      tvec - kernel times (seconds)
      fvec - kernel frequencies
      BlueKernel - Matrix of kernel values
    """
    
    tvec = np.linspace(0,dur,np.size(np.nonzero(t < dur))) #define kernel as length dur
    fvec = f #define frequency span of kernel to match spectrogram
    Kdist = np.zeros((np.size(tvec), np.size(fvec))) #preallocate space for kernel values
    
    for j in range(np.size(tvec)):
        #calculate hat function that is centered on linearly decresing
        #frequency values for each time in tvec
        x = fvec-(f0+(t[j]/dur)*(f1-f0))
        Kval = (1-np.square(x)/(bdwdth*bdwdth))*np.exp(-np.square(x)/(2*(bdwdth*bdwdth)))
        Kdist[j] = Kval #store hat function values in preallocated array
                                                           
    BlueKernel = np.transpose(Kdist) #transpose preallocated array to be plotted vs. tvec and fvec
    
    if plotflag == True:
        plt.figure(PLT_KERNEL)
        plt.pcolormesh(tvec, fvec, BlueKernel) #show plot of kernel
        plt.axis([0, dur, np.min(fvec), np.max(fvec)])
        plt.gca().set_aspect('equal')
        plt.colorbar()
        plt.title('Blue whale B-call kernel')

    return [tvec, fvec, BlueKernel]




def xcorr(t,f,Sxx,tvec,fvec,BlueKernel,plotflag=True,scale_func=defaultScaleFunction):
    """
    Cross-correlate kernel with spectrogram and plot score
    :param np.array f: vector of frequencies returned from plotwav
    :param np.array t: vector of times returned from plotwav
    :param np.array Sxx: 2-D array of spectrogram amplitudes
    :param np.array tvec: vector of times of kernel
    :param np.array fvec: vector of frequencies of kernel
    :param np.array BlueKernel: 2-D array of kernel amplitudes
    :plotflag boolean plotflag: Boolean. If True, plots result. If False, no plot.
    :param function scal_func: Single argument Sxx, returns tuple vmin and vmax floats
    :return numpy.array, numpy.array:
        vector of correlation times, vector of correlation values
    Key variables
        t_scale - correlation times (seconds)
        CorrVal_scale - correlation values
    """

#Cross-correlate kernel with spectrogram
    ind1 = 0
    CorrVal = np.zeros(np.size(t) - (len(tvec)-1)) #preallocate array for correlation values
    corrchunk= np.zeros((np.size(fvec), np.size(tvec))) #preallocate array for element-wise multiplication

    while ind1-1+np.size(tvec) < np.size(t):
        ind2 = ind1 + np.size(tvec) #indices of spectrogram subset to multiply
        for indF in range(np.size(fvec)-1):
            corrchunk[indF] = Sxx[indF][ind1:ind2] #grab spectrogram subset for multiplication
        
        CorrVal[ind1] = np.sum(np.multiply(BlueKernel, corrchunk)) #save cross-correlation value for each frame
        ind1 += 1
   
    CorrVal_scale=CorrVal*1/(np.median(Sxx)*np.size(t))
    t_scale=t[int(len(tvec)/2)-1:-math.ceil(len(tvec)/2)]
#Visualize spectrogram and detection scores of of data
        
    if plotflag==True:
        t1=min(t)
        t2=max(t)
#plot timeseries on upper axis
        plt.figure(PLT_SCORE, figsize=(9, 3))
        plt.subplot(211)
        plt.plot(t_scale,CorrVal_scale) #plot normalized detection scores as a time series.
        plt.axis([t1, t2, 0, np.max(CorrVal_scale)]) #look at only positive values
        plt.xlabel('Seconds')
        plt.ylabel('Detection score')
        plt.title('Spectrogram and detection scores of test data')

#plot spectrogram on lower axis
        cmap = plt.get_cmap('magma')
        vmin, vmax = scale_func(Sxx)
        norm = color.Normalize(vmin=vmin, vmax=vmax)
        plt.subplot(212)
        plt.pcolormesh(t, f, 20*np.log10(Sxx), cmap=cmap,norm=norm)   
        plt.axis([t1, t2, 12, 18]) #look at spectrogram segment between given time boundaries
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [seconds]')
        return  [t_scale, CorrVal_scale] 

    #plt.savefig('Spectrogram_scores.png')
