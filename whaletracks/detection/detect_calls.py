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
import matplotlib.cm as cm
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

def defaultScaleFunction(Sxx):
    vmin=np.median(10*np.log10(Sxx))+0*np.std(10*np.log10(Sxx)) 
    vmax=np.median(10*np.log10(Sxx))+2*np.std(10*np.log10(Sxx)) 
    return vmin, vmax

def defaultKernelLims(f0,f1,bdwdth):
    ker_1=f1-4*bdwdth
    ker_2=f0+4*bdwdth
    ker_min=np.min([ker_1,ker_2])
    ker_max=np.max([ker_1,ker_2])
    return ker_min, ker_max

def plotwav(samp, data, filt_type='bandpass', filt_freqlim=[8, 17], 
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
        plt.pcolormesh(t, f, 10*np.log10(Sxx), cmap=cmap, norm=norm)    
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.ylim(ylim)
    return [f, t, Sxx]


#Define function to build 2-D kernel linear sweep to cross-correlate with spectrograms
def buildkernel(f0, f1, bdwdth, dur, f, t, samp, plotflag=True,kernel_lims=defaultKernelLims):
    
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
    :param tuple kernel_lims: Tuple of minimum kernel range and maximum kernel range
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
    ker_min, ker_max=kernel_lims(f0,f1,bdwdth)
    
    for j in range(np.size(tvec)):
        #calculate hat function that is centered on linearly decresing
        #frequency values for each time in tvec
        x = fvec-(f0+(t[j]/dur)*(f1-f0))
        Kval = (1-np.square(x)/(bdwdth*bdwdth))*np.exp(-np.square(x)/(2*(bdwdth*bdwdth)))
        Kdist[j] = Kval #store hat function values in preallocated array
                                                           
    BlueKernel_full = np.transpose(Kdist) #transpose preallocated array to be plotted vs. tvec and fvec
    freq_inds=np.where(np.logical_and(fvec>=ker_min, fvec<=ker_max))
    
    fvec_sub=fvec[freq_inds]
    BlueKernel=BlueKernel_full[freq_inds,:][0]
    
    
    if plotflag == True:
        plt.figure(PLT_KERNEL)
        plt.pcolormesh(tvec, fvec_sub, BlueKernel) #show plot of kernel
        plt.axis([0, dur, np.min(fvec), np.max(fvec)])
        plt.gca().set_aspect('equal')
        plt.colorbar()
        plt.title('Blue whale B-call kernel')

    return [tvec, fvec_sub, BlueKernel, freq_inds]




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
    
    
    CorrVal_scale=CorrVal*1/(np.median(Sxx)*np.size(tvec)*np.size(fvec))
    #CorrVal_scale=CorrVal*1/(np.median(Sxx))
    #CorrVal_scale=CorrVal
    CorrVal_scale[0]=0
    CorrVal_scale[-1]=0
    neg_ind=CorrVal_scale<0
    CorrVal_scale[neg_ind]=0
    t_scale=t[int(len(tvec)/2)-1:-math.ceil(len(tvec)/2)]
#Visualize spectrogram and detection scores of of data  
    
    if plotflag==True:
        
        t1=min(t)
        t2=max(t)
#plot timeseries on upper axis
        plt.figure(PLT_SCORE, figsize=(9, 3))
        fig, (ax0, ax1) = plt.subplots(nrows=2)
        ax0.plot(t_scale,CorrVal_scale) #plot normalized detection scores as a time series.
        
        ax0.set_xlim([t1, t2]) #look at only positive values
        ax0.set_ylim([0, np.max(CorrVal_scale)])
        ax0.set_xlabel('Seconds')
        ax0.set_ylabel('Detection score')
        ax0.set_title('Spectrogram and detection scores of test data')

#plot spectrogram on lower axis
        cmap = plt.get_cmap('magma')
        vmin, vmax = scale_func(Sxx)
        norm = color.Normalize(vmin=vmin, vmax=vmax)
        #plt.subplot(212)
        im = ax1.pcolormesh(t, f, 10*np.log10(Sxx), cmap=cmap,norm=norm) 
        fig.colorbar(im, ax=ax1,orientation='horizontal')
        ax1.set_xlim([t1, t2]) #look at spectrogram segment between given time boundaries
        ax1.set_ylim([12, 18])
        ax1.set_ylabel('Frequency [Hz]')
        ax1.set_xticks([])
        #ax1.set_xlabel('Time [seconds]')
        fig.tight_layout()
        fig
    return  [t_scale, CorrVal_scale] 

    #plt.savefig('Spectrogram_scores.png')



def xcorr_log(t,f,Sxx,tvec,fvec,BlueKernel,plotflag=True,scale_func=defaultScaleFunction):
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
    Sxx_log1=10*np.log10(Sxx)
    Sxx_log=Sxx_log1-np.min(Sxx_log1)
#Cross-correlate kernel with spectrogram
    ind1 = 0
    CorrVal = np.zeros(np.size(t) - (len(tvec)-1)) #preallocate array for correlation values
    corrchunk= np.zeros((np.size(fvec), np.size(tvec))) #preallocate array for element-wise multiplication

    while ind1-1+np.size(tvec) < np.size(t):
        ind2 = ind1 + np.size(tvec) #indices of spectrogram subset to multiply
        for indF in range(np.size(fvec)-1):
            corrchunk[indF] = Sxx_log[indF][ind1:ind2] #grab spectrogram subset for multiplication
        
        CorrVal[ind1] = np.sum(np.multiply(BlueKernel, corrchunk)) #save cross-correlation value for each frame
        ind1 += 1
    
    
    CorrVal_scale=CorrVal*1/(np.median(Sxx_log)*np.size(tvec))
    #CorrVal_scale=CorrVal*1/(np.median(Sxx))
    #CorrVal_scale=CorrVal
    CorrVal_scale[0]=0
    CorrVal_scale[-1]=0
    neg_ind=CorrVal_scale<0
    CorrVal_scale[neg_ind]=0
    t_scale=t[int(len(tvec)/2)-1:-math.ceil(len(tvec)/2)]
#Visualize spectrogram and detection scores of of data  
    
    if plotflag==True:
        
        t1=min(t)
        t2=max(t)
#plot timeseries on upper axis
        plt.figure(PLT_SCORE, figsize=(9, 3))
        fig, (ax0, ax1) = plt.subplots(nrows=2)
        ax0.plot(t_scale,CorrVal_scale) #plot normalized detection scores as a time series.
        
        ax0.set_xlim([t1, t2]) #look at only positive values
        ax0.set_ylim([0, np.max(CorrVal_scale)])
        ax0.set_xlabel('Seconds')
        ax0.set_ylabel('Detection score')
        ax0.set_title('Spectrogram and detection scores of test data')

#plot spectrogram on lower axis
        cmap = plt.get_cmap('magma')
        vmin=np.median(Sxx_log)+2*np.std(Sxx_log)
        vmax=np.median(Sxx_log)
        #vmin, vmax = scale_func(Sxx_log)
        norm = color.Normalize(vmin=vmin, vmax=vmax)
        #plt.subplot(212)
        im = ax1.pcolormesh(t, f, Sxx_log, cmap=cmap,norm=norm) 
        fig.colorbar(im, ax=ax1,orientation='horizontal')
        ax1.set_xlim([t1, t2]) #look at spectrogram segment between given time boundaries
        ax1.set_ylim([12, 18])
        ax1.set_ylabel('Frequency [Hz]')
        ax1.set_xticks([])
        #ax1.set_xlabel('Time [seconds]')
        fig.tight_layout()
        fig
    return  [t_scale, CorrVal_scale] 

    #plt.savefig('Spectrogram_scores.png')

