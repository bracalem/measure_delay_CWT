#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:41:40 2022

@author: bracalem
"""

import numpy as np
from sem_tools import *
import scipy as sc
from scipy.signal import butter, lfilter
from scipy.fft import fft, fftfreq
from tqdm import tqdm
from scipy import signal


#%% MEASURE OF TIME DELAY WITH WAVELET TANSFORM
# This scipt was made to measure the delay induced by a seismic anomaly on Rayleigh waves. 
# We performed the seismic simulation twice, the first time considering a totally homogeneous medium, then placing a velocity anomaly beneath the seismic array.
# The time delay has been then computed using the cross wavelet transform developed in the work of SHUJUAN MAO and AURELIEN MORDRET.
# https://github.com/shujuanmao/dt-wavelet
# https://github.com/Qhig/cross-wavelet-transform/releases/tag/v1.0.0



#Location of the seismograms of where the medium is homogeneous
hom_loc = ''
#Location of the seismograms of where the medium is heterogeneous
het_loc = ''

#Sampling frequency
fs = 2500
#Number of samples in the seismograms
time_step = 11000

#Samples to cut starting from the end of the seismogram
cut_s_time_f = 0
#Samples to cut starting from the beginning of the seismogram
cut_s_time = 0

#Name of the seismograms to analyze 
#The seismograms need to have the same name, only the location changes
name = ['x']
sta_del_complete = []
for i in tqdm(name):    

    het = np.fromfile(het_loc+ i, dtype=data_type)
    hom = np.fromfile(hom_loc+ + i, dtype=data_type)

    het_split = split_seismograms(het,time_step,cut_s_time,cut_s_time_f)
    hom_split = split_seismograms(hom,time_step,cut_s_time,cut_s_time_f)
    
    #Corner frequency for lowpassing the signal
    freq_max = 15
    sampling = 2500
    
    #Filtering the waveforms, loop on the array receivers
    for p in range(200):
        het_split[p] = butter_bandpass_filter(het_split[p],fs/time_step,freq_max,fs,3)
        hom_split[p] = butter_bandpass_filter(hom_split[p],fs/time_step,freq_max,fs,3)
    
    #Table to store the time delay averaged in different frequency ranges
    v = []
    #Table to store the time delay
    time_distribution = []
    #Table to store the amplitude of the signal
    amp_distribution = []
    for j in tqdm(range(0,len(het_split))):

        len_win = 7700
        #t_rl is a table that identifies the arrival time of the Rayleigh waves, t_rl2[j] identifies the end the phase.
        pad_width = int((len_win-len(hom_split[j][int(t_rl[j]):int(t_rl2[j])]))/2)
        hom_seis = np.pad(hom_split[j][int(t_rl[j]):int(t_rl2[j])], pad_width, mode='linear_ramp')
        het_seis = np.pad(het_split[j][int(t_rl[j]):int(t_rl2[j])], pad_width, mode='linear_ramp')
        
        #Number of frequencies to analyze
        numf = 141
        #Full description of this function, parameters and output can be found in the cited SHUJUAN MAO github repository (README file)
        WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(hom_seis,het_seis, fs, 4, 4, 14, 1, 15, 141) ####PLOT FOR TIME DELAY 4
        #Number of frequencies to average
        range_freq = 10
        k = int((numf-range_freq+1)/10)+1
        for i in range(k):
            time_distribution.append(WXdt[i*10:i*10+range_freq,:])
            amp_distribution.append(WXamp[i*10:i*10+range_freq,:])
     
    #Frequencies to display
    frequencies = np.arange(15,0,-1) #Ellipse settings

    time_distribution_mean = []
    #Average of delay considering the amplitude as a weight
    for l in range(len(time_distribution)):
        ma = np.ma.MaskedArray(time_distribution[l], mask=np.isnan(time_distribution[l]))
        ma2 = np.ma.MaskedArray(amp_distribution[l], mask=np.isnan(time_distribution[l]))
        time_distribution_mean.append(np.ma.average(ma, weights=ma2))


    time_delay_complete = []
    for j in range(k):
        time_delay = []
        for f in range(200):
            time_delay.append(int(np.mean(time_distribution_mean[f*k+j])*1000))
    
#        freq_all = np.linspace(10.2, 0.2, 201)
        v = np.array(time_delay)
        v = v.reshape(10,20)
        v = v.T
        v = v[::-1]
    #Table with computed delays, each component identifies a specific range of frequency. For example time_delay_complete[0] is related to frequencies[0].
    #Each component of time_delay_complete[#] identifies the delay recorded at a specific receiver. For example time_delay_complete[1][10] will be the delay recorded
    #At the receiver number 11 (10+1), in the frequency band given by frequencies[1]
    time_delay_complete.append(v)