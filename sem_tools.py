#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:17:17 2021

@author: bracalem
"""

import numpy as np
import os
import pandas as pd
from scipy.signal import butter, lfilter
from scipy.fft import fft, fftfreq
from scipy import stats

def read_output(location,input_name,time_step):
    os.chdir(location)
    data_type = np.dtype('float32').newbyteorder('<')
    data_in = np.fromfile(input_name, dtype=data_type)
    return data_in

def split_seismograms(data_in,time_step,cutend=0,cut_start=0):
    trks = []
    for sensor_n in range(int(len(data_in)/time_step)):
        trks.append(data_in[sensor_n*time_step+cut_start:(sensor_n+1)*time_step-cutend])
    return trks

def compute_angles(input_file,angle0_xy,angle0_xz,angle0_yz,xc,yc,zc,degree= True):
    data = pd.read_csv(input_file, header = None,delimiter =' ')
    x = np.array(data.iloc[1:,1])
    y = np.array(data.iloc[1:,2])
    z = np.array(data.iloc[1:,0])
    angle_xy = (np.arctan((y-yc)/(x-xc))-angle0_xy)-angle0_xy#/np.pi)*180-angle0_xy
    angle_xz = (np.arctan((z-zc)/(x-xc))-angle0_xz)-angle0_xz#/np.pi)*180-angle0_xz
    angle_yz = (np.arctan((z-zc)/(y-yc))-angle0_yz)-angle0_yz#/np.pi)*180-angle0_yz
    return angle_xy,angle_xz,angle_yz

def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def time_derivative(trk,time_step,out_name):    
    data_type = np.dtype('float32').newbyteorder('<')
    trk = np.fromfile(trk, dtype=data_type)
    trks = split_seismograms(trk,time_step,cutend=0,cut_start=0)
    out = []
    for sensor_n in range(len(trks)):
        deriv = trks[sensor_n][1:] - trks[sensor_n][:-1]
        deriv = np.append(deriv,deriv[-1])
        out.append(deriv)
    out = np.array(out)
    np.ndarray.tofile(out, format='%',file=out_name)
    return out

def fft_(trk,time_step):
    N = len(trk)
    trk = np.pad(trk, 1000, mode='linear_ramp')
    x = np.linspace(0.0, N*time_step, N, endpoint=False)
    yf = fft(trk)
    xf = fftfreq(N, time_step)[:N//2]
    yf = 2.0/N * np.abs(yf[0:N//2])
    return xf,yf