# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 20:23:50 2018

@author: GGPC

Octave and third octaves
"""
import numpy as np
from scipy.signal import butter, lfilter


def octave_freq(low_freq, high_freq):
    f_center = []
    f_low = []
    f_high = []
    for band_num in np.arange(-6, 10):
        CF = 10**3 * 2.**band_num
        HF = CF * (2.**0.5)
        LF = CF / (2.**0.5)
        if CF >= low_freq and CF <= high_freq:
            f_center.append(CF)
            f_high.append(HF)
            f_low.append(LF)
    return f_center, f_low, f_high


def third_oct_freq(low_freq, high_freq):
    f_center = []
    f_low = []
    f_high = []
    for band_num in np.arange(-20, 15):
        CF = 10**3 * 2.**(band_num/3)
        HF = CF * (2.**0.5)
        LF = CF / (2.**0.5)
        if CF >= low_freq and CF <= high_freq:
            f_center.append(CF)
            f_high.append(HF)
            f_low.append(LF)
    return f_center, f_low, f_high


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5*fs
    low = lowcut/nyq
    high = highcut/nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def average_level(pressure_data, P_ref):
    mean_pressure = np.mean(pressure_data)
    return 10*np.log10(mean_pressure/P_ref)

def third_oct_filter(data, low_freq=100, high_freq=8000, fs=51.2*10**3):
    f_center, f_low, f_high = third_oct_freq(low_freq, high_freq)
    output = []
    for f0, f1 in list(zip(f_low, f_high)):
        output.append(butter_bandpass_filter(data=data, lowcut=f0, highcut=f1, fs=fs, order=5))
    return output
