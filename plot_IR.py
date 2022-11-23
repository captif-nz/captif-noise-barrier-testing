import numpy as np
from scipy.signal import correlate
import scipy.fftpack
import sys
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from bokeh.plotting import figure, output_file, show, gridplot
import os
# Import EN1793 specific libraries
import EN_1793_funcs as en1793
import octave_filters as octave

import itertools
import random
import string
import csv


################ Input files to plot ########################
# Impulse response
ff_folder = "BarrierTesting/TeAtatu_Panel_TL_FF"
ff_filename = os.path.join(ff_folder, 'IR')
############################################################

with open(ff_filename, 'rb') as file1:
    raw_data = pickle.load(file1)
avg_IR_raw = raw_data["averagedImpulseResponse"]
t_array = np.arange(len(avg_IR_raw['mic1']))*1./51200

# Flatten curve
IR_avg = {}
mics = list(avg_IR_raw.keys())
for key in mics:
    IR_avg[key] = en1793.remove_drift(data=avg_IR_raw[key], time_array=t_array)

lowest_band = 100
highest_band = 8000
print(IR_avg)
print(t_array)

key = 'mic1'
p = figure(title="Impulse Response", plot_width=800, plot_height=400)
p.line(x=t_array, y=IR_avg[key], legend=key)
show(p)

# en1793.simple_plot(X=t_array, Y=IR_avg["mic1"], labels=["mic1"])