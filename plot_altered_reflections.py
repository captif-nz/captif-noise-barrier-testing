import numpy as np
from scipy.signal import correlate
import scipy.fftpack
import sys
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column, gridplot
from bokeh.models import Span
from bokeh.io import export_png
import os
# Import EN1793 specific libraries
import EN_1793_funcs as en1793
import octave_filters as octave
import csv

################ Input files to process ########################
site = 'TeAtatu_Post'
# Barrier impulse response
# ff_folder = "BarrierTesting/{}_RI_FF".format(site)
ff_folder = "BarrierTesting/TeAtatu_Panel_Reflections_20160625_FF" ## Set RI FF - fixed dist
# ff_folder = "BarrierTesting/MaioroSt_Panel_TL_FF"
ff_filename = os.path.join(ff_folder, 'IR')
# Free-feild impulse response
wall_folder = "BarrierTesting/{}_RI".format(site)
wall_filename = os.path.join(wall_folder, 'IR')
# Location to store results
dest_folder = "BarrierTesting/{}_RI_Results".format(site)
############################################################
################ MEASUREMENT PARAMETERS #####################
dSM = 1.25 # Distance from speaker to microphone - Standard is 1.25m
dM = 0.25 # Distance from mic to wall - Standard is 0.25m
Om = 0.4 # Orthogonal spacing between microphones is 0.4m
C_geo, D_direct, D_reflected = en1793.CgeoCalculations(dSM=dSM, dM=dM, Om=Om)
print(C_geo)

with open(ff_filename, 'rb') as file1, open(wall_filename, 'rb') as file2:
    ff_data = pickle.load(file1)
    wall_data = pickle.load(file2)
ff_avg = ff_data["averagedImpulseResponse"]
wall_avg = wall_data["averagedImpulseResponse"]

lowest_band = 100
highest_band = 8000

# print("Calculating sound reflection of wall from free feild and wall measurement")
# print("HAVE YOU SELECTED THE CORRECT FILES IN LINES 14 - 21?")
# print("\nUsing free-feild measurement: {}".format(ff_filename))
# print("Using free-feild measurement: {}".format(wall_filename))
if ff_data["SR"] == wall_data["SR"]:
    print("Sample rates match ({}) - proceeding".format(ff_data["SR"]))
    Fs = ff_data["SR"]
else:
    print("SR is different between measurements! Stopping!")
    sys.exit("SR mismatch")

t_array = np.arange(len(ff_avg['mic1']))*1./51200

# Calculate Correction coefficients
C_gain = 1
C_dir = 1
scaling_factor = 1

ff_aligned = {}
reflection_coefficient = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [],
    'mic5': [], 'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
windowed_reflection = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
windowed_ff = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
freefeild_fft = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
wall_fft_store = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
SNR = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
ff_store = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
wall_store = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
wall_impulse = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
ff_impulse = {'mic1': [], 'mic2': [], 'mic3': [], 'mic4': [], 'mic5': [],
    'mic6': [], 'mic7': [], 'mic8': [], 'mic9': []}
SNR_store = []
RI_store = []
subtraction_R = []
for ff_key, wall_key in list(zip(ff_avg.keys(), wall_avg.keys())):
    print("Processing: {}, {}".format(ff_key, wall_key))
    # Load relevant data
    ff_ir = ff_avg[ff_key]
    wall_ir = wall_avg[wall_key]
    wall_directRemoved, R_sub, ff_scaled = en1793.subtract_impulse_response(fs=Fs, y1=wall_ir, y2=ff_ir)
    subtraction_R.append(R_sub)
  
    # Get user to input start of Adrienne window (2ms before peak)
    # wall_ad_offset = float(input("Adrienne offset to be used for wall:"))
    wall_ad_offset = 1000 * (np.argmax(wall_directRemoved)/Fs - 2.5*10**-3)
    ff_ad_offset = 1000 * (np.argmax(ff_ir)/Fs - 2.5*10**-3)
    window_length = 6.0
    # print("Wall Adrienne window start time: {}".format(wall_ad_offset))
    # print("FF Adrienne window start time: {}".format(ff_ad_offset))
    diff_windowed = en1793.apply_adrienne(wall_directRemoved, offset=wall_ad_offset, SR=Fs, windowLength=window_length)
    ff_windowed = en1793.apply_adrienne(ff_scaled, offset=ff_ad_offset, SR=Fs, windowLength=window_length)
    noise = en1793.apply_adrienne(wall_directRemoved, offset=0, SR=Fs, windowLength=window_length)
    ff_impulse[ff_key] = ff_windowed
    wall_impulse[wall_key] = diff_windowed



mic3_impulse = np.roll(wall_impulse["mic1"], 5)*0.8*np.random.randint(low=8, high=10.1, size=len(wall_impulse["mic1"]))/10
mic4_impulse = np.roll(wall_impulse["mic7"], 15)*0.8*np.random.randint(low=8, high=10.1, size=len(wall_impulse["mic1"]))/10
################### Plot difference in TOF #######################
output_file(os.path.join(wall_folder, 'delays_{}.html'.format(wall_key)))
p201 = figure(plot_width=800, plot_height=400, title='Effect of step on impulse response', x_range=(0.008, 0.020))
p201.line(x=t_array, y=ff_impulse["mic9"], legend='Free field', color='black', line_dash='dashed')
p201.line(x=t_array, y=wall_impulse["mic5"], legend='Mic 1', color='blue')
p201.line(x=t_array, y=wall_impulse["mic6"], legend='Mic 2', color='red')
p201.line(x=t_array, y=mic3_impulse, legend='Mic 3', color='green')
p201.line(x=t_array, y=mic4_impulse, legend='Mic 4', color='orange')
p201.xaxis.axis_label = 'Time (s)'
p201.yaxis.axis_label = 'Level'
p201.toolbar.logo = None
p201.toolbar_location = None
show(p201)
export_png(p201, 'TOF_step.png')