import numpy as np
from scipy.signal import correlate
import scipy.fftpack
import sys
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.io import export_png
import os
# Import EN1793 specific libraries
import EN_1793_funcs as en1793
import octave_filters as octave

import itertools
import random
import string
import csv


################ Input files to process ########################
# Barrier impulse response
ff_folder = "BarrierTesting/MaioroSt_Post_TL_FF"
ff_filename = os.path.join(ff_folder, 'IR')
# Free-feild impulse response
wall_folder = "BarrierTesting/MaioroSt_Post_TL"
wall_filename = os.path.join(wall_folder, 'IR')
# Location to store results
dest_folder = "BarrierTesting/MaioroSt_Post_TL_Results"
############################################################
plot_title = 'TeAtatu - Post'

with open(ff_filename, 'rb') as file1, open(wall_filename, 'rb') as file2:
    ff_data = pickle.load(file1)
    wall_data = pickle.load(file2)
ff_avg_raw = ff_data["averagedImpulseResponse"]
wall_avg_raw = wall_data["averagedImpulseResponse"]
t_array = np.arange(len(ff_avg_raw['mic1']))*1./51200

# Flatten curves
ff_avg = {}
wall_avg = {}
for ff_key, wall_key in list(zip(ff_avg_raw.keys(), wall_avg_raw.keys())):
    ff_avg[ff_key] = en1793.remove_drift(data=ff_avg_raw[ff_key], time_array=t_array)
    wall_avg[wall_key] = en1793.remove_drift(data=wall_avg_raw[wall_key], time_array=t_array)

lowest_band = 100
highest_band = 8000

####################### Display impulse response files ##################
# import webbrowser
# ff_url = os.path.join(ff_folder, 'ImpulseResponse.html')
# print("Displaying: {}".format(ff_url))
# webbrowser.open(ff_url,new=1)
# wall_url = os.path.join(wall_folder, 'ImpulseResponse.html')
# print("Displaying: {}".format(wall_url))
# webbrowser.open(wall_url,new=1)
#######################################################################

################################ Plot IR for each mic ####################
print('Plotting Measurement Results')
titles = [["Mic 9 Average IR","Mic 8 Average IR","Mic 7 Average IR"],
    ["Mic 6 Average IR","Mic 5 Average IR","Mic 4 Average IR"],
    ["Mic 3 Average IR","Mic 2 Average IR","Mic 1 Average IR"]]
x_data = t_array
y_data_1 = [[ff_avg['mic9'], ff_avg['mic8'], ff_avg['mic7']],
    [ff_avg['mic6'], ff_avg['mic5'], ff_avg['mic4']],
    [ff_avg['mic3'], ff_avg['mic2'], ff_avg['mic1']]]
y_data_2 = [[wall_avg['mic9'], wall_avg['mic8'], wall_avg['mic7']],
    [wall_avg['mic6'], wall_avg['mic5'], wall_avg['mic4']],
    [wall_avg['mic3'], wall_avg['mic2'], wall_avg['mic1']]]
legends = ['Free-feild', 'Wall']
output_file(os.path.join(dest_folder, "allIR.html"))
# p_f = en1793.mic_grid_plot(x_data=x_data, y_data_1=y_data_1, y_data_2=y_data_2, titles=titles,
#     output_filename=os.path.join(dest_folder, "allIR.html"), legends=legends,
#     x_labels='time [s]', y_labels='Amplitude', l_style=['solid','dashed'],
#     y_range=(-1, 2))
# show(p_f)
#########################################################################

print("Calculating sound insulation of wall from free feild and wall measurement")
print("HAVE YOU SELECTED THE CORRECT FILES IN LINES 14 - 21?")
print("\nUsing free-feild measurement: {}".format(ff_filename))
print("Using Wall measurement: {}".format(wall_filename))
# ad_offset_FF = float(input("Adienne offset to be used for FF (default is 8ms):") or "8")
ad_offset_FF = 8
# ad_end_TL = float(input("End of Adienne window TL (default is 9ms):") or "9")
ad_end_TL = 9
# ad_offset_TL = float(input("Start of Adienne window TL (default is 4ms):") or "4")
ad_offset_TL = 4
ad_length_TL = ad_end_TL - ad_offset_TL
noise_length_TL = ad_offset_TL
print("Resulting Adrienne window length: {}ms".format(ad_length_TL))

if ff_data["SR"] == wall_data["SR"]:
    print("Sample rates match ({}) - proceeding".format(ff_data["SR"]))
    Fs = ff_data["SR"]
else:
    print("SR is different between measurements! Stopping!")
    sys.exit("SR mismatch")



# soundInsulation = []
ratio_store = []
SNR_store = []
total_store = []
for ff_key, wall_key in list(zip(ff_avg.keys(), wall_avg.keys())):
    print("Processing: {}, {}".format(ff_key, wall_key))
    # Load relevant data
    ff_ir= ff_avg[ff_key]
    wall_ir= wall_avg[wall_key]
    if len(ff_ir) != len(wall_ir):
        print("Length of impulse responses do not match!")
        exit()
    tt = np.array(range(len(ff_ir)))*1./51200
    # ff_ir = en1793.remove_drift(data=ff_ir_raw, time_array=tt)
    # wall_ir = en1793.remove_drift(data=wall_ir_raw, time_array=tt)
    # en1793.simple_plot(X=[tt, tt], Y=[ff_ir_raw, ff_ir], labels=["Raw", "Flattened"])
    # en1793.simple_plot(X=[tt, tt], Y=[wall_ir_raw, wall_ir], labels=["Raw", "Flattened"])
    # Apply Adrienne window
    ff_windowed = en1793.apply_adrienne(ff_ir, offset=ad_offset_FF, SR=Fs)
    wall_windowed = en1793.apply_adrienne(wall_ir, offset=ad_offset_TL, SR=Fs, windowLength=ad_length_TL)
    noise_windowed = en1793.apply_adrienne(wall_ir, offset=0.0, SR=Fs, windowLength=noise_length_TL)
    # en1793.simple_plot(X=[tt, tt], Y=[ff_windowed, wall_windowed], labels=["Free feild", "Wall"])
    # exit()
    # Get start and end of 1/3rd octave bands
    f_center, f_low, f_high = octave.third_oct_freq(low_freq=100, high_freq=5000)

    # Generate fft of windowed data
    dt = 1.0/51200.0
    ff_fft = abs(scipy.fftpack.fft(ff_windowed))**2
    ff_xf = np.linspace(0.0, 1.0/dt, len(ff_windowed))
    wall_fft = abs(scipy.fftpack.fft(wall_windowed))**2
    wall_xf = np.linspace(0.0, 1.0/dt, len(wall_windowed))
    noise_fft = abs(scipy.fftpack.fft(noise_windowed))**2
    noise_xf = np.linspace(0.0, 1.0/dt, len(noise_windowed))
    df = ff_xf[1]
    # Generate 1/3rd octave bins and bin data into third octaves (this is a "square" window)
    bins = np.concatenate((f_low, f_high[-3:]), axis=None)/df
    ff_sum = stats.binned_statistic(x=np.arange(len(ff_fft)), values=ff_fft, statistic='sum', bins=bins)[0]
    wall_sum = stats.binned_statistic(x=np.arange(len(wall_fft)), values=wall_fft, statistic='sum', bins=bins)[0]
    noise_sum = stats.binned_statistic(x=np.arange(len(noise_fft)), values=noise_fft, statistic='sum', bins=bins)[0]
    # Calculate ratio beteen third octave levels transmitted through wall against the free-feild measurements
    ratio_wall_ff = np.divide(wall_sum, ff_sum)
    ratio_store.append(ratio_wall_ff)
    # SNR calculation
    # SNR_store.append(np.divide(wall_sum, noise_sum))
    # print(SNR_store)
    ratio_sum = np.divide(np.sum(wall_sum), np.sum(ff_sum))
    print(ratio_sum)
    total_store.append(ratio_sum)

soundInsulation = -10.0*np.log10(np.sum(total_store, axis=0)/9.0)
print(soundInsulation)
exit()
# Calculate the sound insulation in each 1/3rd octave band
soundInsulation = -10.0*np.log10(np.sum(ratio_store, axis=0)/9.0)
SNR = 10.0*np.log10(np.sum(SNR_store, axis=0)/9.0)

# Display results
print('FF Sum -----')
print(ff_sum)
print('Wall Sum -----')
print(wall_sum)
print('SI -----')
print(soundInsulation[1:-1])
print('Bins Hz -----')
print(bins*df)
print("Centre Frequencies (Hz)")
print(f_center)
print('SNR -----')
print(SNR)
# Write data to csv to plot against other barriers
with open(os.path.join(dest_folder, "results.csv"), 'w') as writeFile:
    writer = csv.writer(writeFile, delimiter=',')
    writer.writerow(["FF levels (dB)"])
    writer.writerow(ff_sum)
    writer.writerow(["Wall levels (dB)"])
    writer.writerow(wall_sum)
    writer.writerow(["SI (dB)"])
    writer.writerow(soundInsulation)
    writer.writerow(["Frequency Bins (Hz)"])
    writer.writerow(bins*df)
    writer.writerow(["Centre Frequency (Hz)"])
    writer.writerow(f_center)
    writer.writerow(["Signal to Noise (dB)"])
    writer.writerow(SNR)
writeFile.close()



# Calculate single number rating of Sound Insultion

# Plot sound insulation
# output_file(os.path.join(dest_folder, "SoundInsulation.html"))
# p4 = figure(plot_width=800, plot_height=400, title="Sound Insulation", x_axis_type="log", y_range=(0,50))
# p4.line(f_center, soundInsulation[1:-1], line_width=2)
# show(p4)

# Plot sound insulation
CF = [125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000]
output_file(os.path.join(dest_folder, "SNR.html"))
p4 = figure(plot_width=800, plot_height=400, title="SNR: {}".format(plot_title), x_axis_type="log", y_range=(0,50))
p4.line(f_center, SNR[1:-1], line_width=2)
p4.xaxis.axis_label = 'One Third Octace Band Centre Frequency (Hz)'
p4.yaxis.axis_label = 'Signal to Noise Ratio'
p4.xaxis.ticker = CF
p4.xgrid.grid_line_color=None
p4.xaxis.minor_tick_line_color = None
show(p4)
export_png(p4, 'BarrierTesting/MaioroSt_Panel_TL_Results/SNR.png')

