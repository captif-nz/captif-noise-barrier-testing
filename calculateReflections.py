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
ff_folder = "BarrierTesting/TeAtatu_Post_Reflections_20160625_FF" ## Set RI FF - fixed dist
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
    # Check the R_sub value
    # if R_sub < 10:
    #     print("R_sub: {}".format(R_sub))
    #     print("This indicates that the subtraction of the free-feild from the reflected is not acceptable")
    # else:
    #     print("R_sub: {}".format(R_sub))
    #     print("This indicates that the subtraction of the free-feild from the reflected is acceptable")

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

    diff_scaling = max(abs(diff_windowed))
    ff_scaling = max(abs(ff_windowed))
    noise_scaling = max(abs(noise))
    windowed_reflection[wall_key] = diff_windowed
    windowed_ff[ff_key] = ff_windowed

    # Get start and end of 1/3rd octave bands
    f_center, f_low, f_high = octave.third_oct_freq(low_freq=100, high_freq=5000)
    # print(f_center)
    # Generate fft of windowed data
    dt = 1.0/51200.0
    # start = int(0.8*ff_ad_offset*Fs/1000)
    start = 0
    stop =len(ff_windowed) # int((wall_ad_offset*Fs + 2*window_length*Fs)/1000)
    ff_clipped = ff_windowed[start:stop]
    diff_clipped = diff_windowed[start:stop]
    noise_clipped = noise[start:stop]
    ff_fft = abs(scipy.fftpack.fft(ff_clipped) * ff_scaling)**2
    # ff_freqs = scipy.fftpack.fftfreq(len(ff_clipped)) * Fs
    ff_xf = np.linspace(0.0, 1.0/dt, len(ff_clipped))
    wall_fft = abs(scipy.fftpack.fft(diff_clipped) * diff_scaling)**2
    wall_xf = np.linspace(0.0, 1.0/dt, len(diff_clipped))
    noise_fft= abs(scipy.fftpack.fft(noise_clipped) * noise_scaling)**2
    noise_xf = np.linspace(0.0, 1.0/dt, len(noise_clipped))
    df = ff_xf[1]
    ##################

    # Generate 1/3rd octave bins and bin data into third octaves (this is a "square" window)
    bins = np.concatenate((f_low, f_high[-3:]), axis=None)

    # ff_sum = stats.binned_statistic(x=np.arange(len(ff_fft)), values=ff_fft, statistic='sum', bins=bins)[0]
    ff_sum = stats.binned_statistic(x=ff_xf, values=ff_fft, statistic='mean', bins=bins)[0]
    # print('FF Sum: {}'.format(ff_sum))
    # exit()
    wall_sum = stats.binned_statistic(x=wall_xf, values=wall_fft, statistic='mean', bins=bins)[0]
    # print(wall_sum)
    noise_sum = stats.binned_statistic(x=noise_xf, values=noise_fft, statistic='mean', bins=bins)[0]
    print(noise_sum)
    # Calculate ratio between third octave levels transmitted through wall against the free-feild measurements
    # fft_ratio = np.divide(wall_fft, ff_fft)*C_geo[ff_key]*C_dir*C_gain
    # reflection_ratio_fft = (stats.binned_statistic(x=ff_xf, values=fft_ratio, statistic='sum', bins=bins)[0])*df

    # fft_ratio[fft_ratio > 1] = 0
    # print(reflection_ratio_fft)
    raw_ratio = np.divide(wall_sum, ff_sum)
    reflection_ratio = (np.divide(wall_sum, ff_sum))*C_geo[ff_key]*C_dir*C_gain
    # print(reflection_ratio)
    reflection_ratio[reflection_ratio >= 1] = np.nan
    SNR_store.append(np.divide(wall_sum, noise_sum))
    SNR[wall_key] = 10*np.log10(np.divide(wall_sum, noise_sum))
    # print(SNR[wall_key])
    # print(reflection_ratio)
    RI_store.append( reflection_ratio )
    reflection_coefficient[ff_key] = reflection_ratio
    freefeild_fft[ff_key] = ff_fft
    wall_fft_store[wall_key] = wall_fft

    ff_store[wall_key] = ff_sum
    wall_store[wall_key] = wall_sum
    ###################
    SNR_plot = SNR[wall_key]
    CF = [125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000]
    # output_file(os.path.join(wall_folder, 'measurement_{}.html'.format(wall_key)))
    # p1 = figure(plot_width=800, plot_height=400, x_range=(100, 5000), x_axis_type="log", title='FFT of windowed data for {}'.format(wall_key))
    # p1.line(x=ff_xf, y=ff_fft, legend='Free field', color='red')
    # for low, mid, high in zip(f_low, f_center, f_high):
    #     p1.add_layout(Span(location=low, dimension='height', line_color='grey', line_dash='dashed', line_width=0.5))
    #     p1.add_layout(Span(location=mid, dimension='height', line_color='green', line_dash='dashed', line_width=1))
    #     p1.add_layout(Span(location=high, dimension='height', line_color='grey', line_dash='dashed', line_width=0.5))
    # p1.line(x=wall_xf, y=wall_fft, legend='Wall reflections', color='blue')
    # p1.xaxis.ticker = CF
    # p2 = figure(plot_width=800, plot_height=400, x_range=(100, 5000), x_axis_type="log", title='1/3rd Octave of windowed data' )
    # p2.line(x=f_center, y=ff_sum[1:-1], legend='Free field', color='red')
    # p2.line(x=f_center, y=wall_sum[1:-1], legend='Wall reflections', color='blue')
    # # p2.line(x=f_center, y=ff_int[1:-1], legend='Wall reflections', color='green')
    # p2.xaxis.ticker = CF
    # p3 = figure(plot_width=800, plot_height=400, x_range=(100, 5000), y_range=(0, 1.5), x_axis_type="log", title='Reflection ratio' )
    # p3.line(x=f_center, y=reflection_ratio[1:-1], color='red')
    # # p3.line(x=f_center, y=raw_ratio[1:-1], color='green')
    # p4 = figure(plot_width=800, plot_height=400, x_range=(100, 5000), x_axis_type="log", title='SNR' )
    # p4.line(x=f_center, y=SNR_plot[1:-1], color='red')
    # p5 = figure(plot_width=800, plot_height=400, x_range=(0, 0.03), title='Impulse Response - Raw' )
    # p5.line(x=t_array, y=ff_ir, line_color="blue", legend="FF", line_dash='dashed')
    # p5.line(x=t_array, y=wall_ir, line_color="red", legend="Wall")
    # p6 = figure(plot_width=800, plot_height=400, x_range=(0.008, 0.022), title='Impulse Response - Scaled' )
    # p6.line(x=t_array, y=ff_scaled, line_color="blue", legend='Scaled FF IR')
    # p6.line(x=t_array, y=wall_ir, line_color="red", legend="Wall", line_dash='dashed')
    # p7 = figure(plot_width=800, plot_height=400, x_range=(0.008, 0.025), title='Impulse Response - Difference' )
    # p7.line(x=t_array, y=wall_directRemoved, line_color="red", legend='wall')
    # p7.line(x=t_array, y=noise, line_color="blue", legend='noise')
    # p8 = figure(plot_width=800, plot_height=400, x_range=(0.008, 0.022), title='Impulse Response - Windowed' )
    # p8.line(x=t_array, y=diff_windowed, line_color="red", legend="Diff")
    # p8.line(x=t_array, y=ff_windowed, line_color="blue", legend="FF", line_dash='dashed')
    # # p9 = figure(plot_width=800, plot_height=400, x_range=(100, 5000), x_axis_type='log', y_range=(0, 1), title='Reflection ratio - fft' )
    # # p9.line(x=ff_xf, y=fft_ratio, line_color="red")
    # p100 = column([p1, p2, p3, p4, p5, p6, p7, p8])
    # show(p100)
    # export_png(p100, '{}/calculation_{}.png'.format(dest_folder, wall_key))
    # input('Next measurement')
    # exit()
    # show(column([p1, p2, p3, p4, p5, p6, p7]))



################### Plot difference in TOF #######################
# output_file(os.path.join(wall_folder, 'delays_{}.html'.format(wall_key)))
# p201 = figure(plot_width=800, plot_height=400, title='Effect of step on impulse response', x_range=(0.008, 0.020))
# p201.line(x=t_array, y=ff_impulse["mic9"], legend='Free field', color='black', line_dash='dashed')
# p201.line(x=t_array, y=wall_impulse["mic9"], legend='Mic 9 - Left', color='blue')
# p201.line(x=t_array, y=wall_impulse["mic8"], legend='Mic 8 - Middle', color='red')
# p201.line(x=t_array, y=wall_impulse["mic7"], legend='Mic 7 - Right', color='green')
# p201.xaxis.axis_label = 'Time (s)'
# p201.yaxis.axis_label = 'Level'
# p201.toolbar.logo = None
# p201.toolbar_location = None
# show(p201)
# export_png(p201, '{}/TOF.png'.format(dest_folder))

# exit()
##################### PLOT third octave CPB #######################
titles = [["Mic 9 - 1/3rd Octave","Mic 8 - 1/3rd Octave","Mic 7 - 1/3rd Octave"],
    ["Mic 6 - 1/3rd Octave","Mic 5 - 1/3rd Octave","Mic 4 - 1/3rd Octave"],
    ["Mic 3 - 1/3rd Octave","Mic 2 - 1/3rd Octave","Mic 1 - 1/3rd Octave"]]
x_data = [100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000]
y_data_1 = [[ff_store['mic9'], ff_store['mic8'], ff_store['mic7']],
    [ff_store['mic6'], ff_store['mic5'], ff_store['mic4']],
    [ff_store['mic3'], ff_store['mic2'], ff_store['mic1']]]
y_data_2 = [[wall_store['mic9'], wall_store['mic8'], wall_store['mic7']],
    [wall_store['mic6'], wall_store['mic5'], wall_store['mic4']],
    [wall_store['mic3'], wall_store['mic2'], wall_store['mic1']]]
legends = ['Free-feild', 'Reflection']
pg1 = en1793.mic_grid_plot(x_data=x_data, y_data_1=y_data_1, y_data_2=y_data_2, titles=titles,
    output_filename=os.path.join(dest_folder, "FFT.html"), legends=legends,
    x_labels='Frequencies [Hz]', y_labels='Amplitude', l_style=['solid','dashed'],
    y_range=(0, 1000), x_range=(100, 1000), p_height=300, p_width=500)
show(pg1)
export_png(pg1, '{}/ThirdOctave.png'.format(dest_folder))

exit()
##################### Plot FFT ################################
titles = [["Mic 9 - FFT","Mic 8 - FFT","Mic 7 - FFT"],
    ["Mic 6 - FFT","Mic 5 - FFT","Mic 4 - FFT"],
    ["Mic 3 - FFT","Mic 2 - FFT","Mic 1 - FFT"]]

x_data = ff_xf
y_data_1 = [[freefeild_fft['mic9'], freefeild_fft['mic8'], freefeild_fft['mic7']],
    [freefeild_fft['mic6'], freefeild_fft['mic5'], freefeild_fft['mic4']],
    [freefeild_fft['mic3'], freefeild_fft['mic2'], freefeild_fft['mic1']]]
y_data_2 = [[wall_fft_store['mic9'], wall_fft_store['mic8'], wall_fft_store['mic7']],
    [wall_fft_store['mic6'], wall_fft_store['mic5'], wall_fft_store['mic4']],
    [wall_fft_store['mic3'], wall_fft_store['mic2'], wall_fft_store['mic1']]]
legends = ['Free-feild', 'Reflection']
pg1 = en1793.mic_grid_plot(x_data=x_data, y_data_1=y_data_1, y_data_2=y_data_2, titles=titles,
    output_filename=os.path.join(dest_folder, "ThirdOctave.html"), legends=legends,
    x_labels='Frequencies [Hz]', y_labels='Amplitude', l_style=['solid','dashed'],
    y_range=(0, 500), x_range=(100, 10000), p_height=300, p_width=500)
show(pg1)
export_png(pg1, '{}/FFT.png'.format(dest_folder))

exit()
##################### Plot Windowed data ################################
titles = [["Mic 9 - Windowed","Mic 8 - Windowed","Mic 7 - Windowed"],
    ["Mic 6 - Windowed","Mic 5 - Windowed","Mic 4 - Windowed"],
    ["Mic 3 - Windowed","Mic 2 - Windowed","Mic 1 - Windowed"]]
x_data = t_array
y_data_1 = [[windowed_ff['mic9'], windowed_ff['mic8'], windowed_ff['mic7']],
    [windowed_ff['mic6'], windowed_ff['mic5'], windowed_ff['mic4']],
    [windowed_ff['mic3'], windowed_ff['mic2'], windowed_ff['mic1']]]
y_data_2 = [[windowed_reflection['mic9'], windowed_reflection['mic8'], windowed_reflection['mic7']],
    [windowed_reflection['mic6'], windowed_reflection['mic5'], windowed_reflection['mic4']],
    [windowed_reflection['mic3'], windowed_reflection['mic2'], windowed_reflection['mic1']]]
legends = ['Free-feild', 'Reflection']
pg2 = en1793.mic_grid_plot(x_data=x_data, y_data_1=y_data_1, y_data_2=y_data_2, titles=titles,
    output_filename=os.path.join(dest_folder, "allWindowed.html"), legends=legends,
    x_labels='Time [s]', y_labels='Amplitude', l_style=['solid','dashed'],
    y_range=(-1, 2), x_range=(0, 0.025), p_height=300, p_width=500)
show(pg2)
export_png(pg2, '{}/WindowedResponse.png'.format(dest_folder))

# exit()

################################ Plot RI ###########################
titles = [["Reflection Coefficient - Mic 9","Reflection Coefficient - Mic 8","Reflection Coefficient - Mic 7"],
    ["Reflection Coefficient - Mic 6","Reflection Coefficient - Mic 5","Reflection Coefficient - Mic 4"],
    ["Reflection Coefficient - Mic 3","Reflection Coefficient - Mic 2","Reflection Coefficient - Mic 1"]]
y_data_1 = [[reflection_coefficient['mic9'], reflection_coefficient['mic8'], reflection_coefficient['mic7']],
    [reflection_coefficient['mic6'], reflection_coefficient['mic5'], reflection_coefficient['mic4']],
    [reflection_coefficient['mic3'], reflection_coefficient['mic2'], reflection_coefficient['mic1']]]
CF = [125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000]
x_data = CF
plots_temp = [[None, None, None],[None, None, None],[None, None, None]]
r = 0
for row1, row2 in zip(y_data_1, titles):
        c = 0
        for item1, title_val in zip(row1, row2):
                p = figure(plot_width=500, plot_height=300, title=title_val, y_range=(0, 1), x_range=(100, 5000))
                p.xaxis.axis_label = 'One-third octave centre frequency [Hz]'
                p.yaxis.axis_label = 'Reflection Coefficient'
                p.xgrid.minor_grid_line_color = 'grey'
                p.xgrid.minor_grid_line_alpha = 0.1
                p.ygrid.minor_grid_line_color = 'grey'
                p.ygrid.minor_grid_line_alpha = 0.1
                p.line(x_data, item1[1:-1], line_color="blue")
                p.toolbar.logo = None
                p.toolbar_location = None
                plots_temp[r][c] = p
                c+=1
        r+=1
pg3 = gridplot(plots_temp, toolbar_options={'logo': None})
show(pg3)
export_png(pg3, '{}/ReflectionCoeff.png'.format(dest_folder))

# exit()

########################## Plot SNR ######################################
titles = [["Signal to Noise - Mic 9","Signal to Noise - Mic 8","Signal to Noise - Mic 7"],
    ["Signal to Noise - Mic 6","Signal to Noise - Mic 5","Signal to Noise - Mic 4"],
    ["Signal to Noise - Mic 3","Signal to Noise - Mic 2","Signal to Noise - Mic 1"]]
y_data_2 = [[SNR['mic9'], SNR['mic8'], SNR['mic7']],
    [SNR['mic6'], SNR['mic5'], SNR['mic4']],
    [SNR['mic3'], SNR['mic2'], SNR['mic1']]]
print(y_data_1)
CF = [125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000]
x_data = CF
plots_temp = [[None, None, None],[None, None, None],[None, None, None]]
r = 0
for row1, row2 in zip(y_data_2, titles):
        c = 0
        for item1, title_val in zip(row1, row2):
                p = figure(plot_width=500, plot_height=300, title=title_val, y_range=(0, 100), x_range=(100, 5000))
                p.xaxis.axis_label = 'One-third octave centre frequency [Hz]'
                p.yaxis.axis_label = 'Signal to Noise Ratio'
                p.xgrid.minor_grid_line_color = 'grey'
                p.xgrid.minor_grid_line_alpha = 0.1
                p.ygrid.minor_grid_line_color = 'grey'
                p.ygrid.minor_grid_line_alpha = 0.1
                p.line(x_data, item1[1:-1], line_color="blue")
                p.toolbar.logo = None
                p.toolbar_location = None
                plots_temp[r][c] = p
                c+=1
        r+=1
pg4 = gridplot(plots_temp, toolbar_options={'logo': None})
show(pg4)
export_png(pg4, '{}/SNR.png'.format(dest_folder))

# exit()
# print(SNR_store)
SNR_total = 10.0*np.log10(np.nanmean(SNR_store, axis=0))
print(' --------------- SNR -----------------')
print(SNR_total)
RI_total = np.nanmean(RI_store, axis=0)
print(' --------------- RI -----------------')
print(RI_total)
with open(os.path.join(dest_folder, "results.csv"), 'w') as writeFile:
    print('Writing CSV')
    writer = csv.writer(writeFile, delimiter=',')
    writer.writerow(["Centre Frequency (Hz)"])
    writer.writerow(f_center)
    writer.writerow(["Signal to Noise (dB)"])
    writer.writerow(SNR_total)
    writer.writerow(["Reflection Coeff"])
    writer.writerow(RI_total)
    writer.writerow(["Subtractoin Ratio"])
    writer.writerow(subtraction_R)
writeFile.close()

################################ Plot IR for each mic ####################
# print('Plotting Measurement Results')
# titles = [["Impulse Response - Mic 9","Impulse Response - Mic 8","Impulse Response - Mic 7"],
#     ["Impulse Response - Mic 6","Impulse Response - Mic 5","Impulse Response - Mic 4"],
#     ["Impulse Response - Mic 3","Impulse Response - Mic 2","Impulse Response - Mic 1"]]
# x_data = t_array
# y_data_1 = [[ff_avg['mic9'], ff_avg['mic8'], ff_avg['mic7']],
#     [ff_avg['mic6'], ff_avg['mic5'], ff_avg['mic4']],
#     [ff_avg['mic3'], ff_avg['mic2'], ff_avg['mic1']]]
# y_data_2 = [[wall_avg['mic9'], wall_avg['mic8'], wall_avg['mic7']],
#     [wall_avg['mic6'], wall_avg['mic5'], wall_avg['mic4']],
#     [wall_avg['mic3'], wall_avg['mic2'], wall_avg['mic1']]]
# legends = ['Free-feild', 'Reflection']
# p1 = en1793.mic_grid_plot(x_data=x_data, y_data_1=y_data_1, y_data_2=y_data_2, titles=titles,
#     output_filename=os.path.join(dest_folder, "allIR.html"), legends=legends,
#     x_labels='Time [s]', y_labels='Amplitude', l_style=['solid','dashed'],
#     y_range=(-1, 2), x_range=(0, 0.025), p_height=300, p_width=500)
# show(p1)
# export_png(p1, '{}/Impulse_Response.png'.format(dest_folder))

# exit()

#########################################################################