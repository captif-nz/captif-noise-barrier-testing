# -*- coding: utf-8 -*-

# Import python libraries
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from time import sleep
import json
import os
from scipy.signal import correlate
import sys
import time
from bokeh.plotting import figure, output_file, show, gridplot
import pickle

# Import EN1793 specific libraries
import EN_1793_funcs as en1793

# Import NI libraries
import nidaqmx
from nidaqmx.constants import AcquisitionType

def NImeasurement(num_samp, n_runs, test_signal, t_meas):
        # Initialise measurement arrays
        mic0 = []
        mic1 = []
        mic2 = []
        # Run nidaqmx as two tasks
        with nidaqmx.task.Task("OutputTask") as ao_task, nidaqmx.task.Task("InputTask") as ai_task:

                ################## Setup analogue output ################
                ao_task.ao_channels.add_ao_voltage_chan(ao_device_name+"/ao0",
                                                        name_to_assign_to_channel='speaker_output',
                                                        min_val=-4.0,
                                                        max_val=4.0)
                ao_task.timing.cfg_samp_clk_timing(SR,
                                                sample_mode=AcquisitionType.FINITE,
                                                samps_per_chan=num_samp)
                ao_task.write(test_signal,
                        auto_start=False)

                ############ Setup microphone inputs #############
                ai_task.ai_channels.add_ai_microphone_chan(ai_device_name+"/ai0",
                                                        name_to_assign_to_channel='microphone_input0',
                                                        mic_sensitivity=22.4,
                                                        max_snd_press_level=140,
                                                        units=nidaqmx.constants.SoundPressureUnits.PA)
                ai_task.ai_channels.add_ai_microphone_chan(ai_device_name+"/ai1",
                                                        name_to_assign_to_channel='microphone_input1',
                                                        mic_sensitivity=22.4,
                                                        max_snd_press_level=140,
                                                        units=nidaqmx.constants.SoundPressureUnits.PA)
                ai_task.ai_channels.add_ai_microphone_chan(ai_device_name+"/ai2",
                                                        name_to_assign_to_channel='microphone_input2',
                                                        mic_sensitivity=22.4,
                                                        max_snd_press_level=140,
                                                        units=nidaqmx.constants.SoundPressureUnits.PA)
                ai_task.timing.cfg_samp_clk_timing(SR,
                                                sample_mode=AcquisitionType.FINITE,
                                                samps_per_chan=num_samp)

                ########## Perform required number of measurements (typically 5)
                for ndx in range(n_runs):
                        print("Measurement: " + str(ndx))
                        # Start output signal and measurements
                        ao_task.start()
                        ai_task.start()
                        # Wait until both tasks are complete
                        ao_task.wait_until_done(timeout=t_meas+5)
                        ai_task.wait_until_done(timeout=t_meas+5)
                        # Record data from 3 mics
                        data = ai_task.read(number_of_samples_per_channel=num_samp)

                        print("storing data from last run")
                        mic0.append(data[0])
                        mic1.append(data[1])
                        mic2.append(data[2])

                        print("Pausing for 10 seconds before starting next measurement")
                        sleep(5)
                        ao_task.stop()
                        ai_task.stop()
                        sleep(5)

        return mic0, mic1, mic2


# General parameters
dest_folder = "CurrentTest"
SR = 51200
t_meas = 30
n_runs = 5
f0 = 40
f1 = 20000
ao_device_name = "AO1-NI9260"
ai_device_name = "AI1-NI9234"
num_samp = int(SR*t_meas)

test = True

print("Starting measurement of impulse response using ESS")
dir_setup=False
measurement_name = input("Enter the name of this test:")
while dir_setup == False:
        if measurement_name in os.listdir(dest_folder):
                measurement_name = input("Measurement folder already exists try another:")
                dir_setup=False
        else:
                dir_setup=True

print("\nStarting test using the following paramenters")
print("Sampling rate: {}Hz".format(SR))
print("Measurement duration: {} seconds".format(t_meas))
print("Number of runs: {}".format(n_runs))
print("Frequency range of ESS: {}Hz to {}Hz".format(f0, f1))
print("Input device ID: {}".format(ai_device_name))
print("Output device ID: {}".format(ao_device_name))

# Make directory to save data in
os.mkdir(os.path.join(dest_folder, measurement_name))
# Generate ESS signal
test_signal, ess_sweep, test_t, ess_t = en1793.ess_signal(t1=t_meas,
                                                          f0=f0,
                                                          f1=f1,
                                                          SR=SR,
                                                          scale=0.2)


print("\nStarting Measurements ---------------")

mic_data = {
    "mic1": [],
    "mic2": [],
    "mic3": [],
    "mic4": [],
    "mic5": [],
    "mic6": [],
    "mic7": [],
    "mic8": [],
    "mic9": []
}

############# TEST DATA ##########
if test == True:
        print("USING TEST DATA!!!!!!")
        test_data_file= "Results/test20190122/bottom_CAPTIF_ff_20190122.txt"
        with open(test_data_file) as testing:
                test_data = json.load(testing)
                mic0 = test_data['mic0']
                mic1 = test_data['mic1']
                mic2 = test_data['mic2']

for location in ['bottom', 'middle', 'top']:
        # input("Place microphone array at {} location and press enter...".format(location))
        if not test==True:
                mic0, mic1, mic2 = NImeasurement(num_samp, n_runs, test_signal, t_meas)
        else:
                print("THIS IS USING TEST DATA!!!!")
                time.sleep(30)
        print("{} measurement completed - saving data".format(location))
        json_data = {
                "mic0": mic0,
                "mic1":mic1,
                "mic2":mic2,
                "test_signal":test_signal.tolist(),
                "test_t":test_t.tolist(),
                "ess_sweep":ess_sweep.tolist(),
                "ess_t":ess_t.tolist(),
                "f0":f0,
                "f1":f1,
                "SR":SR,
                "t_meas":t_meas,
                "n_runs":n_runs
        }
        filename = os.path.join(dest_folder, measurement_name, location)
        with open(filename, 'wb') as outfile:
                pickle.dump(json_data, outfile)

        if location == 'bottom':
                mic_data['mic1'] = mic0
                mic_data['mic2'] = mic1
                mic_data['mic3'] = mic2
        elif location == 'middle':
                mic_data['mic4'] = mic0
                mic_data['mic5'] = mic1
                mic_data['mic6'] = mic2
        elif location == 'top':
                mic_data['mic7'] = mic0
                mic_data['mic8'] = mic1
                mic_data['mic9'] = mic2

print("Measurements completed")
print("---------------------\n------------------")

###### Calculate impulse response
print("Calculating impulse responses for measurement")

# json_IR_data = {}
micIR = {}
avgIR = {}
for key in mic_data:
    print("Processing: {}".format(key))
    measured_data = mic_data[key]
    # Cycle through the measurements
    ir = []
    ir_formatted = []
    counter1 = 0
    for row in measured_data:
        print("Evaluating measurement: {}".format(counter1))
        counter1+=1
        measurement = row
        # find offset between the input and output signals
        xcorr = correlate(measurement, ess_sweep, 'valid')
        # Use max() to get the sample delay - need to check this will work with the new excitation. Should only have one peak as wormup is not ESS
        i_offset = xcorr.tolist().index(xcorr.max())
        print("Using offset of: {}".format(i_offset))
        # Shift entire measured data set and select only the ESS section
        measurement_clipped = measurement[i_offset: len(ess_sweep)+i_offset]
        # Calculate impulse response
        ir_temp = en1793.impulse_response(ess_sweep, measurement_clipped)
        # Roll to avoid splitting IR at 0
        ir_temp2 = np.roll(ir_temp, int(0.01*51200))
        ir.append(ir_temp2)
        ir_clipped = ir_temp2[0:int(0.05*SR)]
        # print(ir_clipped.mean())
        ir_normalised = ir_clipped-ir_clipped.mean()
        ir_formatted.append(ir_normalised)
        # json_IR_data[key+"_"+str(counter1)] = ir

    micIR[key] = ir
    avgIR[key] = np.mean(ir_formatted, axis=0) - np.mean(np.mean(ir_formatted, axis=0))


json_IR_data = {
        "impulseResponse": micIR,
        "averagedImpulseResponse": avgIR,
        "f0":f0,
        "f1":f1,
        "SR":SR,
        "t_meas":t_meas
}

ir_file = os.path.join(dest_folder, measurement_name, "IR")
with open(ir_file, 'wb') as outfile:
    pickle.dump(json_IR_data, outfile)

print("Impulse response calculation completed")
print("---------------------\n------------------")

print('Plotting Measurement Results')
# output to static HTML file
output_file(os.path.join(dest_folder, measurement_name, "ImpulseResponse.html"))

p1 = figure(plot_width=800, plot_height=400, title="Mic 1 Average IR")
p1.line(np.arange(len(avgIR['mic1']))*1./51200, avgIR['mic1'], line_width=2)

p2 = figure(plot_width=800, plot_height=400, title="Mic 2 Average IR")
p2.line(np.arange(len(avgIR['mic2']))*1./51200, avgIR['mic2'], line_width=2)

p3 = figure(plot_width=800, plot_height=400, title="Mic 3 Average IR")
p3.line(np.arange(len(avgIR['mic3']))*1./51200, avgIR['mic3'], line_width=2)

p4 = figure(plot_width=800, plot_height=400, title="Mic 4 Average IR")
p4.line(np.arange(len(avgIR['mic4']))*1./51200, avgIR['mic4'], line_width=2)

p5 = figure(plot_width=800, plot_height=400, title="Mic 5 Average IR")
p5.line(np.arange(len(avgIR['mic5']))*1./51200, avgIR['mic5'], line_width=2)

p6 = figure(plot_width=800, plot_height=400, title="Mic 6 Average IR")
p6.line(np.arange(len(avgIR['mic6']))*1./51200, avgIR['mic6'], line_width=2)

p7 = figure(plot_width=800, plot_height=400, title="Mic 7 Average IR")
p7.line(np.arange(len(avgIR['mic7']))*1./51200, avgIR['mic7'], line_width=2)

p8 = figure(plot_width=800, plot_height=400, title="Mic 8 Average IR")
p8.line(np.arange(len(avgIR['mic8']))*1./51200, avgIR['mic8'], line_width=2)

p9 = figure(plot_width=800, plot_height=400, title="Mic 9 Average IR")
p9.line(np.arange(len(avgIR['mic9']))*1./51200, avgIR['mic9'], line_width=2)

p = gridplot([[p1, p2, p3], [p4, p5, p6], [p7, p8, p9]])
# show the results
show(p)



