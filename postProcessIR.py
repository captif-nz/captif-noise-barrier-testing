
from scipy import signal
from scipy.fftpack import ifft, fft, fftshift
import numpy as np
import math
import scipy.fftpack

from scipy.signal import correlate
import sys
import pickle
from scipy import stats
from bokeh.plotting import figure, output_file, show, gridplot
import os
# Import EN1793 specific libraries
import EN_1793_funcs as en1793
import octave_filters as octave
import pickle
import sys


measurement_name = "MaioroSt_Post_20190625_FF"
base_dir = "BarrierTesting"

# for measurement_name in measurement_names:
measurement_dir = os.path.join(base_dir, measurement_name)
print("Loading data from {}".format(measurement_dir))
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

with open(os.path.join(measurement_dir, 'bottom'), 'rb') as bottom_data, open(os.path.join(measurement_dir, 'middle'), 'rb') as middle_data, open(os.path.join(measurement_dir, 'top'), 'rb') as top_data:
        bottom = pickle.load(bottom_data)
        print(bottom.keys())
        middle = pickle.load(middle_data)
        top = pickle.load(top_data)
        ess_sweep = bottom['ess_sweep']
        f0 = bottom['f0']
        f1= bottom['f1']
        SR= bottom['SR']
        t_meas = bottom['t_meas']

        mic_data['mic1'] = bottom['mic0']
        mic_data['mic2'] = bottom['mic1']
        mic_data['mic3'] = bottom['mic2']
        mic_data['mic4'] = middle['mic0']
        mic_data['mic5'] = middle['mic1']
        mic_data['mic6'] = middle['mic2']
        mic_data['mic7'] = top['mic0']
        mic_data['mic8'] = top['mic1']
        mic_data['mic9'] = top['mic2']

print(mic_data.keys())
# sys.exit()

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

ir_file = os.path.join(measurement_dir, "IR")
with open(ir_file, 'wb') as outfile:
    pickle.dump(json_IR_data, outfile)

print("Impulse response calculation completed")