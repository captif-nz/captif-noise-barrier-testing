import numpy as np
import json
from scipy.signal import correlate
import sys
import matplotlib.pyplot as plt
import pickle
import json

# Import EN1793 specific libraries
import EN_1793_funcs as en1793

# Helper functions
def check_ess_equivalence(ess_bottom, ess_middle, ess_top):
    bool1 = np.array_equal(ess_bottom, ess_middle)
    bool2 = np.array_equal(ess_bottom, ess_top)
    bool3 = np.array_equal(ess_top, ess_middle)
    bool4 = bool1 & bool2 & bool3
    return bool4

src_folder = "Results/test20190401/"
dest_folder = "Results/test20190401/"
bottom_file = "Results/test20190401/reflectionCAPTIF_ff_2.txt"

ir_file = "reflectionCAPTIF_ff_2_IR"
plotting = False

#bottom_data = np.load(bottom_file)
#middle_data = np.load(middle_file)
#top_data = np.load(top_file)

# Mic array numbering:
# top
# [7 8 9]
# [4 5 6]
# [1 2 3]
# bottom
with open(bottom_file) as bottom: #, open(middle_file) as middle, open(top_file) as top:
    bottom_data = json.load(bottom)
    middle_data = bottom_data
    top_data = bottom_data
# with open(middle_file, 'rb') as middle:
#     middle_data = pickle.loac(middle):
# with open(middle_file, 'rb') as middle:
#     middle_data = pickle.loac(middle):

mic_data = {
    "mic1": bottom_data['mic0'],
    "mic2": bottom_data['mic1'],
    "mic3": bottom_data['mic2'],
    "mic4": middle_data['mic0'],
    "mic5": middle_data['mic1'],
    "mic6": middle_data['mic2'],
    "mic7": top_data['mic0'],
    "mic8": top_data['mic1'],
    "mic9": top_data['mic2']
}

# input_signals = {
#     "bottom_ess": bottom_data['bottom_ess'],
#     "middle_ess": middle_data['middle_ess'],
#     "top_ess": top_data['top_ess']
# }

# Check all parameters match for the measurements
if bottom_data["f0"] == middle_data["f0"] == top_data["f0"]:
    print("f0 matches for all measurements - proceeding")
    f0 = bottom_data["f0"]
else:
    print("f0 is different between measurements! Stopping!")
    sys.exit("f0 mismatch")

if bottom_data["f1"] == middle_data["f1"] == top_data["f1"]:
    print("f1 matches for all measurements - proceeding")
    f1 = bottom_data["f1"]
else:
    print("f1 is different between measurements! Stopping!")
    sys.exit("f1 mismatch")

if bottom_data["t_meas"] == middle_data["t_meas"] == top_data["t_meas"]:
    print("t_meas matches for all measurements - proceeding")
    t_meas = bottom_data["t_meas"]
else:
    print("t_meas is different between measurements! Stopping!")
    sys.exit("t_meas mismatch")

if bottom_data["SR"] == middle_data["SR"] == top_data["SR"]:
    print("SR matches for all measurements - proceeding")
    SR = bottom_data["SR"]
else:
    print("SR is different between measurements! Stopping!")
    sys.exit("SR mismatch")

if check_ess_equivalence(bottom_data['ess_sweep'], bottom_data['ess_sweep'], bottom_data['ess_sweep']):
    print("ESS matches for all measurements - proceeding")
    ess = bottom_data['ess_sweep']
else:
    print("ESS is different between measurements! Stopping!")
    sys.exit("ESS mismatch")

# ess_start = int(round(9.5*SR))
# ess_end = int(round((t_meas - 0.5)*SR))
print("Calculating impulse responses")
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
        xcorr = correlate(measurement, ess, 'valid')
        # Use max() to get the sample delay - need to check this will work with the new excitation. Should only have one peak as wormup is not ESS
        i_offset = xcorr.tolist().index(xcorr.max())
        print("Using offset of: {}".format(i_offset))
        # Shift entire measured data set and select only the ESS section
        measurement_clipped = measurement[i_offset: len(ess)+i_offset]
        # Calculate impulse response
        ir_temp = en1793.impulse_response(ess, measurement_clipped)
        # Roll to avoid splitting IR at 0
        ir_temp2 = np.roll(ir_temp, int(0.01*51200))
        ir.append(ir_temp2)
        ir_clipped = ir_temp2[0:int(0.05*SR)]
        # print(ir_clipped.mean())
        ir_normalised = ir_clipped-ir_clipped.mean()
        ir_formatted.append(ir_normalised)

        if plotting == True:
            plt.subplot(4,1,1)
            plt.plot(np.arange(len(measurement))*1./51200, measurement)
            plt.plot(np.arange(len(bottom_data['test_signal']))*1./51200, bottom_data['test_signal'])
            plt.subplot(4,1,2)
            plt.plot(np.arange(len(xcorr)), xcorr)
            plt.subplot(4,1,3)
            plt.plot(np.arange(len(measurement_clipped))*1./51200, measurement_clipped)
            plt.plot(np.arange(len(bottom_data['ess_sweep']))*1./51200, bottom_data['ess_sweep'])
            plt.subplot(4,1,4)
            plt.plot(np.arange(len(ir_temp))*1./51200, ir_temp)
            plt.xlim([0,0.1])
            plt.ylim
            plt.show()
    micIR[key] = ir
    avgIR[key] = np.mean(ir_formatted, axis=0) - np.mean(np.mean(ir_formatted, axis=0))
    # print( np.mean(np.mean(ir_formatted, axis=0)) )

print("Saving impulse responses as: {}.txt".format(ir_file))
json_data = {
        "impulseResponse": micIR,
        "averagedImpulseResponse": avgIR,
        "f0":f0,
        "f1":f1,
        "SR":SR,
        "t_meas":t_meas,
}
with open(dest_folder + ir_file, 'wb') as outfile:
    pickle.dump(json_data, outfile)
# np.savez(ir_file, impulse_response=IR, f0=f0, f1=f1, SR=SR, t_meas=t_meas)

# tt = np.array(range(len(ess)))*1./SR
# ndx = 1
# for key, value in micIR.items():
#     plt.subplot(9,1,ndx)
#     plt.plot(tt, value[0], label="measurement 1")
#     # plt.xlim([0,1])
#     plt.title(key)
#     ndx+=1

# plt.xlabel("Time")
# plt.ylabel("Impulse response")
# plt.show()
