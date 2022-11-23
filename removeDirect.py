import numpy as np
import EN_1793_funcs as en1793
import octave_filters
import json
import matplotlib.pyplot as plt
import pickle
import time

# Import EN1793 specific libraries
import EN_1793_funcs as en1793
import octave_filters as octave


# bottom
ff_filename = "test20190122/CAPTIF_ff_IR"
wall_filename = "test20190122/CAPTIF_wall_IR"

with open(ff_filename, 'rb') as file1: #, open(middle_file) as middle, open(top_file) as top:
    ff_data = pickle.load(file1)

with open(wall_filename, 'rb') as file1: #, open(middle_file) as middle, open(top_file) as top:
    wall_data = pickle.load(file1)

ff_avg = ff_data["averagedImpulseResponse"]
wall_avg = wall_data["averagedImpulseResponse"]

startOfReflection = int(0.035*51200)

reflected = np.asarray(ff_avg['mic1'])
free_feild = np.copy(reflected)
free_feild[startOfReflection:] = 0


tt = np.array(range(len(reflected)))*1./51200
# plt.subplot(2,1,1)
# plt.plot(tt, reflected)
# plt.xlim([0, 0.05])
# plt.ylim(-2, 2)
# plt.title("Wall")
# plt.subplot(2,1,2)
# plt.plot(tt, free_feild)
# plt.title("Free-feild")
# plt.xlim([0, 0.05])
# plt.ylim(-2, 2)
# plt.tight_layout()
# plt.show()

diff = reflected - free_feild

# plt.subplot(1,1,1)
# plt.plot(tt, diff)
# plt.xlim([0, 0.05])
# plt.ylim(-2, 2)
# plt.title("Difference between free-feild and reflected signals")
# plt.tight_layout()
# plt.show()

diff_windowed = en1793.apply_adrienne(diff, offset=0.035, SR=51200)

plt.subplot(1,1,1)
plt.plot(tt, diff_windowed)
plt.xlim([0, 0.05])
plt.ylim(-2, 2)
plt.title("Windowed difference")
plt.tight_layout()
plt.show()