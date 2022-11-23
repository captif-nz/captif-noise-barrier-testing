
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from scipy.signal import correlate
from scipy.fftpack import fft, ifft


ess = np.load('testing_20191501/input_array_30sec.npy')
ff1 = np.load('testing_20191501/ff_test1.npy')
ff2 = np.load('testing_20191501/ff_test2.npy')
fs = 51200

ff1 = ff1[2, :]  # Limit to ch2 for this example.
ff2 = ff2[2, :]

# There are 2 periods in the recording. We need to process each period
# seperately.

# There will also be an offset at the start, which is made up of time-of-flight
# and delays in the measurement system. It is essential that the chuck of
# recording that is processed is exactly the same chuck as the excitation
# signal. For MLS this is less of an issue as the signal is naturally periodic,
# so you just grab any chunck out of the middle of the recording. But it's not
# as easy with ESS.

# We can find the initial delay by doing a cross-correlation between the two
# signals. This should show 2 peaks (as there are two full periods in the
# recorded signal).

xcorr1 = correlate(ff1, ess, 'valid')
xcorr2 = correlate(ff2, ess, 'valid')
# plt.plot(xcorr)
# plt.show()

# But there is only one peak! This is probably because the first period isn't
# fully captured by the measurement system (i.e. the system starts playing before
# it starts recording). You can get around this by adding a few seconds of
# silence to the start of the excitation signal, or better yet add and extra
# period of the excitation signal and treat the first period as a warm-up phase.

# Since there is only one peak, we can just use max() to get the sample delay:
i_offset1 = xcorr1.tolist().index(xcorr1.max())
i_offset2 = xcorr2.tolist().index(xcorr2.max())
print(i_offset1/fs)
print(i_offset1/fs)

# Grab the second period of the recorded signal:
fff1 = ff1[i_offset1: len(ess)+i_offset1]
fff2 = ff2[i_offset2: len(ess)+i_offset2]

# Take fft of excitation and recorded signals:
fff1_dft = fft(fff1)
fff2_dft = fft(fff2)
ess_dft = fft(ess)

# Deconvolution (in frequency domain):
tf1 = np.divide(fff1_dft, ess_dft)
tf2 = np.divide(fff2_dft, ess_dft)

# Convert back to time domain:
ir1 = ifft(tf1).real  # Also take the real component (I guess this is a division precision thing)
ir2 = ifft(tf2).real

# Circular shift the impulse response 10ms to the right so it's not split across the plot:
ir1 = np.roll(ir1, int(0.01*fs))
ir2 = np.roll(ir2, int(0.01*fs))
ir_avg = np.mean(np.array([ir1, ir2]), axis=0)

# Time array
tt1 = np.array(range(len(ir1)))*1./fs
tt2 = np.array(range(len(ir2)))*1./fs

plt.subplot(3,1,1)
plt.plot(tt1, ir1)
plt.xlim([0, 0.020])
plt.subplot(3,1,2)
plt.plot(tt2, ir2)
plt.xlim([0, 0.020])  # Plot first 20ms
plt.subplot(3,1,3)
plt.plot(tt2, ir_avg)
plt.xlim([0, 0.020])
plt.show()
