from scipy import signal as signal
from scipy.fftpack import ifft, fft, fftshift
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.fftpack


def test_adrienne(SR, windowLength):
    totL = windowLength*10**-3 # Can be either 7.9ms (standard length) or 6 ms
#    SR = SR*10**3
    # Generate sample lengths for each section
    preL = 0.5 *10**-3
    mainL = (windowLength*10**-3 - preL)*(7./10.)
    trailL = (windowLength*10**-3 - preL)*(7./10.)

    ## Generate windows
    preSL = math.ceil(SR*preL)  # Length of rising window
    prewin = signal.blackmanharris(preSL*2) # Generate rising window section

    mainSL = math.ceil(SR*mainL)    # Length of flat window section
    mainwin = signal.boxcar(mainSL) #

    trailSL = math.ceil(SR*trailL)
    trailwin = signal.blackmanharris(trailSL*2)

    totwin = np.empty(preSL+mainSL+trailSL)

    totwin[:preSL] = prewin[:preSL]
    totwin[preSL:(preSL+mainSL)] = mainwin[:mainSL]
    totwin[(preSL+mainSL):] = trailwin[(trailSL):]
    return totwin

def pad_adrienne(window, totLen, SR=51200):
    start = int(0.*SR*10**-3)
    end = int(totLen*SR*10**-3 - len(window))
    padded_window = np.pad(window, (start,end), 'constant', constant_values=0)
    return padded_window

## Enter height and width of sample
height = 4.5
width = 10
diffractionTime = 1000*min([height/2, width/2])/340
print("Time of flight for diffracted signal: {}ms".format(diffractionTime))

Fs = 51200
winLen = diffractionTime
totLen = winLen + 200.

win = test_adrienne(SR=Fs, windowLength=winLen)
padded_win = pad_adrienne(window=win, totLen=totLen, SR=Fs)
filtLen = len(padded_win)/Fs
timeSeries = np.arange(0, filtLen, 1/Fs)

yf = np.abs(np.fft.fft(padded_win))
xf = np.linspace(0.0, 1.0/(2.0/Fs), len(padded_win)/2)

peaks, _ = signal.find_peaks(-1.*yf)
first_min_x = xf[peaks[0]]
first_min_y = yf[peaks[0]]
print("First minimum occurs at: {}Hz".format(first_min_x))

plt.subplot(2,1,1)
plt.plot(timeSeries, padded_win)
plt.title("Adrienne")
plt.xlabel("Time (s)")
plt.ylabel("Filter shape")
plt.subplot(2,1,2)
plt.plot(xf, 20*np.log10(yf[:len(padded_win)//2]))
plt.plot(first_min_x, 20*np.log10(first_min_y), 'o', mfc='none')
plt.title("Adrienne")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Filter response")
plt.ylim([0,100])
plt.xlim([0,1000])
plt.tight_layout()
plt.show()