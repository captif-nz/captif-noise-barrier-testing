from scipy import signal
from scipy.fftpack import ifft, fft, fftshift
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.fftpack
from scipy.signal import correlate
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.palettes import Set1_9 as palette
import itertools
import random
import string

def remove_drift(data, time_array):
    p = np.poly1d(np.polyfit(time_array, data, 3))
    p_fitted = p(time_array)
    return data - p_fitted

def random_string(length):
    return ''.join(random.choice(string.ascii_letters) for m in range(length))

def mic_grid_plot(x_data, y_data_1, y_data_2, titles, output_filename='temp.html', x_labels='',
    y_labels='', legends='', p_height=300, p_width=600, l_width=1, l_style=['solid','solid'],
    y_range=(-2, 2), x_range=(0, 1)):
        n_rows = np.shape(x_data)[0]
        r = np.arange(int(np.ceil(n_rows/2)))
        plots_temp = [[None, None, None],[None, None, None],[None, None, None]]
        r = 0
        for row1, row2, row3 in zip(y_data_1, y_data_2, titles):
                c = 0
                for item1, item2, title_val in zip(row1, row2, row3):
                        p = figure(plot_width=p_width, plot_height=p_height, title=title_val, y_range=y_range, x_range=x_range)
                        p.xaxis.axis_label = x_labels
                        p.yaxis.axis_label = y_labels
                        p.xgrid.minor_grid_line_color = 'grey'
                        p.xgrid.minor_grid_line_alpha = 0.1
                        p.ygrid.minor_grid_line_color = 'grey'
                        p.ygrid.minor_grid_line_alpha = 0.1
                        p.line(x_data, item1, line_width=l_width, line_color="blue", legend=legends[0],  line_dash=l_style[0])
                        p.line(x_data, item2, line_width=l_width, line_color="red", legend=legends[1],  line_dash=l_style[1])
                        p.toolbar.logo = None
                        p.toolbar_location = None
                        plots_temp[r][c] = p
                        c+=1
                r+=1
        p_f = gridplot(plots_temp, toolbar_options={'logo': None})
        # show(p_f)
        return p_f

def simple_plot(X, Y, labels):
    output_file('plots/'+random_string(10)+'.html')
    p = figure(plot_width=800, plot_height=400)
    colors = itertools.cycle(palette)
    for x, y, label, color in zip(X, Y, labels, colors):
        p.line(x=x, y=y, legend=label, color=color)
    show(p)
    return True

def impulse_response(data_in, data_out):
    return ifft(np.divide(fft(data_out), fft(data_in))).real

# Function to generate Arienne window
def adrienne(SR, windowLength=7.9):
    "window = ardienne(windowLength, SR) windowLength = length in ms, SR = sample rate in kHz"
    totL = windowLength*10**-3 # Standard length is 7.9ms
#    SR = SR*10**3
    # Generate sample lengths for each section
    preL = 0.5 *10**-3
    if totL == 7.9*10**-3:
        mainL = 5.18*10**-3
        trailL = 2.22*10**-3
    elif totL == 6.0*10**-3:
        mainL = 3.85*10**-3
        trailL = 1.65*10**-3
    else:
        flat_BH_len = windowLength*10**-3-preL
        mainL = flat_BH_len*(0.7)
        trailL = flat_BH_len*(0.3)

    ###### Print lengths ################
    # print("Pre-window: {}".format(preL))
    # print("Flat section: {}".format(mainL))
    # print("Tail: {}".format(trailL))
    #####################################
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
    ############## Plot window ###################
    # simple_plot(X=[np.arange(0,totL,1./SR),], Y=[totwin,])
    ############################################
    return totwin

def apply_adrienne(data, offset=0, SR=51200, windowLength=7.9):
    offset = offset*10**-3
    window = adrienne(SR=SR, windowLength=windowLength)
    start = int(offset*SR)
    end = len(data)-len(window)-int(offset*SR)
    padded_window = np.pad(window, (start,end), 'constant', constant_values=0)
    ######################## Plot window and data ###############################
    # simple_plot(X=[np.arange(0, len(data)/SR, 1/SR),np.arange(0, len(padded_window)/SR, 1/SR)], Y=[data, padded_window*max(data)], labels=["Raw data","Window"])
    #################################################
    return data*padded_window

def ess_signal(t1=30, f0=10, f1=16000, SR=51200, scale=1):
    """
    Generates a standard ESS signal.
    This signal has a 7 second warm up cycle - this ensures the whole ess is recorded
    Has a 3 second silence at the end of the test to avoid any "clicks"
    There is a short delay between signal starting and the recording starting 
    """
    dt = 1/SR
    ess_time = int(t1 - 10)
    test_t = np.arange(0, t1, dt)

    warmup1 = np.zeros(int(1*SR))
    warmup2 = scale*mlsSignal(t=5, N=8, SR=51.2*10**3)*(np.arange(0,5,dt)/5)
#    warmup1 = signal.chirp(t=warmup1_t, f0=f0, t1=5, f1=f1, method='linear')
    warmup3 = np.zeros(int(2*SR))

    ess_t = np.arange(0, ess_time, dt)
    ess_sweep = scale*signal.chirp(t=ess_t, f0=f0, t1=ess_time, f1=f1, method="logarithmic")

    cooldown_t = np.arange(0,2,dt)
    cooldown = np.zeros(len(cooldown_t))

    test_signal = np.concatenate((warmup1, warmup2, warmup3, ess_sweep, cooldown))
    return test_signal, ess_sweep, test_t, ess_t

def mlsSignal(t, N=8, SR=51200):
    mls_signal = signal.max_len_seq(N)[0] - 0.5
    mls_length = len(mls_signal)
    sample_length = int(SR*t)
    tiles = math.ceil(sample_length/mls_length)
    tiled_mls = np.tile(mls_signal, tiles)
    mls_signal = tiled_mls[0:sample_length]
    return mls_signal

def permutation_matricies(msl, N):
    m = np.zeros((len(msl),len(msl)), dtype = np.int8)
    for i in range(len(msl)):
        z = np.roll(msl, -i)
        m[i][:] = z[:]
    s = m[:N]
    delta = s[:, :N]
    delta_inv = np.linalg.inv(delta)
    delta_inv = np.array(delta_inv, dtype=bool)
    delta_inv = np.array(delta_inv, dtype=int)
    sT = np.transpose(s)
    l = np.matmul(sT, delta_inv)
    l = l==1
    l = np.array(l, dtype=int)
    return s, l

def bool2int(x):
    y = int("".join(str(int(e)) for e in x),2)
    return y

def find_min_shift(fs, y1, y2):
    """ This function performs small (dt/50) shifts to y2 and finds the minimum least squares

    fs = sampling rate \n
    y1 =input signal - held statioary \n
    y2 = input signal - moved"""

    dt = 1.0/fs
    dtau = dt/50
    shift_range = np.arange(-2*dt, 2*dt, dtau)
    # print("Range of shfts: {} to {} with steps {}".format(-2*dt, 2*dt, dtau))
    N = len(y1)
    ndx = 0
    least_squares = np.zeros(len(shift_range))
    y2f = scipy.fftpack.fft(y2)
    xf = np.linspace(0.0, 1.0/dt, N)
    for shift in shift_range:
        s = np.exp(-2.0j*np.pi*xf*shift)
        y2fs = s*y2f
        y2s = scipy.fftpack.ifft(y2fs).real
        least_squares[ndx] = np.sum((y1-y2s)**2)
        ndx += 1
    min_shift = shift_range[np.argmin(least_squares)]
    # min_shift_ndx = np.argmin(least_squares)

    # print("Shift (s): {}, index: {}, value: {}".format(min_shift, min_shift_ndx, shift_range[min_shift_ndx]))
    s = np.exp(-2.0j*np.pi*xf*min_shift)
    y2fs = s*y2f
    y_shifted = scipy.fftpack.ifft(y2fs).real
    return y_shifted

def align_signals(y1, y2):
    """ Uses correlation to align y2 against y1

    y1 = signal held constant \n
    y2 = shifted signal
    """
    xcorr = correlate(y1, y2, 'full')
    # Use max() to get the sample delay - need to check this will work with the new excitation. Should only have one peak as wormup is not ESS
    i_offset = xcorr.tolist().index(xcorr.max())
    rolled_y2 = np.roll(y2, i_offset)
    return rolled_y2

def subtract_impulse_response(fs, y1, y2):
        # y1 Reflected impulse response
        # y2 Free-feild impulse response
    # Subtracts y2 from y1
    y2_aligned = align_signals(y1=y1, y2=y2)
    y2_shifted = find_min_shift(fs=fs, y1=y1, y2=y2_aligned)
    y2_scaled = (max(y1)/max(y2_shifted))*y2_shifted
    difference = y1 - y2_scaled
    t_array = np.arange(len(y1))*1./51200
    dt = 1.0/fs
    delta = int((0.5*10**-3)/dt)
    peak_loc = np.argmax(y2_scaled)
    # print("Location of peak for subtraction: {}".format(peak_loc))
    min_ndx = peak_loc - delta
    max_ndx = peak_loc + delta+1
    # print("Ratio calculated between {} and {}".format(min_ndx, max_ndx))
    pre_subtraction = np.sum(np.square(np.absolute(y2_scaled[min_ndx:max_ndx])))
    # print("Pre: {}".format(pre_subtraction))
    post_subtraction = np.sum(np.square(np.absolute(difference[min_ndx:max_ndx])))
    t_array = np.arange(len(y1))*1./51200

    # print("Post: {}".format(post_subtraction))
    # print("Ratio: {}".format(pre_subtraction/post_subtraction))
    R_sub = 10*np.log10(pre_subtraction/post_subtraction)
    return difference, R_sub, y2_scaled

def CgeoCalculations(dSM=1.25, dM=0.25, Om=0.4):
    """ Calculates Cgeo (Eqn 2 from EN 1793-5:2016)

    dSM = distance from Mic to speaker (1.25m)
    dM = distance from mic to wall (0.25m)
    Om = orthogonal distance between mics (0.4m)"""

    middle_direct = dSM
    sides_direct = math.sqrt(dSM**2 + Om**2)
    corners_direct = math.sqrt(dSM**2 + 2*Om**2)

    middle_reflected = dSM + dM*2
    sides_reflected = math.sqrt((dSM + 2*dM)**2 + Om**2)
    corners_reflected = math.sqrt((dSM + 2*dM)**2 +2*Om**2)

    middle_Cgeo = (middle_reflected/middle_direct)**2
    sides_Cgeo = (sides_reflected/sides_direct)**2
    corners_Cgeo = (corners_reflected/corners_direct)**2

    D_reflected = {'mic1' : [corners_reflected],
            'mic2' : [sides_reflected],
            'mic3' : [corners_reflected],
            'mic4' : [sides_reflected],
            'mic5' : [middle_reflected],
            'mic6' : [sides_reflected],
            'mic7' : [corners_reflected],
            'mic8' : [sides_reflected],
            'mic9' : [corners_reflected]
            }
    D_direct = {'mic1' : [corners_direct],
            'mic2' : [sides_direct],
            'mic3' : [corners_direct],
            'mic4' : [sides_direct],
            'mic5' : [middle_direct],
            'mic6' : [sides_direct],
            'mic7' : [corners_direct],
            'mic8' : [sides_direct],
            'mic9' : [corners_direct]
            }
    C_geo = {'mic1' : [corners_Cgeo],
            'mic2' : [sides_Cgeo],
            'mic3' : [corners_Cgeo],
            'mic4' : [sides_Cgeo],
            'mic5' : [middle_Cgeo],
            'mic6' : [sides_Cgeo],
            'mic7' : [corners_Cgeo],
            'mic8' : [sides_Cgeo],
            'mic9' : [corners_Cgeo]
            }
    return C_geo, D_direct, D_reflected

def referencePathDifference():
    """ Based on Table 3 from EN 1793-5:2016
    """
    dist_reflected = {'mic1' : 0.122,
            'mic2' : 0.062,
            'mic3' : 0.122,
            'mic4' : 0.062,
            'mic5' : 0.000,
            'mic6' : 0.062,
            'mic7' : 0.122,
            'mic8' : 0.062,
            'mic9' : 0.122
            }
    dist_direct = {'mic1' : 0.467,
            'mic2' : 0.483,
            'mic3' : 0.467,
            'mic4' : 0.483,
            'mic5' : 0.500,
            'mic6' : 0.483,
            'mic7' : 0.467,
            'mic8' : 0.483,
            'mic9' : 0.467
            }
    tolerance = 0.025

    return dist_reflected, dist_direct, tolerance

def checkPathDifference(reflected_path):
    """ Based on Table 3 from EN 1793-5:2016
    """
    dist_reflected, dist_direct, tolerance = referencePathDifference()

