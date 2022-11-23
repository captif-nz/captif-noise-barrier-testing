from bokeh.plotting import figure, output_file, show, save, reset_output
from bokeh.models import ColumnDataSource, BasicTickFormatter, DatetimeTickFormatter, BoxAnnotation, Label, NumeralTickFormatter
import bokeh.models as bmo
from bokeh.layouts import row, column
from bokeh.palettes import d3
from bokeh.io import export_png


import numpy as np
import EN_1793_funcs as en1793
import octave_filters
import json
import matplotlib.pyplot as plt
import pickle
import time
from scipy.signal import correlate

# # ESS Signal
# SR = 51200
# t_meas = 30
# n_runs = 5
# f0 = 40
# f1 = 20000
# test_signal, ess_sweep, test_t, ess_t = en1793.ess_signal(t1=t_meas,
#                                                           f0=f0,
#                                                           f1=f1,
#                                                           SR=SR,
#                                                           scale=0.2)

# output_file("test.html")
# p1 = figure(x_axis_label='Time (seconds)', y_axis_label='Generator voltage (V)',
#             plot_width=1000, plot_height=300, title="Signal for ESS input",
#             y_range=[-0.25, 0.25])
# p1.line(x=test_t, y=test_signal, line_width=1)
# p1.circle(x=test_t, y=test_signal, size=0.2)
# warmup_box = BoxAnnotation(left=0, right=8, fill_alpha=0.1, fill_color='red')
# p1.add_layout(warmup_box)
# ess_box = BoxAnnotation(left=8, right=28, fill_alpha=0.1, fill_color='green')
# p1.add_layout(ess_box)
# warmuplabel = Label(x=4, y=0.1, text='Warmup',background_fill_color='white', background_fill_alpha=1.0)
# p1.add_layout(warmuplabel)
# ess_label = Label(x=18, y=0.1, text='ESS', background_fill_color='white', background_fill_alpha=1.0)
# p1.add_layout(ess_label)
# show(p1)

# Plot correlation
# results_file = "test20190122/bottom_CAPTIF_wall_20190122"

# with open(results_file, 'rb') as f: #, open(middle_file) as middle, open(top_file) as top:
#     data = pickle.load(f)

# measurement = data['mic0'][0]

# ess_sweep = data['ess_sweep']

# print("Length of measurement: {}".format(len(measurement)))
# print("Length of ESS sweep: {}".format(len(ess_sweep)))

# # # find offset between the input and output signals
# xcorr = correlate(measurement, ess_sweep, 'valid')
# steps = np.arange(len(xcorr))
# i_offset = xcorr.tolist().index(xcorr.max())
# peak = xcorr[i_offset]
# print("Max: {}, index: {}".format(peak, i_offset))
# output_file("test.html")
# p2 = figure(x_axis_label='Step', y_axis_label='Cross correlation amplitude',
#             plot_width=1000, plot_height=300, title="Cross correlation between measurement and ESS sweep",
#             x_range=[0,len(xcorr)], y_range=[-10000,10000])
# p2.line(x=steps, y=xcorr)
# p2.annulus(x=i_offset, y=peak, inner_radius=1500, outer_radius=2000, color="red")
# p2.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
# p2.yaxis[0].formatter = NumeralTickFormatter(format="0,0")
# show(p2)


## plot free feild and wall data
wall_file = "test20190122/bottom_CAPTIF_wall_20190122"
ff_file = "test20190122/bottom_CAPTIF_ff_20190122"
with open(wall_file, 'rb') as f: #, open(middle_file) as middle, open(top_file) as top:
    wall_data = pickle.load(f)

with open(ff_file, 'rb') as f: #, open(middle_file) as middle, open(top_file) as top:
    ff_data = pickle.load(f)

wall_measurement = wall_data['mic0'][2]
wall_t = wall_data['test_t']
ff_measurement = ff_data['mic0'][2]
ff_t = ff_data['test_t']
print(len(ff_measurement))
print(len(ff_t))

# output_file("test.html")
# p3 = figure(x_axis_label='Time (s)', y_axis_label='Amplitude',
#             plot_width=1000, plot_height=300, title="Free-feild measurement",
#             x_range=[0,30], y_range=[-15,15])
# p3.line(x=ff_t, y=ff_measurement)

# p4 = figure(x_axis_label='Time (s)', y_axis_label='Amplitude',
#             plot_width=1000, plot_height=300, title="Wall measurement",
#             x_range=[0,30], y_range=[-15,15])
# p4.line(x=wall_t, y=wall_measurement)

# show(column(p3, p4))

def returnIR(measurement, data):

    xcorr = correlate(measurement, data['ess_sweep'], 'valid')
    i_offset = xcorr.tolist().index(xcorr.max())
    measurement_clipped = measurement[i_offset: len(data['ess_sweep'])+i_offset]
    # Calculate impulse response
    ir_temp = en1793.impulse_response(data['ess_sweep'], measurement_clipped)
    # Roll to avoid splitting IR at 0
    ir_temp2 = np.roll(ir_temp, int(0.01*51200))
    ir_clipped = ir_temp2[0:int(0.05*data['SR'])]
    # print(ir_clipped.mean())
    return ir_clipped-ir_clipped.mean()

ff_ir = returnIR(measurement=ff_measurement, data=ff_data)
ff_ir_t = np.arange(len(ff_ir))/ff_data['SR']

wall_ir = returnIR(measurement=wall_measurement, data=wall_data)
wall_ir_t = np.arange(len(wall_ir))/wall_data['SR']

# output_file("test.html")
# p5 = figure(x_axis_label='Time (s)', y_axis_label='Amplitude',
#             plot_width=1000, plot_height=500, title="Impulse response",
#             x_range=[0,0.05], y_range=[-2,2])
# p5.line(x=ff_ir_t, y=ff_ir, color="blue", legend='Free-feild')
# p5.line(x=wall_ir_t, y=wall_ir, color="red", legend='CAPTIF shed wall')
# show(p5)

ff_windowed = en1793.apply_adrienne(ff_ir, offset=0.008, SR=ff_data['SR'])
wall_windowed = en1793.apply_adrienne(wall_ir, offset=0.008, SR=wall_data['SR'])
tt = np.array(range(len(ff_windowed)))*1./51200


output_file("test.html")
p5 = figure(x_axis_label='Time (s)', y_axis_label='Amplitude',
            plot_width=1000, plot_height=500, title="Impulse response - Adrienne windowed",
            x_range=[0,0.05], y_range=[-2,2])
p5.line(x=tt, y=ff_windowed, color="blue", legend='Free-feild')
p5.line(x=tt, y=wall_windowed, color="red", legend='CAPTIF shed wall')
show(p5)