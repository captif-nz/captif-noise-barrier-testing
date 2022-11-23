import numpy as np
import EN_1793_funcs as en1793
import octave_filters
import json
import matplotlib.pyplot as plt
import pickle
import time
from bokeh.plotting import figure, output_file, show, gridplot
import os

base_dir = "BarrierTesting"

        # "MaioroSt_Post_20190625_FF",
        # "TeAtatu_Panel_Reflections_20190625",
        #                 "TeAtatu_Panel_Reflections_20190625_FF",
        #                 "TeAtatu_PanelTL_20190625",
# measurement_names = ["TeAtatu_Post_Reflections_20190625",
                        # "TeAtatu_PostTL_20190625"]
measurement_names = ["MaioroSt_Post_20190625_FF"]


for measurement_name in measurement_names:
        measurement_dir = os.path.join(base_dir, measurement_name)
        print("Loading data from {}".format(measurement_dir))
        with open(os.path.join(measurement_dir, 'IR'), 'rb') as f:
                data = pickle.load(f)
                avgIR = data['averagedImpulseResponse']

        print(avgIR)
        print('Plotting Measurement Results')
        # output to static HTML file
        output_file(os.path.join(measurement_dir, "ImpulseResponse.html"))

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
