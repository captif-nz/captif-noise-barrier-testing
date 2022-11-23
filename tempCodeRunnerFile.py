 to static HTML file
# output_file(os.path.join(dest_folder, "ReflectionsFF.html"))
# width = 600
# height = 300
# line_w = 1
# p1_1 = figure(plot_width=width, plot_height=height, title="Mic 1 - windowed reflection")
# p1_1.xgrid.minor_grid_line_color = 'grey'
# p1_1.xgrid.minor_grid_line_alpha = 0.1
# p1_1.line(ff_xf, freefeild_fft['mic1'], line_width=line_w, line_dash='dashed', line_color="blue")
# p1_1.line(ff_xf, wall_fft_store['mic1'], line_width=line_w, line_color="red")
# p2_1 = figure(plot_width=width, plot_height=height, title="Mic 2 - windowed reflection")
# p2_1.xgrid.minor_grid_line_color = 'grey'
# p2_1.xgrid.minor_grid_line_alpha = 0.1
# p2_1.line(ff_xf, freefeild_fft['mic2'], line_width=line_w, line_dash='dashed', line_color="blue")
# p2_1.line(ff_xf, wall_fft_store['mic2'], line_width=line_w, line_color="red")
# p3_1 = figure(plot_width=width, plot_height=height, title="Mic 3 - windowed reflection")
# p3_1.xgrid.minor_grid_line_color = 'grey'
# p3_1.xgrid.minor_grid_line_alpha = 0.1
# p3_1.line(ff_xf, freefeild_fft['mic3'], line_width=line_w, line_dash='dashed', line_color="blue")
# p3_1.line(ff_xf, wall_fft_store['mic3'], line_width=line_w, line_color="red")
# p4_1 = figure(plot_width=width, plot_height=height, title="Mic 4 - windowed reflection")
# p4_1.xgrid.minor_grid_line_color = 'grey'
# p4_1.xgrid.minor_grid_line_alpha = 0.1
# p4_1.line(ff_xf, freefeild_fft['mic4'], line_width=line_w, line_dash='dashed', line_color="blue")
# p4_1.line(ff_xf, wall_fft_store['mic4'], line_width=line_w, line_color="red")
# p5_1 = figure(plot_width=width, plot_height=height, title="Mic 5 - windowed reflection")
# p5_1.xgrid.minor_grid_line_color = 'grey'
# p5_1.xgrid.minor_grid_line_alpha = 0.1
# p5_1.line(ff_xf, freefeild_fft['mic5'], line_width=line_w, line_dash='dashed', line_color="blue")
# p5_1.line(ff_xf, wall_fft_store['mic5'], line_width=line_w, line_color="red")
# p6_1 = figure(plot_width=width, plot_height=height, title="Mic 6 - windowed reflection")
# p6_1.xgrid.minor_grid_line_color = 'grey'
# p6_1.xgrid.minor_grid_line_alpha = 0.1
# p6_1.line(ff_xf, freefeild_fft['mic6'], line_width=line_w, line_dash='dashed', line_color="blue")
# p6_1.line(ff_xf, wall_fft_store['mic6'], line_width=line_w, line_color="red")
# p7_1 = figure(plot_width=width, plot_height=height, title="Mic 7 - windowed reflection")
# p7_1.xgrid.minor_grid_line_color = 'grey'
# p7_1.xgrid.minor_grid_line_alpha = 0.1
# p7_1.line(ff_xf, freefeild_fft['mic7'], line_width=line_w, line_dash='dashed', line_color="blue")
# p7_1.line(ff_xf, wall_fft_store['mic7'], line_width=line_w, line_color="red")
# p8_1 = figure(plot_width=width, plot_height=height, title="Mic 8 - windowed reflection")
# p8_1.xgrid.minor_grid_line_color = 'grey'
# p8_1.xgrid.minor_grid_line_alpha = 0.1
# p8_1.line(ff_xf, freefeild_fft['mic8'], line_width=line_w, line_dash='dashed', line_color="blue")
# p8_1.line(ff_xf, wall_fft_store['mic8'], line_width=line_w, line_color="red")
# p9_1 = figure(plot_width=width, plot_height=height, title="Mic 9 - windowed reflection")
# p9_1.xgrid.minor_grid_line_color = 'grey'
# p9_1.xgrid.minor_grid_line_alpha = 0.1
# p9_1.line(ff_xf, freefeild_fft['mic9'], line_width=line_w, line_dash='dashed', line_color="blue")
# p9_1.line(ff_xf, wall_fft_store['mic9'], line_width=line_w, line_color="red")
# p = gridplot([[p9_1, p8_1, p7_1],
#             [p6_1, p5_1, p4_1],
#             [p3_1, p2_1, p1_1]])
# # show the results
# show(p)