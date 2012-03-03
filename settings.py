# settings.py

#antialias mainCanvas plot (False will speed up rendering)
antialiased = True

#arbitrary scale of wavefunction size (height) on plots
wf_scale = 4.5e-10

#minimum height for wavefunctions (wf_scale * wf_min_height = 0.014)
# affects states that are far above band edge
# setting to zero will show all states
wf_min_height = 3.1e7

#wavefunction height to begin pretty plot cutoff
# wf_scale * pretty_plot_factor = 5e-3
pretty_plot_factor = 1.2e6

#wavefunction height to begin intgration for LO phonon lifetime
phonon_integral_factor = pretty_plot_factor/10

#decimate factor depends on xres
# decimate factor should be an integer
# represents plotting wavefunctions at every xth angstrom
plot_decimate_factor = 1