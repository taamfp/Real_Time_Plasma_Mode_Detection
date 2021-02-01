# Libraries

import time as time_processing
import math
import numpy as np
import pandas as pd
import scipy.io as platform
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import clear_output
from scipy.signal import spectrogram
from scipy.signal import welch
from scipy.signal import find_peaks
from scipy.signal import peak_widths

# Print options
np.set_printoptions(precision=10)

# Data
data_path = 'file'
data = platform.readsav(data_path)

magnetic_field = data['s']
time = data['t']

# Vector

N_points = len(magnetic_field)
print('Lenght of data \n', N_points)

t_min = float('{0:.15f}'.format(53))
t_max = float('{0:.15f}'.format(60))

acquisition_time = (t_max-t_min)/(N_points)

print('Acquisition time \n', acquisition_time, 's')

time = np.arange(53.0, 60.0, acquisition_time)

# Processing time CPU

start = float('{0:.15f}'.format(time_processing.clock()))

print('Processing time: ', float('{0:.15f}'.format(time_processing.clock()-start)))

# Without plot

dt = (time[-1]-time[0])
N = magnetic_field.size

fs = np.round(N/dt)

print('Frequency Sampling \n', fs, 'Hz')

time_plot = []
magnetic_plot = []

f_mode = []
t_mode = []

nfft = 16384
n_points = 4096
noverlap_x = n_points // 8

start = [1, 2, 3]

print('Starting acquisition...')
time_processing.sleep(10)

for i in range(1, len(start)+1):
    print(i)
    time_processing.sleep(1)
    
for i in range(9000000, len(time)):
    start = float('{0:.15f}'.format(time_processing.clock()))
    clear_output(wait=True)
    time_plot.append(time[i])
    magnetic_plot.append(magnetic_field[i])
    time_processing.sleep(0.0005)
    print('Processing time: ', float('{0:.15f}'.format(time_processing.clock()-start)))
    
    
    # Spectrogram and Power Spectral Density
    if i % n_points == 0 and i>0:
        
        dt = float(time_plot[-1]-time_plot[0])
        N = float(len(magnetic_plot))
        
        fs = N/dt
        
        f, Pxx_spec = welch(magnetic_plot, fs, nperseg=n_points, noverlap=noverlap_x, nfft=nfft, return_onesided=True)
        peaks, _ = find_peaks(Pxx_spec)
        
        fwhm = peak_widths(Pxx_spec, peaks, rel_height=0.5)
        
        left_freq = np.rint(fwhm[2:3][0])
        right_freq = np.rint(fwhm[3:][0])
        widths = fwhm[1:2][0]
        
        left_frequency = []
        right_frequency = []
        width = []
        
        for j in range(len(right_freq)):
            left_frequency.append(f[int(left_freq[j])])
            right_frequency.append(f[int(right_freq[j])])
            width.append(float(widths[j]))
            
        
        ind = list(Pxx_spec[peaks]).index(max(Pxx_spec[peaks]))
        f_mode.append(f[peaks][ind])
        
        print('Maximum Mode Frequency Amplitude: \n', max(Pxx_spec[peaks]))
        print('Maximum Mode Frequency', f[peaks][ind], 'Hz')
        print('Maximum Left Mode Frequency', left_frequency[ind], 'Hz')
        print('Maximum Right Mode Frequency', right_frequency[ind], 'Hz')
        time_processing.sleep(2)
	
# With plot

dt = (time[-1]-time[0])
N = magnetic_field.size

fs = np.round(N/dt)

print('Frequency Sampling \n', fs, 'Hz')


plt.rcParams["figure.figsize"] = (10,3)
plt.rcParams.update({'font.size': 10})


time_plot = []
magnetic_plot = []

f_mode = []
t_mode = []

nfft = 16384
n_points = 4096
noverlap_x = n_points // 8

start = [1, 2, 3]

print('Starting acquisition...')
time_processing.sleep(10)

for i in range(1, len(start)+1):
    print(i)
    time_processing.sleep(1)
    
for i in range(9000000, len(time)):
    start = float('{0:.15f}'.format(time_processing.clock()))
    clear_output(wait=True)
    time_plot.append(time[i])
    magnetic_plot.append(magnetic_field[i])
    
    
    # Plot data
    plt.title('JET #86469')
    plt.xlabel('time (s)')
    plt.ylabel('M')
    plt.plot(time_plot, magnetic_plot)
    plt.pause(0.0005)
    plt.clf()
    print('Processing time: ', float('{0:.15f}'.format(time_processing.clock()-start)))
    print(i)
    
    # Spectrogram and Power Spectral Density
    if i % n_points == 0 and i>0:
        
        dt = float(time_plot[-1]-time_plot[0])
        N = float(len(magnetic_plot))
        
        fs = N/dt
        
        f, Pxx_spec = welch(magnetic_plot, fs, nperseg=n_points, noverlap=noverlap_x, nfft=nfft, return_onesided=True)
        peaks, _ = find_peaks(Pxx_spec)
        
        fwhm = peak_widths(Pxx_spec, peaks, rel_height=0.5)
        
        left_freq = np.rint(fwhm[2:3][0])
        right_freq = np.rint(fwhm[3:][0])
        widths = fwhm[1:2][0]
        
        left_frequency = []
        right_frequency = []
        width = []
        
        for j in range(len(right_freq)):
            left_frequency.append(f[int(left_freq[j])])
            right_frequency.append(f[int(right_freq[j])])
            width.append(float(widths[j]))
            
        
        plt.figure()
        plt.semilogx(f, Pxx_spec)
        plt.plot(f[peaks], Pxx_spec[peaks], 'v')
        plt.hlines(widths, left_frequency, right_frequency, color="C2")
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PDS')
        plt.grid('True')
        plt.grid(which='minor')
        plt.show()
        
        ind = list(Pxx_spec[peaks]).index(max(Pxx_spec[peaks]))
        f_mode.append(f[peaks][ind])
        
        
        
        plt.specgram(magnetic_plot, Fs=fs, NFFT=n_points, noverlap=n_points//8, cmap='jet')
        plt.xlabel('t (s)')
        plt.ylabel('f (Hz)')
        plt.ylim(0, 40000)
        plt.axhline(left_frequency[ind], color='red')
        plt.axhline(f[peaks][ind], color='green')
        plt.axhline(right_frequency[ind], color='blue')
        plt.pause(2)
        plt.clf()