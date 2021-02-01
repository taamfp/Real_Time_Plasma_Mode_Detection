# Libraries

import time as time_processing
import math
import numpy as np
import pandas as pd
import scipy.io as platform
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.signal import welch
from scipy.signal import find_peaks
from scipy.signal import peak_widths


np.set_printoptions(precision=10)

# Data
data_path = 'file'
data = platform.readsav(data_path)


magnetic_field = data['i803']
time = data['t_i803']

# Vector

N_points = len(magnetic_field)
print('Lenght of data \n', N_points)

t_min = float('{0:.15f}'.format(53))
t_max = float('{0:.15f}'.format(60))

acquisition_time = (t_max-t_min)/(N_points)

print('Acquisition time \n', acquisition_time, 's')

time = np.arange(53.0, 60.0, acquisition_time)

# Plot t vs M

plt.style.use('classic')

plt.rcParams["figure.figsize"] = (20,15)

plt.rcParams.update({'font.size': 15})


plt.plot(time, magnetic_field)
plt.title('JET #86469')
plt.xlabel('time (s)')
plt.ylabel('M')
plt.show()


dt = (time[-1]-time[0])
N = magnetic_field.size

fs = np.round(N/dt)
print('Frequency Sampling \n', fs, 'Hz')

# Parameters
nfft = 16384
nperseg_x = 4096
noverlap_x = nperseg_x // 8

return_onesided_m = True


f, t, S = spectrogram(magnetic_field, fs=fs, window='hann', nperseg=nperseg_x, noverlap=2082, nfft=nfft, return_onesided=return_onesided_m)
if (return_onesided_m - 1):
        f = np.fft.fftshift(f)
        S = np.fft.fftshift(S, axes=0)
        
plt.pcolormesh(time[np.round(t*fs).astype(int)], f, np.log(S), shading='auto', cmap='jet')
plt.xlabel('t (s)')
plt.ylabel('f (Hz)')
plt.ylim(0, 40000)
plt.show()

plt.pcolormesh(time[np.round(t*fs).astype(int)], f, np.log(S), shading='auto', cmap='jet')
plt.xlabel('t (s)')
plt.ylabel('f (Hz)')
plt.xlim(57.49, 57.51)
plt.ylim(0, 20000)
plt.hlines(f_mode, xmin=57.49, xmax=57.5, color='green', linewidth=5)
plt.show()


f, Pxx_spec = welch(magnetic_field, fs, nperseg=nperseg_x, noverlap=noverlap_x, nfft=nfft, return_onesided=return_onesided_m)

# Example 0.000000025

height = input('Enter peak height:')

peaks, _ = find_peaks(Pxx_spec, height = float(height))


plt.figure()
plt.semilogx(f, Pxx_spec)
plt.plot(f[peaks], Pxx_spec[peaks], 'v')
plt.xlabel('frequency [Hz]')
plt.ylabel('PDS')
plt.grid('True')
plt.grid(which='minor')
plt.show()
time_processing.sleep(5)
plt.clf()


# Removing the baseline according to expected frequencies

# Example 5000

threshold = input('Enter threshold frequency in Hz:')

for i in range(len(f)):
    if f[i] <= int(threshold):
        Pxx_spec[i] = 0
    else:
        None

plt.figure()
plt.semilogx(f, Pxx_spec)
plt.xlabel('frequency [Hz]')
plt.ylabel('PDS')
plt.grid('True')
plt.grid(which='minor')
plt.savefig('P_Spectrum_Baseline2.png')
plt.show()
time_processing.sleep(5)
plt.clf()


height = input('Enter peak height:')

peaks, _ = find_peaks(Pxx_spec, height = float(height))


plt.figure()
plt.semilogx(f, Pxx_spec)
plt.plot(f[peaks], Pxx_spec[peaks], 'v', markersize=10)
plt.xlabel('frequency [Hz]')
plt.ylabel('PDS')
plt.grid('True')
plt.grid(which='minor')
plt.savefig('Signal_Peak_finder2.png')
plt.show()
time_processing.sleep(5)
plt.clf()


fwhm = peak_widths(Pxx_spec, peaks, rel_height=0.3)

left_freq = np.rint(fwhm[2:3][0])
right_freq = np.rint(fwhm[3:][0])
widths = fwhm[1:2][0]

left_frequency = []
right_frequency = []
width = []

for i in range(len(right_freq)):
    left_frequency.append(f[int(left_freq[i])])
    right_frequency.append(f[int(right_freq[i])])
    width.append(float(widths[i]))
    
plt.figure()
plt.semilogx(f, Pxx_spec)
plt.plot(f[peaks], Pxx_spec[peaks], 'v', markersize=10)
plt.xlabel('frequency [Hz]')
plt.ylabel('PDS')
plt.hlines(widths, left_frequency, right_frequency, color="C2")
plt.xlim(5000, 50000)
plt.savefig('Signal_Peak_FWHM2.png')
plt.grid('True')
plt.grid(which='minor')
plt.show()