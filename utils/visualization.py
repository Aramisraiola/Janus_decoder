import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from utils.stats import *


# ====================================================
# Coincidence plotting / Autocorrelation plotting
# ====================================================

def plot_y_projection(time, coincidences, ax=None, t_max=None, t_min=None, fit_histogram=False, plot_fit_results=False,
                      bins=50, title='Counts histogram (y projection)'):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot()

    if (t_max is not None) * (t_min is not None):
        mask = (time <= t_max) * (time >= t_min)
        coincidences = coincidences[mask]
        time = time[mask]

    n, bins, _ = ax.hist(coincidences, bins=bins, label='Counts histogram ({} bins)'.format(bins))

    if fit_histogram:
        bin_centers = (bins[:-1] + bins[1:]) / 2
        popt, pcov = curve_fit(gaussian_function, bin_centers, n,
                               p0=[np.mean(coincidences), np.std(coincidences), np.amax(n)])

        if plot_fit_results:
            x = np.arange(bin_centers.min(), bin_centers.max(), 1)
            plt.plot(x, gaussian_function(x, *popt), label='Gaussian fit')
            plt.plot([], [], label='$\mu$ = {:.2f}$\pm${:.2f}'.format(popt[0], np.sqrt(pcov[0, 0])), color='none')
            plt.plot([], [], label='$\sigma$ = {:.2f}$\pm${:.2f}'.format(popt[1], np.sqrt(pcov[1, 1])), color='none')
            plt.plot([], [], label='sqrt($\mu$)/$\sigma$ = {:.2f}'.format(np.sqrt(popt[0]) / popt[1]), color='none')

    ax.set_xlabel('Counts [ ]')
    ax.set_ylabel('Frequency [ ]')
    ax.set_title(title)
    ax.legend(loc='best')
    if fit_histogram:
        return ax, popt, pcov
    else:
        return ax, 0, 0


def plot_coincidences(time, coincidences, rebin=1, ax=None, t_min=None, t_max=None, g2=False, delay=0,
                      title='Coincidences over time', label=None):
    time = np.array(time) + delay
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot()

    if (t_max is not None) * (t_min is not None):
        mask = (time <= t_max) * (time >= t_min)
        coincidences = coincidences[mask]
        time = time[mask]
    if g2:
        coincidences_g2 = coincidences / np.mean(coincidences)
        ax.plot(time, coincidences_g2,
                label='Rebinned data (factor = {})'.format(rebin),
                color='#007E64')
        ax.set_ylabel('g$^2(\\tau)$ [ ]')
        ax.legend()



    else:
        ax.plot(time, coincidences,
                label='Rebinned data (factor = {})'.format(rebin),
                color='#007E64')
        ax.set_ylabel('Coincidences [ ]')
    ax.set_title(title)
    ax.set_xlabel('Time lag $\\tau$ [ps]')
    ax.set_ylabel('Counts [ ]')

    return ax


def plot_coincidences_error(time, coincidences, rebin=1, ax=None, t_min=None, t_max=None, g2=False, delay=0,
                            title='Coincidences over time'):
    time = np.array(time) + delay
    error = np.sqrt(coincidences)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot()

    if (t_max is not None) * (t_min is not None):
        mask = (time <= t_max) * (time >= t_min)
        coincidences = coincidences[mask]
        time = time[mask]
        error = error[mask]
    if g2:
        coincidences_g2 = coincidences / np.mean(coincidences)
        ax.errorbar(time, coincidences_g2, yerr=error / np.mean(coincidences), fmt='o', color='black', ecolor='grey',
                    elinewidth=1,
                    label='Rebinned data (factor = {})'.format(rebin))
        ax.set_ylabel('g$^2(\\tau)$ [ ]')
        ax.legend()

    else:
        ax.errorbar(time, coincidences, yerr=error, fmt='o', color='black',
                    ecolor='grey', elinewidth=1,
                    label='Rebinned data (factor = {})'.format(rebin), zorder=0)
        ax.set_ylabel('Coincidences [ ]')
    ax.set_title(title)
    ax.set_xlabel('Time lag $\\tau$ [ps]')
    ax.set_ylabel('Counts [ ]')

    return ax


# ====================================================
# Signal/Coincidences Fast Fourirer Transform (FFT)
# ====================================================

def filter_fft(freq, fft_result, min, max):
    fft_result[(freq > min) & (freq < max)] = 0
    fft_result[(freq < -min) & (freq > -max)] = 0
    return fft_result


def plot_fft(time, coincidences, filter_min=None, filter_max=None):
    fft_result = np.fft.fft(coincidences)
    fft_freq = np.fft.fftfreq(len(time), np.mean(np.diff(time))) * 1000

    title_2 = 'Signal FFT (no frequency filter)'

    if (filter_min is not None) * (filter_max is not None):
        fft_result = filter_fft(fft_freq, fft_result, filter_min, filter_max)
        title_2 = 'Signal FFT (frequency cut from {} to {} GHz)'.format(filter_min, filter_max)

    ifft_result = np.fft.ifft(fft_result)

    plt.figure(figsize=(12, 10))
    plt.subplot(3, 1, 1)
    plt.plot(time, coincidences)
    plt.title('Original Signal')
    plt.xlabel('Time [ps]')
    plt.ylabel('Coincidences [ ]')

    plt.subplot(3, 1, 2)
    plt.stem(fft_freq[fft_freq != 0], np.abs(fft_result)[fft_freq != 0], markerfmt=" ")
    plt.title(title_2)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Amplitude [ ]')

    plt.subplot(3, 1, 3)
    plt.plot(time, np.real(ifft_result))
    plt.title('IFFT (Reconstructed Signal)')
    plt.xlabel('Time [ps]')
    plt.ylabel('Coincidecnes')

    plt.tight_layout()
    return time, np.real(ifft_result)
