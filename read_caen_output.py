import pandas as pd
from utils.visualization import *
from utils.data_initialization import *
from utils.stats import *
import sys


def chi_squared(observed, expected, errors):
    observed = np.array(observed)
    expected = np.array(expected)
    errors = np.array(errors)

    # Calculate chi-squared
    residuals = (observed - expected) / errors
    chi_squared_value = np.sum(residuals ** 2)

    return chi_squared_value


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Must run the scrpit as: python read_caen_output.py <filename>")
        sys.exit(1)

    file = sys.argv[1]

    rebin = 5
    min = 0
    max = 8000

    time, coincidences = open_caen_output(file, rebin=rebin, caen_filter=True)

    mask = (time <= max) * (time >= min)





    ax = plot_coincidences_error(time, coincidences, rebin=rebin, t_min=0, t_max=4000, g2=True,
                                 title='Cross-correlation between channels 11 and 6 for Na-lamp light through monomode fiber \n and with polarizer filter. The acquisition was performed with Caen TDC (A5203) during 10 minutes.')

    coincidences_baseline = np.mean(coincidences[time < 1400])
    coincidences = coincidences / coincidences_baseline
    max_fit = 4000
    min_fit = 0

    lower_bound_1 = 0
    upper_bound_1 = 1400
    lower_bound_2 = 2100
    upper_bound_2 = 4000

    # Create masks for each interval
    mask_1 = (time >= lower_bound_1) & (time <= upper_bound_1)
    mask_2 = (time >= lower_bound_2) & (time <= upper_bound_2)

    # Combine masks using logical OR
    final_mask = mask_1 | mask_2

    # Apply the final mask to coincidences
    coincidences_baseline = coincidences[final_mask]
    # coincidences = coincidences / np.mean(coincidences_baseline)

    range_fit = (time <= max_fit) * (time >= min_fit)
    mu = time[range_fit][np.argmax(coincidences[range_fit])]
    x_plot = np.arange(np.amin(time[range_fit]), np.amax(time[range_fit]), 0.01)
    sigma = 0.64 / (np.amax(coincidences[range_fit]) - 1)

    popt, pcov = curve_fit(lorentzian_bkg, time[range_fit], coincidences[range_fit],
                           p0=[mu, sigma, np.amax(coincidences[range_fit]) - 1,
                               1], maxfev=8000000, sigma=np.sqrt(coincidences[range_fit]))

    confidence_sigma = 5

    max_error = lorentzian_bkg(x_plot, popt[0] + confidence_sigma * np.sqrt(pcov[0, 0]),
                               popt[1] + confidence_sigma * np.sqrt(pcov[1, 1]),
                               popt[2] + confidence_sigma * np.sqrt(pcov[2, 2]),
                               popt[3] + confidence_sigma * np.sqrt(pcov[3, 3]))

    min_error = lorentzian_bkg(x_plot, popt[0] - confidence_sigma * np.sqrt(pcov[0, 0]),
                               popt[1] - confidence_sigma * np.sqrt(pcov[1, 1]),
                               popt[2] - confidence_sigma * np.sqrt(pcov[2, 2]),
                               popt[3] - confidence_sigma * np.sqrt(pcov[3, 3]))

    ax.plot(x_plot, lorentzian_bkg(x_plot, *popt), label='Lorentzian fit', color='red',linewidth=2)

    ax.fill_between(x_plot, y1=max_error, y2=min_error, label=r'{}$\sigma$ Confidence band'.format(confidence_sigma),
                    color='red', alpha=0.2,
                    edgecolor='none')

    ax.plot([], [], label='$\mu$ = {:.2f}$\pm${:.2f} ps'.format(popt[0], np.sqrt(pcov[0, 0])), color='none')
    ax.plot([], [], label='$\Gamma$ = {:.2f}$\pm${:.2f} ps'.format(popt[1], np.sqrt(pcov[1, 1])), color='none')
    ax.plot([], [], label='$A$ = {:.2f}$\pm${:.2f}'.format(popt[2], np.sqrt(pcov[2, 2])), color='none')
    """ax.plot([], [], label=r'$\chi^2$/dof = {:.2f}'.format(
        chi_squared(coincidences[range_fit], lorentzian_bkg(time[range_fit], *popt),
                    np.sqrt(coincidences[range_fit])) / (coincidences[range_fit].shape[0] - len(popt))), color='none')"""
    print(coincidences[range_fit], lorentzian_bkg(time[range_fit], *popt))
    ax.legend()

    _, popt, _ = plot_y_projection(time, coincidences, t_min=0, t_max=1.2e4, fit_histogram=True,
                                   plot_fit_results=True, bins=100)

    time, coincidences = plot_fft(time[mask], coincidences[mask], 10, 200)

    plot_coincidences(time, coincidences, rebin=10, t_min=0, t_max=1.2e4, g2=True,
                      title='Cross-correlation ch 11 and 6 after fft filtering')

    _, popt, _ = plot_y_projection(time, coincidences, t_min=0, t_max=1.2e4, fit_histogram=True,
                                   plot_fit_results=True, bins=100, title='Counts y-projection after FFT filtering')

    plt.show()
