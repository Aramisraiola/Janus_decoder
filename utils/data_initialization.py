from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import csv
import os

# ====================================================
# DATA OPENING AND INITIALIZATION
# ====================================================

def open_data_files(directory_path, file_number, **kwargs):
    file_list = sorted(os.listdir(directory_path))
    df = pd.read_csv('{}/{}'.format(directory_path, file_list[file_number]), **kwargs)
    return df


def prepare_files_list(path, n_files=None):
    file_list = sorted(os.listdir(path))
    file_list_dict = [{'name': path + '/' + file_name, 'delay': 0} for file_name in file_list]
    if n_files is not None:
        file_list_dict = [{'name': path + '/' + file_name, 'delay': 0} for file_name in file_list[:n_files]]
    return file_list_dict


def data_files_summing(file_list, min_ps, max_ps,hist=1):
    x_values = []
    y_values = []

    for file in file_list:
        with open(file['name'], newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
            xx_values = []
            yy_values = []

            for row in spamreader:
                try:
                    float_value = float(row[0])
                except:
                    continue

                if float(row[0]) < min_ps or float(row[0]) > max_ps:
                    continue
                xx_values.append(float(row[0]))
                yy_values.append(float(row[hist]))

        for i in range(len(xx_values)):
            if len(x_values) < i + 1:
                x_values.append(xx_values[i])
                y_values.append(yy_values[i])

            else:
                if x_values[i] != xx_values[i]:
                    raise ("Timestamps not identical in " + file['name'])
                target_bin = i + file['delay']
                if target_bin >= 0 and target_bin < len(x_values):
                    y_values[target_bin] += yy_values[i]
    return x_values, y_values

def open_caen_output(file, rebin, caen_filter=False, bin_res = 3.125):

    df = pd.read_csv(file, delimiter='\t')
    bins = df['Time']
    height = df['Coincidences']


    if caen_filter:
        mask_800 = (bins % 800 != 0)
        bins = bins[mask_800]
        height = height[mask_800]


    bins = bins.tolist()
    height = height.tolist()

    _, _, time, coincidences = rebin_data_etienne(bins, height, rebin_factor=rebin, rebin_offset=0,
                                                  binning_res_ps=bin_res)

    return time,coincidences
def open_data_files_series(directory_path, files_to_omit=None, **kwargs):
    file_list = sorted(os.listdir(directory_path))

    if files_to_omit is not None:
        file_list.pop(files_to_omit)

    file_list = np.array(file_list)

    frames = []
    for file in file_list:
        frames.append(pd.read_csv('{}/{}'.format(directory_path, file), **kwargs))

    columns_to_sum = [df['Histo01'] for df in frames]

    sum_series = sum(columns_to_sum)
    summed_dataframe = pd.DataFrame({'Time [ps]': frames[0]['BinCenter (ps)'], 'Coincidences []': sum_series})

    return frames, summed_dataframe


# ====================================================
# DATA PROCESSING AND CORRECTION
# ====================================================

def flip_flop_corrector_g2(x, y, num_averages=64, excess_location_ps=5440, binning_res_ps=13):
    # correct for the flip-flop effect on the loaded values
    averages = np.zeros(num_averages)

    while len(y) % num_averages != 0:
        y.pop()
        x.pop()

    band_to_exclude = int((excess_location_ps - x[0]) / binning_res_ps / num_averages)
    for count, value in enumerate(y):
        if int(count / num_averages) != band_to_exclude:
            averages[count % num_averages] += value

    average = 0
    for i in range(num_averages):
        averages[i] /= (len(y) / num_averages) - 1
        average += averages[i]

    average /= num_averages

    corrected_y = []
    centered_y = []
    err_y = np.sqrt(y)
    err_y_corr = []

    for count, value in enumerate(y):
        corrected_y.append((value - averages[count % num_averages] + average) / average)
        centered_y.append(value / average)

        err_y_corr.append(np.abs((err_y[count] - averages[count % num_averages] + average) / average))

    return x, y, centered_y, corrected_y, average, err_y_corr


def flip_flop_corrector_counts(x, y, num_averages=64, excess_location_ps=5440, binning_res_ps=13):
    # correct for the flip-flop effect on the loaded values
    averages = np.zeros(num_averages)

    while len(y) % num_averages != 0:
        y.pop()
        x.pop()

    band_to_exclude = int((excess_location_ps - x[0]) / binning_res_ps / num_averages)
    for count, value in enumerate(y):
        if int(count / num_averages) != band_to_exclude:
            averages[count % num_averages] += value

    average = 0
    for i in range(num_averages):
        averages[i] /= (len(y) / num_averages) - 1
        average += averages[i]

    average /= num_averages

    corrected_y = []
    centered_y = []
    err_y = np.sqrt(y)
    err_y_corr = []

    for count, value in enumerate(y):
        corrected_y.append(value)
        centered_y.append(value)

        err_y_corr.append(np.abs(err_y[count]))

    return x, y, centered_y, corrected_y, average, err_y_corr


def rebin_data(x_data, y_data, multiplication_factor=5):
    x_new = []
    y_new = []

    for i in range(0, len(x_data), multiplication_factor):

        if i + multiplication_factor <= len(x_data):
            x_bin = x_data[..., i:i + multiplication_factor]
            y_bin = y_data[..., i:i + multiplication_factor]
        else:
            x_bin = x_data[..., i:]
            y_bin = y_data[..., i:]

        x_rebinned = np.mean(x_bin)
        y_rebinned = np.mean(y_bin)

        x_new.append(x_rebinned)
        y_new.append(y_rebinned)

    x_new = np.array(x_new)
    y_new = np.array(y_new)

    return x_new[:], y_new[:]


def rebin_data_etienne(x, y, rebin_factor, rebin_offset, binning_res_ps=3.125):
    rebinned_y = np.zeros(int(len(y) / rebin_factor))
    rebinned_x = np.zeros(int(len(x) / rebin_factor))

    for count, value in enumerate(y):
        target_bin = int((count + rebin_offset) / rebin_factor)
        if target_bin < len(rebinned_y):
            rebinned_y[target_bin] += value

    for count, value in enumerate(x):
        target_bin = int((count + rebin_offset) / rebin_factor)
        if target_bin < len(rebinned_x):
            rebinned_x[target_bin] += value

    for i in range(len(rebinned_x)):
        rebinned_x[i] /= rebin_factor
        rebinned_y[i] /= rebin_factor

    for i in range(int(rebin_offset / rebin_factor) + 1):
        rebinned_x = np.delete(rebinned_x, 0)
        rebinned_y = np.delete(rebinned_y, 0)

    histo_x = []
    histo_y = []

    bin_width = binning_res_ps * rebin_factor

    for val in rebinned_x:
        histo_x.append(val - (bin_width / 2))
        histo_x.append(val + (bin_width / 2))

    for val in rebinned_y:
        histo_y.append(val)
        histo_y.append(val)

    return histo_x, histo_y, rebinned_x, rebinned_y


def interactive_data_filtering(x_values, y_coincidences, xlabel='Time [ns]', ylabel='Coincidences'):
    plt.ion()

    fig, ax = plt.subplots()
    ax.plot(x_values, y_coincidences)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.draw()

    print("Click on the plot to place the data threshold to eliminate IDQ correlator drops. Press enter to confirm.")
    thresholds = []
    threshold_line = None

    while True:
        click_data = plt.ginput(1, timeout=-1)

        if not click_data:
            if threshold_line:
                thresholds.append(threshold_line.get_ydata()[0])
            plt.close()  # Close the plot
            break

        x_click, y_click = click_data[0]

        if threshold_line:
            threshold_line.set_visible(False)

        threshold_line = ax.axhline(y=y_click, color='r', linestyle='--')
        plt.draw()

    user_input = input(f"Press Enter to keep {thresholds[-1]}): ")

    # Check if the user provided a new y-value, and update if necessary
    if user_input:
        thresholds[-1] = float(user_input)

    # Further processing with the selected thresholds
    print("Selected thresholds:", thresholds)

    # Turn off interactive mode
    plt.ioff()
    return y_click