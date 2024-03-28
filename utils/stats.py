import numpy as np


# ====================================================
# FUNCTIONS/PDFs
# ====================================================

def gaussian_function(x, mu, sigma, amplitude):
    return amplitude * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    return y


def gaussian_function_bkg(x, mu, sigma, amplitude, c):
    return amplitude * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) + c

    return y


def lorentzian_bkg(x, mu, gamma, A, c):
    return A * gamma*0.5 / ((x - mu) ** 2 + (np.abs(gamma)/2) ** 2) + c


def lorentzian(x, mu, gamma, A):
    return A * (gamma * 0.5 / ((x - mu) ** 2 + (np.abs(gamma) / 2) ** 2))


def sqrt_bkg(x, A, c):
    return (A * np.sqrt(x)) + c


def neg_exponential(t, t_0, t_c, A):
    return A * np.exp(-np.abs(2*(t - t_0) / t_c))

def neg_exponential_bkg(t, t_0, t_c, A, c):
    return A * np.exp(-np.abs(2*(t - t_0) / t_c))+c


# ====================================================
# STATISTICS
# ====================================================

def FWHM(x, y):
    diff = x[1] - x[0]
    half_max = max(y) / 2.
    l = np.where(y > half_max, 1, 0)

    return np.sum(l) * diff

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

