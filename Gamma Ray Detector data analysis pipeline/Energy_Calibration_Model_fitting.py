
import numpy as np
from scipy.optimize import curve_fit
from numpy import inf as INF
import matplotlib.pyplot as plt


TWO_PI = np.pi * 2


# function to display the region of interest using axvspan
def axes_axvspan(ax, region):
    """
    Generate a vertical span on a given axes object
    :param ax: axes object
    :param region: tuple containing bounds of vspan
    :return: Vertical span within given bounds
    """
    return ax.axvspan(*region, alpha=0.3, facecolor='r', edgecolor='black')


# plot the spectrum between channels 0 and slice_index
def plot_spectrum(ax, channels, count_rate, slice_index, source, detector, roi, xlabel='Channels',
                  ylabel='Count Rate', **kwargs):
    """
    Plots the channels vs count_rate up to slice index
    - Adapted from 'plot_spectrum' function in the curve fitting notebook provided with the lab
    :param ax: axes object
    :param channels: list of channels
    :param count_rate: list of count_rate
    :param slice_index: only plot data at indices between 0 and this value
    :param source: name of the source
    :param detector: name of the detector
    :param roi: region of interest
    :param xlabel: x axis label
    :param ylabel: y axis label
    :param kwargs: keyword arguments
    :return: scatter plot of channels vs count_rate
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    axes_axvspan(ax, roi)
    ax.grid(True)
    ax.set_title(f"{detector} spectrum for {source}")
    return ax.scatter(channels[:slice_index], count_rate[:slice_index], **kwargs)


# plot a zoomed in version of the spectrum
def zoom_in(ax, channels, count_rate, zoom_indices, source, detector, roi,
            xlabel='Channels', ylabel='Count Rate', **kwargs):
    """
    Zoomed in plot on region of interest
    :param ax: axes object
    :param channels: list of channels
    :param count_rate: list of count_rate
    :param zoom_indices: indices within which to plot data
    :param source: name of the source
    :param detector: name of the detector
    :param roi: region of interest
    :param xlabel: x axis label
    :param ylabel: y axis label
    :param kwargs: keyword arguments
    :return: scatter plot of channels vs count_rate
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    axes_axvspan(ax, roi)
    ax.grid(True)
    ax.set_title(f"{detector} spectrum for {source}")
    return ax.scatter(channels[zoom_indices[0]:zoom_indices[1]], count_rate[zoom_indices[0]:zoom_indices[1]], **kwargs)


# gaussian function
def gaussian_func(x, mu, sigma, area):
    """
    Gaussian function
    - taken from 'gaussian' function in the curve fitting notebook provided with the lab
    :param x: list of values
    :param mu: mean
    :param sigma: std deviation
    :param area: area
    :return: Gaussian distribution
    """
    return area * np.exp(-0.5 * (x-mu)**2 / sigma**2) / np.sqrt(TWO_PI * sigma**2)


# mask with True for x in [xmin, xmax)
def roi_mask(x, xmin, xmax):
    """
    Boolean mask with value True for x in [xmin, xmax)
    - Taken from 'in_interval' function in the curve fitting notebook provided with the lab
    :param x: list of values
    :param xmin: lower region of interest bound
    :param xmax: upper region of interest bound
    :return: Boolean mask
    """
    _x = np.asarray(x)
    return np.logical_and(xmin <= _x, _x < xmax)


# return elements of x and y for which x falls within [xmin, xmax)
def filter_with_mask(x, y, xmin, xmax):
    """
    Filters values in x and y where xmin <= x < xmax
    - Taken from 'filter_in_interval' function in the curve fitting notebook provided with the lab
    :param x: list of x values
    :param y: list of y values
    :param xmin: lower region of interest value
    :param xmax: upper region of interest value
    :return: filtered arrays
    """
    try:
        assert len(x) == len(y)
    except AssertionError:
        print("Assertion Exception Raised.")
    else:
        _mask = roi_mask(x, xmin, xmax)
        return [np.asarray(x)[_mask] for x in (x, y)]


# estimate centroid
def first_moment(x, y):
    """
    Estimate the centroid of the gaussian peak
    - Taken from 'first_moment' function in the curve fitting notebook provided with the lab
    :param x: list of x values
    :param y: list of y values
    :return: centroid estimate
    """
    return np.sum(x * y) / np.sum(y)


# estimate square of the width of the peak
def second_moment(x, y):
    """
    Estimate the width of the gaussian peak squared
    - Taken from 'second_moment' function in the curve fitting notebook provided with the lab
    :param x: list of x values
    :param y: list of y values
    :return: FWHM^2 estimate
    """
    x0 = first_moment(x, y)
    return np.sum((x - x0)**2 * y) / np.sum(y)


# return estimates of parameters of gaussian function
def initial_gaussian_estimates(channels, count_rate):
    """
    Calculates and returns gaussian estimates
    - Taken from 'gaussian_initial_estimates' function in the curve fitting notebook provided in the lab
    :param channels: array of channels
    :param count_rate: array of count_rate
    :return: tuple of gaussian estimates
    """
    mu0 = first_moment(channels, count_rate)
    sigma0 = np.sqrt(second_moment(channels, count_rate))
    area0 = np.sum(count_rate)

    return (mu0, sigma0, area0)


# calculate parameters for function fit to data
def fit_curve(fit_func, channels, count_rate, **kwargs):
    """
    Calculate the fit parameters and covariance matrix for a given fucntion using scipy.optimize.curve_fit
    :param fit_func: function that will be fitted to channel/count data
    :param channels: list of channels
    :param count_rate: list of count_rate
    :param kwargs: keyword arguments (should include initial guess parameters)
    :return: parameters and covariances
    """
    popt, pcov = curve_fit(fit_func, channels, count_rate, **kwargs)
    return popt, pcov


def line(x, c0, c1):
    """
    Defines a line
    :param x: list of x values
    :param c0: slope
    :param c1: y intercept
    :return: line
    """
    return c0 + c1 * x


# returns list of arrays of components of functions
def gaussian_plus_line_components(x, mu, sigma, area, c0, c1):
    """
    Return list of arrays of gaussian and line functions evaluated at channels
    - Taken from 'gaussian_plus_line_components' function in the curve fitting notebook provided in the lab
    :param x: array of channels
    :param mu: mean
    :param sigma: std dev
    :param area: area
    :param c0: slope
    :param c1: y intercept
    :return: list of arrays
    """
    components = [
        gaussian_func(x, mu, sigma, area),
        line(x, c0, c1),
    ]
    return components


def gaussian_plus_line(x, mu, sigma, area, c0, c1):
    """
    Returns a gaussian on a linear background
    - Taken from 'gaussian_plus_line' function in the curve fitting notebook provided in the lab
    :param x: array of channels
    :param mu: mean
    :param sigma: std dev
    :param area: area
    :param c0: y intercept
    :param c1: slope
    :return: Gaussian on linear background
    """
    _components = gaussian_plus_line_components(x, mu, sigma, area, c0, c1)
    return sum(_components)


def double_gaussian_components(x, mu0, sigma0, area0, mu1, sigma1, area1):
    """
    Returns lits of arrays of two gaussian functions evaluated at channels
    :param x: array of channels
    :param mu0: mean of first gaussian
    :param sigma0: std dev of first gaussian
    :param area0: area of first gaussian
    :param mu1: mean of second gaussian
    :param sigma1: std dev of second gaussian
    :param area1: area of second gaussian
    :return: list of arrays
    """
    components = [
        gaussian_func(x, mu0, sigma0, area0),
        gaussian_func(x, mu1, sigma1, area1)
    ]

    return components


def double_gaussian(x, mu0, sigma0, area0, mu1, sigma1, area1):
    """
    Returns a double gaussian function
    :param x: array of channels
    :param mu0: mean of first gaussian
    :param sigma0: std dev of first gaussian
    :param area0: area of first gaussian
    :param mu1: mean of second gaussian
    :param sigma1: std dev of second gaussian
    :param area1: area of second gaussian
    :return: double gaussian
    """
    _components = double_gaussian_components(x, mu0, sigma0, area0, mu1, sigma1, area1)
    return sum(_components)


def double_gaussian_plus_line_components(x, mu0, sigma0, area0, mu1, sigma1, area1, c0, c1):
    """
    Returns list of arrays of two gaussian functions and a line evaluated at channels
    :param x: array of channels
    :param mu0: mean of first gaussian
    :param sigma0: std dev of first gaussian
    :param area0: area of first gaussian
    :param mu1: mean of second gaussian
    :param sigma1: std dev of second gaussian
    :param area1: area of second gaussian
    :param c0: y-intercept
    :param c1: slope
    :return: list of arrays
    """
    components = [
        gaussian_func(x, mu0, sigma0, area0),
        gaussian_func(x, mu1, sigma1, area1),
        line(x, c0, c1)
    ]

    return components


def double_gaussian_plus_line(x, mu0, sigma0, area0, mu1, sigma1, area1, c0, c1):
    """
    Returns double gaussian on linear background
    :param x: array of channels
    :param mu0: mean of first gaussian
    :param sigma0: std dev of first gaussian
    :param area0: area of first gaussian
    :param mu1: mean of second gaussian
    :param sigma1: std dev of second gaussian
    :param area1: area of second gaussian
    :param c0: y-intercept
    :param c1: slope
    :return: Double gaussian on linear background
    """
    _components = double_gaussian_plus_line_components(x, mu0, sigma0, area0, mu1, sigma1, area1, c0, c1)
    return sum(_components)


# gaussian bounds and params
GAUSSIAN_PARAMS = ('mu', 'sigma', 'area')
GAUSSIAN_BOUNDS = [
        np.array([0, 0, 0]),
        np.array([INF, INF, INF]),
]

# data structure to hold singular gaussian model tools
gaussian_model = {
        "model function": gaussian_func,
        "estimates": initial_gaussian_estimates,
        "params": GAUSSIAN_PARAMS,
        "bounds": GAUSSIAN_BOUNDS,
        "components": gaussian_func,
}

# params and bounds of gaussian plus line
GAUSSIAN_PLUS_LINE_PARAMS = ('mu', 'sigma', 'area', 'c0', 'c1')
GAUSSIAN_PLUS_LINE_BOUNDS = [
    np.array([0, 0, 0, -1e9, -1e9]),
    np.array([INF, INF, INF, INF, INF])
    ]


# data structure to hold gaussian plus line compound model tools
gaussian_line_model = {
        "model function": gaussian_plus_line,
        "estimates": initial_gaussian_estimates,
        "params": GAUSSIAN_PLUS_LINE_PARAMS,
        "bounds": GAUSSIAN_PLUS_LINE_BOUNDS,
        "components": gaussian_plus_line_components,
}

# params and bounds of double gaussian
DOUBLE_GAUSS_PARAMS = ('mu0', 'sigma0', 'area0', 'mu1', 'sigma1', 'area1')
DOUBLE_GAUSS_BOUNDS = [
    np.array([0, 0, 0, 0, 0, 0]),
    np.array([INF, INF, INF, INF, INF, INF]),
    ]

# data structure to hold double gaussian compound model tools
double_gaussian_model = {
    "model function": double_gaussian,
    "estimates": initial_gaussian_estimates,
    "params": DOUBLE_GAUSS_PARAMS,
    "bounds": DOUBLE_GAUSS_BOUNDS,
    "components": double_gaussian_components,
}

# params and bounds of double gaussian and line
DOUBLE_GAUSS_PLUS_LINE_PARAMS = ('mu0', 'sigma0', 'area0', 'mu1', 'sigma1', 'area1', 'c0', 'c1')
DOUBLE_GAUSS_PLUS_LINE_BOUNDS = [
    np.array([0, 0, 0, 0, 0, 0, -1e9, -1e9]),
    np.array([INF, INF, INF, INF, INF, INF, INF, INF]),
    ]

# data structure to hold double gaussian and line compound model tools
double_gaussian_plus_line_model = {
    "model function": double_gaussian_plus_line,
    "estimates": initial_gaussian_estimates,
    "params": DOUBLE_GAUSS_PLUS_LINE_PARAMS,
    "bounds": DOUBLE_GAUSS_PLUS_LINE_BOUNDS,
    "components": double_gaussian_plus_line_components,
}


def format_initial_params(params, guesses):
    """
    Display initial guesses
    - Adapted from 'format_result' function in the curve fitting notebook provided in the lab
    :param params: parameter names
    :param guesses: initial guesses
    :return: lines to be printed
    """
    _lines = (f"{p} = {o}" for p, o in zip(params, guesses))
    return "\n".join(_lines)


def format_fit_params(params, popt, pcov):
    """
    Display final fit parameters
    - Adapted from 'format_result' function in the curve fitting notebook provided in the lab
    :param params: parameter names
    :param popt: final parameters
    :param pcov: covariances
    :return: lines to be printed
    """
    uncertainties = np.sqrt(np.diag(pcov))
    _lines = (f"{p} = {o} +/- {e}" for p, o, e in zip(params, popt, uncertainties))
    return "\n".join(_lines)


def check_model_and_give_estimates(model, _channels, _count_rate, added_function_params):
    """
    Depending on the model, the parameters fed to curve_fit will be different
    :param model: model data structure
    :param _channels: selected channels
    :param _count_rate: selected count rate
    :param added_function_params: parameters of any additional functions other than gaussians
    :return: parameters to be fed to curve_fit
    """
    additional_params = tuple(added_function_params)
    # if the model involves a double gaussian, split the selected channels and count rates in half
    if model == double_gaussian_model or model == double_gaussian_plus_line_model:
        _p00 = model["estimates"](_channels[:len(_channels)//2], _count_rate[:len(_channels)//2])
        _p01 = model["estimates"](_channels[len(_channels)//2:], _count_rate[len(_channels)//2:])
        # if the model also contains a line, add the line parameters
        if model == double_gaussian_plus_line_model:
            _p0 = _p00 + _p01 + additional_params
        else:
            _p0 = _p00 + _p01
    elif model == gaussian_line_model:
        _p0 = model["estimates"](_channels, _count_rate) + additional_params
    else:
        _p0 = model["estimates"](_channels, _count_rate)

    return _p0


def model_fit(model, source, channels, count_rate, roi, added_function_params, **kwargs):
    """
    API that takes model structure and returns estimates and covariance matrix
    - Adapted from suggestions in the curve fitting notebook provided in the lab
    :param source: name of the source
    :param model: the model data structure
    :param channels: array of channels
    :param count_rate: array of count_rate
    :param roi: region of interest
    :param added_function_params: parameters of additional functions other than gaussians
    :param kwargs: keyword arguments
    :return: estimates and covariances
    """
    # select channels and count_rate
    _channels, _count_rate = filter_with_mask(channels, count_rate, *roi)
    _p0 = check_model_and_give_estimates(model, _channels, _count_rate, added_function_params)

    print(f"#-- The initial estimates for {source} --#\n")
    print(f"{format_initial_params(model['params'], _p0)}\n")

    # do fit
    popt, pcov = fit_curve(model["model function"], _channels, _count_rate, p0=_p0, bounds=model["bounds"],
                           maxfev=80000, **kwargs)

    print(f"#-- The fitted estimates for {source}--#\n")
    print(f"{format_fit_params(model['params'], popt, pcov)}\n")

    return popt, pcov


def plot_model(ax, model_function, x_range, model_params, num_pts=5000, **kwargs):
    """
    Plots the fitted curve from the model
    - Taken from 'plot_model' function in the curve fitting notebook provided in the lab
    :param ax: axes object
    :param model_function: the model function
    :param x_range: range of values within which to plot the model
    :param model_params: calculated parameters of the model
    :param num_pts: number of sampling points
    :param kwargs: keyword arguments
    :return: plot of the model
    """
    _channels = np.linspace(*x_range, num_pts)
    _count_rate = model_function(_channels, *model_params)
    return ax.plot(_channels, _count_rate, **kwargs)


def plot_model_components(ax, components, x_range, model_params, labels, num_pts=5000, **kwargs):
    """
    Plot the individual components of the model
    :param ax: axes object
    :param components: components of the model
    :param x_range: range of values within which to plot the model
    :param model_params: calculated parameters of the model
    :param labels: labels for the individual components
    :param num_pts: number of sampling points
    :param kwargs: keyword arguments
    :return: list of calls to plot
    """
    _channels = np.linspace(*x_range, num_pts)
    _count_rate_array_list = components(_channels, *model_params)
    plots = []
    colors = ['r', 'orange', 'green']
    for index, count_rate_arr in enumerate(_count_rate_array_list):
        plots.append(ax.plot(_channels, count_rate_arr, '--', c=colors[index], label=labels[index], **kwargs))
    return plots


def get_model(model):
    """
    Choose the model relevant plot labels
    :param model: name of the model to be used
    :return: model name, model data structure and label list
    """
    if model == "Gaussian":
        chosen_model = gaussian_model
        labels = None
    elif model == "Gaussian and Line":
        chosen_model = gaussian_line_model
        labels = ["Gaussian Fit", "Line Fit"]
    elif model == "Double Gaussian":
        chosen_model = double_gaussian_model
        labels = ["Gaussian Fit 1", "Gaussian Fit 2"]
    elif model == "Double Gaussian and Line":
        chosen_model = double_gaussian_plus_line_model
        labels = ["Gaussian Fit 1", "Gaussian Fit 2", "Line Fit"]
    else:
        raise Exception("The model is invalid!")

    return model, chosen_model, labels


def count_rate_uncertainty(counts, bkg_counts, real_time, bkg_rt, roi):
    """
    Calculate uncertainty in count rate
    :param counts: count array
    :param bkg_counts: background count array
    :param real_time: observation time of detector
    :param bkg_rt: background time
    :param roi: region of interest
    :return: uncertainty calculated in quadrature
    """
    # select counts and background counts within region of interest
    # set any zero values to 1
    counts = counts[roi[0]:roi[1]]
    counts[counts <= 0] = 1
    bkg_counts = bkg_counts[roi[0]:roi[1]]
    bkg_counts[bkg_counts <= 0] = 1
    delta_count_rate = np.sqrt((np.sqrt(counts)/real_time)**2 + (np.sqrt(bkg_counts)/real_time)**2)
    return delta_count_rate


def calculate_and_plot_model(source, detector, spectrum, background, roi, zoomed_params, model, slice_index,
                             added_function_params, **kwargs):
    """
    Put everthing together - run model API and plot model
    :param source: name of the source
    :param detector: name of the detector
    :param spectrum: spectrum data read from the spectrum file
    :param background: background counts
    :param roi: region of interest
    :param zoomed_params: indices within which to plot a closer figure
    :param model: name of the model from the config file
    :param slice_index: index up to which to plot the spectrum
    :param added_function_params: parameters of additional functions other than gaussians
    :param kwargs: keyword arguments
    :return: None
    """
    # counts from spectrum
    counts = np.array(spectrum[0])
    # background counts
    bkg = np.array(background[0])
    # real time from spectrum data
    real_time = spectrum[1][2]
    bkg_real_time = background[1][0]
    # subtract the background
    counts_bkg_subtraction = counts/real_time - bkg/bkg_real_time
    # if any counts rates end up less than zero, set them to very small
    counts_bkg_subtraction[counts_bkg_subtraction < 0] = 1 / real_time
    # count rate
    count_rate = counts_bkg_subtraction

    # calculate the uncertainty
    count_rate_err = count_rate_uncertainty(counts, bkg, real_time, bkg_real_time, roi)
    # generate an array of channels
    channels = np.array(list(range(len(count_rate))))

    # plot the spectrum
    fig, ax = plt.subplots()
    plot_spectrum(ax, channels, count_rate, slice_index, source, detector, roi=roi, **kwargs)
    # plot a zoomed in spectrum
    zoom_fig, zoom_ax = plt.subplots()
    zoom_in(zoom_ax, channels, count_rate, zoomed_params, source, detector, roi=roi, **kwargs)

    # select the model and get the necessary labels
    model_name, chosen_model, labels = get_model(model)
    # fit the models including the uncertainty
    print(f"DETECTOR: {detector}")
    model_params, covariances = model_fit(chosen_model, source, channels, count_rate, roi, added_function_params,
                                          sigma=count_rate_err, absolute_sigma=True)

    # overplot the models on the spectrum
    plot_model(ax, chosen_model["model function"], (0, slice_index), model_params, c='black', label=f"{model_name} fit")
    # overplot the models on the zoomed spectrum
    plot_model(zoom_ax, chosen_model["model function"], zoomed_params, model_params, c='black',
               label=f"{model_name} fit")

    # if the model contains more than just a gaussian, overplot the individual functions
    if chosen_model != gaussian_model:
        plot_model_components(ax, chosen_model["components"], (0, slice_index), model_params, labels)
        plot_model_components(zoom_ax, chosen_model["components"], zoomed_params, model_params, labels)

    fig.legend()
    zoom_fig.legend()

    return model_params, covariances, detector, model_name, count_rate


def params_and_ucerts(source, detector, spectrum, background, roi, zoomed_params, model, slice_index,
                           added_function_params, **kwargs):
    """
    Return centroids and sigmas witht their uncertainties
    :param source: source name
    :param detector: detector name
    :param spectrum: spectrum data
    :param background: background data
    :param roi: region of interest
    :param zoomed_params: indices within which to plot a closer figure
    :param model: name of the model from the config file
    :param slice_index: index up to which to plot the spectrum
    :param added_function_params: initial guesses for line models
    :param kwargs: parameters of additional functions other than gaussians
    :return: centroids, sigmas and their uncertainties
    """
    # calculate parameters and plot model and components
    model_params, covariances, detector, model_name, count_rate = calculate_and_plot_model(source, detector, spectrum,
                                background, roi, zoomed_params, model, slice_index, added_function_params, **kwargs)
    # calculate the uncertainties
    uncertainties = np.sqrt(np.diag(covariances))
    # return the desired model parameters depending on the type of model
    if model_name == "Double Gaussian" or model_name == "Double Gaussian and Line":
        centroids = [model_params[0], model_params[3]]
        sigmas = [model_params[1], model_params[4]]
        peak_areas = [model_params[2], model_params[5]]
        centroid_ucerts = [uncertainties[0], uncertainties[3]]
        sigma_ucerts = [uncertainties[1], uncertainties[4]]
        peak_area_ucerts = [uncertainties[2], uncertainties[5]]
    else:
        centroids = [model_params[0]]
        sigmas = [model_params[1]]
        peak_areas = [model_params[2]]
        centroid_ucerts = [uncertainties[0]]
        sigma_ucerts = [uncertainties[1]]
        peak_area_ucerts = [uncertainties[2]]

    return centroids, centroid_ucerts, sigmas, sigma_ucerts, peak_areas, peak_area_ucerts, count_rate

