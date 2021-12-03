
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from Read_files import extract_spectrum_and_model_input, get_background
from Energy_Calibration_Model_fitting import params_and_ucerts


def collect_model_data(spec_file_path, bkg_file, config_file_path, detector, source, **kwargs):
    """
    Run models and collect relevant data
    :param spec_file_path: path to spectrum file
    :param bkg_file: background file name
    :param config_file_path: path to config file
    :param detector: detector name
    :param source: source name
    :param kwargs: keyword arguments
    :return: relevant model data
    """
    spectrum_and_metadata, ROI, zoomed_params, model, slice_index, added_function_params =\
        extract_spectrum_and_model_input(spec_file_path, config_file_path, source)
    bkg_spectrum_and_meta = get_background(spec_file_path, bkg_file)

    centroids, centroid_uncertainties, sigmas, sigma_uncertainties, peak_areas, peak_area_ucerts, count_rate = \
        params_and_ucerts(source, detector, spectrum_and_metadata, bkg_spectrum_and_meta, ROI, zoomed_params, model,
                          slice_index, added_function_params, **kwargs)
    return centroids, centroid_uncertainties, sigmas, sigma_uncertainties, peak_areas, peak_area_ucerts, count_rate


def choose_sources_and_energies(detector):
    """
    Choose the sources given the detector, select the true enrgies and their uncertainties
    :param detector: name of the detector
    :return: sources, energies and energy uncertainties
    """
    sources = ["AMERICIUM", "BARIUM", "CAESIUM", "COBALT"]
    # separate the sources so peaks can be fitted separately, HPGe has too high a resolution for double peak fitting
    sources_HPGE = ["AMERICIUM_1", "AMERICIUM_2", "BARIUM_1", "BARIUM_2", "CAESIUM", "COBALT_1", "COBALT_2"]
    energies = [26.3446, 59.5409, 302.8508, 356.0129, 661.657, 1173.228, 1332.501]
    energy_ucerts = [0.0002, 0.0001, 0.0005, 0.0007, 0.003, 0.003, 0.005]
    if detector == "CdTe":
        sources = sources[:2]
    elif detector == "NaI(Tl)":
        sources = sources[:3]
        energies = energies[:5]
        energy_ucerts = energy_ucerts[:5]
    elif detector == "BGO":
        sources = sources[:]
        energies = energies[1:]
        energy_ucerts = energy_ucerts[1:]
    elif detector == "HPGe":
        sources = sources_HPGE
        energies = energies
        energy_ucerts = energy_ucerts
    return sources, energies, energy_ucerts


def line(x, a, b):
    """Return a line"""
    return a * x + b


def channel_energy_arrays(detector, detector_file_path, bkg_file, config_file):
    """
    Collect data for each source and populate arrays of relevant data to pass on down the pipeline
    :param detector: detector name
    :param detector_file_path: path to detector folder
    :param bkg_file: background file
    :param config_file: config file name
    :return: centroid, energy, sigma arrays and uncertainties
    """
    centroid_list = []
    peak_uncertainties_list = []
    sigma_list = []
    sigma_uncertainties_list = []
    peak_area_list = []
    peak_area_ucerts_list = []
    sum_count_rates = []

    # choose the energies and sources for the detector
    sources, energies, energy_ucerts = choose_sources_and_energies(detector)

    # for each source, collect the model data
    for source in sources:
        centroids, peak_uncertainties, sigmas, sigma_uncertainties, peak_areas, peak_area_ucerts, count_rate = \
            collect_model_data(detector_file_path, bkg_file,
                                                    config_file, detector, source, marker='+', s=10)
        centroid_list += centroids
        peak_uncertainties_list += peak_uncertainties
        sigma_list += sigmas
        sigma_uncertainties_list += sigma_uncertainties
        peak_area_list += peak_areas
        peak_area_ucerts_list += peak_area_ucerts
        sum_count_rates.append(np.sum(count_rate))

    centroid_arr = np.array(centroid_list)
    peak_uncertainties_arr = np.array(peak_uncertainties_list)
    sigma_arr = np.array(sigma_list)
    sigma_uncertainties_arr = np.array(sigma_uncertainties_list)
    peak_area_arr = np.array(peak_area_list)
    peak_area_ucerts_arr = np.array([peak_area_ucerts_list])
    energies_arr = np.array(energies)
    energy_ucerts_arr = np.array(energy_ucerts)
    sum_count_rate_arr = np.array(sum_count_rates)

    return centroid_arr, peak_uncertainties_arr, energies_arr, energy_ucerts_arr, sigma_arr, sigma_uncertainties_arr,\
           peak_area_arr, peak_area_ucerts_arr, sum_count_rate_arr


def fit_line(centroids, energies, **kwargs):
    """
    Determine the line fit parameters for the calibration
    :param centroids: peak centorids from model fitting
    :param energies: true energies
    :param kwargs: kwargs
    :return: fit params and uncertainties
    """
    # fit centroids as y variable to allow uncertainties to be included
    line_popt, line_pcov = curve_fit(line, energies, centroids, **kwargs)
    # param uncertainties
    line_ucerts = np.sqrt(np.diag(line_pcov))

    # flip the parameters so that energy can be plotted on the y axis
    channel_c1 = 1 / line_popt[0]
    # uncertainty in the first parameter
    channel_c1_ucert = channel_c1 * (line_ucerts[0] / line_popt[0])
    channel_c2 = -line_popt[1] / line_popt[0]
    # uncertainty in the second parameter
    channel_c2_ucert = abs(channel_c2) * uncertainty_mult_or_div(line_popt[0], line_popt[1], line_ucerts[0],
                                                                 line_ucerts[1])

    line_consts = (channel_c1, channel_c2)
    line_const_ucerts = (channel_c1_ucert, channel_c2_ucert)

    return line_consts, line_const_ucerts


def linear_ucert(x, x_ucert, a, b, a_ucert, b_ucert):
    """Uncertainty in using the line fit"""
    term = x**2 * a_ucert**2 + b_ucert**2 + a**2 * x_ucert**2
    return np.sqrt(term)


def channel_energy_fit(detector, detector_file_path, bkg_file, config_file, e_range, **kwargs):
    """
    Calculate the fit and plot the channel energy relationship
    :param detector: detector name
    :param detector_file_path: path to detector folder
    :param bkg_file: background file
    :param config_file: config file name
    :param e_range: array of energies over the range
    :param kwargs: keyword arguments
    :return: array of calculated energies, quadratic params and covariances, plots, sigmas and sigma uncertainties
    """

    centroid_arr, peak_uncertainties_arr, energies_arr, energy_ucerts_arr, sigma_arr,\
    sigma_uncertainties_arr, peak_area_arr, peak_area_ucerts_arr, sum_count_rates_arr = channel_energy_arrays(detector,
                                                                            detector_file_path, bkg_file, config_file)

    e_cal_fit, e_cal_ax = plt.subplots()
    e_cal_ax.set_title(f"Energy Calibration Plot for {detector}")
    e_cal_ax.set_xlabel("Peak Centroids")
    e_cal_ax.set_ylabel("Peak Energies (keV)")

    e_cal_plot = e_cal_ax.errorbar(centroid_arr, energies_arr, xerr=peak_uncertainties_arr, yerr=energy_ucerts_arr,
                                   label="Peak Centroid Data", fmt='.', c='r', **kwargs)

    absolute_sigma = True
    sigma = peak_uncertainties_arr

    # fit the line
    line_consts, line_const_ucerts = fit_line(centroid_arr, energies_arr, absolute_sigma=absolute_sigma, sigma=sigma,
                                              **kwargs)

    # calculate the model over the range
    line_model_e_arr = line(e_range, *line_consts)
    # calculate the model at the specific peaks
    line_fit = line(centroid_arr, *line_consts)
    # calculate uncertainty of energies calculated from fit function
    line_fit_ucerts = linear_ucert(centroid_arr, peak_uncertainties_arr, *line_consts, *line_const_ucerts)
    # label for the fit
    fit_label = f"{line_consts[0]:.5} $\\pm$ {line_const_ucerts[0]:.3} $\\times$ channel + {line_consts[1]:.5}" \
                f" $\\pm$ {line_const_ucerts[1]:.3}"

    plot_e_cal_fit = e_cal_ax.plot(e_range, line_model_e_arr, label=fit_label)
    e_cal_fit.legend()

    return line_fit, line_fit_ucerts, line_consts, line_const_ucerts, centroid_arr, peak_uncertainties_arr, sigma_arr,\
           sigma_uncertainties_arr, peak_area_arr, peak_area_ucerts_arr, sum_count_rates_arr


# energy resolution model as stated in lab guide
def resolution_model_func(energies, a, b, c):
    """
    Model to fit the resolution vs energy data
    :param energies: array of energies
    :param a: const
    :param b: const
    :param c: const
    :return: resolution model array
    """
    _result = (a / energies**2) + (b / energies) + c
    return np.sqrt(_result)


def fit_resolution_model(energies, resolutions, **kwargs):
    """
    Calculate resolution model parameters
    :param energies: energies of peaks
    :param resolutions: resolutions of peaks
    :return: parameters and covariances
    """
    parameter_guesses = (10000, 10000, 100)
    res_popt, res_pcov = curve_fit(resolution_model_func, energies, resolutions, p0=parameter_guesses, **kwargs)
    return res_popt, res_pcov


def uncertainty_mult_or_div(x, y, deltax, deltay):
    """
    Calculates the uncertainity in multiplying or dividing two quantities
    :param x: any quantity
    :param y: any quantity
    :param deltax: x uncertainty
    :param deltay: y uncertainty
    :return: uncertainity in x*y or x/y
    """
    term = (deltax / x)**2 + (deltay / y)**2
    return np.sqrt(term)


def energy_resolution(detector, e_range, line_fit, line_fit_ucerts, line_consts, line_const_ucerts,
                      centroids, centroid_ucerts, sigma_arr,
                      sigma_uncertainties_arr, **kwargs):
    """
    Calculate the energy resolution
    :param detector: detector name
    :param e_range: array of energies over the range
    :param line_fit: array of calculated energies from channel energy relationship
    :param line_fit_ucerts: uncertainty in line fit
    :param line_consts: line params
    :param line_const_ucerts: uncertainties in the line consts
    :param centroids: peak centroids
    :param centroid_ucerts: centroid uncertainties
    :param sigma_arr: array of sigmas from model
    :param sigma_uncertainties_arr: array of sigma uncertainties
    :return: energy resolution plot
    """

    # calculate the FWHM in terms of channels
    fwhm_channel_arr = 2.355 * sigma_arr
    # FWHM in keV
    fwhm_energies_ucerts = linear_ucert(fwhm_channel_arr, 2.355 * sigma_uncertainties_arr, *line_consts,
                                           *line_const_ucerts)

    # energies from line model
    energies = line_fit
    energy_ucerts = line_fit_ucerts
    # resolutions as percentage
    resolution = (fwhm_channel_arr / centroids) * 100
    # resolution uncertainties
    res_ucert = resolution * uncertainty_mult_or_div(centroids, fwhm_channel_arr, centroid_ucerts,
                                                     2.355 * sigma_uncertainties_arr)

    res_fig, res_ax = plt.subplots()
    res_ax.set_title(f"Energy Resolution for {detector}")
    res_ax.set_xlabel("Energy (KeV)")
    res_ax.set_ylabel("Resolution (%)")

    res_ax.set_xscale('log')
    res_ax.set_yscale('log')

    # remove Barium lower energy peak from resolution
    if detector == "NaI(Tl)":
        energies = np.delete(energies, 2)
        resolution = np.delete(resolution, 2)
        energy_ucerts = np.delete(energy_ucerts, 2)
        res_ucert = np.delete(res_ucert, 2)

    # removes the outlying americium point for BGO
    if detector == "BGO":
        energies = energies[1:]
        resolution = resolution[1:]
        energy_ucerts = energy_ucerts[1:]
        res_ucert = res_ucert[1:]
        e_range = np.logspace(2.4, 3.15, 5000)

    if detector == "HPGe":
        energies = energies[1:]
        resolution = resolution[1:]
        energy_ucerts = energy_ucerts[1:]
        res_ucert = res_ucert[1:]
        e_range = np.logspace(1.7, 3.17, 5000)

    res_ax.errorbar(energies, resolution, xerr=energy_ucerts, yerr=res_ucert,
                    label="Resolution at Peak", fmt='.', c='r', **kwargs)

    absolute_sigma = True
    sigma = res_ucert

    # The HPGe resolution fit does not work well at all despite showing a straight line trend in log space
    if detector == "HPGe":
        absolute_sigma = False
        sigma = None

    # calculate fit params
    res_popt, res_pcov = fit_resolution_model(energies, resolution, absolute_sigma=absolute_sigma, sigma=sigma,
                                              maxfev=10000)

    # fit uncertainties
    res_fit_ucerts = np.sqrt(np.diag(res_pcov))
    # calculate model over range
    res_model_arr = resolution_model_func(e_range, *res_popt)

    fit_label = f"({res_popt[0]:.3} $\\pm$ {res_fit_ucerts[0]:.3} $\\times$ E$^{{{-2}}}$ + {res_popt[1]:.3}" \
                f" $\\pm$ {res_fit_ucerts[1]:.3} $\\times$ E$^{{{-1}}}$" \
                f" + {res_popt[2]:.3} $\\pm$ {res_fit_ucerts[2]:.3})$^{{{0.5}}}$"
    plot_res_model = res_ax.plot(e_range, res_model_arr, label=fit_label, **kwargs)

    res_fig.legend()

    return


def calibration_resolution_analysis(detector, data_file_path, bkg_file, config_file, e_range, e_range_log):
    """Call the calibration and resolution routines"""
    line_fit, line_fit_ucerts, line_consts, line_const_ucerts, centroid_arr, peak_uncertainties_arr, sigma_arr,\
    sigma_uncertainties_arr, peak_area_arr, peak_area_ucerts_arr, sum_count_rates_arr = \
        channel_energy_fit(detector, data_file_path, bkg_file, config_file, e_range)
    energy_resolution(detector, e_range_log, line_fit, line_fit_ucerts, line_consts, line_const_ucerts, centroid_arr,
                      peak_uncertainties_arr, sigma_arr, sigma_uncertainties_arr)

    return line_fit
