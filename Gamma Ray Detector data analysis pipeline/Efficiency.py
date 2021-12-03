
import numpy as np
from datetime import datetime
from Read_files import json_energy_calibration_reader, file_reader, get_background
from Energy_Calibration_Model_fitting import params_and_ucerts
from Energy_Calibration_and_Resolution import uncertainty_mult_or_div
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def remove_bad_fits(detector, source, angle):
    """Skip files where peak fitting models fail -  this only happens for the two barium files for NaI(Tl)"""
    skip_file = False
    if detector == "NaI(Tl)":
        if source == "Barium":
            if angle == 100 or angle == 120:
                skip_file = True

    return skip_file


def read_all_angle_files(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list, **kwargs):
    """
    Read every file for every source for the detector and create a dictionary of date, time and peak areas +
     uncertainties
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    :param kwargs: kwargs
    :return: dictionary {source: {angle:[[date, time], peak areas, uncertainties]}, ....}
    """
    source_dict = {}
    # get the background for the detector
    bkg_spectrum_and_meta = get_background(spectrum_file_path, bkg_file)

    # the HPGe config file is layed out differently so unfortunately there is a separate mechanism that must be used to
    # obtain a dictionary at the end witht the same format as for the other detectors
    # list of 'source' strings for HPGe
    json_hpge_sources = ["AMERICIUM_1", "AMERICIUM_2", "BARIUM_1", "BARIUM_2", "CAESIUM", "COBALT_1", "COBALT_2"]

    # for each source in the sources used with the detector
    for source in sources:
        source_dict[source] = {}

        # HPGe collection of config data - spectra with multiple peaks were not fitted together
        # therefore this must be done
        if detector == "HPGe":
            if source == "Americium":
                _, ROI, zoomed_params, model, slice_index, added_function_params \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[0])
                _, ROI1, zoomed_params1, model1, slice_index1, added_function_params1 \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[1])
            elif source == "Barium":
                _, ROI, zoomed_params, model, slice_index, added_function_params \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[2])
                _, ROI1, zoomed_params1, model1, slice_index1, added_function_params1 \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[3])
            elif source == "Caesium":
                _, ROI, zoomed_params, model, slice_index, added_function_params \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[4])
            else:
                _, ROI, zoomed_params, model, slice_index, added_function_params \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[5])
                _, ROI1, zoomed_params1, model1, slice_index1, added_function_params1 \
                    = json_energy_calibration_reader(config_file_path, json_hpge_sources[6])

        else:
            json_source = source.upper()
            # read the config file for the current source
            _, ROI, zoomed_params, model, slice_index, added_function_params \
                = json_energy_calibration_reader(config_file_path, json_source)

        # for each file in the list of possible files
        for file_name in file_list:

            source_and_file = "/" + source + file_name

            # parse the file otherwise skip iteration
            try:
                spectrum_and_metadata = file_reader(spectrum_file_path, source_and_file)
            except FileNotFoundError:
                continue

            angle, date, time = spectrum_and_metadata[1][1], spectrum_and_metadata[1][3], spectrum_and_metadata[1][4]

            # skip poorly fitted spectra
            if remove_bad_fits(detector, source, angle):
                continue

            # fit the peak models to the data
            _, _, _, _, peak_areas, peak_area_ucerts, _ = \
                params_and_ucerts(source, detector, spectrum_and_metadata, bkg_spectrum_and_meta, ROI, zoomed_params,
                                  model, slice_index, added_function_params, **kwargs)

            # for the HPGe detector, each source except Caesium has to be run through again to fit the second peak
            if detector == "HPGe" and source != "Caesium":
                _, _, _, _, peak_areas1, peak_area_ucerts1, _ = \
                    params_and_ucerts(source, detector, spectrum_and_metadata, bkg_spectrum_and_meta, ROI1,
                                      zoomed_params1,
                                      model1, slice_index1, added_function_params1, **kwargs)
                peak_areas = peak_areas + peak_areas1
                peak_area_ucerts = peak_area_ucerts + peak_area_ucerts1

            # populate the dictionary
            source_dict[source][angle] = [[date, time], peak_areas, peak_area_ucerts]
    plt.close('all')

    return source_dict


def activity_func(A_0, decay_const, time):
    """Calculate the activity"""
    return A_0 * np.exp(-decay_const * time)


def decay_const_ucert(half_lives, hl_ucerts):
    """
    Find the uncertainty in calculating the decay constant
    :param half_lives: half life array
    :param hl_ucerts: half life uncertainties array
    :return:
    """
    term = (np.log(2) * half_lives**-2)**2 * hl_ucerts**2
    return np.sqrt(term)


def extend_arrays(array):
    """Extend arrays - for use when quantities need to be used twice for the same source"""
    array = np.insert(array, 1, array[0])
    array = np.insert(array, 2, array[2])
    array = np.insert(array, -1, array[-1])
    return array


def calculate_decay_const():
    """
    Calculate the decay constants of the sources
    :return: decay constants and their uncertainties
    """
    half_lives = np.array([432.20, 10.51, 30.07, 5.2714])
    half_lives = half_lives * 3.154e7
    hl_ucerts = np.array([0.7, 0.05, 0.03, 0.0005])
    hl_ucerts = hl_ucerts * 3.154e7
    dc_ucerts = decay_const_ucert(half_lives, hl_ucerts)

    log2 = np.log(2)
    decay_consts = log2 / half_lives

    # extend the decay constant array
    decay_consts = extend_arrays(decay_consts)
    # extend the uncertainties array
    dc_ucerts = extend_arrays(dc_ucerts)

    return decay_consts, dc_ucerts


def initial_activity_in_bq():
    """
    Calculate the initial activities and their uncertainites in Bq (upstairs and downstairs sources)
    :return: downstairs activities, downstairs activity uncertainties, upstairs activities, upstairs activity
     uncertainties
    """

    activities_in_mc_down = np.array([11.16, 11.42, 11.16, 12.28])  # micro curies
    activities_in_mc_up = np.array([11.92, 10.82, 12.41, 11.32])

    accuracies = np.array([0.05, 0.048, 0.037, 0.019]) # accuracies from lab manual

    activities_in_bq_down = 37000 * activities_in_mc_down # becquerels
    activities_in_bq_up = 37000 * activities_in_mc_up

    # uncertainties
    activity_ucerts_bq_down = accuracies * activities_in_bq_down
    activity_ucerts_bq_up = accuracies * activities_in_bq_up

    # extend arrays
    activities_in_bq_down = extend_arrays(activities_in_bq_down)
    activities_in_bq_up = extend_arrays(activities_in_bq_up)

    activity_ucerts_bq_down = extend_arrays(activity_ucerts_bq_down)
    activity_ucerts_bq_up = extend_arrays(activity_ucerts_bq_up)

    return activities_in_bq_down, activity_ucerts_bq_down, activities_in_bq_up, activity_ucerts_bq_up


def quantities_by_detector(detector):
    """Select the activities and decay constant and their uncertainties for each detetcor - quite ugly"""

    decay_fracs = np.array([0.024, 0.359, 0.1833, 0.6205, 0.851, 0.999, 0.9998]) # decay fractions

    # decay constants
    decay_consts, dc_ucerts = calculate_decay_const()
    if detector == "HPGe":
        # initial activities
        _, _, activities, activity_ucerts = initial_activity_in_bq()
        # multiply by decay fractions
        activities = activities * decay_fracs
    elif detector == "NaI(Tl)":
        # initial activities
        activities, activity_ucerts, _, _ = initial_activity_in_bq()
        # choose activities and decay constants
        activities = activities[:-2]
        activity_ucerts = activity_ucerts[:-2]
        decay_consts = decay_consts[:-2]
        dc_ucerts = dc_ucerts[:-2]
        decay_fracs = decay_fracs[:-2]
        # multiply by decay fractions
        activities = activities * decay_fracs
    elif detector == "BGO":
        # initial activities
        activities, activity_ucerts, _, _ = initial_activity_in_bq()
        # choose activities and decay constants
        activities = activities[1:]
        activity_ucerts = activity_ucerts[1:]
        decay_consts = decay_consts[1:]
        dc_ucerts = dc_ucerts[1:]
        decay_fracs = decay_fracs[1:]
        # multiply by decay fractions
        activities = activities * decay_fracs

    return activities, activity_ucerts, decay_consts, dc_ucerts



def split_date_string(date):
    """Return the month, day and year from the file metadata string"""
    month = date.split("/", 1)[0]
    day = (date.split("/", 1)[1]).split("/", 1)[0]
    year = (date.split("/", 1)[1]).split("/", 1)[1]
    return int(month), int(day), int(year)


def split_time_string(time):
    """Return the hour minute and second from the file metadata string"""
    hour = time.split(":", 1)[0]
    minute = (time.split(":", 1)[1]).split(":", 1)[0]
    sec = (time.split(":", 1)[1]).split(":", 1)[1]
    return int(hour), int(minute), int(sec)


def calculate_time_since_source_calibration(date, time, down=True):
    """Calculate the time in seconds from when sources were calibrated depending on whether they are from upstairs or
     downstairs"""
    calibration_date_down = "12/1/1979"
    calibration_date_up = "2/1/1979"
    calibration_time = "12:00:00"

    # split calibration time string
    cal_hour, cal_minute, cal_sec = split_time_string(calibration_time)

    # split strings from files
    month, day, year = split_date_string(date)
    hour, minute, sec = split_time_string(time)

    if not down:
        # for upstairs sources
        cal_month, cal_day, cal_year = split_date_string(calibration_date_up)
        start = datetime(year=cal_year, month=cal_month, day=cal_day, hour=cal_hour, minute=cal_minute, second=cal_sec)

    else:
        # for downstairs sources
        cal_month, cal_day, cal_year = split_date_string(calibration_date_down)
        start = datetime(year=cal_year, month=cal_month, day=cal_day, hour=cal_hour, minute=cal_minute, second=cal_sec)

    end = datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=sec)

    diff = end - start
    time_in_secs = diff.total_seconds()

    return time_in_secs


def choose_activities_and_decay_consts(detector, source):
    """Choose the activities and decay constants + uncertainties for each source for each detector - ugly"""

    activities, activity_ucerts, decay_consts, dc_ucerts = quantities_by_detector(detector)
    if detector == "NaI(Tl)":
        if source == "Americium":
            activities = activities[:2]
            activity_ucerts = activity_ucerts[:2]
            decay_consts = decay_consts[:2]
            dc_ucerts = dc_ucerts[:2]
        elif source == "Barium":
            activities = activities[2:4]
            activity_ucerts = activity_ucerts[2:4]
            decay_consts = decay_consts[2:4]
            dc_ucerts = dc_ucerts[2:4]
        else:
            activities = activities[-1]
            activity_ucerts = activity_ucerts[-1]
            decay_consts = decay_consts[-1]
            dc_ucerts = dc_ucerts[-1]
    elif detector == "BGO":
        if source == "Americium":
            activities = activities[0]
            activity_ucerts = activity_ucerts[0]
            decay_consts = decay_consts[0]
            dc_ucerts = dc_ucerts[0]
        elif source == "Barium":
            activities = activities[1:3]
            activity_ucerts = activity_ucerts[1:3]
            decay_consts = decay_consts[1:3]
            dc_ucerts = dc_ucerts[1:3]
        elif source == "Caesium":
            activities = activities[3]
            activity_ucerts = activity_ucerts[3]
            decay_consts = decay_consts[3]
            dc_ucerts = dc_ucerts[3]
        else:
            activities = activities[4:]
            activity_ucerts = activity_ucerts[4:]
            decay_consts = decay_consts[4:]
            dc_ucerts = dc_ucerts[4:]
    else:
        if source == "Americium":
            activities = activities[:2]
            activity_ucerts = activity_ucerts[:2]
            decay_consts = decay_consts[:2]
            dc_ucerts = dc_ucerts[:2]
        elif source == "Barium":
            activities = activities[2:4]
            activity_ucerts = activity_ucerts[2:4]
            decay_consts = decay_consts[2:4]
            dc_ucerts = dc_ucerts[2:4]
        elif source == "Caesium":
            activities = activities[4]
            activity_ucerts = activity_ucerts[4]
            decay_consts = decay_consts[4]
            dc_ucerts = dc_ucerts[4]
        else:
            activities = activities[5:]
            activity_ucerts = activity_ucerts[5:]
            decay_consts = decay_consts[5:]
            dc_ucerts = dc_ucerts[5:]

    return activities, decay_consts, activity_ucerts, dc_ucerts


def current_activity_ucerts(initial_ac, ac_ucert, time, delta_t, decay_const, dc_ucert):
    """
    Calculate the uncertainty in the activity calculation -  this should be done using a monte carlo simulation but
    time constaints left no opportunity to perform that -  this just uses simpler error propagation rules
    :param initial_ac: initial activitiy
    :param ac_ucert: activity uncertainty
    :param time: time since calibration
    :param delta_t: uncertainty in time
    :param decay_const: decay constant
    :param dc_ucert: uncertainty in decay constant
    :return:
    """
    term = np.exp(-2 * decay_const * time) * (ac_ucert**2 + (initial_ac * time * dc_ucert)**2 + (initial_ac *
                                decay_const * delta_t)**2)
    return np.sqrt(term)


def calculate_abs_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list):
    """
    Create a dictionary of absolute efficiencies for all data for the detector
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    :return: dictionary {source:{angle:[absolute efficiency, uncertainties]}}
    """

    # generate the data dictionary
    source_dict = read_all_angle_files(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list)

    # upstairs or downstairs sources
    down=True
    if detector == "HPGe":
        down=False

    abs_eff_dict = {}
    for source in source_dict:
        abs_eff_dict[source] = {}
        for angle in source_dict[source]:

            date = source_dict[source][angle][0][0] # get the date of measurement
            time = source_dict[source][angle][0][1] # get the time of measurement
            # time since calibration
            time_since_cal = calculate_time_since_source_calibration(date, time, down)

            # choose the quantities for the source
            initial_ac, decay_consts, activity_ucerts, dc_ucerts = choose_activities_and_decay_consts(detector, source)

            activity_list = []
            activity_ucert_list = []

            # if there is only one peak for the source, make the quantities into lists so zip() can be used generally
            if type(initial_ac) is not np.ndarray:
                initial_ac = [initial_ac]
                decay_consts = [decay_consts]
                activity_ucerts = [activity_ucerts]
                dc_ucerts = [dc_ucerts]

            # calculate the current activity and uncertainties
            for init_ac, dc_const, ac_ucert, dc_ucert in zip(initial_ac, decay_consts, activity_ucerts, dc_ucerts):
                current_act = activity_func(init_ac, dc_const, time_since_cal)
                current_act_ucert = current_activity_ucerts(init_ac, ac_ucert, time_since_cal, 1, dc_const, dc_ucert)

                activity_list.append(current_act)
                activity_ucert_list.append(current_act_ucert)

            peak_areas = source_dict[source][angle][1] # get the peak areas
            peak_areas = np.array(peak_areas)

            peak_area_ucerts = source_dict[source][angle][2] # get the peak area uncertainties
            peak_area_ucerts = np.array(peak_area_ucerts)

            activity_arr = np.array(activity_list)
            activity_ucert_arr = np.array(activity_ucert_list)

            # absolute efficiency and uncertainty calculation
            absolute_efficiency = peak_areas/activity_arr
            absolute_efficiency_ucerts = absolute_efficiency * uncertainty_mult_or_div(peak_areas, activity_arr,
                                                                        peak_area_ucerts, activity_ucert_arr)

            # populate the dictionary
            abs_eff_dict[source][angle] = [absolute_efficiency, absolute_efficiency_ucerts]

    return abs_eff_dict


def calculate_solid_angle(detector, angle):
    """Given detector and angle - calculate the solid angle"""

    # convert to radians
    angle = angle * (np.pi / 180)

    if detector == "NaI(Tl)":
        distance = 15
        diameter = 5.8
        length = 6.5
    elif detector == "BGO":
        distance = 15
        diameter = 5.5
        length = 4
    else:
        distance = 10
        diameter = 4.8
        length = 3.7

    return ((np.pi * (diameter / 2)**2) / distance**2) * (np.cos(angle))**2 + ((2 * (diameter / 2) * length) /
                                                                               distance**2) * (np.sin(angle))**2


def calculate_int_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list):
    """
    Create a dictionary of intrinsic efficiencies for all data for the detector
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    :return: dictionary {source:{angle:[intrinsic efficiency, uncertainties]}}
    """

    # generate absolute efficiency dict
    abs_eff_dict = calculate_abs_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list)

    int_eff_dict = {}
    for source in sources:
        int_eff_dict[source] = {}
        for angle in abs_eff_dict[source]:
            # calculate the solid angle
            solid_angle = calculate_solid_angle(detector, angle)
            # absolute efficiencies from dictionary
            abs_eff_arr = abs_eff_dict[source][angle][0]
            # absolute efficiency uncertainties from dictionary
            abs_eff_ucerts_arr = abs_eff_dict[source][angle][1]
            # multiply by geometric factor to calculate the intrinsic efficiency
            int_eff_arr = ((np.pi * 4) / solid_angle) * abs_eff_arr
            # crudely assuming no uncertainty in geometric factor, calculate the uncertainties in intrinsic efficiency
            int_eff_ucerts_arr = ((np.pi * 4) / solid_angle) * abs_eff_ucerts_arr
            # populate the dictionary
            int_eff_dict[source][angle] = [int_eff_arr, int_eff_ucerts_arr]

    return int_eff_dict


def efficiency_model(e_array, a, b, c):
    """Function to fit to the efficiency vs energy data"""
    return np.exp(a + b * np.log(e_array) + c * (np.log(e_array))**2)


def fit_efficiency(energies, efficiencies, **kwargs):
    """Calculate the efficiency fit params"""
    eff_popt, eff_pcov = curve_fit(efficiency_model, energies, efficiencies, **kwargs)
    return eff_popt, eff_pcov


def abs_efficiency_plot(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list, energies,
                        e_range):
    """
    Plot absolute efficiency vs energy and fit model
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    :param energies: energies from calibration
    :param e_range: samples over the range of energies
    """

    # generate absolute efficiency dict
    abs_eff_dict = calculate_abs_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list)

    abs_eff_arr = []
    abs_eff_ucert_arr = []

    for source in sources:
        # get the 0 angle absolute efficiencies and uncertainties for each source
        current_eff_arr = abs_eff_dict[source][0][0]
        current_eff_ucerts_arr = abs_eff_dict[source][0][1]
        abs_eff_arr += current_eff_arr.tolist()
        abs_eff_ucert_arr += current_eff_ucerts_arr.tolist()

    abs_eff_arr = np.asarray(abs_eff_arr)
    abs_eff_ucert_arr = np.asarray(abs_eff_ucert_arr)

    # the data from the first americium peak is excluded because the uncertainty in efficiency is larger than the
    # efficiency itself
    if detector == "NaI(Tl)":
        energies = energies[1:]
        abs_eff_arr = abs_eff_arr[1:]
        abs_eff_ucert_arr = abs_eff_ucert_arr[1:]
        e_range = np.logspace(1.7, 2.9, 5000)

    # the efficiency for the BGO detector at the Americium energy was far to low to allow for a reasonable fit, it is
    # excluded here but this results in the energy range of efficiencies being reduced for this detector
    if detector == "BGO":
        energies = energies[1:]
        abs_eff_arr = abs_eff_arr[1:]
        abs_eff_ucert_arr = abs_eff_ucert_arr[1:]
        e_range = np.logspace(2.4, 3.15, 5000)

    # the data from the first americium peak is excluded because the uncertainty in efficiency is larger than the
    # efficiency itself
    if detector == "HPGe":
        energies = energies[1:]
        abs_eff_arr = abs_eff_arr[1:]
        abs_eff_ucert_arr = abs_eff_ucert_arr[1:]
        e_range = np.logspace(1.7, 3.17, 5000)

    abs_eff_fig, abs_eff_ax = plt.subplots()
    abs_eff_ax.set_xlabel("Energy (keV)")
    abs_eff_ax.set_ylabel("Absolute Efficiency (%)")
    abs_eff_ax.set_title(f"Absolute Efficiency for {detector}")
    abs_eff_ax.set_xscale('log')
    abs_eff_ax.set_yscale('log')

    # plot
    abs_eff_ax.errorbar(energies, abs_eff_arr * 100, yerr=abs_eff_ucert_arr * 100, label="Absolute Peak Efficiencies",
                        fmt='.', c='r')

    # calculate model params
    eff_popt, eff_pcov = fit_efficiency(energies, abs_eff_arr * 100, absolute_sigma=True, sigma=abs_eff_ucert_arr)

    # uncertainties in model paramas
    eff_ucerts = np.sqrt(np.diag(eff_pcov))

    # absolute efficiencies using fit function
    abs_eff_model_arr = efficiency_model(e_range, *eff_popt)

    fit_label = f"ln($\\epsilon$) = {eff_popt[0]:.3} $\\pm$ {eff_ucerts[0]:.3} + {eff_popt[1]:.3}" \
                f" $\\pm$ {eff_ucerts[1]:.3} $\\times$ ln(E) + {eff_popt[2]:.3} $\\pm$ {eff_ucerts[2]:.3} $\\times$" \
                f"(ln(E))$^{{{2}}}$"

    # model plot
    abs_eff_ax.plot(e_range, abs_eff_model_arr, label=fit_label)

    abs_eff_fig.legend()

    return


def int_efficiency_plot(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list, energies,
                        e_range):
    """
    Plot intrinsic efficiency vs energy and fit model
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    :param energies: energies from the calibration
    :param e_range: samples over the range of energies
    """

    # generate intrinsic efficiency dict
    int_eff_dict = calculate_int_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list)

    int_eff_arr = []
    int_eff_ucert_arr = []

    for source in sources:
        # get the 0 angle intrinsic efficiencies and uncertainties for each source
        current_int_eff_arr = int_eff_dict[source][0][0]
        current_int_eff_ucerts_arr = int_eff_dict[source][0][1]
        int_eff_arr += current_int_eff_arr.tolist()
        int_eff_ucert_arr += current_int_eff_ucerts_arr.tolist()

    int_eff_arr = np.asarray(int_eff_arr)
    int_eff_ucert_arr = np.asarray(int_eff_ucert_arr)

    # the data from the first americium peak is excluded because the uncertainty in efficiency is larger than the
    # efficiency itself
    if detector == "NaI(Tl)":
        energies = energies[1:]
        int_eff_arr = int_eff_arr[1:]
        int_eff_ucert_arr = int_eff_ucert_arr[1:]
        e_range = np.logspace(1.7, 2.9, 5000)

    # the efficiency for the BGO detector at the Americium energy was far to low to allow for a reasonable fit, it is
    # excluded here but this results in the energy range of efficiencies being reduced for this detector
    if detector == "BGO":
        energies = energies[1:]
        int_eff_arr = int_eff_arr[1:]
        int_eff_ucert_arr = int_eff_ucert_arr[1:]
        e_range = np.logspace(2.4, 3.15, 5000)

    # the data from the first americium peak is excluded because the uncertainty in efficiency is larger than the
    # efficiency itself
    if detector == "HPGe":
        energies = energies[1:]
        int_eff_arr = int_eff_arr[1:]
        int_eff_ucert_arr = int_eff_ucert_arr[1:]
        e_range = np.logspace(1.7, 3.17, 5000)

    int_eff_fig, int_eff_ax = plt.subplots()
    int_eff_ax.set_xlabel("Energy (keV)")
    int_eff_ax.set_ylabel("Intrinsic Efficiency (%)")
    int_eff_ax.set_title(f"Intrinsic Efficiency for {detector}")
    int_eff_ax.set_xscale('log')
    int_eff_ax.set_yscale('log')
    # plot

    int_eff_ax.errorbar(energies, int_eff_arr * 100, yerr=int_eff_ucert_arr * 100, label="Intrinsic Efficiencies",
                        fmt='.', c='r')

    # calculate model params
    eff_popt, eff_pcov = fit_efficiency(energies, int_eff_arr * 100, absolute_sigma=True, sigma=int_eff_ucert_arr)
    # uncertainties in model paramas
    eff_ucerts = np.sqrt(np.diag(eff_pcov))

    # intrinsic efficiencies using fit function
    abs_eff_model_arr = efficiency_model(e_range, *eff_popt)

    fit_label = f"ln($\\epsilon$) = {eff_popt[0]:.3} $\\pm$ {eff_ucerts[0]:.3} + {eff_popt[1]:.3}" \
                f" $\\pm$ {eff_ucerts[1]:.3} $\\times$ ln(E) + {eff_popt[2]:.3} $\\pm$ {eff_ucerts[2]:.3} $\\times$" \
                f"(ln(E))$^{{{2}}}$"

    # model plot
    int_eff_ax.plot(e_range, abs_eff_model_arr, label=fit_label)

    int_eff_fig.legend()

    return


def off_axis_int_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list):
    """
    Plot intrinsic efficiency vs angle for each source
    :param spectrum_file_path: path to source folders
    :param config_file_path: config file path
    :param bkg_file: background file
    :param detector: detector name
    :param sources: list of sources for that detector
    :param file_list: list of possible source files
    """

    # generate intrinsic efficiency dict
    int_eff_dict = calculate_int_eff(spectrum_file_path, config_file_path, bkg_file, detector, sources, file_list)

    off_axis_fig, off_axis_ax = plt.subplots()

    colors = ["red", "blue", "green", "orange"]

    # plot the intrinsic efficiencies as a function of angle
    for ind, source in enumerate(sources):
        if detector == "BGO":
            # remove americium from the bgo data
            if source == "Americium":
                continue
        for angle in int_eff_dict[source]:
            int_effs = int_eff_dict[source][angle][0]
            if len(int_effs) > 1:
                int_effs = int_effs[1]
            else:
                int_effs = int_effs[0]
            if angle == 0:
                label = f"Off axis response for {source}"
            else:
                label = None
            # plot
            off_axis_ax.scatter(angle, int_effs * 100, marker='.', c=colors[ind], label=label)

    off_axis_fig.legend()

    off_axis_ax.set_xlabel("Angle (Deg)")
    off_axis_ax.set_ylabel("Intrinsic Efficiency (%)")
    off_axis_ax.set_title(f"Intrinsic Efficiency vs Angle for {detector}")
    off_axis_ax.set_yscale('log')

    return
