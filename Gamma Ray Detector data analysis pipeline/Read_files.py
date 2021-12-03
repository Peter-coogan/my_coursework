
import json


def mca_parser(file_object):
    """Parse .mca files"""
    spectrum_data = []

    add_data = False
    for line in file_object:
        line = line.rstrip()
        if line.startswith("<<END>>"):
            add_data = False

        if add_data:
            spectrum_data.append(line)

        if line.startswith("<<DATA>>"):
            add_data = True

        if line.startswith("REAL_TIME"):
            real_time = line.split("- ", 1)[1]

    return spectrum_data, real_time


def spe_parser(file_object):
    """Parse .spe files"""

    # lis to hold spectrum data
    spectrum_data = []

    line_counter = 0

    add_data = False
    data_start = False
    data_start_line = None
    date_start = False
    date_start_line = None
    meas_start = False
    meas_start_line = None

    # go line by line and use the strings in the files to pull data from lines
    for line in file_object:
        line_counter += 1
        line = line.strip()
        if line.startswith("$ROI"):
            add_data = False

        if add_data:
            spectrum_data.append(line)

        if line.startswith("$DATA"):
            data_start = True
            data_start_line = line_counter

        if data_start:
            if line_counter == data_start_line + 1:
                add_data = True

        if line.startswith("$DATE_MEA"):
            date_start = True
            date_start_line = line_counter

        if date_start:
            if line_counter == date_start_line + 1:
                date_and_time = line.split(" ")
                date, time = date_and_time[0], date_and_time[1]

        if line.startswith("$MEAS_TIM"):
            meas_start = True
            meas_start_line = line_counter

        if meas_start:
            if line_counter == meas_start_line + 1:
                real_time = line.split(" ")[0]

    return spectrum_data, real_time, date, time


def file_reader(file_path, file_name):
    """
    Use the file parsers above to extract the data
    :param file_path: path to folders for each source
    :param file_name: source and file - eg. /Americium/angle_0.spe
    :return: list of lists containing spectrum data and metadata
    """

    # extract the source and angle
    source = file_name.split("/", 1)[1].split("/", 1)[0]
    angle = file_name.split("_", 1)[1].split(".", 1)[0]
    angle = int(angle)

    # get the file type
    file_type = file_name.split(".", 1)[-1]

    with open(file_path + file_name) as f:

        if file_type == "mca":
            spectrum_data, real_time = mca_parser(f)

        elif file_type == "spe":
            spectrum_data, real_time, date, time = spe_parser(f)

        else:
            raise Exception("Invalid file type!")

    # cast to int
    spectrum_data = list(map(int, spectrum_data))
    real_time = float(real_time)
    metadata = [source, angle, real_time, date, time]
    # list of spectrum data and relevant metadata
    spectrum_and_metadata = [spectrum_data, metadata]

    return spectrum_and_metadata


def json_energy_calibration_reader(config_file_path, source):
    """
    Read the json config files
    :param config_file_path: path to file
    :param source: name of radioactive source
    :return: data from config file
    """

    with open(config_file_path) as f:
        data = json.load(f)
        data = data[source]

        source_and_file = data["source_and_file"]
        ROI = data["ROI"]
        zoomed_plot_params = data["Zoomed_plot_params"]
        model = data["Model"]
        slice_index = data["slice_index"]
        added_function_params = data["added function params"]

    return source_and_file, ROI, zoomed_plot_params, model, slice_index, added_function_params


def extract_spectrum_and_model_input(spectrum_file_path, config_file_path, source):
    """Call the functions above to obtain spectrum data and config file data"""

    source_and_file, ROI, zoomed_params, model, slice_index, added_function_params\
        = json_energy_calibration_reader(config_file_path, source)
    spectrum_and_metadata = file_reader(spectrum_file_path, source_and_file)
    return spectrum_and_metadata, ROI, zoomed_params, model, slice_index, added_function_params


def get_background(bkg_file_path, bkg_file):
    """Read the background file and extract the data"""

    file_type = bkg_file.split(".", 1)[-1]

    with open(bkg_file_path + bkg_file) as f:

        if file_type == "mca":
            spectrum_data, real_time = mca_parser(f)

        elif file_type == "spe":
            spectrum_data, real_time, date, time = spe_parser(f)

        else:
            raise Exception("Invalid file type!")

    # cast to int
    spectrum_data = list(map(int, spectrum_data))
    real_time = float(real_time)
    metadata = [real_time, date, time]
    # list of spectrum data and relevant metadata
    bkg_spectrum_and_meta = [spectrum_data, metadata]

    return bkg_spectrum_and_meta
