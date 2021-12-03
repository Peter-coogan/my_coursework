
import matplotlib.pyplot as plt

plt.rc('legend', fontsize=5)


NAI_FILE_PATH = "C:/Users/pfbas/Desktop/MSc Space Science and Technology/Space Detector Lab/" \
                         "Gamma Ray Lab/Lab Data/NaI(Tl)"
NAI_BKG_FILE = "/bkground.spe"
NAI_ENERGY_CAL_CONFIG_FILE = NAI_FILE_PATH + "/NAI_Energy_Calibration_config.json"


BGO_FILE_PATH = "C:/Users/pfbas/Desktop/MSc Space Science and Technology/Space Detector Lab/" \
                         "Gamma Ray Lab/Lab Data/BGO"
BGO_BKG_FILE = "/background.spe"
BGO_ENERGY_CAL_CONFIG_FILE = BGO_FILE_PATH + "/BGO_Energy_Calibration_config.json"


HPGe_FILE_PATH = "C:/Users/pfbas/Desktop/MSc Space Science and Technology/Space Detector Lab/" \
                         "Gamma Ray Lab/Lab Data/HPGe"
HPGe_BKG_FILE = "/background.spe"
HPGe_ENERGY_CAL_CONFIG_FILE = HPGe_FILE_PATH + "/HPGe_Energy_Calibration_config.json"

def main():

    # UNCOMMENT EACH BLOCK TO RUN ENTIRE PIPELINE FOR EACH DETECTOR

    # This is not very readable and could certainly be condensed into a run_pipeline() function that takes the detector
    # as an arguement

    '''
    ####################################################################################################################
    detector = "NaI(Tl)"
    e_range_nai = np.linspace(15, 700, 5000)
    e_range_log_nai = np.logspace(1.3, 2.9, 5000)

    line_fit = calibration_resolution_analysis(detector, NAI_FILE_PATH, NAI_BKG_FILE, NAI_ENERGY_CAL_CONFIG_FILE,
     e_range_nai, e_range_log_nai)
    plt.show()

    NAI_folders = ["Americium", "Barium", "Caesium"]
    NAI_files = ["/angle_0.spe", "/angle_20.spe", "/angle_40.spe", "/angle_60.spe", "/angle_80.spe",
                     "/angle_100.spe",
                     "/angle_120.spe", "/angle_-20.spe", "/angle_-40.spe", "/angle_-60.spe", "/angle_-80.spe",
                     "/angle_-100.spe", "/angle_-120.spe"]

    abs_efficiency_plot(NAI_FILE_PATH, NAI_ENERGY_CAL_CONFIG_FILE, NAI_BKG_FILE, "NaI(Tl)", NAI_folders, NAI_files,
                        line_fit, e_range_log_nai)
    plt.show()

    int_efficiency_plot(NAI_FILE_PATH, NAI_ENERGY_CAL_CONFIG_FILE, NAI_BKG_FILE, "NaI(Tl)", NAI_folders, NAI_files,
                        line_fit, e_range_log_nai)
    plt.show()

    off_axis_int_eff(NAI_FILE_PATH, NAI_ENERGY_CAL_CONFIG_FILE, NAI_BKG_FILE, "NaI(Tl)", NAI_folders, NAI_files)
    plt.show()
    ####################################################################################################################
    '''
    '''
    ####################################################################################################################
    detector = "BGO"
    e_range_bgo = np.linspace(25, 560, 5000)
    e_range_log_bgo = np.logspace(1.7, 3.15, 5000)
    line_fit = calibration_resolution_analysis(detector, BGO_FILE_PATH, BGO_BKG_FILE, BGO_ENERGY_CAL_CONFIG_FILE,
                                               e_range_bgo, e_range_log_bgo)
    plt.show()

    BGO_folders = ["Americium", "Barium", "Caesium", "Cobalt"]
    BGO_files = ["/angle_0.spe", "/angle_20.spe", "/angle_30.spe", "/angle_40.spe", "/angle_60.spe", "/angle_80.spe",
                 "/angle_90.spe", "/angle_100.spe", "/angle_120.spe", "/angle_-20.spe", "/angle_-30.spe",
                 "/angle_-40.spe", "/angle_-60.spe", "/angle_-80.spe", "/angle_-90.spe", "/angle_-100.spe",
                 "/angle_-120.spe"]

    abs_efficiency_plot(BGO_FILE_PATH, BGO_ENERGY_CAL_CONFIG_FILE, BGO_BKG_FILE, "BGO", BGO_folders, BGO_files,
                        line_fit, e_range_log_bgo)
    plt.show()

    int_efficiency_plot(BGO_FILE_PATH, BGO_ENERGY_CAL_CONFIG_FILE, BGO_BKG_FILE, "BGO", BGO_folders, BGO_files,
                        line_fit, e_range_log_bgo)
    plt.show()

    off_axis_int_eff(BGO_FILE_PATH, BGO_ENERGY_CAL_CONFIG_FILE, BGO_BKG_FILE, "BGO", BGO_folders, BGO_files)
    plt.show()
    ####################################################################################################################
    '''
    '''
    ####################################################################################################################
    detector = "HPGe"
    e_range_hpge = np.linspace(40, 5500, 5000)
    e_range_log_hpge = np.logspace(1.25, 3.17, 5000)
    line_fit = calibration_resolution_analysis(detector, HPGe_FILE_PATH, HPGe_BKG_FILE, HPGe_ENERGY_CAL_CONFIG_FILE,
     e_range_hpge, e_range_log_hpge)
    plt.show()

    HPGe_folders = ["Americium", "Barium", "Caesium", "Cobalt"]
    HPGe_files = ["/angle_0.spe", "/angle_30.spe", "/angle_60.spe", "/angle_90.spe", "/angle_-45.spe", "/angle_-90.spe"]

    abs_efficiency_plot(HPGe_FILE_PATH, HPGe_ENERGY_CAL_CONFIG_FILE, HPGe_BKG_FILE, "HPGe", HPGe_folders,
                                       HPGe_files, line_fit, e_range_log_hpge)
    plt.show()

    int_efficiency_plot(HPGe_FILE_PATH, HPGe_ENERGY_CAL_CONFIG_FILE, HPGe_BKG_FILE, "HPGe", HPGe_folders,
                                       HPGe_files, line_fit, e_range_log_hpge)
    plt.show()

    off_axis_int_eff(HPGe_FILE_PATH, HPGe_ENERGY_CAL_CONFIG_FILE, HPGe_BKG_FILE, "HPGe", HPGe_folders, HPGe_files)
    plt.show()
    ####################################################################################################################
    '''


if __name__ == "__main__":
    main()
