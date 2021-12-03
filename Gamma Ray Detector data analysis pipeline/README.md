
# Source Code Explanation
## Read_files.py
This file contains the functions which read the spectrum and background files. They output the channel and count data as well as the metadata. This file also reads the json configuration files which contain data on the regions of interest and other useful parameters for each source used for each detector.
## Energy_Calibration_Model_fitting.py
This file is responsible for fitting models to the peaks in the spectra. It plots the spectra, uses the regions of interest to located the peaks and depending on the detector and source, it fits the corresponding model. The peak model parameters and their uncertainties are calculated by functions in this file.
## Energy_Calibration_and_Resolution.py
This file uses the peak parameters from the previous file to plot and fit models for energy calibration and energy resolution.
## Efficiency.py
This file calculates the on axis absolute and intrinsic efficiencies as well as the off axis intrinsic efficiencies. It plots these and fits models to the first two. 
## Driver_script.py
This file contains the paths to the respective files and folders for each detector. It contains a single main() function in which sections for each detector can be uncommented and run individually. It imports the necessary functions from the other files above.
