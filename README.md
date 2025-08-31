## Observations of Reduced Ventilation in Meridional Overturning Circulation: Evidence from Physical and Biogeochemical Changes in Repeat Observations Along 110°E
This project provides the data and codes to reproduce the figures of a paper submitted to JGR-Oceans
The raw data are saved in MATLAB data format. Please download all files into a folder and run each script under MATLAB. The code and data to reproduce key results are shown below. 


## Data

### 'd63.mat'
d63 is the interpolated voyage data from DM63
### 'v19.mat'
v19 is the interpolated voyage hydrology data from IN19
### 'ctd19.mat'
v19 is the interpolated voyage hydrology data from IN19
### 'DM63deep_station_position.mat' and 'IN19deep_station_position.mat'
These two mat files store the position information for the deep stations in two voyages
### 'Figure1.mat'
Figure1.mat includes the data for the run script 'Figure1.m' to create Figure1 in the manuscript.
This data contains dynamic height data from CARS and salinity, PV from CARS
### 'WOA23_oxy.mat'
WOA23_oxy is the oxygen data converted from WOA23 for the script 'Figure1.m'
### 'argo_interp_noise.mat' and 'oras4_interp_noise_04to17.mat'
They are the noise data calculated from Argo monthly and ORAS4, respectively
### 'bio_rms.mat'
This data includes the biogeochemical RMS error calculated from voyage data.

## script
### 'ra_projvar.m'
This is the script for interpolating vertical profiles
### 'oxy_auo_pho_nox.m'
This is the code used to generate Figure 7 in the manuscript.
### 'cmocean.m'
This is the code to add a colorbar to the figures
### 'Thickness_change.m'
This is the code for calculating the thickness change and crease figure 6 in the manuscript.


 
