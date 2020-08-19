# TIP: press Ctrl+Shift+S to run the whole script
# TIP: to see only the column names and first few rows of a dataframe, use head(DFname), where DFname is the name of the dataframe

# load the needed libraries; install them first if needed
library(plotly) # dynamic plotting
library(deSolve) # ODE/PDE solver
library(ggplot2) # plotting
library(ggpubr) # package for combining multiple ggplots

# enter simulation parameters
date2day = '20200819' # todays date
molecule = 'SrF' # molecule to use: CaH, BaH, SrF, CaF.
scan_type = 'position' # 'velocity'
save_results = F # save the data for acceleration vs position or velocity
plot_results = T # plot the results
save_plots = T
plot_populations = F # make sure that the final value for the scan range is set to something suitable for proper population dynamics to occur
laser_power = 200 # [mW] TOTAL laser power for all beams
beam_waist=7e-3 # [m] 1/e2 RADIUS of the trapping beams used in the experiment
freq_number = 4 # number of laser frequencies addressing molecules

# load molecular parameters and used functions
source('Load_molecular_params.R')
source('MOT_auxiliary_functions.R')
# set simulation parameters 
t_step=0.01/Gamma_n
t_steps_num=100 # number of time steps

# specify the polarization vector for the laser frequencies

#### correct polarizations here for config I
polarization_vector=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1))
pol_encode = 'configI_pppm' # mnemonic for encoding laser polarization to use in the name of the saved file; does not effect the simulation results
rel_detun=c(-1*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4])

#### polarizations here for config II
# polarization_vector=rbind(c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1))
# pol_encode = 'configII_pmmmm' 
# rel_detun=c(-1*Gamma_n+Xstate_split[1],-2*Gamma_n+Xstate_split[1],-1.2*Gamma_n+Xstate_split[2],-1.2*Gamma_n+Xstate_split[3],-1.2*Gamma_n+Xstate_split[4]) 

#### polarizations here for config III
# polarization_vector=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0))
# pol_encode = 'configIII_pppmp' # mnemonic for encoding laser polarization to use in the name of the saved file; does not effect the simulation results
# rel_detun=c(-1.2*Gamma_n+Xstate_split[1],-1.2*Gamma_n+Xstate_split[2],-1.2*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4],-2*Gamma_n+Xstate_split[4]) 

freq_number = dim(polarization_vector)[1] # extract the number of frequencies
# delta_lup = configure_all_detunings(Xstate_split,rel_detun)
delta_lup = configure_detunings(freq_number,Xstate_split,rel_detun)

# polarization_vector = rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1)) # has to be the same length as freq_num
# set_detun = c(-1,-1,-1,2,2,2,2,-1) # set detunings in units of natural linewidth for the transition
# EOMfreq=52*(2*pi*1e6) # [MHz]

B_field_grad = 7.5 #[Gauss/cm] gradient of the magnetic field in the xy plane

if (scan_type == 'position'){
  max_val = 10e-3 # [meters] max value of position
  step_size = 1e-3 # [meters]
  source('DC_MOT_simulation_position_scan.R')
}

if (scan_type == 'velocity'){
  max_val = 20 # [m/s] max value of position
  step_size = 1 # [m/s]
  source('DC_MOT_simulation_velocity_scan.R')
}
