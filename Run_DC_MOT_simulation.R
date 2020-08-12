# load the needed libraries

library(plotly) # dynamic plotting
library(deSolve) # ODE/PDE solver
library(ggplot2) # plotting
# run MOT simulation

date2day = '20200812' # todays date
molecule = 'CaH'
scan_type = 'position' # 'velocity'
save_results = F
plot_results = T
plot_populations = T # make sure that the final value for the scan range is set to something suitable for proper population dynamics to occur
laser_power = 250 # [mW] total laser power for all beams
beam_waist=12e-3 # [m] 1/e2 RADIUS of the trapping beams used in the experiment
freq_number = 8 # number of laser frequencies
# specify the polarization vector for the laser frequencies
polarization_vector = rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1)) # has to be the same length as freq_num
set_detun = c(-1,-1,-1,2,2,2,2,-1) # set detunings in units of natural linewidth for the transition
EOMfreq=52*(2*pi*1e6) # [MHz]


B_field_grad = 15 #[Gauss/cm] gradient of the magnetic field in the xy plane

if (scan_type == 'position'){
  max_val = 5e-3 # [meters] max value of position
  step_size = 1e-3 # [meters]
  source('DC_MOT_simulation_position_scan.R')
}

if (scan_type == 'velocity'){
  max_val = 20 # [m/s] max value of position
  step_size = 1 # [m/s]
  source('DC_MOT_simulation_velocity_scan.R')
}



###### set the polarizations now for different configurations refer to Fig. 8
# correct polarizations here for config I
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0))
# polarizations here for config II
# pol_vec=rbind(c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#               c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1))
# polarizations here for config III
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),
#             c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0))
# wrong polarization here for checking
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1))
# # for CaH with 7 beams
# pol_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(1,0,0),
#               c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(1,0,0))

# for CaH with 7 beams: use EOM to adress both J''=1/2 and J''=3/2; add one extra beam to the upper manifold
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1))

# # for CaH with 8 beams: use EOM to adress both J''=1/2 and J''=3/2; add two extra beams to the upper manifold
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1))
#                c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1))

# # for CaH with 8 beams: use EOM to adress both J''=1/2 and J''=3/2; add two extra beams to the upper manifold
# # using Tarbutt DC MOT config
# pol_vec=rbind(c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(1,0,0),c(1,0,0))

# using Tarbutt DC MOT config added to J=1/2, F=1
#pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1))
  
  