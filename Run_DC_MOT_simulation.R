# run MOT simulation

date2day = '20200806' # todays date
molecule = 'CaH'
scan_type = 'position' # 'velocity'
save_results = T
plot_results = T
laser_power = 250 # [mW] total laser power for all beams
beam_waist=12e-3 # [m] 1/e2 RADIUS of the trapping beams used in the experiment
freq_number = 8 # number of laser frequencies
# specify the polarization vector for the laser frequencies
polarization_vector = rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1)) # has to be the same length as freq_num
# encode relative detunings now for all beams
rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,2*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state


B_field_grad = 15 #[Gauss/cm] gradient of the magnetic field in the xy plane


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
  
  