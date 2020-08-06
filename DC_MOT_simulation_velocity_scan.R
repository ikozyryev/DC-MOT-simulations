# use rate equations to model magneto-optical trapping forces for BaH
# following Tarbutt, NJP (2015)

setwd("C:/Users/ikozy/Dropbox/work/calculations/BaH_NJP")
options(digits=12) #set digits

####### define experimental parameters ##########
date2day = "20200525"

molecule2use = "CaH"

# molecule2use = "SrF"
# molecule2use = "CaF"

res2save = T # indicator whether to save the results to file

source('20200121_BaH_MOT_functions.R')

if (molecule2use == "BaH"){
  # read in relevant functions and constants
  source('20200121_BaH_MOT_constants.R')
}

if (molecule2use == "SrF"){
  # read in relevant functions and constants
  source('20200121_SrF_MOT_constants.R')
}


if (molecule2use == "CaF"){
  # read in relevant functions and constants
  source('20200121_CaF_MOT_constants.R')
}

if (molecule2use == "CaH"){
  # read in relevant functions and constants
  source('20200121_CaH_MOT_constants.R')
}

file2save = paste(date2day,'_',molecule2use,'_DC_250mW_tot_10Gpercm_1by2meppp3by2pppme_vel_scan.txt',sep='')

# file2save = paste(molecule2use,'_DC_MOT_pppm_tot200mW_15Gpercm_12mm_zero_gu.txt',sep='')
# file2save = paste(molecule2use,'_DC_dual_MOT_pmmpmp_inverted_vel_scan_tot200mW_15Gpercm_12mm_true_gl_true_gu_20200306_.txt',sep='')

# load the libraries
library(plotly)
library(deSolve)
library(ggplot2)

w=12e-3 # [m] 1/e2 radius of the trapping beams used in the experiment
t_step=0.01/Gamma_n

# six laser beams from each direction because of the EOM 20 MHz splittings

# ignore angle effects
#k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))
# for 4 beams from each direction
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#            c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))
# for 5 beams from each direction for congs II and III
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#             c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))

# for 3 beams from each direction
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),
#             c(0,0,-1),c(0,0,-1),c(0,0,-1))

# for BaH dual
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#             c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))

# for BaH dual
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#             c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))

# for 8 beams from each direction
k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
            c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))

# for 7 beams from each direction
# k_vec=rbind(c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#             c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1),c(0,0,-1))


Ptot=250e-3 # [W] total laser power in the APi light

beam_nums=dim(k_vec)[1] # number of beams/freq from one direction
P_laser=rep(Ptot/beam_nums*2,dim(k_vec)[1]) # [W] laser power per frequency component, note that beam are rertoreflected (i.e. x2 factor)

# beam_pols=c(1/2,1/sqrt(2),1/2) # polarization vector at 45 degree angle relative to the B-field
# pol_vec=beam_pols
# for (cnt in 2:beam_nums){
#   pol_vec=rbind(pol_vec,beam_pols)
# }

###### set the polarizations now for different configurations refer to Fig. 8
# correct polarizations here for config I
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1))
# correct polarizations here for config I for BaH
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(0,0,1))
# polarizations here for config II
# pol_vec=rbind(c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1),
#               c(1,0,0),c(0,0,1),c(0,0,1),c(0,0,1),c(0,0,1))
# polarizations here for config III
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),
#             c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0))
# wrong polarization here for checking
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1))
# for BaH dual
# pol_vec=rbind(c(0,0,1),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1),
#               c(0,0,1),c(1,0,0),c(1,0,0),c(0,0,1),c(1,0,0),c(0,0,1))

# for BaH dual
# pol_vec=rbind(c(1,0,0),c(0,0,1),c(0,0,1),c(1,0,0),c(0,0,1),c(1,0,0),
#               c(1,0,0),c(0,0,1),c(0,0,1),c(1,0,0),c(0,0,1),c(1,0,0))

# for CaH with 8 beams: use EOM to adress both J''=1/2 and J''=3/2; add to extra beam to the upper manifold
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1))

# for CaH with 7 beams: use EOM to adress both J''=1/2 and J''=3/2; add one extra beam to the upper manifold
# pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),
#               c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1))

# using Tarbutt DC MOT config added to J=1/2, F=1
pol_vec=rbind(c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1),
              c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(1,0,0),c(0,0,1),c(0,0,1))

# # for config I
# rel_detun=c(-1*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4])

# for config I for BaH with 3 beams only
# rel_detun=c(-1*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4])

# for config II
# rel_detun=c(-1*Gamma_n+Xstate_split[1],-2*Gamma_n+Xstate_split[1],-1.2*Gamma_n+Xstate_split[2],-1.2*Gamma_n+Xstate_split[3],-1.2*Gamma_n+Xstate_split[4]) 

# for config III
# rel_detun=c(-1.2*Gamma_n+Xstate_split[1],-1.2*Gamma_n+Xstate_split[2],-1.2*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4],-2*Gamma_n+Xstate_split[4]) 

# for BaH dual freq
#rel_detun=c(2*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[1],2*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[4],-1*Gamma_n+Xstate_split[4]) 

EOMfreq=52*2*pi*1e6 # [MHz]
# 7 laser beams total; refer to Tarbutt NJP paper for the diagram
# for 8 beams
# rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,-1*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state
# for 7 beams
# rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state

rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,2*Gamma_n+Xstate_split[1],-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state


#plot(seq(1,10),rep(rel_detunMHz[1],10),ylim=c(-50,50),type='l',lwd=2,col='black')
#plot(seq(1,10),rep(rel_detunMHz[1],10),ylim=-1*c(8642.8-50,8642.8+50),type='l',lwd=2,col='black')
#for (i in 2:length(rel_detunMHz)){
#abline(h=rel_detunMHz[i],lwd=2)
#}
# rel_detun=rel_detunMHz*2*pi*1e6
# delta_lup=laser_beam_detun12(Xstate_split,rel_detun,EOMfreq) # import all the relative detunings
# delta_lup=rep(0,6) # make sure the units are correct
# initialize polarization vecotors for each laser beam with multiple frequencies
# gu=c(-0.51,-0.51,-0.51,0)
#gu=rep(0,4)
# delta_lup=laser_beam_detun14EOM(Xstate_split,rel_detun)# import all the relative detunings
delta_lup=laser_beam_detun16EOM(Xstate_split,rel_detun)# import all the relative detunings

# 6 beams for config I
#delta_lup=laser_beam_detun6(Xstate_split,rel_detun)# import all the relative detunings
# 10 beams for configs II and III
#delta_lup=laser_beam_detun12(Xstate_split,rel_detun)# import all the relative detunings
#delta_lup=laser_beam_detun_10iii(Xstate_split,rel_detun)# import all the relative detunings

# EOMfreq=20
# rel_detunMHz=c(EOMfreq,0,-EOMfreq,-1*8642.8,-1*(8642.8+EOMfreq),-1*(8642.8+2*EOMfreq))# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state
#plot(seq(1,10),rep(rel_detunMHz[1],10),ylim=c(-50,50),type='l',lwd=2,col='black')
#plot(seq(1,10),rep(rel_detunMHz[1],10),ylim=-1*c(8642.8-50,8642.8+50),type='l',lwd=2,col='black')
#for (i in 2:length(rel_detunMHz)){
  #abline(h=rel_detunMHz[i],lwd=2)
#}
# rel_detun=rel_detunMHz*2*pi*1e6
# delta_lup=laser_beam_detun12(Xstate_split,rel_detun,EOMfreq) # import all the relative detunings
# delta_lup=rep(0,6) # make sure the units are correct
# initialize polarization vecotors for each laser beam with multiple frequencies
# gu=c(-0.51,-0.51,-0.51,0)
#gu=rep(0,4)

Afield=10e-2 # [T/m] 
#Afield=0 #1e-4 # [T] B field constant

t_steps_num=5000 # number of time steps

#vel=c(0,0,0) # [m/s] velocity vector

#pol_vec_local=rep(NA,3) # for storing local laser polarization

##### Generic loop for calculating the rate equations
parmax=20 # max value of the parameter
#parvals=seq(-parmax,parmax,by=1) # parameter values for the loop iterations
parvals=seq(0,parmax,by=0.5)
parsteps=length(parvals) # number of steps to take
# accelz_store=matrix(NA, nrow = t_steps_num, ncol = parsteps)
### store acceleration value and photons scattered
accelz_store=rep(NA,parsteps)
#gamma_store=rep(NA,parsteps)
#pos_final=rep(NA,parsteps)
#vel_final=rep(NA,parsteps)
# gamma_sp=matrix(0, nrow = t_steps_num, ncol = parsteps) # record the number of scattered photons
time_seq = (1:t_steps_num)*t_step


pos=c(0,0,0) # the spatial potition to use
# calculate the laser power at the specific position
Ip=2*P_laser[1]/(pi*w^2)*exp(-2*(pos[1]^2+pos[2]^2)/w^2)
Rabi_ij=kij*trans_dipole*sqrt(2*Ip/(hbar^2*c_light*eps0)) # calculate the Rabi frequency at that position
#flup_satp=4*kij^2*Ip*trans_dipole^2/(hbar^2*c_light*eps0*Gamma_n^2)
#sat_p=Ip/Isat2level
#vel = c(0,0,1)
#setwd("C:/Users/ikozy/Dropbox/work/calculations/BaH_NJP/velocity_data9")

# velz_seq = seq(50,320,by = 15)
# vel_seq = seq(-1.75,1.75,by = 0.25)
#vel_seq = 

#for (velcnt in 1:length(vel_seq)) {
 # print(c('velcnt ',velcnt))
  for (parcnt in 1:parsteps){ # loop over the parameter counter now which is position here
    #parcnt=2
    print(parcnt)
    #rel_detun=(rel_detunMHz+parvals[parcnt])*2*pi*1e6
    #delta_lup=laser_beam_detun12EOM(Xstate_split,rel_detun) # import all the relative detunings
    deltalup_all=array(0,dim=c(12,4,beam_nums)) # create array for storing scattering rates from all beams
    # here you can account for the hyperfine splitting in the excited state
    for (q in 1:(beam_nums)){
      deltalup_all[,,q] = cbind(delta_lup[,q],delta_lup[,q],delta_lup[,q],delta_lup[,q])
    }
    # Rlup_all=scattering_rate_matrix_2OMvec(pos,vel,Afield,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij) # calculate the t_steps_numing array  
    #(pos,vel,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij)
    #Rlup_all=scattering_rate_matrix_2OMvec(pos,vel,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij)
    
    #### this part to be replaced by the RK4 functions; scan the velocity now at the center of the trap (almost to have a well defined quantization axis)
    init_vals = c(1e-4,parvals[parcnt],as.matrix(num_X_init),as.matrix(num_A_init),0)
    #init_vals = c(1e-3,1,rep(1/16,16),0)
    #init_vals=sol[2,2:20]
    prefac = (h_Planck/(mass*lambda))
    # params = list(Rlup_all,r_lu,prefac) # for derivs
    pos = c(0,0,1e-4)
    Bfield=Afield*c(pos[1],pos[2],-2*pos[3]) # magnetic field
    
    # print(angle)
   
    Rlup_all=scattering_rate_matrix_2OMvec(pos,c(180,0,parvals[parcnt]),k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij,Bfield)
    params = list(Rlup_all,r_lu,prefac) # for derivs where Rlup_all is calculated outside the loop
    # sol <- rk(init_vals,time_seq,derivs,params,method="rk45dp7") # implements dynamic time stepping 
    sol <- rk(init_vals,time_seq,derivs,params,method="rk4") # implements dynamic time stepping 
    max_ind = dim(sol)[1]
    Nl = sol[max_ind,4:15]
    Nu = sol[max_ind,16:19]
    accelz_store[parcnt] = calc_accel(Nl,Nu,Rlup_all,prefac)
    #pos_final[parcnt] = sol[max_ind,2] # final position
    #vel_final[parcnt] = sol[max_ind,3] # final position
    #gamma_store[parcnt] = sol[max_ind,20] # final scattered photons
    #plot(sol[,1],sol[,2])
    
    #plot(sol[,1],sol[,3])
    #dftot = data.frame('time_lab'=sol[,1],sol[,4:19])
    #dftotsum = data.frame('time_lab'=sol[,1],'tot_sum'=apply(sol[,4:19],MARGIN=1,sum))
    # dfexcited = data.frame('time_lab'=sol[,1],sol[,16:19])
    
    # ggplot(dftot,aes(time_lab)) + geom_line(aes(y=X3),size=2) + geom_line(aes(y=X4),size=2) + geom_line(aes(y=X5),size=2) + geom_line(aes(y=X6),size=2) +
    #   geom_line(aes(y=X7),size=2) + geom_line(aes(y=X8),size=2) + geom_line(aes(y=X9),size=2) + geom_line(aes(y=X10),size=2) + geom_line(aes(y=X11),size=2) + geom_line(aes(y=X12),size=2) +
    #   geom_line(aes(y=X13),size=2) + geom_line(aes(y=X14),size=2) + geom_line(aes(y=X15),col='blue',size=2) + geom_line(aes(y=X16),col='blue',size=2) + geom_line(aes(y=X17),col='blue',size=2) + geom_line(aes(y=X18),col='blue',size=2)
    # 
    # ggplot(dftotsum,aes(time_lab)) + geom_line(aes(y=tot_sum))
    
    # for (i in 1:t_steps_num){
    #   Nldot=rep(0,dim(Xstates)[1])
    #   Nudot=rep(0,dim(Astates)[1])
    #   # print(i)
    #   for (u in 1:dim(Astates)[1]) {
    #     Nldot=Nldot+(Rlup_all[,u,1]+Rlup_all[,u,2])*(num_A[i,u]-num_X[i,])+Gamma_n*r_lu[,u]*num_A[i,u]
    #      if (i<t_steps_num) gamma_sp[i+1,parcnt]=gamma_sp[i,parcnt]+num_A[i,u] # count photons
    #   }
    #   for (l in 1:dim(Xstates)[1]) {
    #     Nudot=Nudot+(Rlup_all[l,,1]+Rlup_all[l,,2])*(num_X[i,l]-num_A[i,])
    #   }
    #   Nudot=Nudot-Gamma_n*num_A[i,]
    #   if (i<t_steps_num){
    #    num_X[i+1,]=num_X[i,]+Nldot*t_step
    #    num_A[i+1,]=num_A[i,]+Nudot*t_step
    #   }
    # }
    # for (i in 1:t_steps_num){
    #  #print(i)
    #  # accel=c(0,0,0)
    #   accel=0
    #  # delta_accel
    #  for (l in 1:dim(Xstates)[1]) {
    #    for (u in 1:dim(Astates)[1]) {
    #      # accel=accel+(h_Planck/(mass*lambda))*k_vec[p,]*R_lup*(num_X[i,l]-num_A[i,u])
    #      accel=accel+(h_Planck/(mass*lambda))*(Rlup_all[l,u,1]-Rlup_all[l,u,2])*(num_X[i,l]-num_A[i,u])
    #      #     # print(c(delta_accel*1e-3,Xstates[l,2],Astates[u,2]))
    #     }
    #   }
    #  accelz_store[i,parcnt]=accel
    # }
    ###########
  }

if (res2save == T){
  data2save = data.frame(parvals,accelz_store*1e-3)
  write.table(data2save,file=file2save,row.names = F)
}
  
#}

#plot(parvals,gamma_store,col='blue',lwd=2)
plot(parvals,accelz_store*1e-3,col='blue',lwd=2,xlab="Velocity (m/s)",ylab="Acceleration (km/s^2)",xlim=c(0,20))
abline(h=0,lty=2,lwd=2)

plot(parvals,pos_final,col='blue',lwd=2)
plot(parvals,vel_final,col='green',type='l',lwd=2)
