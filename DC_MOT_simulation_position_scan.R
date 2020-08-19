# use rate equations to model magneto-optical trapping forces for molecules
# following Tarbutt, NJP 17, 015007 (2015)

options(digits=12) #set the number of digits
####### define experimental parameters ##########

file2save = paste(date2day,'_',molecule,'_',scan_type,'_',laser_power,'mW','_',B_field_grad,'Gpercm_',pol_encode,'.txt',sep='')

# initialize the k-vector for laser beams
# one can encode here angle dependence, etc.

# beams are coming from +z and -z directions
k_vec = c(0,0,1) 

# encode k-vector for +z
for (i in 2:freq_number){ # starts from 2 since we bind one entry already
  k_vec = rbind(k_vec,c(0,0,1)) # initialize
}
# encode k-vectors for -z
for (i in 1:freq_number){ # counter starts from 1 since we are appending to the existing list
  k_vec = rbind(k_vec,c(0,0,-1)) # initialize
}

beam_nums=2*freq_number # total number of laser beams involved (freq_numberx2 since beams are from both directions)
P_laser=rep(laser_power*1e-3/freq_number,freq_number) # [W] laser power per frequency component; notice conversion from mW to W

pol_vec = rbind(polarization_vector,polarization_vector) # to account for the fact that power comes from all sides

# encode relative detunings now for all beams
# beams are encoded by the way they are generated; refer to the associated writeup for further details
## TIP: use Ctrl+Shift+C to comment or uncomment a selected fraction
# rel_detun=c(set_detun[1]*Gamma_n+Xstate_split[2]+EOMfreq,
#             set_detun[2]*Gamma_n+Xstate_split[2],
#             set_detun[3]*Gamma_n+Xstate_split[2]-EOMfreq,
#             set_detun[4]*Gamma_n+Xstate_split[3],
#             set_detun[5]*Gamma_n+Xstate_split[3]-EOMfreq,
#             set_detun[6]*Gamma_n+Xstate_split[3]-2*EOMfreq,
#             set_detun[7]*Gamma_n+Xstate_split[1],
#             set_detun[8]*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state

# 7 laser beams total; refer to Tarbutt NJP paper for the diagram
# relative detunings for 8 beams
# rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,-1*Gamma_n+Xstate_split[3],-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state



# relative detunings for 7 beams
# rel_detun=c(-1*Gamma_n+Xstate_split[2]+EOMfreq,-1*Gamma_n+Xstate_split[2],-1*Gamma_n+Xstate_split[2]-EOMfreq,2*Gamma_n+Xstate_split[3],2*Gamma_n+Xstate_split[3]-EOMfreq,2*Gamma_n+Xstate_split[3]-2*EOMfreq,-1*Gamma_n+Xstate_split[4])# #rep(0,beam_nums) #-1*Gamma_n # relative detuning for each state


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

#delta_lup=laser_beam_detun16EOM(Xstate_split,rel_detun)# import all the relative detunings

Afield=B_field_grad*1e-2 # [T/m] 
#Afield=0 #1e-4 # [T] B field constant

#vel=c(0,0,0) # [m/s] velocity vector

#pol_vec_local=rep(NA,3) # for storing local laser polarization

##### Generic loop for calculating the rate equations
parmax=max_val # max value of the parameter
#parvals=seq(-parmax,parmax,by=1) # parameter values for the loop iterations
parvals=seq(0,parmax,by=step_size)
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
Ip=2*P_laser[1]/(pi*beam_waist^2)*exp(-2*(pos[1]^2+pos[2]^2)/beam_waist^2)
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
    # NOTICE: deltalup_all is the 3D array of size 12x4xbeam_nums
    # here you can account for the hyperfine splitting in the excited state
    for (q in 1:(beam_nums)){
      ### TIP: this is the place to introduce hyperfine structure in the excited state!
      deltalup_all[,,q] = cbind(delta_lup[,q],delta_lup[,q],delta_lup[,q],delta_lup[,q])
    }
    # Rlup_all=scattering_rate_matrix_2OMvec(pos,vel,Afield,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij) # calculate the t_steps_numing array  
    #(pos,vel,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij)
    #Rlup_all=scattering_rate_matrix_2OMvec(pos,vel,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij)
    
    init_vals = c(parvals[parcnt],0,as.matrix(num_X_init),as.matrix(num_A_init),0)
    #init_vals = c(1e-3,1,rep(1/16,16),0)
    #init_vals=sol[2,2:20]
    prefac = (h_Planck/(mass*lambda))
    # params = list(Rlup_all,r_lu,prefac) # for derivs
    pos = c(0,0,parvals[parcnt])
    Bfield=Afield*c(pos[1],pos[2],-2*pos[3]) # magnetic field
    
    # print(angle)
   
    Rlup_all=scattering_rate_matrix_2OMvec(c(0,0,parvals[parcnt]),c(180,0,0),k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij,Bfield)
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

if (save_results == T){
  data2save = data.frame(parvals*1e3,accelz_store*1e-3)
  write.table(data2save,file=file2save,row.names = F)
}

# plots the populations
#Xstate_names = c('JhalfF1Mneg1','JhalfF1M0','JhalfF1M')
#Xstates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0),c(1.5,1,-1),c(1.5,1,0),c(1.5,1,1),c(1.5,2,-2),c(1.5,2,-1),c(1.5,2,0),c(1.5,2,1),c(1.5,2,2))
#Astates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0))

if (plot_populations == T){
  dftot = data.frame('time_lab'=sol[,1]*1e6,sol[,4:19])
  names(dftot) <- c('time_lab','GS1','GS2','GS3','GS4','GS5','GS6','GS7','GS8','GS9','GS10','GS11','GS12','ES1','ES2','ES3','ES4')
  
  dftotsum = data.frame('time_lab'=sol[,1]*1e6,'tot_sum'=apply(sol[,4:19],MARGIN=1,sum))
  # dfexcited = data.frame('time_lab'=sol[,1],sol[,16:19])
  lwd4plts = 1 # set the plot linewidth
  pops <- ggplot(dftot,aes(time_lab)) + geom_line(aes(y=GS1),size=lwd4plts) + geom_line(aes(y=GS2),size=lwd4plts) + geom_line(aes(y=GS3),size=lwd4plts) + geom_line(aes(y=GS4),size=lwd4plts) +
    geom_line(aes(y=GS5),size=lwd4plts) + geom_line(aes(y=GS6),size=lwd4plts) + geom_line(aes(y=GS7),size=lwd4plts) + geom_line(aes(y=GS8),size=lwd4plts) + geom_line(aes(y=GS9),size=lwd4plts) + geom_line(aes(y=GS10),size=lwd4plts) +
    geom_line(aes(y=GS11),size=lwd4plts) + geom_line(aes(y=GS12),size=lwd4plts) + geom_line(aes(y=ES1),col='blue',size=lwd4plts) + geom_line(aes(y=ES2),col='blue',size=lwd4plts) + geom_line(aes(y=ES3),col='blue',size=lwd4plts) + geom_line(aes(y=ES4),col='blue',size=lwd4plts) +
    xlab("Time (microsec)") + ylab("Population fraction")
  
  # plot scattered photon number
  gamma_store = data.frame('time_lab'=sol[,1]*1e6,'photons' = sol[,20]) # scattered photons
  
  scat_plt <- ggplot(gamma_store,aes(time_lab,photons)) + geom_line(size = 1) + xlab("Time (microsec)") + ylab("Scattered photons (#)")
  
  pop_sum <- ggplot(dftotsum,aes(time_lab)) + geom_line(aes(y=tot_sum)) + xlab("Time (microsec)") + ylab("Population sum")
  
  dfaccel <- data.frame('pos'=parvals*1e3,'accel'=accelz_store*1e-3)
  accel_plt <- ggplot(dfaccel,aes(pos,accel)) + geom_point(size = 2) + xlab("Position (mm)") + ylab("Acceleration (km/s^2)")
  # plot(parvals*1e3,accelz_store*1e-3,col='blue',lwd=2,xlab="Position (mm)",ylab="Acceleration (km/s^2)")
  # abline(h=0,lty=2,lwd=2)
  figure <- ggarrange(pops,scat_plt,pop_sum,accel_plt,labels = c("A","B","C","D"),ncol=2, nrow=2)
  figure
}else{
  dfaccel <- data.frame('pos'=parvals*1e3,'accel'=accelz_store*1e-3)
  accel_plt <- ggplot(dfaccel,aes(pos,accel)) + geom_point(size = 2) + xlab("Position (mm)") + ylab("Acceleration (km/s^2)")
}

if (save_plots == T){
  file2save4plts = paste(date2day,'_',molecule,'_',scan_type,'_',laser_power,'mW','_',B_field_grad,'Gpercm_',pol_encode,'.pdf',sep='')
  ggsave(file2save4plts)
}

# plot(parvals*1e3,accelz_store*1e-3,col='blue',lwd=2,xlab="Position (mm)",ylab="Acceleration (km/s^2)")
# abline(h=0,lty=2,lwd=2)


#}

#plot(parvals,gamma_store,col='blue',lwd=2)


#plot(parvals,pos_final,col='blue',lwd=2)
#plot(parvals,vel_final,col='green',type='l',lwd=2)
