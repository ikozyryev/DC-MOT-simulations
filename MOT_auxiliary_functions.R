# auxilliary function for running the simulations of MOTs
configure_detunings <- function(freq_number,Xstate_split,rel_detun){
  if (freq_number == 4){
    return(laser_beam_detun8(Xstate_split,rel_detun))
  }else if (freq_number == 5){
    return(laser_beam_detun10(Xstate_split,rel_detun))
  }else if (freq_number == 6){
    return(laser_beam_detun12EOM(Xstate_split,rel_detun))
  }else if (freq_number == 7){
    return(laser_beam_detun14EOM(Xstate_split,rel_detun))
  }else if (freq_number == 8){
    return(laser_beam_detun16EOM(Xstate_split,rel_detun))
  }else{
    return(F)
  }
}

# functions for SrF calculations
laser_beam_detun <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun+Xstate_split[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun+Xstate_split[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun+Xstate_split[3]-Xstate_detun
  delta_lup[,4]=rel_detun+Xstate_split[4]-Xstate_detun
  delta_lup[,5]=rel_detun+Xstate_split[1]-Xstate_detun
  delta_lup[,6]=rel_detun+Xstate_split[2]-Xstate_detun
  delta_lup[,7]=rel_detun+Xstate_split[3]-Xstate_detun
  delta_lup[,8]=rel_detun+Xstate_split[4]-Xstate_detun
  return(delta_lup)
}

laser_beam_detun_var <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun[1]+Xstate_split[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]+Xstate_split[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]+Xstate_split[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]+Xstate_split[4]-Xstate_detun
  delta_lup[,5]=rel_detun[1]+Xstate_split[1]-Xstate_detun
  delta_lup[,6]=rel_detun[2]+Xstate_split[2]-Xstate_detun
  delta_lup[,7]=rel_detun[3]+Xstate_split[3]-Xstate_detun
  delta_lup[,8]=rel_detun[4]+Xstate_split[4]-Xstate_detun
  return(delta_lup)
}

laser_beam_detun_var12 <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun[1]+Xstate_split[1]-Xstate_detun 
  delta_lup[,2]=rel_detun[2]+Xstate_split[1]-Xstate_detun 
  delta_lup[,3]=rel_detun[3]+Xstate_split[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]+Xstate_split[3]-Xstate_detun
  delta_lup[,5]=rel_detun[5]+Xstate_split[4]-Xstate_detun
  delta_lup[,6]=rel_detun[6]+Xstate_split[4]-Xstate_detun
  
  delta_lup[,7]=rel_detun[1]+Xstate_split[1]-Xstate_detun 
  delta_lup[,8]=rel_detun[2]+Xstate_split[1]-Xstate_detun 
  delta_lup[,9]=rel_detun[3]+Xstate_split[3]-Xstate_detun
  delta_lup[,10]=rel_detun[4]+Xstate_split[3]-Xstate_detun
  delta_lup[,11]=rel_detun[5]+Xstate_split[4]-Xstate_detun
  delta_lup[,12]=rel_detun[6]+Xstate_split[4]-Xstate_detun
  
  return(delta_lup)
}

## for BaH
laser_beam_detun12 <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun[1]-Xstate_detun 
  delta_lup[,2]=rel_detun[2]-Xstate_detun 
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  delta_lup[,5]=rel_detun[5]-Xstate_detun
  delta_lup[,6]=rel_detun[6]-Xstate_detun
  
  delta_lup[,7]=rel_detun[1]-Xstate_detun 
  delta_lup[,8]=rel_detun[2]-Xstate_detun 
  delta_lup[,9]=rel_detun[3]-Xstate_detun
  delta_lup[,10]=rel_detun[4]-Xstate_detun
  delta_lup[,11]=rel_detun[5]-Xstate_detun
  delta_lup[,12]=rel_detun[6]-Xstate_detun
  
  return(delta_lup)
}

laser_beam_detun8 <- function(Xstate_split,rel_detun){
  # there are a total of 12 ground MF states
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun[1]-Xstate_detun 
  delta_lup[,2]=rel_detun[2]-Xstate_detun 
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  
  delta_lup[,5]=rel_detun[1]-Xstate_detun 
  delta_lup[,6]=rel_detun[2]-Xstate_detun 
  delta_lup[,7]=rel_detun[3]-Xstate_detun
  delta_lup[,8]=rel_detun[4]-Xstate_detun
  
  return(delta_lup)
}

# 6 beams only for BaH
laser_beam_detun6 <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  delta_lup[,1]=rel_detun[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[1]-Xstate_detun
  delta_lup[,5]=rel_detun[2]-Xstate_detun
  delta_lup[,6]=rel_detun[3]-Xstate_detun
  #delta_lup[,7]=rel_detun+Xstate_split[3]-Xstate_detun
  #delta_lup[,8]=rel_detun+Xstate_split[4]-Xstate_detun
  return(delta_lup)
}

laser_beam_detun10 <- function(Xstate_split,rel_detun){
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12))
  # the Xstate detuning in included in rel_detun
  delta_lup[,1]=rel_detun[1]-Xstate_detun #c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun #(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  delta_lup[,5]=rel_detun[5]-Xstate_detun
  
  delta_lup[,6]=rel_detun[1]-Xstate_detun
  delta_lup[,7]=rel_detun[2]-Xstate_detun
  delta_lup[,8]=rel_detun[3]-Xstate_detun
  delta_lup[,9]=rel_detun[4]-Xstate_detun
  delta_lup[,10]=rel_detun[5]-Xstate_detun
  return(delta_lup)
}

laser_beam_detun12EOM <- function(Xstate_split,rel_detun){ # assume EOM driven sidebands
  # hyperfine and SR splittings for the X state
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # total number of beams
  # columns are for each laser beam and rows are for each hyperfine state
  delta_lup[,1]=rel_detun[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  delta_lup[,5]=rel_detun[5]-Xstate_detun
  delta_lup[,6]=rel_detun[6]-Xstate_detun
  # repeat now
  delta_lup[,7]=rel_detun[1]-Xstate_detun
  delta_lup[,8]=rel_detun[2]-Xstate_detun
  delta_lup[,9]=rel_detun[3]-Xstate_detun
  delta_lup[,10]=rel_detun[4]-Xstate_detun
  delta_lup[,11]=rel_detun[5]-Xstate_detun
  delta_lup[,12]=rel_detun[6]-Xstate_detun
  # delta_lup/(2*pi*1e6)
  return(delta_lup)
}

laser_beam_detun14EOM <- function(Xstate_split,rel_detun){ # assume EOM driven sidebands
  # hyperfine and SR splittings for the X state
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # total number of beams
  # columns are for each laser beam and rows are for each hyperfine state
  delta_lup[,1]=rel_detun[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  delta_lup[,5]=rel_detun[5]-Xstate_detun
  delta_lup[,6]=rel_detun[6]-Xstate_detun
  delta_lup[,7]=rel_detun[7]-Xstate_detun
  # repeat now
  delta_lup[,8]=rel_detun[1]-Xstate_detun
  delta_lup[,9]=rel_detun[2]-Xstate_detun
  delta_lup[,10]=rel_detun[3]-Xstate_detun
  delta_lup[,11]=rel_detun[4]-Xstate_detun
  delta_lup[,12]=rel_detun[5]-Xstate_detun
  delta_lup[,13]=rel_detun[6]-Xstate_detun
  delta_lup[,14]=rel_detun[7]-Xstate_detun
  # delta_lup/(2*pi*1e6)
  return(delta_lup)
}

laser_beam_detun16EOM <- function(Xstate_split,rel_detun){ # assume EOM driven sidebands
  # hyperfine and SR splittings for the X state
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # total number of beams
  # columns are for each laser beam and rows are for each hyperfine state
  delta_lup[,1]=rel_detun[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  delta_lup[,5]=rel_detun[5]-Xstate_detun
  delta_lup[,6]=rel_detun[6]-Xstate_detun
  delta_lup[,7]=rel_detun[7]-Xstate_detun
  delta_lup[,8]=rel_detun[8]-Xstate_detun
  # repeat now
  delta_lup[,9]=rel_detun[1]-Xstate_detun
  delta_lup[,10]=rel_detun[2]-Xstate_detun
  delta_lup[,11]=rel_detun[3]-Xstate_detun
  delta_lup[,12]=rel_detun[4]-Xstate_detun
  delta_lup[,13]=rel_detun[5]-Xstate_detun
  delta_lup[,14]=rel_detun[6]-Xstate_detun
  delta_lup[,15]=rel_detun[7]-Xstate_detun
  delta_lup[,16]=rel_detun[8]-Xstate_detun
  # delta_lup/(2*pi*1e6)
  return(delta_lup)
}

laser_beam_detun8AOM <- function(Xstate_split,rel_detun){ # assume AOM driven sidebands
  # hyperfine and SR splittings for the X state
  Xstate_detun=c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])#
  delta_lup=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # total number of beams
  # columns are for each laser beam and rows are for each hyperfine state
  delta_lup[,1]=rel_detun[1]-Xstate_detun#c(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])# -1.2*Gamma_n
  delta_lup[,2]=rel_detun[2]-Xstate_detun#(Xstate_split[1],Xstate_split[1],Xstate_split[1],Xstate_split[2],Xstate_split[3],Xstate_split[3],Xstate_split[3],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4],Xstate_split[4])-1.2*Gamma_n
  delta_lup[,3]=rel_detun[3]-Xstate_detun
  delta_lup[,4]=rel_detun[4]-Xstate_detun
  #delta_lup[,5]=rel_detun[5]-Xstate_detun
  #delta_lup[,6]=rel_detun[6]-Xstate_detun
  # repeat now
  delta_lup[,5]=rel_detun[1]-Xstate_detun
  delta_lup[,6]=rel_detun[2]-Xstate_detun
  delta_lup[,7]=rel_detun[3]-Xstate_detun
  delta_lup[,8]=rel_detun[4]-Xstate_detun
  #delta_lup[,11]=rel_detun[5]-Xstate_detun
  #delta_lup[,12]=rel_detun[6]-Xstate_detun
  delta_lup/(2*pi*1e6)
  return(delta_lup)
}

scattering_rate_matrix <- function(pos,vel,Afield,k_vec,pol_vec,gl,gu,delta_lup,sat_p){
  Rlup_init=cbind(rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]))
  Rlup_all=array(Rlup_init,dim=c(dim(Rlup_init)[1],dim(Rlup_init)[2],2)) # create array for storing scattering rates from all beams
  Bfield=Afield*c(pos[1],pos[2],-2*pos[3]) # magnetic field
  Bfield_mag=sqrt(Bfield[1]^2+Bfield[2]^2+Bfield[3]^2) # calculate B field magnitude
  for (l in 1:dim(Xstates)[1]){
    for (u in 1:dim(Astates)[1]){
      delta_omega=(gu[u]*Astates[u,3]-gl[l]*Xstates[l,3])*Bfield_mag*Bohr_magneton/hbar
      for (p in 1:dim(pol_vec)[1]){# loop over all laser beams
        # if ((pos[1]^2+pos[2]^2)<(rcut)^2){
        #   Ip=2*P_laser[p]/(pi*w^2)*exp(-2*(pos[1]^2+pos[2]^2)/w^2)
        # #   #sat_p=1
        # }else{
        #   Ip=0
        # #   #sat_p=0
        # }
        # sat_p=Ip/Isat2level
        #sat_p=1 # assume all transitions saturated
        #sat_param[i]=sat_p
        if(Bfield_mag!=0) {
          angle=acos((Bfield[1]*k_vec[p,1]+Bfield[2]*k_vec[p,2]+Bfield[3]*k_vec[p,3])/Bfield_mag) # angle between magnetic field and k-vector of laser beam
        }else{
          angle=0
        }
        # print(angle)
        pol_vec_local[1]=0.5*((1+cos(angle))*pol_vec[p,1]-sqrt(2)*sin(angle)*pol_vec[p,2]+(1-cos(angle))*pol_vec[p,3])
        pol_vec_local[2]=0.5*(sqrt(2)*sin(angle)*pol_vec[p,1]+2*cos(angle)*pol_vec[p,2]-sqrt(2)*sin(angle)*pol_vec[p,3])
        pol_vec_local[3]=0.5*((1-cos(angle))*pol_vec[p,1]+sqrt(2)*sin(angle)*pol_vec[p,2]+(1+cos(angle))*pol_vec[p,3])
        #print(pol_vec[p,])
        #print(pol_vec_local)
        # pol_vec_local=pol_vec[p,] # don't change the polarization
        if (abs(pol_vec_local[1])>1e-10 & (Xstates[l,3]+1==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=r_lu[l,u]*pol_vec_local[1]^2
          #print(c('sigma+',p,l,u))
        }else if(abs(pol_vec_local[2])>1e-10 & (Xstates[l,3]==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=r_lu[l,u]*pol_vec_local[2]^2
          #print(c('pi',p,l,u))
        }else if(abs(pol_vec_local[3])>1e-10 & (Xstates[l,3]-1==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=r_lu[l,u]*pol_vec_local[3]^2
          #print(c('sigma-',p,l,u))
        }else{
          f_lup=0
        }
        R_lup=0.5*Gamma_n*f_lup*sat_p/(1+sat_p+4*(delta_lup[l,p]-2*pi/lambda*(k_vec[p,1]*vel[1]+k_vec[p,2]*vel[2]+k_vec[p,3]*vel[3])-delta_omega)^2/Gamma_n^2)
        # R_lup=0.5*Gamma_n*f_lup*sat_p/(1+4*(delta_lup[l,p]-2*pi/lambda*(k_vec[p,1]*vel[1]+k_vec[p,2]*vel[2]+k_vec[p,3]*vel[3])-delta_omega)^2/Gamma_n^2)
        if (k_vec[p,3]==1) Rlup_all[l,u,1]=Rlup_all[l,u,1]+R_lup
        if (k_vec[p,3]==(-1)) Rlup_all[l,u,2]=Rlup_all[l,u,2]+R_lup
      }# sum over p ends here
    }# sum over u ends here
  }# sum over l ends here
  return(Rlup_all) # return the scattering array at this positions
} # function ends here

# no saturation parameter in the denominator this is the correct one to use for optica molasses
scattering_rate_matrix_2OM <- function(pos,vel,Afield,k_vec,pol_vec,gl,gu,delta_lup,Rabi_ij){
# scattering_rate_matrix_2OM <- function(pos,vel,Afield,k_vec,pol_vec,gl,gu,delta_lup,sat_p){
  Rlup_init=cbind(rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]))
  Rlup_all=array(Rlup_init,dim=c(dim(Rlup_init)[1],dim(Rlup_init)[2],2)) # create array for storing scattering rates from all beams
  Bfield=Afield*c(1,0,0)# uniform B field is in the x-direction *c(pos[1],pos[2],-2*pos[3]) # magnetic field
  Bfield_mag=sqrt(Bfield[1]^2+Bfield[2]^2+Bfield[3]^2) # calculate B field magnitude
  for (l in 1:dim(Xstates)[1]){
    for (u in 1:dim(Astates)[1]){
      delta_omega=(gu[u]*Astates[u,3]-gl[l]*Xstates[l,3])*Bfield_mag*Bohr_magneton/hbar
      for (p in 1:dim(pol_vec)[1]){# loop over all laser beams
        # if ((pos[1]^2+pos[2]^2)<(rcut)^2){
        #   Ip=2*P_laser[p]/(pi*w^2)*exp(-2*(pos[1]^2+pos[2]^2)/w^2)
        # #   #sat_p=1
        # }else{
        #   Ip=0
        # #   #sat_p=0
        # }
        # sat_p=Ip/Isat2level
        #sat_p=1 # assume all transitions saturated
        #sat_param[i]=sat_p
        if(Bfield_mag!=0) {
          angle=acos((Bfield[1]*k_vec[p,1]+Bfield[2]*k_vec[p,2]+Bfield[3]*k_vec[p,3])/Bfield_mag) # angle between magnetic field and k-vector of laser beam
        }else{
          angle=0
        }
        # print(angle)
        pol_vec_local[1]=0.5*((1+cos(angle))*pol_vec[p,1]-sqrt(2)*sin(angle)*pol_vec[p,2]+(1-cos(angle))*pol_vec[p,3])
        pol_vec_local[2]=0.5*(sqrt(2)*sin(angle)*pol_vec[p,1]+2*cos(angle)*pol_vec[p,2]-sqrt(2)*sin(angle)*pol_vec[p,3])
        pol_vec_local[3]=0.5*((1-cos(angle))*pol_vec[p,1]+sqrt(2)*sin(angle)*pol_vec[p,2]+(1+cos(angle))*pol_vec[p,3])
        #print(pol_vec[p,])
        #print(pol_vec_local)
        # pol_vec_local=pol_vec[p,] # don't change the polarization
        if (abs(pol_vec_local[1])>1e-10 & (Xstates[l,3]+1==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=pol_vec_local[1]^2
          #print(c('sigma+',p,l,u))
        }else if(abs(pol_vec_local[2])>1e-10 & (Xstates[l,3]==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=pol_vec_local[2]^2
          #print(c('pi',p,l,u))
        }else if(abs(pol_vec_local[3])>1e-10 & (Xstates[l,3]-1==Astates[u,3])){ # this will ensure that they are coupled
          f_lup=pol_vec_local[3]^2
          #print(c('sigma-',p,l,u))
        }else{
          f_lup=0
        }
        # here we use the correct expression in terms of Rabi frequency instead of saturation intensity
        R_lup=0.5*Gamma_n*f_lup*2*Rabi_ij[l,u]^2/Gamma_n^2/(1+4*(delta_lup[l,p]-2*pi/lambda*(k_vec[p,1]*vel[1]+k_vec[p,2]*vel[2]+k_vec[p,3]*vel[3])-delta_omega)^2/Gamma_n^2)
        # R_lup=0.5*Gamma_n*f_lup*kij[l,u]^2*sat_p/(1+4*(delta_lup[l,p]-2*pi/lambda*(k_vec[p,1]*vel[1]+k_vec[p,2]*vel[2]+k_vec[p,3]*vel[3])-delta_omega)^2/Gamma_n^2)
        if (k_vec[p,3]>0) Rlup_all[l,u,1]=Rlup_all[l,u,1]+R_lup
        if (k_vec[p,3]<0) Rlup_all[l,u,2]=Rlup_all[l,u,2]+R_lup
        
        #if (k_vec[p,3]==1) Rlup_all[l,u,1]=Rlup_all[l,u,1]+R_lup
        #if (k_vec[p,3]==(-1)) Rlup_all[l,u,2]=Rlup_all[l,u,2]+R_lup
      }# sum over p ends here
    }# sum over u ends here
  }# sum over l ends here
  return(Rlup_all) # return the scattering array at this positions
} # function ends here

pol_rotation <- function(angle,pol_vec){
  # rotate the polarization vector by a specific angle
  pol_vec_local = rep(NA,3)
  # refer to Hanley et al, J. Modern Optics 65:5-6, 667-676 (2017)
  pol_vec_local[1]=0.5*((1+cos(angle))*pol_vec[1]-sqrt(2)*sin(angle)*pol_vec[2]+(1-cos(angle))*pol_vec[3])
  pol_vec_local[2]=0.5*(sqrt(2)*sin(angle)*pol_vec[1]+2*cos(angle)*pol_vec[2]-sqrt(2)*sin(angle)*pol_vec[3])
  pol_vec_local[3]=0.5*((1-cos(angle))*pol_vec[1]+sqrt(2)*sin(angle)*pol_vec[2]+(1+cos(angle))*pol_vec[3])
  return(pol_vec_local)
}

# no saturation parameter in the denominator this is the correct one to use for optica molasses
## will try to vectorize this function now
scattering_rate_matrix_2OMvec <- function(pos,vel,k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij,Bfield){

  Rlup_all=array(0,dim=c(12,4,2)) # create array for storing scattering rates from all beams
  dotprod = apply(k_vec*matrix(rep(vel,dim(k_vec)[1]),byrow = T,nrow=dim(k_vec)[1]),MARGIN = 1,sum)*2*pi/lambda
  # matrix(dotprod,nrow = 12, ncol=4)
  # Bfield_mag = 1e-4 # 1 Gauss
  Bfield_mag=sqrt(Bfield[1]^2+Bfield[2]^2+Bfield[3]^2) # calculate B field magnitude
  for (p in 1:dim(pol_vec)[1]) { # loop over all the laser beams
    if(Bfield_mag!=0) {
      angle=acos((Bfield[1]*k_vec[p,1]+Bfield[2]*k_vec[p,2]+Bfield[3]*k_vec[p,3])/Bfield_mag) # angle between magnetic field and k-vector of laser beam
    }else{
      angle=0
    }
    pol_vec_local = pol_rotation(angle,pol_vec[p,])
    Rlup_init = Rabi_ij^2*(pi_mat*pol_vec_local[2]^2+sigmap_mat*pol_vec_local[1]^2+sigmam_mat*pol_vec_local[3]^2)/Gamma_n #cbind(rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]),rep(0,dim(Xstates)[1]))
    
    R_lup = Rlup_init/(1+4*(deltalup_all[,,p] - dotprod[p] - delta_omega*Bfield_mag)^2/Gamma_n^2)
      if (k_vec[p,3]>0) Rlup_all[,,1]=Rlup_all[,,1]+R_lup
      if (k_vec[p,3]<0) Rlup_all[,,2]=Rlup_all[,,2]+R_lup
  }
     
  return(Rlup_all) # return the scattering array at this positions
} # function ends here

derivs <- function(t,init_vals,params){
  # unbundle the lists of variables and parameters
  r = init_vals[1]
  v = init_vals[2]
  Nl = as.matrix(init_vals[3:14])
  Nu = as.matrix(init_vals[15:18])
  #Nl=Nl/(sum(Nl)+sum(Nu))
  #Nu=Nu/(sum(Nl)+sum(Nu))
  gam = init_vals[19]
  Rlup_all = params[[1]]
  r_lu = params[[2]]
  prefac = params[[3]]
  Nl_mat=matrix(rep(Nl,4),nrow=12)
  Nu_mat=matrix(rep(Nu,each=12),nrow=12)
  
  #Rlup_all=scattering_rate_matrix_2OMvec(c(0,0,r),c(0,0,v),k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij)
  
  dr <- v
  dv <- prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nl_mat,MARGIN = 2,sum))-
    prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nu_mat,MARGIN = 1,sum)) #+9.8 # added g at the end
  
  #Rtot=Rlup_all[,,1]+Rlup_all[,,2]
  
  dN_l <- apply((Rlup_all[,,1]+Rlup_all[,,2])*Nu_mat,MARGIN = 1,sum)-
    Nl*apply((Rlup_all[,,1]+Rlup_all[,,2]),MARGIN = 1,sum)+
    Gamma_n*apply(r_lu*Nu_mat,MARGIN = 1,sum)

  
  dN_u <- -Gamma_n*Nu+apply((Rlup_all[,,1]+Rlup_all[,,2])*Nl_mat,MARGIN = 2,sum)-
    Nu*apply((Rlup_all[,,1]+Rlup_all[,,2]),MARGIN = 2,sum)
  
  #sum(Gamma_n*apply(r_lu*Nu_mat,MARGIN = 1,sum))

  #print(sum(Nl)+sum(Nu))

  #sum(-Gamma_n*Nu)

  #print(sum(dN_l)+sum(dN_u))
  
  dgam <- Gamma_n*sum(Nu)
  
  return(list(c(dr,dv,dN_l,dN_u,dgam)))
}

derivs2 <- function(t,init_vals,params){
  # params includes r_lu and prefac only; Rlup_all is calculated here
  # assume position and velocity is changing now so the scattering rate is effected
  # unbundle the lists of variables and parameters
  r = init_vals[1]
  v = init_vals[2]
  Nl = as.matrix(init_vals[3:14])
  Nu = as.matrix(init_vals[15:18])
  #Nl=Nl/(sum(Nl)+sum(Nu))
  #Nu=Nu/(sum(Nl)+sum(Nu))
  gam = init_vals[19]
  #Rlup_all = params[[1]]
  r_lu = params[[1]]
  prefac = params[[2]]
  vz = params[[3]]
  Bfield = params[[4]]
  Nl_mat=matrix(rep(Nl,4),nrow=12)
  Nu_mat=matrix(rep(Nu,each=12),nrow=12)
  
  Rlup_all=scattering_rate_matrix_2OMvec(c(0,0,r),c(vz,0,v),k_vec,pol_vec,delta_omega,deltalup_all,Rabi_ij,Bfield)
  
  dr <- v
  dv <- prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nl_mat,MARGIN = 2,sum))-
    prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nu_mat,MARGIN = 1,sum)) #+9.8 # added g at the end
  
  #Rtot=Rlup_all[,,1]+Rlup_all[,,2]
  
  dN_l <- apply((Rlup_all[,,1]+Rlup_all[,,2])*Nu_mat,MARGIN = 1,sum)-
    Nl*apply((Rlup_all[,,1]+Rlup_all[,,2]),MARGIN = 1,sum)+
    Gamma_n*apply(r_lu*Nu_mat,MARGIN = 1,sum)
  
  
  dN_u <- -Gamma_n*Nu+apply((Rlup_all[,,1]+Rlup_all[,,2])*Nl_mat,MARGIN = 2,sum)-
    Nu*apply((Rlup_all[,,1]+Rlup_all[,,2]),MARGIN = 2,sum)
  
  #sum(Gamma_n*apply(r_lu*Nu_mat,MARGIN = 1,sum))
  
  #print(sum(Nl)+sum(Nu))
  
  #sum(-Gamma_n*Nu)
  
  #print(sum(dN_l)+sum(dN_u))
  
  dgam <- Gamma_n*sum(Nu)
  
  return(list(c(dr,dv,dN_l,dN_u,dgam)))
}

calc_accel <- function(Nl,Nu,Rlup_all,prefac){
  Nl_mat=matrix(rep(Nl,4),nrow=12)
  Nu_mat=matrix(rep(Nu,each=12),nrow=12)
  accelz <- prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nl_mat,MARGIN = 2,sum))-
    prefac*sum(apply((Rlup_all[,,1]-Rlup_all[,,2])*Nu_mat,MARGIN = 1,sum)) #+9.8 # added g at the end
  return(accelz)
}


