# initialize constants and parameters for BaH transverse laser cooling simulations
amu=1.66e-27 # [kg]
pi=3.14159
h_Planck=6.626e-34 # [J s]
hbar=1.0546e-34 # [J s]
c_light= 299792458 # [m/s] light speed
e_charge=1.6022e-19 # [C] electron charge
# Bohr_magneton=1.40e6 # [MHz/G]
Bohr_magneton=9.27e-24 # [J/T] in SI units
eps0=8.854e-12 # permittivity of free space constant

#lambda=635.1e-9 # [m]
#tau=58e-9 # [s]
lambda=695.1e-9 # [m] X-A
tau=33e-9 # [s]
FCF00=0.95
Gamma_n=1/tau*FCF00# [1/s]
mass=41*amu # CaH mass
Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3) # saturation intensity for a corresponding two level system
# Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3)



# lambda=1060e-9 # [m]
# tau=137e-9 # [s]
# FCF00=0.987
# Gamma_n=1/tau*FCF00# [1/s]
# mass=139*amu # BaH
# #Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3) # saturation intensity for a corresponding two level system
# Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3)

trans_dipole=sqrt(3*eps0*hbar*lambda^3/(8*tau*pi^2)) # calculate the value of the transition dipole moment

# define the states to use
#Xstates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0),c(1.5,1,-1),c(1.5,1,0),c(1.5,1,1),c(1.5,2,-2),c(1.5,2,-1),c(1.5,2,0),c(1.5,2,1),c(1.5,2,2))
#Astates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0))
# notice that for J=1/2 states we specify M_J projection and not M_F
Xstates=rbind(c(0.5,1,-0.5),c(0.5,1,-0.5),c(0.5,1,0.5),c(0.5,0,0.5),c(1.5,1,-1),c(1.5,1,0),c(1.5,1,1),c(1.5,2,-2),c(1.5,2,-1),c(1.5,2,0),c(1.5,2,1),c(1.5,2,2))

# again here specify the M_J projection and not M_F
Astates=rbind(c(0.5,1,-0.5),c(0.5,1,-0.5),c(0.5,1,0.5),c(0.5,0,0.5))

Is=pi*h_Planck*c_light*Gamma_n/(3*lambda^3)*0.1 # [mW/cm^2] saturation intensity for a corresponding 2-level system
# print(Is) # for a corresponding two-level system
ng=12 # number of ground state mf sublevel
ne=4 # number of excited state mf sublevels
seff=2*ng^2/(ng+ne) # scaling for the effective saturation parameter in a multi-level system
Is_eff=seff*Isat2level # assume rotationally closed transition
Is_eff_mWpercm2=seff*Is
scale_factor=1 # match the actual scattering rate measured in the experiment
gamma_n_eff=scale_factor*2*Gamma_n*ne/(ng+ne) # effective linewidth

# load all the states in |J,F,mf> format
pop_init=1/dim(Xstates)[1]
num_X_init=rep(pop_init,dim(Xstates)[1]) # initial population for the ground states; assume equally populated for simplicity

num_A_init=rep(0,dim(Astates)[1]) # initial population for the excited states

# Xstate_split=-2*pi*c(0,4,8642.8,40+8642.8)*1e6 # spin-rotation and hyperfine splittings
Xstate_split=-2*pi*c(0,53,1858+53,1857+53+102)*1e6 # CaH

# use the fitted g_J and g_F values 
#gl=c(-0.65,-0.65,-0.65,-0.65,0.87,0.87,0.7851,0.5058,0.4615,0.70,0.5545,0.50058) # BaH X state emprical g-factors
# gl=rep(NA,12)
# for (i in 1:length(gl)){
#   gl[i]=2*(Xstates[i,1]*(Xstates[i,1]+1)-1*2+0.5*1.5)*(Xstates[i,2]*(Xstates[i,2]+1)-0.5*1.5+Xstates[i,1]*(Xstates[i,1]+1))/(2*Xstates[i,1]*(Xstates[i,1]+1))/(2*Xstates[i,2]*(Xstates[i,2]+1))
# }
# gl[4]=0
gl=c(-0.755,-0.0877,-0.718,0.048,0.8800,0.8800,0.8567,0.50058,0.4803,0.522,0.522,0.50058) # CaH X state emprical g-factors
#gu=c(-0.088,-0.088,-0.088,0) # upper hyperfine sublevels g-factors
# gu=c(-0.51,-0.51,-0.51,-0.51) # BaH upper state g-factors
gu=c(-0.021,-0.021,-0.021,0) # upper hyperfine sublevels g-factors for CaH calculated from data in Chen...Steimle, Brown PRA 73, 012502 (2006)

#gu=rep(0,4) # zero out the excited state g-factors
#gu=c(-0.051,-0.051,-0.051,0) 0.1xg-factor
# gl=c(-0.47,-0.47,-0.47,0,0.97,0.97,0.97,0.5,0.5,0.5,0.5,0.5) # lower hyperfine sublevels g-factors
# gu=c(-0.088,-0.088,-0.088,0) # upper hyperfine sublevels g-factors

Zeeman_det = matrix(NA, nrow = 12, ncol = 4)

# delta_omega=(gu[u]*Astates[u,3]-gl[l]*Xstates[l,3])*Bfield_mag*Bohr_magneton/hbar
# set up the matrix with Zeeman detunings
for (i in 1:12){
  for (j in 1:4) {
    Zeeman_det[i,j] = gu[j]*Astates[j,3]-gl[i]*Xstates[i,3]
  }
}
delta_omega = Zeeman_det*Bohr_magneton/hbar # relative Zeeman shift per Tesla

kij=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # initialize the relative electric dipole matrix
###### feel the kij matrix with values taken from Leland Aldridge's thesis for SrOH (p. 105)
kij[1,1]=-sqrt(2)/3
kij[1,2]=-sqrt(2)/3
kij[1,4]=-sqrt(2)/3
##
kij[2,1]=sqrt(2)/3
kij[2,3]=-sqrt(2)/3
kij[2,4]=sqrt(2)/3
##
kij[3,2]=sqrt(2)/3
kij[3,3]=sqrt(2)/3
kij[3,4]=-sqrt(2)/3
##
kij[4,1]=sqrt(2)/3
kij[4,2]=sqrt(2)/3
kij[4,3]=sqrt(2)/3
##
kij[5,1]=1/6
kij[5,2]=1/6
kij[5,4]=-1/3
##
kij[6,1]=-1/6
kij[6,3]=1/6
kij[6,4]=1/3
##
kij[7,2]=-1/6
kij[7,3]=-1/6
kij[7,4]=-1/3
##
kij[8,1]=-1/sqrt(6)
##
kij[9,1]=sqrt(3)/6
kij[9,2]=-sqrt(3)/6
##
kij[10,1]=-1/6
kij[10,2]=1/3
kij[10,3]=-1/6
##
kij[11,2]=-sqrt(3)/6
kij[11,3]=sqrt(3)/6
##
kij[12,3]=-1/sqrt(6)

##### 
r_lu=kij^2 # spontaneous branching ratios; you can check that all columns sum to one

pi_mat0 = matrix(0, nrow = 12, ncol = 4) # matrix of pi transitions
pi_mat0[1,1] = 1
pi_mat0[2,2] = 1
pi_mat0[3,3] = 1
pi_mat0[4,4] = 1
pi_mat0[2,4] = 1
pi_mat0[4,2] = 1
pi_mat0[5,1] = 1
pi_mat0[6,2] = 1
pi_mat0[6,4] = 1
pi_mat0[7,3] = 1
pi_mat0[9,1] = 1
pi_mat0[10,2] = 1
pi_mat0[10,4] = 1
pi_mat0[11,3] = 1

pi_mat = pi_mat0*(r_lu>0)

sigmap_mat0 = matrix(0, nrow = 12, ncol = 4) # matrix of sigma+ transitions
sigmap_mat0[1,2] = 1
sigmap_mat0[1,4] = 1
sigmap_mat0[2,3] = 1
sigmap_mat0[4,3] = 1
sigmap_mat0[5,2] = 1
sigmap_mat0[5,4] = 1
sigmap_mat0[6,3] = 1
sigmap_mat0[8,1] = 1
sigmap_mat0[9,2] = 1
sigmap_mat0[9,4] = 1
sigmap_mat0[10,3] = 1

sigmap_mat = sigmap_mat0*(r_lu>0)

sigmam_mat0 = matrix(0, nrow = 12, ncol = 4) # matrix of sigma- transitions

sigmam_mat0[2,1] = 1
sigmam_mat0[3,2] = 1
sigmam_mat0[3,4] = 1
sigmam_mat0[4,1] = 1
sigmam_mat0[6,1] = 1
sigmam_mat0[7,2] = 1
sigmam_mat0[7,4] = 1
sigmam_mat0[10,1] = 1
sigmam_mat0[11,2] = 1
sigmam_mat0[11,4] = 1
sigmam_mat0[12,3] = 1

sigmam_mat = sigmam_mat0*(r_lu>0)
