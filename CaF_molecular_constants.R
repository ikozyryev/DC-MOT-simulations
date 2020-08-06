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

lambda=606e-9 # [m]
tau=19e-9 # [s]
Gamma_n=1/tau# [1/s]
mass=56*amu # BaH
#Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3) # saturation intensity for a corresponding two level system
Isat2level=pi*h_Planck*c_light*Gamma_n/(3*lambda^3)

trans_dipole=sqrt(3*eps0*hbar*lambda^3/(8*tau*pi^2)) # calculate the value of the transition dipole moment

# define the states to use
Xstates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0),c(1.5,1,-1),c(1.5,1,0),c(1.5,1,1),c(1.5,2,-2),c(1.5,2,-1),c(1.5,2,0),c(1.5,2,1),c(1.5,2,2))
Astates=rbind(c(0.5,1,-1),c(0.5,1,0),c(0.5,1,1),c(0.5,0,0))

Is=pi*h_Planck*c_light*Gamma_n/(3*lambda^3)*0.1 # [mW/cm^2] saturation intensity for a corresponding 2-level system
# print(Is) # for a corresponding two-level system
ng=12 # number of ground state mf sublevel
ne=4 # number of excited state mf sublevels
seff=2*ng^2/(ng+ne) # scaling for the effective saturation parameter in a multi-level system
Is_eff=seff*Isat2level # assume rotationally closed transition
Is_eff_mWpercm2=seff*Is
scale_factor=1 # match the actual scattering rate measured in the experiment
FCF00=0.987
gamma_n_eff=scale_factor*2*Gamma_n*ne/(ng+ne) # effective linewidth

# load all the states in |J,F,mf> format
pop_init=1/dim(Xstates)[1]
num_X_init=rep(pop_init,dim(Xstates)[1]) # initial population for the ground states; assume equally populated for simplicity

num_A_init=rep(0,dim(Astates)[1]) # initial population for the excited states

# Xstate_split=-2*pi*c(0,4,8642.8,40+8642.8)*1e6 # spin-rotation and hyperfine splittings
# Xstate_split=-2*pi*c(0,55,129,170)*1e6 # relative detunings of SrF hyperfine ground state sublevels
Xstate_split=-2*pi*c(0,76,123,148)*1e6 # CaF

gl=c(-0.295,-0.295,-0.295,0,0.795,0.795,0.795,0.5,0.5,0.5,0.5,0.5) # lower hyperfine sublevels g-factors taken from Tarbutt and Steimle PRA 92, 053401 (2015)

gu=c(-0.0211,-0.0211,-0.0211,0) # upper hyperfine sublevels g-factors for CaF taken from Tarbutt NJP (2015)


#gl=c(-0.65,-0.65,-0.65,0,0.97,0.97,0.97,0.56,0.56,0.56,0.56,0.56) # BaH X state g-factors
#gu=c(-0.51,-0.51,-0.51,0) # BaH upper state g-factors
#gl=c(-0.47,-0.47,-0.47,0,0.97,0.97,0.97,0.5,0.5,0.5,0.5,0.5) # lower hyperfine sublevels g-factors
#gu=c(-0.088,-0.088,-0.088,0) # upper hyperfine sublevels g-factors

#gu=c(-0.051,-0.051,-0.051,0) 0.1xg-factor
# gl=c(-0.47,-0.47,-0.47,0,0.97,0.97,0.97,0.5,0.5,0.5,0.5,0.5) # lower hyperfine sublevels g-factors
# gu=c(-0.088,-0.088,-0.088,0) # upper hyperfine sublevels g-factors

# for SrF
#gl=c(-0.47,-0.47,-0.47,0,0.97,0.97,0.97,0.5,0.5,0.5,0.5,0.5) # lower hyperfine sublevels g-factors
#gu=c(-0.088,-0.088,-0.088,0) # upper hyperfine sublevels g-factors


Zeeman_det = matrix(NA, nrow = 12, ncol = 4)

# delta_omega=(gu[u]*Astates[u,3]-gl[l]*Xstates[l,3])*Bfield_mag*Bohr_magneton/hbar
# set up the matrix with Zeeman detunings
for (i in 1:12){
  for (j in 1:4) {
    Zeeman_det[i,j] = gu[j]*Astates[j,3]-gl[i]*Xstates[i,3]
  }
}
delta_omega = Zeeman_det*Bohr_magneton/hbar # relative Zeeman shift per Tesla

r_lu=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # initialize the relative electric dipole matrix
# kij=cbind(rep(0,12),rep(0,12),rep(0,12),rep(0,12)) # initialize the relative electric dipole matrix
# ###### directly feel the r_lu matrix for CaF using values from Eunmi Chae thesis p. 37
# kij[1,1]=-sqrt(2)/3
# kij[1,2]=-sqrt(2)/3
# kij[1,4]=-sqrt(2)/3
# ##
# kij[2,1]=sqrt(2)/3
# kij[2,3]=-sqrt(2)/3
# kij[2,4]=sqrt(2)/3
# ##
# kij[3,2]=sqrt(2)/3
# kij[3,3]=sqrt(2)/3
# kij[3,4]=-sqrt(2)/3
# ##
# kij[4,1]=sqrt(2)/3
# kij[4,2]=sqrt(2)/3
# kij[4,3]=sqrt(2)/3
# ##
# kij[5,1]=1/6
# kij[5,2]=1/6
# kij[5,4]=-1/3
# ##
# kij[6,1]=-1/6
# kij[6,3]=1/6
# kij[6,4]=1/3
# ##
# kij[7,2]=-1/6
# kij[7,3]=-1/6
# kij[7,4]=-1/3
# ##
# kij[8,1]=-1/sqrt(6)
# ##
# kij[9,1]=sqrt(3)/6
# kij[9,2]=-sqrt(3)/6
# ##
# kij[10,1]=-1/6
# kij[10,2]=1/3
# kij[10,3]=-1/6
# ##
# kij[11,2]=-sqrt(3)/6
# kij[11,3]=sqrt(3)/6
# ##
# kij[12,3]=-1/sqrt(6)



r_lu[1,1]=0.0667298
r_lu[1,2]=0.0667298
r_lu[1,4]=0.331582

r_lu[2,1]=0.0667298
r_lu[2,3]=0.0667298
r_lu[2,4]=0.331582

r_lu[3,2]=0.0667298
r_lu[3,3]=0.0667298
r_lu[3,4]=0.331582

r_lu[4,1]=2/9
r_lu[4,2]=2/9
r_lu[4,3]=2/9

r_lu[5,1]=0.18327
r_lu[5,2]=0.18327
r_lu[5,4]=0.00175166

r_lu[6,1]=0.18327
r_lu[6,3]=0.18327
r_lu[6,4]=0.00175166

r_lu[7,2]=0.18327
r_lu[7,3]=0.18327
r_lu[7,4]=0.00175166

r_lu[8,1]=1/6

r_lu[9,1]=1/12
r_lu[9,2]=1/12

r_lu[10,1]=1/36
r_lu[10,2]=1/9
r_lu[10,3]=1/36

r_lu[11,2]=1/12
r_lu[11,3]=1/12

r_lu[12,3]=1/6
#spontaneous branching ratios; you can check that all columns sum to one
##### 
kij=sqrt(r_lu) # 

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
