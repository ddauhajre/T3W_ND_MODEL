###########################
run_ID = 'TEST'
###########################

##########################
# VERTICAL GRID PARAMETERS
'''
THESE DIMENSIONAL UNITS
ARE HERE JUST TO RELATE
BACK TO DIMENSIONAL MODEL

IN ALL DISCRETIZATIONS,
THE N.D. EQNS USE
N.D. VERTICAL AND TEMPORAL
SPACING

i.e., dz = z_{k+1}/H - z_{k}/H

z_grd_N --> leave True to determine vertical spacing by vertical levels
             as opposed to setting a dz

N_levs  --> number of vertical levels

theta_b --> bottom stretching parameter (degrees)

'''

##########################
z_grd_N=True
N_levs=300
theta_b = 10
##########################
# TIME STEPPING CONSTANTS
'''
AGAIN DIMENSIONAL UNITS
ARE JUST TO RELATE BACK TO
DIMENSIONAL SYSTEM AND FOR
THE CONSTRUCTUION OF THE 
VERTICAL MIXING TEMPORAL
SHAPE FUNCTION


dt        --> time step (seconds)
tend_days --> simulation length (days)
T_d       --> length of day in seconds to normalize to non-dimensonal time
'''
##########################
dt = 300 
tend_days = 10 
T_d = 86400.

####################################
# SET N.D. PARAMETERS
'''
N.D. PARAMETERS DERIVED
AND DESCRIBED IN WRITEUP
T3W_REDUCED_ANA.PDF

DEFINED W/ DIMENSIONAL PARAMETERS AS:

EK0   = K0 / (f_0 * H^2)
Gamma = (T * Delta_K) / (2 * H^2)
Omega = T * f_0
k     = 2K0 / Delta_K 


Ek_sign_conv --> sets sign of dubar / dz = 1 or -1
                 vertical boundary condition
               'neg' --> dubar / dz = -1 at z=0
                'pos' --> dubar / dz = 1 at z = 0

'''
###################################
Ek0   = 0.42 
Gamma = 1
Omega = 1
k     = 1.01 


Ek_sign_conv = 'neg'
##########################
#TEMPORAL SHAPE FUNCTION

'''
K_shape --> 'comp' for composite sigmoid and sinusoidal functional shape of vertical diffusivity
             K(t) = 0.5 * ( alpha * P(t) + beta * cos(t))
                    where P(t) is a sigmoid shape of a diurnal cycle (see T3W_ND_funcs.py)

alpha, beta --> parameters for K(t), 0<= alpha,beta <= 2
hour_shift --> sets symmetry of P(t) diurnal cycle (=4, equal day and night periods)
const_days --> number of diurnal periods to leave K(t) constant after initial condition
'''

#########################
K_shape = 'comp'
alpha =0
beta = 2
hour_shift = 4
const_days = 2
######################################
# VERTICAL BOUNDARY CONDITION CHOICE
'''
'free' --> free-slip at surface and bottom
'no_slip_bot' ---> u,v = 0 at z=-h
'''
######################################
BC_choice = 'free'


######################################
#TIME STEPPING SCHEME
'''
AB2 --> ADAMS BASHFORTH
LF  ---> LEAP FROG
'''
######################################
tstep_scheme = 'AB3'







