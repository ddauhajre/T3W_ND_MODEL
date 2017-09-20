######################################
__title__          = "T3W_ND_funcs.py"
__author__         = "Daniel Dauhajre"
__date__           = "September 2017"
__email__          = "ddauhajre@atmos.ucla.edu"
__python_version__ = "2.7.9"

'''
SUITE OF FUNCTIONS USED FOR EVOLVING
THE NON-DIMENSIONAL MODEL

CONTAINS SETUP/GRID CREATION AS WELL
AS TIME-STEPPING FUNCTIONS TO SOLVE VARIOUS
N.D. EQNS

'''
######################################

#########################################
#IMPORT MODULES
import os
import sys
import numpy as np
import scipy as sp
import pickle as pickle
from scipy.sparse import *
from scipy.sparse.linalg import spsolve
from scipy import integrate
from copy import copy
#########################################




      ############################################################
      #             VARIABLE DICTIONARY SETUP
'''
                    ALL DATA RELEVANT TO MODEL IS STORED IN
                    DICTIONARY THAT WILL BE SAVED AS A PICKLE FILE
                    
                    NEED FUNCTIONALITY TO SET UP THIS DICTIONARY
'''
      ############################################################

def add_keys_var_dict(keys_add,keys_type,var_dict={}):
    '''
    var_dict --> dictionary to add keys to (default is empty {} or
                 but can be an actual pre-existing dictionary)

    keys_add ---> list of string of keys to set up dictionary with

    keys_typ e--> list of data types corresponding to each key (i.e., empty 
                       list or string or numpy array
    '''
    for ke in range(len(keys_add)):
        var_dict[keys_add[ke]] = keys_type[ke]

    return var_dict
    ##################################################

def create_empty_list_key(keys_add,var_dict={}):
    for ke in range(len(keys_add)):
        var_dict[keys_add[ke]] = []
    return var_dict

      
def setup_sol_vars_for_dict(nt,N):
    """
    CREATE LIST OF KEYS / TYPES FOR SOLUTION
    DICTIONARY

    # NECESSARY DIMENSIONS FOR SOLN ARRAYS
    nt --> number of time points
    N ---> vertical levels
    """
    keys = ['ubar', 'vbar', 'u_p', 'v_p','ustar','vstar','u_pp', 'v_pp']
    key_types = [np.zeros(N), np.zeros(N), np.zeros([nt,N]), np.zeros([nt,N]),np.zeros([nt,N]),np.zeros([nt,N]),np.zeros([nt,N]),np.zeros([nt,N])]
    return [keys, key_types]


      ############################################################
      #             SAVE SOLUTION DICTIONARY WITH PICKLE
      ############################################################

def change_dir_py(dir_name):
     if not os.path.exists(dir_name):
        print 'Creating directory: ' + dir_name
        os.makedirs(dir_name)
     print 'Moving into directory: ' + dir_name
     os.chdir(dir_name)
     ##########################################

def save_to_pickle(var_dict,out_name):
    pickle.dump(var_dict,open(out_name +".p", "wb"))
    print 'Saved output as: ' + out_name + '.p'
    #####################################################

      ############################################################
      #             TIME-SPACE GRID SET-UP
      ############################################################




def make_zt_grid_N_levs(N,dt,tend_days,theta_b=0.):
    """
    CREATE VERTICAL GRID OF rho-
    and w-levels based on number of 
    total rho-levels
    """
    print ''
    print 'Creating (z,t)-grid with following dimensions'
    print '#######################'
    print 'N= ' + str(N) +  ' levels'
    print 'dt = ' + str(dt) + 's'
    print 'tend_days = ' + str(tend_days) + ' days'
    print '#######################'

    z_r_temp = np.zeros([N])
    z_w_temp = np.zeros([N+1])
    for k in range(1,N):
        z_r_temp[k-1] = (k - N - 0.5) / N
        z_w_temp[k] = (k - float(N)) / N
    
    z_r_temp[-1] = -0.5 / N 
    z_w_temp[0] = -1.
    

    ############################################
    # INCREASE RESOLUTION AT BOTTOM WITH THETA_B
    ############################################
    if theta_b >0:
       z_w = ( np.exp(theta_b * z_w_temp) - 1) / ( 1 - np.exp(-theta_b))
       z_r = ( np.exp(theta_b * z_r_temp) - 1) / ( 1 - np.exp(-theta_b))
    else:
       z_w = z_w_temp
       z_r = z_r_temp 

    #########################
    # TIME VECTOR
    #########################
    tend_sec = tend_days * 86400
    tvec_sec = np.arange(0,tend_sec,dt)


    ##########################   
    #STORE IN GRID DICTIONARY
    ##########################
    grid_dict = {}
    grid_dict['z_r']       = z_r
    grid_dict['z_w']       = z_w
    grid_dict['dt']        = dt
    grid_dict['tend_days'] = tend_days
    grid_dict['tvec_sec']  = tvec_sec


    
    return grid_dict  
    #####################################################




def make_zt_grid(Lz,dz,dt,tend_days):
    """
    CREATE STAGGERED VERTICAL GRID OF
    rho- and w-levels based on input parameters
    
    AS WELL AS TIME-GRID
    """
    print ''
    print 'Creating (z,t)-grid with following dimensions'
    print '#######################'
    print 'Lz = ' + str(Lz) +  'm'
    print 'dz = ' + str(dz) +  'm'
    print 'dt = ' + str(dt) + 's'
    print 'tend_days = ' + str(tend_days) + ' days'
    print '#######################'

    z_r = np.arange(-Lz + (dz/2.), dz/2.,dz)
    z_w = np.arange(-Lz,0+dz,dz)
   
    
    #########################
    # TIME VECTOR
    #########################
    tend_sec = tend_days * 86400
    tvec_sec = np.arange(0,tend_sec,dt)


    ##########################   
    #STORE IN GRID DICTIONARY
    ##########################
    grid_dict = {}
    grid_dict['z_r']       = z_r
    grid_dict['z_w']       = z_w
    grid_dict['dt']        = dt
    grid_dict['tend_days'] = tend_days
    grid_dict['tvec_sec']  = tvec_sec


    
    return grid_dict  
    #####################################################


    ############################################################
    #            CALCULATING GRADIENT FUNCTIONS 
    ############################################################

def calc_vert_shears(var,grid_dict):
    """
    CALCULATE VERTICAL SHEARS
    OF A FIELD

    CALCULATE FIRST AND SECOND 
    DERIVATIVES
    """     
    N = len(var)
    var_z = np.zeros([N-1])
    var_zz = np.zeros([N-2])
    
    H = abs(grid_dict['z_w'][0])
    z_r = grid_dict['z_r'] / H
    z_w = grid_dict['z_w'] / H
   
    for k in range(N-1):
        var_z[k] = (var[k+1] - var[k]) / (z_r[k+1] - z_r[k])
    
    for k in range(N-2):
        var_zz[k] = (var_z[k+1] - var_z[k]) / (z_w[k+1+1] - z_w[k+1]) 

    return var_z, var_zz 
    ##############################################################

def calc_u_v_zz_flux(u,v,grid_dict,BC_opt,bar_BCs=True, BC_signs='neg'):
    """
    Calculate 2nd derivative of a component of flow
    in flux form using staggered vertical grid

    BC_opt --> same BC choices as in rest of module

    bar_BCs --> if calculting 2nd derivatives of 
                ubar,vbar need to enforce correct BCs

    BC_signs
      'neg' u_z = -1
       'pos' u_z = 1
    """
    N = len(u)
    u_zz = np.zeros(N)
    v_zz = np.zeros(N)
    H = abs(grid_dict['z_w'][0])
    z_w = grid_dict['z_w'] / H
    z_r = grid_dict['z_r'] / H
    dz = z_r[1:] - z_r[:-1]
    
    #######################
    #   INTERIOR
    #######################
    for k in range(1,N-1):
        u_zz[k] = ( (u[k+1] - u[k]) / dz[k])  - ( (u[k] - u[k-1]) / dz[k-1])
        v_zz[k] = ( (v[k+1] - v[k]) / dz[k])  - ( (v[k] - v[k-1]) / dz[k-1])
    
    #######################
    #   BOUNDARIES
    #######################
 
    #SURFACE
    if bar_BCs:
       if BC_signs == 'neg':
          u_zz[-1] = -1 - ( (u[-1] - u[-2]) / dz[-1])
       if BC_signs == 'pos':
          u_zz[-1] = 1 - ( (u[-1] - u[-2]) / dz[-1])
       v_zz[-1] = -( (v[-1] - v[-2]) / dz[-1])
    else:
        u_zz[-1] = - ( (u[-1] - u[-2]) / dz[-1])
        v_zz[-1] = - ( (v[-1] - v[-2]) / dz[-1])
     
    #BOTTOM     
    if bar_BCs:
       if BC_opt == 'free':
          if BC_signs == 'neg':
             u_zz[0] = ( (u[1] - u[0]) / dz[0]) +1
          if BC_signs == 'pos':
             u_zz[0] = ( (u[1] - u[0]) / dz[0]) - 1
          v_zz[0] = ( (v[1] - v[0]) / dz[0]) 
       if BC_opt == 'no_slip_bot':
          u_zz[0] = ( (u[1] - u[0]) / dz[0]) + (2*u[0]/dz[0])
          v_zz[0] = ( (v[1] - v[0]) / dz[0]) + (2*v[0]/dz[0])
            
    else:
        if BC_opt == 'free':
           u_zz[0] = ( (u[1] - u[0]) / dz[0]) 
           v_zz[0] = ( (v[1] - v[0]) / dz[0]) 
        if BC_opt == 'no_slip_bot':
           u_zz[0] = ( (u[1] - u[0]) / dz[0]) + (2*u[0]/dz[0])
           v_zz[0] = ( (v[1] - v[0]) / dz[0]) + (2*v[0]/dz[0])
        
    return u_zz,v_zz
    ###############################################################





      ############################################################
      #             CREATING OF FORCING TEMPORAL SHAPE FUNCTIONS
      ############################################################

def make_sigmoid_tseries(tvec_sec, shift_time, twidth):
    nt_curve = len(tvec_sec)
    F = np.zeros([nt_curve])
    for t in range(nt_curve):
        F[t] = np.tanh( (shift_time - tvec_sec[t]) / (0.5 * twidth))
    return F

def make_K_day_composite(dt, tvec_sec, hour_shift,alpha,beta):
    '''
    CREATE A COMPOSITED VERSION OF A DIURNAL CYCLE
    FOR A SINGLE DAY USING THE FOLLOWING FORMULATION

    K(t) = 0.5 * (alpha * P(t,t_w) + beta * cos(t/T))

    Where P(t,t_w) is a simgoid mirroring time-series and alpha 
    and beta are tuning factors for the shape and speed of the composite
    time series

    alpha = 0 is the default case which is also the slowest transition
    '''
    ###################################
    # CREATE P(t,t_w)
    ####################################     
    
    # SET UP PARAMETERS FOR P(t,t_w)
    hrs_day = hour_shift / (4./24.)
    len_day_tsteps_hr_shift = len(np.arange(0,hrs_day*60*60,dt))
    if len_day_tsteps_hr_shift %2!=0:
       len_day_tsteps_hr_shift+=1

    len_day_sec = hrs_day * 60 * 60
    delta = 2
    shift_time = (hrs_day * 60 * 60) / (2 * delta)
    T  = len(np.arange(0,24*60*60,dt))
    hour_diff = 4 - hour_shift #4 hour maxmimum shift
    T_d = hour_diff * 60 * 60 / dt
     
    #CREATE HALF-DAY WITH SIGMOID FUNCTION
    twidth = hour_shift * 60 * 60
    S_1 = make_sigmoid_tseries(tvec_sec[0:len_day_tsteps_hr_shift/delta +1], shift_time,twidth)
    #NORMALIZE S SO THAT INITIAL VALUE = 1
    S = S_1 / S_1[0]     

    P_day = np.zeros([T])
    P_day[0:T_d] = S[0]
    T_s = len(S)
    P_day[T_d:T_d + T_s] = S
     
    if T_d + T_s<=T/2:
       P_day[T_d + T_s:T - (T_d + T_s)] = S[-1]
       P_day[T - (T_d + T_s):T-T_d] = S[::-1]
       P_day[T - T_d::] = S[0]
    else:
       P_day[T_s + T_d::] = S[::-1][0:T-(T_s + T_d)]
       
    ##########################################
    # CREATE cos(t/T)
    ##########################################
    K_cos = np.zeros([T])
    for t in range(T):
        K_cos[t] = np.cos(t / float(T) * (2 * np.pi))
   
    ################################################
    # CREATE COMPOSITE
    ################################################
    K_comp = 0.5 * (alpha * P_day + K_cos * beta)

    #RETURN NORMALIZED VERSION THAT HAS A TIME-MEAN = 0
    return K_comp - np.mean(K_comp)
    ########################################################
    

def make_K_ND_comp(dt,tvec_sec,tend_days,hour_shift,const_days,alpha=0,beta=2):
    """
    CREATE DIURNAL CYCLE OF 
    NON-DIMENSIONAL MIXING
    BASED ON COMPOSITE FORMULATION 
    OF SIGMOID BASED AND COSINE BASED
    FUNCTIONS 

    TUNING OF THIS IS SET BY hour_shift, alpha, beta


    const_days is number of days to hold N.D mixing constant
    for a spinup or similar type of period
    """

    ############################################
    # CREATE DIURNAL CYCLE FOR 1-DAY
    ############################################
    K_day = make_K_day_composite(dt, tvec_sec, hour_shift,alpha,beta)

    #####################################################
    # CREATE FULL VERSION W/ CONSTANT DAYS IF NECESSARY
    ####################################################
    K_full = np.zeros([len(tvec_sec)])
    #FILL W/ CONSTANT
    T = len(np.arange(0,24*60*60,dt))
    len_const_days = T * const_days
    K_full[0:len_const_days] = K_day[0]
    
    #FILL DIURNALLY CYCLING DAYS
    tind = len_const_days 
    for d in range(tend_days - const_days):
        K_full[tind:tind+T] = K_day
        tind = tind + T
    return K_full
    ###########################################################







      ############################################################
      #          TIME-STEPPING FUNCTION
      ############################################################
def tstep_system(var_dict):
    """
    WITH ubar, vbar calculated

    TIME-STEP THE SYSTEM TO CALCULATE u_p, v_p at
    every time-point
    """
    ##########################
    # DECLARE SOLUTION ARRAYS
    #########################
    u_out = copy(var_dict['u_p'])
    v_out = copy(var_dict['v_p'])
    ustar_out = copy(var_dict['ustar'])
    vstar_out = copy(var_dict['vstar'])


    nt = u_out.shape[0]    

    ################################
    # CALCULATE 2ND DERIVATIVES
    # OF ubar, vbar for R.H.S FORCING
    ################################
    ubar_zz,vbar_zz = calc_u_v_zz_flux(var_dict['ubar'],var_dict['vbar'],var_dict['grid'],var_dict['BC_choice'],bar_BCs=True, BC_signs = var_dict['Ek_sign_conv'])
   
    #################################
    # GET COMPUTATIONAL INITIAL 
    # CONDITION AT t=1 w/ implicit
    # time-stepping
    #################################
   
    ###############################################
    #INITIAL CONDITION FOR u_p, v_p based on Ek_max
    ################################################
    Ek_max = var_dict['Ek0'] +  (var_dict['Gamma'] / var_dict['Omega'])
    u_max,v_max = solve_ubar_vbar(var_dict['grid'],Ek_max, var_dict['BC_choice'], sign_conv = var_dict['Ek_sign_conv'])
    u_out[0,:] = u_max - var_dict['ubar']
    v_out[0,:] = v_max - var_dict['vbar']
    

    #CALCULAT R.H.S. TERMS  
    Urhs, Vrhs = calc_Urhs_Vrhs(var_dict['grid'], u_out[0,:], v_out[0,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][0],var_dict['T_d'])
    Gamma_0 = var_dict['Gamma'] * (var_dict['K_ND'][0]  + var_dict['k'])
    print 'Calculating u_p, v_p at t=1'
    u_out[1,:],v_out[1,:] = tstep_up_vp(var_dict['grid'],u_out[0,:],v_out[0,:],var_dict['grid']['dt'],var_dict['T_d'],Gamma_0,Urhs,Vrhs,var_dict['BC_choice'],scheme='implicit')    
    print '     ustar,vstar at t=' + str(0)
    ustar_out[0,:],vstar_out[0,:]=solve_ustar_vstar(var_dict['grid'],var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][0],var_dict['k'],ubar_zz,vbar_zz,var_dict['BC_choice'])

    ###############################
    # CALCUATE AT n=2 w/ LF scheme
    ##############################
    Urhs, Vrhs = calc_Urhs_Vrhs(var_dict['grid'], u_out[1,:], v_out[1,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][1],var_dict['T_d'])
    Gamma_0 = var_dict['Gamma'] * (var_dict['K_ND'][1]  + var_dict['k'])
    print 'Calculating u_p, v_p at t=2'
    u_out[2,:],v_out[2,:] = tstep_up_vp(var_dict['grid'],u_out[0,:],v_out[0,:],var_dict['grid']['dt'],var_dict['T_d'],Gamma_0,Urhs,Vrhs,var_dict['BC_choice'],scheme='LF')    
    print '     ustar,vstar at t=' + str(1)
    ustar_out[1,:],vstar_out[1,:]=solve_ustar_vstar(var_dict['grid'],var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][1],var_dict['k'],ubar_zz,vbar_zz,var_dict['BC_choice'])

    if var_dict['tstep_scheme'] == 'AB2' or var_dict['tstep_scheme'] == 'AB3':
       ###############################################
       # SOLVER INTERIOR POINTS W/ AB2 or AB3 METHOD
       ##############################################
      for n in range(2,nt-1):
           if var_dict['tstep_scheme'] == 'AB2' or n==2:
              #URHS, VHRS --> 3/2 urhs^{n} - 1/2 urhs^{n}
              Urhs_n, Vrhs_n = calc_Urhs_Vrhs(var_dict['grid'],u_out[n,:],v_out[n,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['T_d'])
              Urhs_n_1, Vrhs_n_1 = calc_Urhs_Vrhs(var_dict['grid'],u_out[n-1,:],v_out[n-1,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n-1],var_dict['T_d'])
              Urhs = (3/2. * Urhs_n) - 0.5 * Urhs_n_1
              Vrhs = (3/2. * Vrhs_n) - 0.5 * Vrhs_n_1
           if var_dict['tstep_scheme'] == 'AB3':
              #URHS, VHRS --> 23/12 urhs^{n}  4/3 urhs^{n-1} + 5/12 urhs^{n-2}
              Urhs_n, Vrhs_n = calc_Urhs_Vrhs(var_dict['grid'],u_out[n,:],v_out[n,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['T_d'])
              Urhs_n_1, Vrhs_n_1 = calc_Urhs_Vrhs(var_dict['grid'],u_out[n-1,:],v_out[n-1,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n-1],var_dict['T_d'])
              Urhs_n_2, Vrhs_n_2 = calc_Urhs_Vrhs(var_dict['grid'],u_out[n-2,:],v_out[n-2,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n-1],var_dict['T_d'])
              Urhs = (23./12. * Urhs_n) - (4./3. * Urhs_n_1)  + (5./12. * Urhs_n_2)
              Vrhs = (23./12. * Vrhs_n) - (4./3. * Vrhs_n_1)  + (5./12. * Vrhs_n_2)

           
           Gamma_n = var_dict['Gamma'] * (var_dict['K_ND'][n] + var_dict['k'])
           print 'Calculating u_p,v_p at t=' + str(n+1)
           u_out[n+1,:],v_out[n+1,:] = tstep_up_vp(var_dict['grid'],u_out[n,:],v_out[n,:],var_dict['grid']['dt'],var_dict['T_d'],Gamma_n,Urhs,Vrhs,var_dict['BC_choice'],scheme=var_dict['tstep_scheme'])
           print '     ustar,vstar at t=' + str(n)
           ustar_out[n,:],vstar_out[n,:]=solve_ustar_vstar(var_dict['grid'],var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['k'],ubar_zz,vbar_zz,var_dict['BC_choice'])
           ##################################################################################    



    if var_dict['tstep_scheme'] == 'LF':
       #################################
       # SOLVE INTERIOR POINTS w/ LF method
       ################################
       for n in range(2,nt-1):
           Urhs,Vrhs = calc_Urhs_Vrhs(var_dict['grid'],u_out[n,:],v_out[n,:],ubar_zz,vbar_zz,var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['T_d'])
           Gamma_n = var_dict['Gamma'] * (var_dict['K_ND'][n] + var_dict['k'])
           print 'Calculating u_p,v_p at t=' + str(n+1)
           u_out[n+1,:],v_out[n+1,:] = tstep_up_vp(var_dict['grid'],u_out[n-1,:],v_out[n-1,:],var_dict['grid']['dt'],var_dict['T_d'],Gamma_n,Urhs,Vrhs,var_dict['BC_choice'],scheme='LF')
           print '     ustar,vstar at t=' + str(n)
           ustar_out[n,:],vstar_out[n,:]=solve_ustar_vstar(var_dict['grid'],var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['k'],ubar_zz,vbar_zz,var_dict['BC_choice'])
           ##################################################################################    


    n+=1
    ustar_out[n,:],vstar_out[n,:]=solve_ustar_vstar(var_dict['grid'],var_dict['Omega'],var_dict['Gamma'],var_dict['K_ND'][n],var_dict['k'],ubar_zz,vbar_zz,var_dict['BC_choice'])


 
    
    return u_out, v_out,ustar_out,vstar_out,ubar_zz,vbar_zz


      ############################################################
      #             FUNCTIONS TO SOLVE N.D. EQNS
      ############################################################

def calc_Urhs_Vrhs(grid_dict,u_p_n, v_p_n, ubar_zz, vbar_zz, Omega, Gamma, K_ND_n,T_d):
    """
    CALCULATE R.H.S TERMS TO BE USED IN
    tstep_up_vp

    u_p_n --> u_p at t=n
    v_p_n --> v_p at t=n
    K_ND_n --> K(t) ND shape function at t=n
    """
    H = abs(grid_dict['z_w'][0])
    Hz = (grid_dict['z_w'][1:]/H) - (grid_dict['z_w'][:-1]/H)
    Beta = Gamma * K_ND_n
    Urhs = ((-Omega * v_p_n)*Hz) - (Beta * ubar_zz) 
    Vrhs = ((Omega * u_p_n)*Hz) - (Beta * vbar_zz)
    return Urhs, Vrhs


def tstep_up_vp(grid_dict,u_pre,v_pre,dt,T_d,Gamma_n,Urhs,Vrhs,BC_opt,scheme='LF'):
    """
    SOLVE TRANSIENT, PERIODIC
    N.D. EQN

    Urhs, Vrhs obtained in separate function

    THIS FUNCTION ADVANCES u_p, v_p to the next time-step
    USING u_pre, v_pre to be either u^{n}, v^{n} 
                                         or 
                                    u^{n-1}, v^{n-1}
    
    DEPENDING ON WHETHER IF scheme =='LF'
                          or scheme == 'implicit'


    IMPLICIT --> u_pre == u^n
    EXPLICIT --> u_pre == u^{n-1}
 
 
    FOR IMPLICIT SCHEME --> DIFFUSIVE IMPLICIT METHOD
                           dt terms has factor of 1 
    FOR EXPLICIT SCHEME --> LEAP FROG METHOD AND
                            dt terms will have a factor of 2 multiplitication


    Gamma_n = Gamma * (K_ND(t=n) + k) 
    """

    ##########################
    # SET UP dt-coefficient
    # based on scheme
    #########################
    if scheme == 'AB2' or 'AB3':
       dt_coeff = dt / T_d
    if scheme == 'LF':
       dt_coeff = 2 * dt / T_d
    if scheme == 'implicit':
       dt_coeff = dt/T_d
    #################################
    # DIMENSIONS
    ################################## 
    H = abs(grid_dict['z_w'][0])
    #Hz = grid_dict['z_w'][1:] - grid_dict['z_w'][:-1]
    #dz = grid_dict['z_r'][1:] - grid_dict['z_r'][:-1]
    Hz = (grid_dict['z_w'][1:]/H) - (grid_dict['z_w'][:-1]/H)
    dz = (grid_dict['z_r'][1:]/H) - (grid_dict['z_r'][:-1]/H)
    
    N = len(grid_dict['z_r'])

    ############################
    # DECLARE SOLUTION ARRAYS
    ############################
    u_out = np.zeros(N)
    v_out = np.zeros(N)
    ################################
    # SET UP MATRICES
    ###############################
    ndim = 2 * N
    kv   = 2

    A = lil_matrix((ndim, ndim))
    R = np.zeros(ndim)


    #######################
    # BOTTOM
    ########################
    if BC_opt == 'free':
       D1 = (1*Hz[0] + (dt_coeff * Gamma_n / dz[0]))
       D2 = (dt_coeff * Gamma_n / dz[0])
    if BC_opt == 'no_slip_bot':
       D1 = (1*Hz[0] + (dt_coeff * Gamma_n / dz[0]) + (2*Gamma_n*dt_coeff/dz[0]))
    
    D2 = (dt_coeff * Gamma_n/dz[0])
   
    idx = 0   
    A[idx,idx+kv] = -D2
    A[idx,idx] = D1
   
    A[idx+1,idx+1+kv] = -D2
    A[idx+1,idx+1]    = D1

    ##################################
    # INTERIOR
    ##################################
    for k in range(1,N-1):
        D_k       = (1*Hz[k] + (dt_coeff * Gamma_n / dz[k-1]) + (dt_coeff * Gamma_n / dz[k]))
        D_k_minus = dt_coeff * Gamma_n / dz[k-1]
        D_k_plus  = dt_coeff * Gamma_n / dz[k]

        idx = 2 * k
 
        A[idx,idx+kv] = -D_k_plus
        A[idx,idx]    = D_k
        A[idx,idx-kv] = -D_k_minus

        A[idx+1,idx+1+kv] = -D_k_plus
        A[idx+1,idx+1]    = D_k
        A[idx+1,idx+1-kv]  = -D_k_minus

    ##########################
    #   TOP
    ##########################
    D_N = (1*Hz[-1] + (Gamma_n * dt_coeff /dz[-1]))
    D_N_minus = Gamma_n * dt_coeff /dz[-1] 

    idx = (2*N) - kv
    
    A[idx,idx-kv] = -D_N_minus
    A[idx,idx]    = D_N
    
    A[idx+1,idx+1-kv] = -D_N_minus
    A[idx+1,idx+1]    = D_N
    ##########################################################
    ##########################################################

    ##############################
    # FILL R.H.S MATRIX
    ##############################     
        
    for k in range(N):
        idx = 2*k
        R[idx] = u_pre[k]*Hz[k] - dt_coeff * Urhs[k]
        R[idx+1] = v_pre[k]*Hz[k] - dt_coeff * Vrhs[k]
    
    ##########################################################
    ##########################################################

    ##############################
    # SOLVE AND REORDER RESULTS
    ##############################     
    A = A.tocsr()
    X = spsolve(A,R)

    for k in range(N):
        idx = 2*k
        u_out[k] = X[idx]
        v_out[k] = X[idx+1]

    return u_out, v_out
    ################################################################# 



def solve_ubar_vbar(grid_dict, Ek0, BC_opt, sign_conv = 'neg'):
    """
    SOLVE FOLLOWING N.D. SYSTEM
    u = Ek0 * v_zz
    -v = Ek0 * u_zz

    where u,v --> ubar, vbar

    BC_opt 
     'free'
         u_z = -1 at  z = 0, -1
         v_z = 0 at   z = 0, -1 
     
     'no_slip_bot'
         u_z = -1 at z = 0
         v_z = 0 at  z = 0
         u,v = 0 at  z = -1
    sign_conv  --> controls positive/negative surface and bottom u,v (e.g.,
                   choose to have positive v at surface with 'neg' or negative
                   v at surface with 'pos')
       'neg'
          u_z = -1 at boundary where necessary
        'pos'
           u_z = 1 at boundary where necessary
    """
    #VERTICAL SPACING
    H = abs((grid_dict['z_w'][0]))
    dz = (grid_dict['z_r'][1:]/H) - (grid_dict['z_r'][:-1]/H)
    Hz = (grid_dict['z_w'][1:]/H) - (grid_dict['z_w'][:-1]/H)
        
    ###################
    # SET UP OUTPUT
    ###################
    N = len(grid_dict['z_r'])
    u_out = np.zeros([N])
    v_out = np.zeros([N])
    

    #############################
    # CREATE MATRICES TO SOLVE
    ############################
    ndim = 2 *N
    kv   = 2

    #COEFFICIENT MATRIX
    A = lil_matrix((ndim,ndim))
    #R.H.S MATRIX
    R = np.zeros(ndim)
   
    ###################
    # FILL BOTTOM
    ###################
    idx = 0
    C0 = Ek0 / dz[0]
    A[idx,idx+1]  = -1*Hz[0]
    A[idx,idx+kv] = -C0
   
    A[idx+1,idx] = 1*Hz[0]
    A[idx+1,idx+1+kv] = -C0    

    if BC_opt =='free':
       A[idx,idx] = C0
       A[idx+1,idx+1] = C0
    if BC_opt =='no_slip_bot':
       A[idx,idx] = C0 + (2*Ek0 / dz[0])
       A[idx+1,idx+1] = C0 + (2*Ek0 / dz[0])
       
    
    #FILL INTERIOR
    for k in range(1,N-1):
        C_k = Ek0 / dz[k-1]
        C_k1 = Ek0 / dz[k]
        C_k2 = C_k + C_k1 

        idx = 2*k
        #FIRST ROW
        A[idx,idx+kv] = -C_k1
        A[idx,idx]    = C_k2
        A[idx,idx-kv] = -C_k
        A[idx,idx+1]  = -1*Hz[k]

        #SECOND ROW
        A[idx+1,idx+1+kv] = -C_k1
        A[idx+1,idx+1]    = C_k2
        A[idx+1,idx+1-kv] = -C_k
        A[idx+1,idx]      = 1 *Hz[k]
    
    ###############
    # TOP
    ################
    idx = (2*N) - kv
    C_N = Ek0 / dz[-1]
    
    #FIRST ROW
    A[idx,idx-kv] = -C_N
    A[idx,idx]    = C_N
    A[idx,idx+1]  = -1*Hz[-1]

    #SECOND ROW
    A[idx+1,idx+1-kv] = -C_N
    A[idx+1,idx+1]    = C_N
    A[idx+1,idx]      = 1*Hz[-1]

    ##############################
    # FILL R.H.S MATRIX
    ##############################
    ##############
    # SURFACE
    ###############
    if sign_conv == 'neg':
       R[-2] = -Ek0
    if sign_conv == 'pos':
       R[-2] = Ek0
    
    ##############
    # BOTTOM
    ###############
    if BC_opt =='free':
       if sign_conv == 'neg':
          R[0] = Ek0
       if sign_conv == 'pos':
          R[0]  = -Ek0
    if BC_opt == 'no_slip_bot':
       R[0] = 0

    ######################
    # SOLVE SYSTEM
    #######################
    A = A.tocsr()
    X = spsolve(A,R)

    #REORDER RESULTS
    for k in range(N):
        idx = 2*k
        u_out[k] = X[idx]
        v_out[k] = X[idx+1]

    return u_out, v_out
    ###############################################


def solve_ustar_vstar(grid_dict,Omega,Gamma,K_ND_n,k,ubar_zz,vbar_zz,BC_opt):
    """
     SOLVE FOLLOWING DIAGNOSTIC SYSTEM AT A TIME-STEP
  
     Omega * (-vstar,ustar) - Gamma * (K_ND[n] + k) * (ustar_zz,vstar_zz) 
                          = Gamma * K_ND[n] *(ubar_zz,vbar_zz)
     
     K_ND_n --> K_ND at t=n
     
     DEFINE THESE AS CONSTANTS
     C = (Gamma * (K_ND[n] + k)) / Omega
     B = Gamma * K_ND[n] / Omega


     B.C.S 
     ustar_z,vstar_z = 0 at z=0,-1
     
     AS IN THE N.D. PROBLEM, VERTICAL STRUCTURE IS PASSED ON
     THROUGH THE STEADY, ubar,vbar SOLN IN THE R.HS. TERM ABOVE
    """

    ######################################
    # CALCULATE CONSTANTS FOR MATRICES
    #####################################
    C      = (Gamma * (K_ND_n + k)) / Omega
    B      = (Gamma * K_ND_n) / Omega  

    #VERTICAL SPACING
    H = abs((grid_dict['z_w'][0]))
    dz = (grid_dict['z_r'][1:]/H) - (grid_dict['z_r'][:-1]/H)
    Hz = (grid_dict['z_w'][1:]/H) - (grid_dict['z_w'][:-1]/H)
    
    
    ###################
    # SET UP OUTPUT
    ###################
    N = len(grid_dict['z_r'])
    u_out = np.zeros([N])
    v_out = np.zeros([N])
    

    #############################
    # CREATE MATRICES TO SOLVE
    ############################
    ndim = 2 *N
    kv   = 2

    #COEFFICIENT MATRIX
    A = lil_matrix((ndim,ndim))
    #R.H.S MATRIX
    R = np.zeros(ndim)
  

    ###################
    # FILL BOTTOM
    ###################
    idx = 0
    C0 = C / dz[0]
    A[idx,idx+1]  = -1*Hz[0]
    A[idx,idx+kv] = -C0
   
    A[idx+1,idx] = 1*Hz[0]
    A[idx+1,idx+1+kv] = -C0    

    if BC_opt =='free':
       A[idx,idx] = C0
       A[idx+1,idx+1] = C0
    if BC_opt =='no_slip_bot':
       A[idx,idx] = C0 + (2*C / dz[0])
       A[idx+1,idx+1] = C0 + (2*C / dz[0])
    
    #FILL INTERIOR
    for k in range(1,N-1):
        C_k = C / dz[k-1]
        C_k1 = C / dz[k]
        C_k2 = C_k + C_k1 

        idx = 2*k
        #FIRST ROW
        A[idx,idx+kv] = -C_k1
        A[idx,idx]    = C_k2
        A[idx,idx-kv] = -C_k
        A[idx,idx+1]  = -1*Hz[k]

        #SECOND ROW
        A[idx+1,idx+1+kv] = -C_k1
        A[idx+1,idx+1]    = C_k2
        A[idx+1,idx+1-kv] = -C_k
        A[idx+1,idx]      = 1*Hz[k]
    
    ###############
    # TOP
    ################
    idx = (2*N) - kv
    C_N = C / dz[-1]
    
    #FIRST ROW
    A[idx,idx-kv] = -C_N
    A[idx,idx]    = C_N
    A[idx,idx+1]  = -1*Hz[-1]

    #SECOND ROW
    A[idx+1,idx+1-kv] = -C_N
    A[idx+1,idx+1]    = C_N
    A[idx+1,idx]      = 1*Hz[-1]


    ##############################
    # FILL R.H.S MATRIX
    ##############################
    for k in range(N):
        idx = 2*k
        R[idx] = B * ubar_zz[k]
        R[idx+1] = B * vbar_zz[k]

    
    ######################
    # SOLVE SYSTEM
    #######################
    A = A.tocsr()
    X = spsolve(A,R)

    #REORDER RESULTS
    for k in range(N):
        idx = 2*k
        u_out[k] = X[idx]
        v_out[k] = X[idx+1]

    return u_out, v_out
    ###############################################


