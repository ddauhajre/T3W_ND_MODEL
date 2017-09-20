######################################
__title__          = "T3W_ND_main.py"
__author__         = "Daniel Dauhajre"
__date__           = "September 2017"
__email__          = "ddauhajre@atmos.ucla.edu"
__python_version__ = "2.7.9"

'''
MAIN---> RUN N.D. MODEL AND SAVE OUTPUT
'''
######################################

######################################
# IMPORT NECESSARY MODULES
import os
import sys
import pickle as pickle
import numpy as np
import T3W_ND_funcs as T3W_ND_func
from pylab import *
import matplotlib.pyplot as plt
plt.ion()
######################################

######################################
# SET CONSTANTS
#####################################
code_path = './src_code/'
execfile(code_path + 'T3W_ND_params.py')

print '         ################################################'
print '                         RUN ID: ' + run_ID
print '         ################################################'


print '                 #############################    '
print '                 SIMULATION PARAMATERS/CHOICES    '
print '                 #############################    '
print '                 EK0 = ' + str(Ek0) 
print '                 Gamma = ' + str(Gamma) 
print '                 Omega = ' + str(Omega) 
print '                 k  = ' + str(k) 
print '                 alpha = ' + str(alpha)
print '                 beta = ' + str(beta)
print '                 hour_shift = ' + str(hour_shift)
print '                 BC_choice --> ' + BC_choice
print '                 tstep_scheme -->' + tstep_scheme
print ''


######################################
# SOLUTION DICTIONARY SETUP
######################################
if z_grd_N:
   ################################################
   # MAKE GRID BASED OFF NUMBER OF VERTICAL LEVELS
   ################################################
   grid_dict = T3W_ND_func.make_zt_grid_N_levs(N_levs,dt,tend_days,theta_b=theta_b)
else:
   grid_dict = T3W_ND_func.make_zt_grid(Lz_m, dz_m, dt, tend_days)

#u,v decomposed fields
[var_keys, var_key_types] = T3W_ND_func.setup_sol_vars_for_dict(len(grid_dict['tvec_sec']),len(grid_dict['z_r']))
var_sol_dict = T3W_ND_func.add_keys_var_dict(var_keys, var_key_types)

#grid variables
var_sol_dict = T3W_ND_func.add_keys_var_dict(['grid'], [grid_dict],var_dict=var_sol_dict)

#PARAMETERS
param_keys = ['Ek0', 'Gamma', 'Omega', 'k','K_shape', 'BC_choice','T_d','tstep_scheme','theta_b','hour_shift','alpha', 'beta', 'Ek_sign_conv']
param_key_types = [Ek0, Gamma, Omega, k, K_shape, BC_choice,T_d,tstep_scheme,theta_b,hour_shift,alpha,beta, Ek_sign_conv]
var_sol_dict = T3W_ND_func.add_keys_var_dict(param_keys, param_key_types, var_dict=var_sol_dict)

print ''
print 'SOLUTION DICTIONARY CREATED'
################################################################################################


print ''
print ''
print '         #################### '
print ''
print ''
print '         SOLVING N.D. EQNS'
print ''
print ''
print '         #################### '
##################################
# SOLVE FOR UBAR, VBAR
##################################
print 'Solving for ubar, vbar...'
var_sol_dict['ubar'], var_sol_dict['vbar'] = T3W_ND_func.solve_ubar_vbar(var_sol_dict['grid'],var_sol_dict['Ek0'], var_sol_dict['BC_choice'], sign_conv = var_sol_dict['Ek_sign_conv'])



###################################
# CREATE N.D. TEMPORAL SHAPE 
# FUNCTION FOR TEMPORAL FORCING 
# (i.e., mixing)
##################################
var_sol_dict['K_ND'] = T3W_ND_func.make_K_ND_comp(var_sol_dict['grid']['dt'], var_sol_dict['grid']['tvec_sec'], tend_days, var_sol_dict['hour_shift'], const_days,alpha=var_sol_dict['alpha'],beta=var_sol_dict['beta'])

###############################
# TIME STEP SYSTEM
##############################
print ''
print ''
print '         ########################## '
print '              TIME STEPPING               '
print '         ########################## '
var_sol_dict['ubar_zz'] = np.zeros(len(var_sol_dict['ubar']))
var_sol_dict['vbar_zz'] = np.zeros(len(var_sol_dict['ubar']))


var_sol_dict['u_p'], var_sol_dict['v_p'],var_sol_dict['ustar'],var_sol_dict['vstar'],var_sol_dict['ubar_zz'],var_sol_dict['vbar_zz'] = T3W_ND_func.tstep_system(var_sol_dict)

#####################################################
# CALCULATE DEVIATION of u_p, v_p from ustar,vstar
'''
i.e., periodic deviation from diagnostic steady
solution at each point in time
'''
##################################################
var_sol_dict['u_pp'] = var_sol_dict['u_p'] - var_sol_dict['ustar']
var_sol_dict['v_pp'] = var_sol_dict['u_p'] - var_sol_dict['vstar']

#####################################
# CALCULATE FULL FLOW
'''
u = u' + ubar
'''
####################################
var_sol_dict['u'] = var_sol_dict['u_p'] + var_sol_dict['ubar']
var_sol_dict['v'] = var_sol_dict['v_p'] + var_sol_dict['vbar']



#############################
#   SAVE OUTPUT
############################
T3W_ND_func.change_dir_py('output')
T3W_ND_func.save_to_pickle(var_sol_dict,run_ID)
os.chdir('..')






