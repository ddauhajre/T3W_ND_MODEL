# T3W_ND_MODEL
TO RUN THE MODEL

1) Declare parameters in /src_code/T3W_ND_params.py
2) From main directory, 
    >> sh run_T3W_ND.sh
    or 
    >> python src_code/T3W_ND_main.py


Once model runs, output will be saved in /output directory that is created as model runs

OVERVIEW OF CODE (in /src_code directory)
T3W_ND_params.py --> parameters file that is to be edited by user to specify model run
T3W_ND_funcs.py  --> library of functions used for time-stepping, solution setup, etc.
T3W_ND_main.py   --> main script that initiliazes parameters and calls time-stepping routines
