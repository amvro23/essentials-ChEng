from nist_data import \
    (A, B, C, D, E, F, G, H, DELTA_A, DELTA_B, DELTA_C, DELTA_D, DELTA_E, DELTA_F, DELTA_G, DELTA_H, DELTA_HR298, DELTA_SR298)
from nist_thermodynamics import (get_Cp, calc_delta_hr, calc_delta_sr, calc_delta_gibbs, calc_equilibrium_constant)
#---------------------------------------------------------------------------------------------------------------------------------------
T = 873.15

#---------------------------------------------------------------------------------------------------------------------------------------
Cp = get_Cp(T, A, B, C, D, E)
# J/mol/K
Cp

#---------------------------------------------------------------------------------------------------------------------------------------
delta_hr = calc_delta_hr(T, DELTA_HR298, DELTA_A, DELTA_B, DELTA_C, DELTA_D, DELTA_E, DELTA_F, DELTA_G, DELTA_H)
# J/mol
delta_hr*1000

#---------------------------------------------------------------------------------------------------------------------------------------
delta_sr = calc_delta_sr(T, DELTA_SR298, DELTA_A, DELTA_B, DELTA_C, DELTA_D, DELTA_E, DELTA_F, DELTA_G)
# J/mol/K
delta_sr

#---------------------------------------------------------------------------------------------------------------------------------------
delta_gibbs = calc_delta_gibbs(T, delta_hr, delta_sr, multireaction=False)
# J/mol
delta_gibbs

#---------------------------------------------------------------------------------------------------------------------------------------
Keq = calc_equilibrium_constant(T, delta_gibbs, R = 8.314)

Keq