import numpy as np

# ------------------------------------------------------------------------------------------------
# COMPONENTS
# ------------------------------------------------------------------------------------------------

CH4 = 0
CO2 = 1
H2O = 2
H2 = 3
CO = 4
AR = 5

# Nist data
A = np.array([-7.03029e-1,  2.499735e1,  3.0092e1, 3.3066178e1,  2.556759e1,  2.0786e1])
B = np.array([ 1.084773e2,  5.518696e1,  6.832514,-1.1363417e+01,  6.09613,  2.825911e-7])
C = np.array([-4.252157e1, -3.36913700e1,  6.793435, 1.1432816e+01,  4.054656, -1.464191e-7])
D = np.array([ 5.862788,  7.948387, -2.53448, -2.772874, -2.671301,  1.092131e-8])
E = np.array([ 6.78565e-1, -1.36638e-1,  8.2139e-2, -1.58558e-1,  1.31021e-1, -3.66137100e-8])
F = np.array([-7.684376e1, -4.036075e2, -2.50881e2, -9.980797, -1.180089e2, -6.19735])
G = np.array([ 1.587163e2,  2.282431e2,  2.233967e2, 1.72707974e2,  2.273665e2,  1.79999e2])
H = np.array([-7.48731e1, -3.935224e2, -2.418264e2, 0.0, -1.105271e2,  0.0])

HF298 = np.array([-74.9, -393.8, -242.0, 0.0, -110.6, 0.0]) # kJ/mol
GF298 = np.array([-50.9, -394.6, -228.8, 0.0, -137.4, 0.0]) # kJ/mol
SF298 = (np.array(HF298) - np.array(GF298)) / 298.15 # kJ/mol/K

def get_delta_coefficient(coef, positions, y):
    
    coef = np.array(coef)
    positions = np.array(positions)
    components = np.zeros(coef.shape, dtype=int) 
    components[positions] = y
    
    return coef.dot(components)

def get_delta_multiple_coefficient(coefs, positions, y):
    
    return np.array([get_delta_coefficient(coef, positions, y) for coef in coefs])

DELTA_DRM = get_delta_multiple_coefficient([A, B, C, D, E, F, G, H, HF298, GF298, SF298], [CO2, CH4, CO, H2], [-1, -1, 2, 2])
DELTA_SRM1 = get_delta_multiple_coefficient([A, B, C, D, E, F, G, H, HF298, GF298, SF298], [CH4, H2O, CO, H2], [-1, -1, 1, 3])
DELTA_SRM2 = get_delta_multiple_coefficient([A, B, C, D, E, F, G, H, HF298, GF298, SF298], [CH4, H2O, CO2, H2], [-1, -2, 1, 4])
DELTA_WGS = get_delta_multiple_coefficient([A, B, C, D, E, F, G, H, HF298, GF298, SF298], [CO, H2O, CO2, H2], [-1, -1,  1,  1])

def get_deltas(*args):
    
    mat = np.vstack(args)
    
    return mat.T

DELTA_A, DELTA_B, DELTA_C, DELTA_D, DELTA_E, DELTA_F, DELTA_G, DELTA_H, DELTA_HR298, DELTA_GF298, DELTA_SR298 = \
    get_deltas(DELTA_DRM, DELTA_SRM1, DELTA_SRM2, DELTA_WGS)