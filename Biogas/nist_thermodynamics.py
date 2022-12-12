import numpy as np

def heat_capacity(a, b, c, d, e, T):
    """
    Returns the heat capacity of one component in J/mol at temperature T in K using Nist equation.
    Parameters
    ----------
    a : float
        Coefficient from Nist.
    b : float
        Coefficient from Nist.
    c : float
        Coefficient from Nist.
    d : float
        Coefficient from Nist.
    e : float
        Coefficient from Nist.
    T : float or int
        Temperature in K.
    Returns
    -------
    Cp: float
        Heat capacity of one component in J/mol.
    """
    T = T/1000
    
    return a + b*T + c*T**2 + d*T**3 + e/T

def get_Cp(T, a, b, c, d, e):
    """
    Returns the heat capacities of species at temperature T in K
    Parameters
    ----------
    T : float or int
        Temperature in K.
    a : 1d array or float
        Coefficients for heat capacity from Nist equation.
    b : 1d array or float
        Coefficients for heat capacity from Nist equation.
    c : 1d array or float
        Coefficients for heat capacity from Nist equation.
    d : 1d array or float
        Coefficients for heat capacity from Nist equation.
    e : 1d array or float
        Coefficients for heat capacity from Nist equation.
    Returns
    -------
    res : 1d array or float
        Heat capacities in J/mol/K.
    """
    
    return heat_capacity(np.array(a), np.array(b), np.array(c), np.array(d), np.array(e), T)

# ------------------------------------------------------------------------------------------------
# DELTAS OF REACTION
# ------------------------------------------------------------------------------------------------

def calc_delta_hr(T, delta_hr298, va, vb, vc, vd, ve, vf, vg, vh):
    """
    Returns the heat of reaction in kJ/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    delta_hr298 : 1d array, float or int
        Contains heat of reaction at 298.15K.
    va : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vb : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vc : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vd : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    ve : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vf : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vg : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vh : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    Returns
    -------
    1d array, float or int
        Contains heat of reaction at T in kJ/mol.
    """
    
    T = T/1000
    
    return np.array(delta_hr298) + np.array(va)*(T) + np.array(vb)/2*(T**2)\
        + np.array(vc)/3*(T**3) + np.array(vd)/4*(T**4) - np.array(ve)/(T)\
        + 1*np.array(vf) + 0*np.array(vg) - 1*np.array(vh)


def calc_hf_temp(T, hf_298, a, b, c, d, e, f, g, h):
    """
    Returns the heat of formation in kJ/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    hf_298 : 1d array, float or int
        Contains heat of formation at 298.15K.
    a : 1d array, float or int
        Contains a coefficient of the heat capacity.
    b : 1d array, float or int
        Contains a coefficient of the heat capacity.
    c : 1d array, float or int
        Contains a coefficient of the heat capacity.
    d : 1d array, float or int
        Contains a coefficient of the heat capacity.
    e : 1d array, float or int
        Contains a coefficient of the heat capacity.
    f : 1d array, float or int
        Contains a coefficient of the heat capacity.
    g : 1d array, float or int
        Contains a coefficient of the heat capacity.
    h : 1d array, float or int
        Contains a coefficient of the heat capacity.
    Returns
    -------
    1d array, float or int
        Contains heat of formation at T in kJ/mol.
    """
    T = T/1000
    
    return np.array(hf_298) + np.array(a)*(T) + np.array(b)/2*(T**2)\
        + np.array(c)/3*(T**3) + np.array(d)/4*(T**4) - np.array(e)/(T)\
        + 1*np.array(f) + 0*np.array(g) - 1*np.array(h)


def calc_delta_sr(T, delta_sr298, va, vb, vc, vd, ve, vf, vg):
    """
    Returns the entropy of reaction in J/mol/K.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    delta_sr298 : 1d array, float or int
        Contains entropy of reaction at 298.15K.
    va : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vb : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vc : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vd : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    ve : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vf : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    vg : 1d array, float or int
        Contains delta coefficient of the reaction(s).
    Returns
    -------
    1d array, float or int
        Contains entropy of reaction at T in J/mol/K.
    """
    T = T/1000
    
    return np.array(delta_sr298) + np.array(va)*np.log(T) + np.array(vb)*(T)\
        + np.array(vc)/2*(T**2) + np.array(vd)/3*(T**3) - np.array(ve)/(2*T**2)\
        + 0*vf + 1*vg

def calc_delta_gibbs(T, delta_hr, delta_sr, multireaction=False):
    """
    Returns the gibbs energy of reaction in J/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    delta_hr : float, int, or 1d array
        Contains heat of reaction at T in kJ/mol.
    delta_sr : float, int, or 1d array
        Contains entropy of reaction at T in J/mol/K.
    Returns
    -------
    float, int, or 1d array
        Contains gibbs energy of reaction at T in J/mol.
    """
    return np.array(delta_hr)*1000 - T*np.array(delta_sr)


def calc_equilibrium_constant(T, delta_gibbs, R = 8.314):
    """
    Returns the dimensionless equilibrium constant of reaction.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    delta_gibbs : float, int, or 1d array
        Contains heat of reaction at T in J/mol.
    R : float
        Gas constant is 8.314 J/mol/K.
    Returns
    -------
    float, int, or 1d array
        Contains dimensionless equilibrium constant of reaction at T.
    """
    return np.exp(-np.array(delta_gibbs)/R/T)