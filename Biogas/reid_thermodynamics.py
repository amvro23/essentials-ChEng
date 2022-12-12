import numpy as np

def heat_capacity(a, b, c, d, T):
    """
    Returns the heat capacity of one component in J/mol at temperature T in K using Reid (1987) equation.
    Parameters
    ----------
    a : float
        Coefficient from Reid.
    b : float
        Coefficient from Reid.
    c : float
        Coefficient from Reid.
    d : float
        Coefficient from Reid.
    T : float or int
        Temperature in K.
    Returns
    -------
    Cp: float
        Heat capacity of one component in J/mol.
    """
    return a + b*T + c*T**2 + d*T**3


def get_Cp(T, a, b, c, d):
    """
    Returns the heat capacities of species at temperature T in K
    Parameters
    ----------
    T : float or int
        Temperature in K.
    a : 1d array or float
        Coefficients for heat capacity from Reid equation.
    b : 1d array or float
        Coefficients for heat capacity from Reid equation.
    c : 1d array or float
        Coefficients for heat capacity from Reid equation.
    d : 1d array or float
        Coefficients for heat capacity from Reid equation.
    Returns
    -------
    res : 1d array or float
        Heat capacities in J/mol/K.
    """
    return heat_capacity(np.array(a), np.array(b), np.array(c), np.array(d), T)

# ------------------------------------------------------------------------------------------------
# DELTAS OF REACTION
# ------------------------------------------------------------------------------------------------

def calc_delta_hr(T, delta_hr298, va, vb, vc, vd):
    """
    Returns the heat of reaction in J/mol.
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
    Returns
    -------
    1d array, float or int
        Contains heat of reaction at T in J/mol.
    """
    return np.array(delta_hr298) + np.array(va)*(T-298.15) + np.array(vb)/2*(T**2-298.15**2)\
        + np.array(vc)/3*(T**3-298.15**3) + np.array(vd)/4*(T**4-298.15**4)


def calc_hf_temp(T, hf_298, a, b, c, d):
    """
    Returns the heat of formation in J/mol.
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
    Returns
    -------
    1d array, float or int
        Contains heat of formation at T in J/mol.
    """
    return np.array(hf_298) + np.array(a)*(T-298.15) + np.array(b)/2*(T**2-298.15**2)\
        + np.array(c)/3*(T**3-298.15**3) + np.array(d)/4*(T**4-298.15**4)


def calc_delta_sr(T, delta_sr298, va, vb, vc, vd):
    """
    Returns the entropy of reaction in J/mol.
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
    Returns
    -------
    1d array, float or int
        Contains entropy of reaction at T in J/mol/K.
    """
    return np.array(delta_sr298) + np.array(va)*np.log(T/298.15) + np.array(vb)*(T-298.15)\
        + np.array(vc)/2*(T**2-298.15**2) + np.array(vd)/3*(T**3-298.15**3)


def calc_delta_gibbs(T, delta_hr, delta_sr, multireaction=False):
    """
    Returns the gibbs energy of reaction in J/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    delta_hr : float, int, or 1d array
        Contains heat of reaction at T in J/mol.
    delta_sr : float, int, or 1d array
        Contains entropy of reaction at T in J/mol/K.
    Returns
    -------
    float, int, or 1d array
        Contains gibbs energy of reaction at T in J/mol.
    """
    return np.array(delta_hr) - T*np.array(delta_sr)


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
        Gas constant at J/mol/K.
    Returns
    -------
    float, int, or 1d array
        Contains dimensionless equilibrium constant of reaction at T.
    """
    return np.exp(-np.array(delta_gibbs)/R/T)