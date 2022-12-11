import numpy as np

# CLASSES USED
class REACTION:
    """
    Store reaction info here
    """
    def __init__(self, coeffs, HF298, SF298, v):
        """
        Pass parameters desribing reactions
        """
        # NIST coefficients of species participating in the reaction
        self.coeffs = coeffs
        # enthalpy of formation in kJ/mol
        self.HF298 = HF298
        # entropy of formation in J/mol/K
        self.SF298 =  SF298
        # stoichiometry
        self.v = v
        
    def print_params(self):
        """
        Print molecule parameters.
        """
        print("""coeffs: {}
        \t HF298 = {} in kJ/mol
        \t SF298 = {} in J/mol/K
        \t stoichiometry = {} """.format(self.coeffs, self.HF298, self.SF298, self.v))
    
class MIXTURE:
    """
    Store reaction info here
    """
    def __init__(self, coeffs, HF298, y, Mr):
        """
        Pass parameters desribing species in the mixture
        """
        # NIST coefficients of species in the mixture
        self.coeffs = coeffs
        # enthalpy of formation in kJ/mol
        self.HF298 = HF298
        # molar ratio of species in the mixture dimensionless
        self.y = y
        # molar mass of species in the mixture in kg/mol
        self.Mr =  Mr
    
    def print_params(self):
        """
        Print mixture parameters.
        """
        print("""coeffs: {}
        \t HF298 = {} in kJ/mol
        \t molar ratio = {}
        \t molar mass = {} in kg/mol""".format(self.coeffs, self.HF298, self.y, self.Mr))

# EQUATIONS USED
def enthalpy(T):
    """
    - param T: float temperature in Kelvin
    
    Returns an array of enthalpy shomate equation in kJ/mol.
    """
    t = T/1000
    
    return  np.array([t,  t**2 / 2.0, t**3 / 3.0, t**4 / 4.0, -1.0 / t, 1.0, 0.0, -1.0], dtype = "object")

def entropy(T):
    """
    - param T: float temperature in Kelvin
    
    Returns an array of entropy shomate equation in J/mol/K.
    """    
    t = T/1000
    
    return np.array([np.log(t), t,  t**2 / 2.0,  t**3 / 3.0, -1.0 / (2.0 * t**2), 0.0, 1.0, 0.0], dtype = "object")

def specific_heat(T):
    """
    - param T: float temperature in Kelvin
    
    Returns an array of specific shomate equation in J/mol/K.
    """        
    t = T/1000
    
    return np.array([1,  t, t**2, t**3, 1.0/t], dtype = "object")

def enthalpy_rxn_kj_mol(rxn, T): 
    """
    - param rxn: reaction of interest
    - param T: float temperature in Kelvin
    
    Returns an array of reaction enthalpy in kJ/mol at this T and pressure 1 bar.
    """
    Hrxn298 = rxn.HF298.dot(rxn.v)

    Hrxn = Hrxn298 + (rxn.coeffs.dot(enthalpy(T))).dot(rxn.v)

    return np.array([Hrxn])

def entropy_rxn_kj_mol(rxn, T):  
    """
    - param rxn: reaction of interest
    - param T: float temperature in Kelvin
    
    Returns an array of reaction entropy in kJ/mol/K at this T and pressure 1 bar.
    """
    Srxn = (rxn.coeffs.dot(entropy(T)/1000.0)).dot(rxn.v)

    return np.array([Srxn])

def gibbs_rxn_kj_mol(rxn, T):  
    """
    - param rxn: reaction of interest
    - param T: float temperature in Kelvin
    
    Returns an array of reaction gibbs energy in kJ/mol at this T and pressure 1 bar.
    """
    Grxn = enthalpy_rxn_kj_mol(rxn, T) - T*entropy_rxn_kj_mol(rxn, T)

    return Grxn

def equilibrium_constant(rxn, T):
    """
    - param rxn: reaction of interest
    - param T: float temperature in Kelvin
    
    Returns an array of reaction Keq at this T and pressure 1 bar.
    """
    return np.exp(-gibbs_rxn_kj_mol(rxn, T)/8.314e-3/T)

def partial_enthalpies_kj_mol(species, T):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    
    Returns an array of partial enthalpies in kJ/mol at this T and pressure 1 bar.
    """
    return species.coeffs.dot(enthalpy(T)) + species.HF298


def cpi_j_mol_K(species, T):  
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    
    Returns an array of partial specific heats in J/mol/K at this T and pressure 1 bar.
    """
    return np.array([species.coeffs.dot(specific_heat(T))])


def cp_mix_j_mol_K(species, T):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    
    Returns an array of mixture's specific heat in J/mol/K at this T and pressure 1 bar.
    """
    
    return np.array([species.coeffs.dot(specific_heat(T)).dot(species.y)])

def cp_mix_kg_mol_K(species, T):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    
    Returns an array of mixture's specific heat in J/kg/K at this T and pressure 1 bar.
    """
    num = (species.coeffs.dot(specific_heat(T))).dot(species.y)
    den = species.Mr.dot(species.y)
    
    return np.array([num/den])

def ci_mol_m3(species, T, P):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    - param P: float pressure in bar
    
    Returns an array of partial concentrations of species mol/m3 at specific T and P.
    """
    return np.array([species.y*P/(8.314e-05*T)])

def c_total_mol_m3(species, T, P):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    - param P: float pressure in bar
    
    Returns an array of the mixture's total concentration at specific T and P.
    """
    return np.array([np.sum(ci_mol_m3(species, T, P))])

def d_mix_kg_m3(species, T, P):
    """
    - param species: species of interest
    - param T: float temperature in Kelvin
    - param P: float pressure in bar
    
    Returns an array of the mixture's density in kg/m3 at specific T and P.
    """
    return c_total_mol_m3(species, T, P)*species.Mr.dot(species.y)

# molar mass of the mixture in kg/mol
def mr_mix_kg_mol(species):
    """
    - param species: species of interest
    
    Returns an array of the mixture's molar mass in kg/mol
    """
    return np.array([species.Mr.dot(species.y)])
