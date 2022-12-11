import thermodynamics
import numpy as np

######################################################################################################################
# INPUT DATA
# coefficients of shomate equation - nist data from https://webbook.nist.gov/

co2 = [24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224]
ch4 = [-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, 158.7163, -74.87310]
h2o = [30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264]
h2 = [33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0]
co = [25.56759, 6.096130, 4.054656,  -2.671301,  0.131021, -118.0089, 227.3665, -110.5271]
ar = [20.78600, 2.825911e-07, -1.464191e-07,  1.092131e-08, -3.661371e-08, -6.197350, 179.9990, 0.00000]

######################################################################################################################
# enthalpy of formation in kJ/mol
HF298_ch4 = -74.85
HF298_co2 = -393.51
HF298_h2o = -241.826
HF298_h2 = 0.0
HF298_co = -110.53
HF298_ar = 0.0

######################################################################################################################
# entropy of formation in J/mol/K
SF298_ch4 = 186.25
SF298_co2 = 197.66
SF298_h2o = 188.84
SF298_h2 = 130.68
SF298_co = 213.79
SF298_ar = 154.84

######################################################################################################################
# molar masses of the species containing in the mixture in kg/mol
mr_ch4 = 16.04e-03
mr_co2 = 44.01e-03
mr_h2o = 18.01528e-03
mr_h2 = 2.01588e-03
mr_co = 28.01e-03
mr_ar = 39.948e-03

######################################################################################################################
# determine temperature and pressure for simulation
T = 873.15 # K
P = 1 # bar

######################################################################################################################
# CO2 + CH4 ⇌ 2CO + 2H2 - DRM reaction
coeffs_DRM = np.array([co2, ch4, co, h2])

# heats of formation at 298.15 K for CO2, CH4, CO, H2
HF298_DRM = np.array([HF298_co2, HF298_ch4, HF298_co, HF298_h2]) # kJ/mol
SF298_DRM = np.array([SF298_co2, SF298_ch4, SF298_co, SF298_h2]) # J/mol/K

# stoichiometry of reactants and products
v_DRM = np.array([-1, -1, 2, 2])

# create a reaction object and print its stored parameters:
reaction_DRM = thermodynamics.REACTION(coeffs_DRM, HF298_DRM, SF298_DRM, v_DRM)

reaction_DRM.print_params()

print("DRM enthalpy at {} K is {}".format(T, thermodynamics.enthalpy_rxn_kj_mol(reaction_DRM, T)[0]), "kJ/mol")
print("DRM entropy at {} K is {}".format(T, thermodynamics.entropy_rxn_kj_mol(reaction_DRM, T)[0]), "kJ/mol/K")
print("DRM gibbs energy at {} K is {}".format(T, thermodynamics.gibbs_rxn_kj_mol(reaction_DRM, T)[0]), "kJ/mol")
print("DRM Keq at {} K is {}".format(T, np.exp(-thermodynamics.equilibrium_constant(reaction_DRM, T)[0])))

######################################################################################################################
# CH4 + H2O ⇌ CO + 3H2 - SRM1 reaction
coeffs_SRM1 = np.array([ch4, h2o, co, h2])

# heats of formation at 298.15 for CH4, H2O, CO, H2
HF298_SRM1 = np.array([HF298_ch4, HF298_h2o, HF298_co, HF298_h2]) # kJ/mol
SF298_SRM1 = np.array([SF298_ch4, SF298_h2o, SF298_co, SF298_h2]) # J/mol/K

# stoichiometry of reactants and products
v_SRM1 = np.array([-1, -1, 1, 3])

# create a reaction object and print its stored parameters:
reaction_SRM1 = thermodynamics.REACTION(coeffs_SRM1, HF298_SRM1, SF298_SRM1, v_SRM1)

reaction_SRM1.print_params()

print("SRM1 enthalpy at {} K is {}".format(T, thermodynamics.enthalpy_rxn_kj_mol(reaction_SRM1, T)[0]), "kJ/mol")
print("SRM1 entropy at {} K is {}".format(T, thermodynamics.entropy_rxn_kj_mol(reaction_SRM1, T)[0]), "kJ/mol/K")
print("SRM1 gibbs energy at {} K is {}".format(T, thermodynamics.gibbs_rxn_kj_mol(reaction_SRM1, T)[0]), "kJ/mol")
print("SRM1 Keq at {} K is {}".format(T, np.exp(-thermodynamics.equilibrium_constant(reaction_SRM1, T)[0])))

######################################################################################################################
# CH4 + 2H2O ⇌ CO2 + 4H2 - SRM2 reaction
coeffs_SRM2 = np.array([ch4, h2o, co2, h2])

# heats of formation at 298.15 K for CH4, H2O, CO2, H2
HF298_SRM2 = np.array([HF298_ch4, HF298_h2o, HF298_co2, HF298_h2]) # kJ/mol
SF298_SRM2 = np.array([SF298_ch4, SF298_h2o, SF298_co2, SF298_h2]) # J/mol/K

# stoichiometry of reactants and products
v_SRM2 = np.array([-1, -2, 1, 4])

# create a reaction object and print its stored parameters:
reaction_SRM2 = thermodynamics.REACTION(coeffs_SRM2, HF298_SRM2, SF298_SRM2, v_SRM2)

reaction_SRM2.print_params()

print("SRM2 enthalpy at {} K is {}".format(T, thermodynamics.enthalpy_rxn_kj_mol(reaction_SRM2, T)[0]), "kJ/mol")
print("SRM2 entropy at {} K is {}".format(T, thermodynamics.entropy_rxn_kj_mol(reaction_SRM2, T)[0]), "kJ/mol/K")
print("SRM2 gibbs energy at {} K is {}".format(T, thermodynamics.gibbs_rxn_kj_mol(reaction_SRM2, T)[0]), "kJ/mol")
print("SRM2 Keq at {} K is {}".format(T, np.exp(-thermodynamics.equilibrium_constant(reaction_SRM2, T)[0])))

######################################################################################################################
# CO + H2O ⇌ CO2 + H2 - WGS reaction
coeffs_WGS = np.array([co, h2o, co2, h2])

# heats of formation at 298.15 K for CO, H2O, CO2, H2
HF298_WGS = np.array([HF298_co, HF298_h2o, HF298_co2, HF298_h2]) # kJ/mol
SF298_WGS = np.array([SF298_co, SF298_h2o, SF298_co2, SF298_h2]) # J/mol/K

# stoichiometry of reactants and products
v_WGS = np.array([-1, -1,  1,  1])

# create a reaction object and print its stored parameters:
reaction_WGS = thermodynamics.REACTION(coeffs_WGS, HF298_WGS, SF298_WGS, v_WGS)

reaction_WGS.print_params()

print("WGS enthalpy at {} K is {}".format(T, thermodynamics.enthalpy_rxn_kj_mol(reaction_WGS, T)[0]), "kJ/mol")
print("WGS entropy at {} K is {}".format(T, thermodynamics.entropy_rxn_kj_mol(reaction_WGS, T)[0]), "kJ/mol/K")
print("WGS gibbs energy at {} K is {}".format(T, thermodynamics.gibbs_rxn_kj_mol(reaction_WGS, T)[0]), "kJ/mol")
print("WGS Keq at {} K is {}".format(T, np.exp(-thermodynamics.equilibrium_constant(reaction_WGS, T)[0])))

######################################################################################################################
# determine species in the mixture
coeffs_species = np.array([ch4, co2, h2o, h2, co, ar])

# heats of formation at 298.15 K for CH4, CO2, H2O, H2, CO, Ar in kJ/mol
HF298_species = np.array([HF298_ch4, HF298_co2, HF298_h2o,  HF298_h2, HF298_co, HF298_ar])

# molar mass of species in the mixture in kg/mol
Mr_species = np.array([mr_ch4, mr_co2, mr_h2o, mr_h2, mr_co, mr_ar])

# molar ratio of the species in the mixture dimensionless 50% biogas - 50% argon
y_species = np.array([0.3, 0.2, 0.0, 0.0, 0.0, 0.5])

######################################################################################################################
# create an object with mixture parameters
species_enthalpy = thermodynamics.MIXTURE(coeffs_species, HF298_species, y_species, Mr_species)

species_enthalpy.print_params()

T = 873.15 # K

print("partial enthalpies of species in the mixture at {} K is \n{}"
.format(T, thermodynamics.partial_enthalpies_kj_mol(species_enthalpy, 873.15)[0]), "\n in kJ/mol")

######################################################################################################################
# create an object with mixture parameters
species_Cp = thermodynamics.MIXTURE(coeffs_species[:, :-3], HF298_species, y_species, Mr_species)

species_Cp.print_params()

print("partial specific heats of species in the mixture at {} K is \n{}"
      .format(T, thermodynamics.cpi_j_mol_K(species_Cp, T)[0]), "\n in kJ/mol")

print("concentrations of species in the mixture at {} K and {} bar is \n{}"
      .format(T, P, thermodynamics.ci_mol_m3(species_Cp, T, P)[0]), "in mol/m3 ")

print("total concentration of mixture at {} K and {} bar is \n{}"
      .format(T, P, thermodynamics.c_total_mol_m3(species_Cp, T, P)[0]), "in mol/m3 ")

print("density of the mixture at {} K and {} bar is \n{}"
      .format(T, P, thermodynamics.d_mix_kg_m3(species_Cp, T, P)[0]), "in kg/m3 ")

print("molar mass of the mixture at {} K and {} bar is \n{}"
      .format(T, P, thermodynamics.mr_mix_kg_mol(species_Cp)[0]), "in kg/m3 ")

######################################################################################################################