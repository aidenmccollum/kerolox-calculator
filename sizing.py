#Initial Sizing Calculator for a Rapidly Iterable Kerolox Engine

import math
from pint import UnitRegistry
ureg = UnitRegistry()
#For finding available units: ureg.get_compatible_units('[mass]')

# ASSUMPTIONS
P0 = 300 * ureg.pound_force_per_square_inch #Chamber pressure
Pa = 14.7 * ureg.pound_force_per_square_inch #atmospheric pressure
Pe = 10.2 * ureg.pound_force_per_square_inch #Exit Pressure (optimized for 10,000 ft)
of_ratio = 2
mdot_ox = 2 * ureg.kilogram / ureg.second
mdot_fuel = mdot_ox/of_ratio
mdot_total = mdot_fuel+mdot_ox

Lstar = 45 * ureg.inch #assumed for a kerolox engine
print(f'\n\nTotal Mass Flow Rate: {mdot_total}\n')

#CEA OUTPUTS
gamma = 1.1549
Cp = 5035.9 * ureg.joule / ureg.kilogram / ureg.kelvin   #specific heat
R = (1-(1/gamma))*Cp #specific gas constant
print(f'Specific Gas Constant: {R}')
T0 = 3378.67 * ureg.kelvin #chamber temperature

#Sizing the nozzle Throat
At = mdot_total/(P0.to(ureg.pascal)*gamma*(math.sqrt(((2/(gamma+1))**((gamma+1)/(gamma-1))))/math.sqrt(gamma*R*T0 / (ureg.joule/ureg.kilogram)))*ureg.second/ureg.meter) #throat area, mass flow relation

#finding the mach number
M_chamber = math.sqrt((((P0/Pe)**((gamma-1)/gamma))-1)*2/(gamma-1))

#sizing the area ratios 
Ae_At = (((gamma+1)/2)**-((gamma+1)/(2*(gamma-1))))*(((1+((gamma-1)/2)*M_chamber**2)**((gamma+1)/(2*(gamma-1))))/M_chamber)
Dt = 2*math.sqrt(At/math.pi/(ureg.inch**2)) * ureg.inch
Ae = (Ae_At*At).to(ureg.inch**2)
print(f'Nozzle Characteristics:\n   Exit Mach: {M_chamber}\n   Throat Diameter: {Dt}\n   Area Ratio: {Ae_At}\n   Exit Area: {Ae}\n')

#Sizing the chamber 
V_chamber = Lstar*At
L_chamber = 8 * ureg.inch #from throat diameter to chamber length chart
A_chamber = V_chamber/(L_chamber)
print(f'Chamber Characteristics:\n   Chamber Length: {L_chamber.to(ureg.inch)}\n   Chamber Diameter: {2*math.sqrt(A_chamber/math.pi/(ureg.inch**2))*ureg.inch}\n')

#finding the thrust 
M_exit = 2.712
T_exit = T0/(1+0.5*(gamma-1)*M_exit**2)
U_exit = M_exit*math.sqrt(gamma*R*T_exit / (ureg.joule/ureg.kilogram)) * ureg.meter / ureg.second
thrust = mdot_total*U_exit + (Pe-Pa)*Ae
print(f'Thrust: {thrust.to(ureg.pound_force)}')

