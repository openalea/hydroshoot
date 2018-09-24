from math import *
from scipy import exp

R = 8.314510 # L kPa mol-1 K-1


def saturated_air_vapor_pressure(Tac):
    """ saturated vapor pressure in the ambiant air (kPa)"""
    return 0.611*exp(17.27*Tac/(237.3+Tac))
 
def relative_humidity(ea, es_a):
    """  compute Relative humidity (pourcent)"""
    return (ea/es_a)*100.

def air_vapor_pressure_deficit(Tac, hs):
    """ compute air Vapor pressure deficit at air temperature"""
    es_a = saturated_air_vapor_pressure(Tac) #% saturated vapor pressure in the ambiant air (kPa)#saturated vapor pressure in the ambiant air (kPa)
    ea = es_a*hs/100 #% vapor pressure in the ambiant air (kPa)
    return es_a-ea #% Vapor pressure deficit at air temperature (Kpa)

def solar_declination (DOY):
    """ Declinaison (rad) du soleil en fonction du jour de l'annee """
    alpha=2*3.14*(DOY-1)/365
    return (0.006918-0.399912*cos(alpha)+0.070257*sin(alpha))

def extraterrestrial_solar_irradiance (DOY,HU,latitude):
    """ rayonnement extraterrestre horarire """
    hrad=2*3.14/24*(HU-12)
    lat=radians(latitude)
    dec=solar_declination (DOY)
    costheta=sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(hrad)
    Io=1370*(1+0.033*cos(2*3.14*(DOY-4)/366))#eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
    So=Io*costheta #eclairement dans un plan parallele a la surface du sol
    return So
    #extraterrestrial_solar_irradiance (100,11,44)

def clear_sky_global_radiation(Ra, z):
    """ compute clear sky global radiation according to extraterestrial radiation and altitude (m)
    eq 37 - FA056, p 51"""
    return Ra*(3600./1e6)*(0.75 + 2e-5*z) #formule pour Ra en MJ.m-2.h-1
    #clear_sky_global_radiation(extraterrestrial_solar_irradiance (100,11,44),10.)


def canopy_net_longwave_radiation_loss(Tac, Rs, Rs0):
    """ compute LongwaveNetRadiation = f(Tk_a,Tac_x, ea, Rs,Rs0) after FA0 56 eq #(39)
    Tak  : absolute air temperature (K)
    es_a : saturated vapor pressure in the ambiant air (kPa)
    Rs : measured solar radiation (MJ m2 hour)
    Rs0 : calculated clear-sky solar radiation (MJ m2 hour)"""

    ## Declaration des constantes
    sigma = 2.0412*1e-10 # Stefan-Boltzmann constant per surface area (MJ m-2 K-4 H-1)

    ## Compute longwave net radiation
    es_a = saturated_air_vapor_pressure(Tac) #% saturated vapor pressure in the ambiant air (kPa)
    Tak = Tac+273.16
    ratio = min(1., Rs/Rs0)
    Rnl= sigma*(Tak**4)*(0.34-0.14*(es_a**0.5))*(1.35*ratio-0.35)
    return Rnl*1e6/3600. #(en w.m-2)
    #canopy_net_longwave_radiation_loss(25., 600.*3600./1e6, clear_sky_global_radiation(extraterrestrial_solar_irradiance(100,11,44),10.))

def net_absorbed_radiation(Rg, Tac, DOY, HU, latitude = 0.44, altitude = 0., albedo=0.2):
    """ """
    Rextra = extraterrestrial_solar_irradiance(DOY,HU,latitude)
    Rs0_ = clear_sky_global_radiation(Rextra, altitude)
    Rnl = canopy_net_longwave_radiation_loss(Tac, Rg*3600./1e6, Rs0_)
    Rabs=Rg*(1-albedo)-Rnl

    return Rabs
    #net_absorbed_radiation(600., 25., 100, 11)


def boundary_layer_conductance(u, w=0.1):
    """ compute Boundary layer conductance for CO2= f(u, w) after Kim and Lieth (2003)
    u : wind speed m s-1
    w : leaf Characteristic dimension in relation to wind speed (m)"""

    d = 0.72*w # leaf dimension (m)
    gb = 0.147*(u/d)**0.5  # boundary layer conductance (mol m2 s-1)
    # gb = 3.33 #  pour les simulations a niveau de feuille, il faut introduire la valeur de gb du LcPro (rb =0.33 mol m-2 s-1)
    return gb
    #boundary_layer_conductance(2.)

def Kelvin(T):
    """
    Converts from Celsius to absolute temperature.
    """
    Tak = T + 273.
    return Tak
    
def VPD_leaf_air(Tac, Tlc, hs):
    """
    Returns leaf to air vapour pressure deficit [kPa].
    
    Parameters:
    -`Tac`: air temperature [degreeC]
    -`Tlc`: leaf temperature [degreeC]
    - `hs`: air relative humidity (%)
    """
    es_l = saturated_air_vapor_pressure(Tlc) #% saturated vapor pressure in the leaf (kPa)
    es_a = saturated_air_vapor_pressure(Tac) #% saturated vapor pressure in the ambiant air (kPa)
    ea = es_a*hs/100 #% vapor pressure in the ambiant air (kPa)
    return es_l-ea

def cmol2cpa(temp,conc=400.):
    """
    Returns CO2 partial pressure [ubar] as a function of temperature [degrees]
    """
    V = 1.e6 * R * (temp+273) / 101.3 # Volume of 1000000 mols of air [L]
    P = conc* R * (temp+273) / V # CO2 partial pressure [kPa]
    return P*1.e4

def cpa2cmol(temp,partpress):
    """
    Returns CO2 concentration [umol mol-1] as a function of temperature [degrees]
    """
    V = 1.e6 * R * (temp+273) / 101.3 # Volume of 1000000 mols of air [L]
    conc = partpress* V / (R * (temp+273)) # CO2 partial pressure [kPa]
    return conc*1.e-4