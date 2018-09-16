# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Gas-exchange module of HydroShoot.

This module computes net photosynthesis and stomatal conductance rates.

"""
from copy import deepcopy
from scipy import exp, arccos, sqrt, cos, log

from hydroshoot import meteo_utils as mutils

# Constants
O = 210 # Oxygen partial pressure [mmol mol-1]
R = 0.0083144598 # Ideal gaz constant [kJ K-1 mol-1]


#TODO: Update the `meteo_utils` module.
#==============================================================================
# Farquhar Parameters
#==============================================================================

# Initial values (before accountinf for Na)
def par_photo_default(Vcm25= 89.0, Jm25 = 143.0, cRd = 0.008, TPU25 = 10.0,
       Kc25 = 404.9, Ko25 = 278.4, Tx25 = 42.75,
       alpha =  [0.2, 0.2, 0.19, 0.19, 0.14, 0.12], alpha_T_limit = [15, 20, 25, 30, 34, 50],
       a1 = 0.98, a2 = 0.98, a3 = 0.98, ds = 0.635, dHd = 200,
       RespT_Kc = {'model':'Ahrenius', 'c':38.05, 'deltaHa':79.43},
       RespT_Ko = {'model':'Ahrenius', 'c':20.30, 'deltaHa':36.38},
       RespT_Vcm = {'model':'Ahrenius', 'c':26.35, 'deltaHa':65.33},
       RespT_Jm = {'model':'Ahrenius', 'c':17.57, 'deltaHa':43.54},
       RespT_TPU = {'model':'Ahrenius', 'c':21.46, 'deltaHa':53.1},
       RespT_Rd = {'model':'Ahrenius', 'c':18.72, 'deltaHa':46.39},
       RespT_Tx = {'model':'Ahrenius', 'c':19.02, 'deltaHa':37.83}):
    """
    Generates a dictionary containing default **25 degrees** paramter values of the Farquhar's model for Vitis vinifera cv. Syrah.

   :Parameters:
   - **Vcm25**: Maximum RuBP-saturated rate of carboxylation [umol m-2 s-1]
   - **Jm25**: Maximum of electron transport [umol m-2 s-1]
   - **cRd**: the coefficient of mitochondrial respiration to `Vcm25`
   - **TPU25**: the rate of triose phosphate transport [umol m-2 s-1]
   - **Kc25**: Michaelis-Menten constant for the carboxylase [umol mol-1]
   - **Ko25**: Michaelis-Menten constant for the oxygenase [mmol mol-1]
   - **Tx25**: CO2 compensation point in the absence of mitochondrial respiration [umol mol-1]
   - **alpha**, **alpha_T_limit**, **a1**, **a2**, **a3**: parameters regulating electron transport for sunlit and shaded leaves (**not used**)
   - **c** and **deltaHa** are empirical parameters defining the temperature response curves of each of Kc, Ko, Vcm, Jm, TPU, Rd and Tx
   - **ds**: float, enthalpie of activation [KJ mol-1]
   - **dHd**: float, enthalpie of deactivation [KJ mol-1]
   """

    par_photodef = {}
    par_photodef['Vcm25']= Vcm25
    par_photodef['Jm25'] = Jm25
    par_photodef['Rd'] = cRd*Vcm25
    par_photodef['TPU25'] = TPU25
    par_photodef['Kc25'] = Kc25
    par_photodef['Ko25'] = Ko25
    par_photodef['Tx25'] = Tx25
    par_photodef['alpha'] = alpha
    par_photodef['alpha_T_limit'] = alpha_T_limit
    par_photodef['a1'] = a1
    par_photodef['a2'] = a2
    par_photodef['a3'] = a3
    par_photodef['ds'] = ds
    par_photodef['dHd'] = dHd
    par_photodef['RespT_Kc'] = RespT_Kc
    par_photodef['RespT_Ko'] = RespT_Ko
    par_photodef['RespT_Vcm'] = RespT_Vcm
    par_photodef['RespT_Jm'] = RespT_Jm
    par_photodef['RespT_TPU'] = RespT_TPU
    par_photodef['RespT_Rd'] = RespT_Rd
    par_photodef['RespT_Tx'] = RespT_Tx

    return par_photodef

def par_25_N_dict(Vcm25_N = (34.02, -3.13), Jm25_N = (78.27, -17.3),
             Rd_N = (0.42, -0.01), TPU25_N = (6.24, -1.92)):
    """
    Generates a dictionary containing the (slope, intercept) values of
    the linear relationship between photosynthetic capacity parameters
    (Vcmax, Jmax, TPU, Rd) and surface-based leaf Nitrogen content,
    according to Prieto et al. (2012, doi: 10.1111/j.1365-3040.2012.02491.x)

    :Parameters:
    - **The first value of each tuple**: float, the **slope** (parameter/Na) [umol g-1 s-1]
    - **The second value of each tuple**: float, the **intercept** [umol m-2 s-1]
    """
    par_dict = {
    'Vcm25_N' : Vcm25_N,
    'Jm25_N' : Jm25_N,
    'Rd_N' : Rd_N,
    'TPU25_N' : TPU25_N}

    return par_dict


def leaf_Na(ageTT, PPFD10, aN=-0.0008, bN=3.3, aM=6.471, bM=56.635):
    """
    Rerturns Nitrogen content per area (Na) [g m-2] according to Prieto et al. 2012 (doi: 10.1111/j.1365-3040.2012.02491.x).
    Deflaut parameter values are given for Montpellier Experiments on Syrah.

   :Parameters:
   - **ageTT**: float, cumulative degree-day temperature from budburst
   - **PPFD10**: float, cumulative intercepted irradiance (PPFD) over the past 10 days prior to simulation period
   - **aN** and **bN**: float, shape parameters for the linear relationship Nm=f(TT), (cf. Eq.5 in Prieto et al. 2012)
   - **aM** and **bM**: float, shape parameters for the linear relationship LMA=f(PPFD10), (cf. Eq.6 in Prieto et al. 2012)
    """

    LMA = aM * log(max(1.e-3,PPFD10)) + bM
    Nmass = aN*ageTT + bN
    Na = LMA*Nmass/100. #a cause d'erreur d'unite ds le papier

    return max(Na, 0.)
    # petit difference d'arrondi sur aN et bN?


#==============================================================================
# compute An
#==============================================================================

def Arrhenius_1(param,Tlc,par_photo):
    """
    Estimates the effect of temperature on the photosynthetic parameters `Tx`, `Kc`, `Ko` as described in Bernacchi et al. (2003)
    
    :Parameters:
    - **param**: string, one of 'Tx', 'Kc', 'Ko'
    - **Tlc**: float, leaf temperature [degreeC]
    - **par_photo**: a dictionart containing paramter default values of the Farquhar's model, at 25 degrees temperature
    """

    param_list = {
    'Tx':('Tx25','RespT_Tx'),
    'Kc':('Kc25','RespT_Kc'),
    'Ko':('Ko25','RespT_Ko')}

    Tak = mutils.Kelvin(Tlc)
#    indice1 = param_list[param][0]
    indice2 = param_list[param][1]
#    p25 = par_photo[indice1]
    c, deltaHa = [par_photo[indice2][x] for x in ('c','deltaHa')]
    p = exp(c - (deltaHa /(R*Tak)))

    return p


def Arrhenius_2(param,Tlc,par_photo):
    """
    Estimates the effect of temperature on the photosynthetic parameters `Vcmax`, `Jmax`, `TPUmax`, `Rdmax` as described in Bernacchi et al. (2003)

    :Parameters:
    - **param**: string, one of the following 'Vcmax', 'Jmax', 'TPUmax', 'Rdmax'
    - **Tlc**: float, leaf temperature [degreeC]
    - **par_photo**: a dictionart containing paramter values of the Farquhar's model
    """
    
    ds = par_photo['ds']
    dHd = par_photo['dHd']
    
    param_list = {
    'Vcmax':('Vcm25','RespT_Vcm'),
    'Jmax':('Jm25','RespT_Jm'),
    'TPUmax':('TPU25','RespT_TPU'),
    'Rdmax': ('Rd','RespT_Rd')}

    Tak = mutils.Kelvin(Tlc)
    indice1 = param_list[param][0]
    indice2 = param_list[param][1]
    p25 = par_photo[indice1]
    c, deltaHa = [par_photo[indice2][x] for x in ('c','deltaHa')]
    p = p25*(exp(c - (deltaHa /(R*Tak)))) / (1. + exp((ds*Tak-dHd)/(R*Tak)))

    return p


def dHd_sensibility(psi, Tleaf, dHd_max=200., dHd_min1=195.,dHd_min2=190.,
                    psi_max=-.75, psi_min=-2., Tleaf_1=35, Tleaf_2=40):
    """
    """
    dHd_min = dHd_min1 - (dHd_min1 - dHd_min2)*min(1,max(0,(Tleaf - Tleaf_1))/float(Tleaf_2-Tleaf_1))
    dHd = dHd_max - max(0, (dHd_max - dHd_min)*min(1, (psi - psi_max)/float(psi_min-psi_max)))


    return dHd


#def compute_an_2(par_photo, PPFD, Tlc, Ci):
#    """
#    Estimates net CO2 assimilation according to Farquhar et al. (1980).
#    Grapevine data are collected from Schultz FPB 2003 & Nikolov 1995
#
#    :Parameters:
#    - **par_photo**: a dictionary of Farquhar's model parameters
#    - **PPFD**: float, ?incident? PPFD [umol m-2 s-1]
#    - **Tlc**: float, leaf temperature [degreeC]
#    - **Ci**: float, intercellular CO2 concentration [umol mol-1]
#
#    :Returns:
#    - **An**: net CO2 assimilation rate [umol m-2 s-1]
#    """
#
#    T = Arrhenius_1('Tx',Tlc,par_photo) # CO2 compensation point in the absence of mitochondrial respiration [umol mol-1]
#    Kc = Arrhenius_1('Kc',Tlc,par_photo) # Michaelis-Menten constant for the carboxylase [umol mol-1]
#    Ko = Arrhenius_1('Ko',Tlc,par_photo) # Michaelis-Menten constant for the oxygenase [mmol mol-1]
#
#    Vcmax = Arrhenius_2('Vcmax',Tlc,par_photo) # Maximum RuBP-saturated rate of carboxylation [umol m-2 s-1]
#    Jmax = Arrhenius_2('Jmax',Tlc,par_photo) # Maximum of electron transport [umol m-2 s-1]
#    TPU = Arrhenius_2('TPUmax',Tlc,par_photo) # The rate of triose phosphate transport [umol m-2 s-1]
#    Rd = Arrhenius_2('Rdmax',Tlc,par_photo) # Mitochondrial respiration rate in the light [umol m-2 s-1]
#
#    alpha = par_photo['alpha'][0] # Leaf absorptance to photosynthetic photon flux [-]
#    for i in range(1,len(par_photo['alpha'])):
#        if par_photo['alpha_T_limit'][i-1]< Tlc < par_photo['alpha_T_limit'][i]:
#            alpha = par_photo['alpha'][i]
#
#    J = (alpha*PPFD)/((1+((alpha**2*PPFD**2)/(Jmax**2)))**0.5)  # Flux d'electrons en fonction de la temperature de la feuille et du niveau d'eclairement
#
#    # TODO: check ci*0.1013 formula correctness   
#    # The Ci values must be expressed in Pa, so multiplied by 0.1
#    Ac = (Vcmax*(Ci-T))/(Ci+Kc*(1+O/Ko))         # Rubisco-limited An (see Eq. 2 in Schultz 2003 for values)
#    Aj = (J*(Ci-T))/(4*(Ci+2.*T))                # RuBP-limited An (see Eq. 2 in Schultz 2003 for values)
#    Ap = (TPU*3)                                 # Export-limited An (see Sharkey, 2007 for values)
#
#    An = min([Ac, Aj, Ap]) - Rd
#    #An = Vc*(1-T/(Ci))-Rd
#
#    return An

def compute_an_2par(par_photo, PPFD, Tlc):
    """
    Getting photosynthetic variables for the anaytical solution of Evers et al. (2010).

    :Parameters:
    - **par_photo**: a dictionart containing paramter values of the Farquhar's model
    - **PPFD**: float, photosynthetic photon flux density [umol m-2 s-1]
    - **Tlc**: float, leaf temperature [degreeC]
    """

    T = Arrhenius_1('Tx',Tlc,par_photo)
    Kc = Arrhenius_1('Kc',Tlc,par_photo)
    Ko = Arrhenius_1('Ko',Tlc,par_photo)

    Vcmax = Arrhenius_2('Vcmax',Tlc,par_photo)
    Jmax = Arrhenius_2('Jmax',Tlc,par_photo)
    TPU = Arrhenius_2('TPUmax',Tlc,par_photo)
    Rd = Arrhenius_2('Rdmax',Tlc,par_photo)

#    alpha = par_photo['alpha'][0] # Leaf absorptance to photosynthetic photon flux [-]
#    for i in range(1,len(par_photo['alpha'])):
#        if par_photo['alpha_T_limit'][i-1]< Tlc < par_photo['alpha_T_limit'][i]:
#            alpha = par_photo['alpha'][i]
    alpha = .24

    J = (alpha*PPFD)/((1+((alpha**2*PPFD**2)/(Jmax**2)))**0.5)  # Flux d'electrons en fonction de la temperature de la feuille et du niveau d'eclairement

    x1c = Vcmax
    x2c = Kc*(1+O/Ko)

    x1j = J/4.
    x2j = 2.*T

    x1t = 3*TPU
    x2t = -T

    return x1c,x2c,x1j,x2j,x1t,x2t,Rd


#==============================================================================
# Numerical solution of the An - gs - ci variables
#==============================================================================

#def incrementCi(Ca, An, gs, gb): # 22/06/2011: sacar el meteo_data, ea, es_l del parentesis y del parentesis de la linea 44
#    """
#    Calculates the chloroplast CO2 partial pressure issued from Fick's first law.
#    
#    :Parameters:
#    - **Ca**: air CO2 concentration [umol mol-1]
#    - **An**: net CO2 assimilation rate [umol m-2 s-1]
#    - **gs**: stomatal conductance for water [mol m-2 s-1]
#    - **gb**: boundary layer conductance for water [mol m-2 s-1].
#    
#    :Note:
#    gs and gb are converted to CO2 conductances via division by 1.6 and 1.37, respectively.
#    
#    """
#    return Ca-An*(1.6/gs+1.37/gb)


#def coupling_Anci_iter(par_photo, meteo_leaf, psi, Ci, Tlc, w=0.1, iterCi=50, deltaci=0.0001):
#    """
#    Estimates simultaneously the values of net CO2 assimilation rate (An) and
#    intercellular CO2 concentration (Ci), by an iterative scheme.
#    
#    :Warning:
#    This function does not work under water deficit.
#    Use instead :func:`compute_an_analytic`.
#    
#    :Parameters:
#    - **par_photo**: a dictionary of Farquhar's model parameters
#    - **meteo_leaf**: list of meteorological data
#    - **w**: float, leaf Characteristic dimension in relation to wind speed [m]
#    - **iterCi**: integer, the maximum number of desired iterations
#    - **deltaci**: float, the tolerance in error between two consecutive values of `Ci` [umol mol-1]
#    - **Tlc**: float, leaf temperature [degreeC].
#    
#    :Returns:
#    - **An**: net CO2 assimilation rate [umol m-2 s-1]
#    - **Cinew**: intercellular CO2 concentration [umol mol-1]
#    - **gs**, **gb**: floats, resp. stomatal and boundary layer conductances for water [mol m2 s-1]
#    """
#    Tac = meteo_leaf['Tac']
#    PPFD = meteo_leaf['PPFD'] 
#    hs = meteo_leaf['hs'] 
#    Ca = meteo_leaf['Ca']
#    
#    VPD = mutils.VPD_leaf_air(Tac,Tlc,hs)
#    gb = mutils.computeBoundaryLayerConductance(meteo_leaf['u'], w)
#    
#    i = 0
#    while i<iterCi :
#        An = compute_an_2(par_photo, PPFD, Tlc, Ci)
#
#        Cs = Ca-An*(1.37/gb)
#        x1c,x2c,x1j,x2j,x1t,x2t,Rd = compute_an_2par(par_photo, PPFD, Tlc)
#        Tx = x2j/2.
#        gs = BWB_gs.gs_Leuning(An,VPD,Cs,Tx,psi)
#
#        Cinew = Ca-An*(1.6/gs+1.37/gb)
#        
#        if abs(Cinew-Ci) < deltaci:
#            #print 'nb iteration Ci_'+ str(i)
#            #Ci = Cinew
#            break 
#        else:
#            Ci = Cinew
#            if i>iterCi-2:
#                print 'warning ! Ci calculation does not converge to a solution'   
#        i=i+1
#
#    return An, Cinew, gs, gb


#==============================================================================
# Analytical solution of the An - gs - ci processes
#==============================================================================
def fvpd_3(model, vpd, psi, psi_crit=-0.37, c1=5.278, c2=1.85, D0=30.):
    """
    """
    if model =='misson':
        m = 1. / (1. + (psi/psi_crit)** c2)
    elif model == 'tuzet':
        m = (1. + exp(c2 * psi_crit)) / (1. + exp(float(c2 * (psi_crit - psi))))
    elif model == 'linear':
        m = 1 - min(1,float(psi/psi_crit))
    elif model == 'vpd':
        m = 1. / (1. + vpd / float(D0))
    else:
        raise ValueError ("The 'model' argument must be one of the following ('misson','tuzet', 'linear').")
    return c1 * m


def g0_sensibility(psi, psi_crit=-0.1, n=3):
    """
    """
    return 1. / (1. + (psi/psi_crit)**n)


def gm(temp, gm25=0.1025, EA_gm=49600., DEA_gm=437400., Sgm=1400.):
    """
    Returns gm, the mesophyll conductance to CO2 diffusion [umol m-2 s-1 ubar-1].

    :Parameters:
    - **temp**: float, leaf temperature [degreeC]
    - **gm25**: Mesophyll conductance for CO2 diffusion [mol m-2 leaf s-1 bar-1]. Default value from (Bernacci 2002 for tobaco)
    - **EA_gm**: the activation energy for parameter X [J mol-1]; Yin 2009; Bernacci 2002 plant physiology
    - **DEA_gm**: the deactivation energy, resource?
    - **Sgm**: an entropy term [J K-1 mol -1]
    """

    R2 = R * 1.e3 #Ideal gaz constant in [J K-1 mol-1]
    gm = gm25*exp((1./298.-1./(temp+273.))*EA_gm/R2)*(1.+exp(Sgm/R2-DEA_gm/298./R2))/(1.+exp(Sgm/R2-1./(temp+273.)*DEA_gm/R2))

    return gm


def compute_amono_analytic(x1, x2, temp, vpd, gammax, Rd, psi, model='misson',
                           g0=0.019, rbt=2./3., ca=400., m0=5.278, psi0=-0.1,
                           D0=30., n=1.85):
    """
    Returns Amono, the gross CO2 assimilation rate [umolCO2 m-2 s-1] according to the analytical solution of Evers et al. (2010) jxb.
    
    A is called 'mono' since it is calculated for either of RuBisco- or e transport-limited phases. The net A (An) is calculated using the function 'compute_an_analytic'
    
    :Parameters:
    - **x1** and **x2**: float, are Farquhar's An parameters (x1= Vcmax or J/4.; x2= KmC(1+O/KmO) or 2 gamma_asterisk)
    - **temp**: float, leaf temperature [degreeC]
    - **vpd**: float, vapor pressure deficit [kPa]
    - **gammax**: float, CO2 compensation point [ubar]
    - **Rd**: float, the mitochondrial respiration rate in the light [umol m-2 s-1]
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters
    """

    f_VPD = fvpd_3(model, vpd, psi, psi_crit=psi0, c1=m0, c2=n, D0=D0)
#    if psi < -1.5: f_VPD = 0

    cube_a = g0*(x2+gammax)+(g0/gm(temp)+f_VPD)*(x1-Rd)
    cube_b = mutils.cmol2cpa(temp,ca)*(x1-Rd)-gammax*x1-Rd*x2
    cube_c = mutils.cmol2cpa(temp,ca)+x2+(1./gm(temp)+rbt)*(x1-Rd)
    cube_d = x2+gammax+(x1-Rd)/gm(temp)
    cube_e = 1./gm(temp)+(g0/gm(temp)+f_VPD)*(1./gm(temp)+rbt)
    cube_p = -(cube_d+(x1-Rd)/gm(temp)+cube_a*(1./gm(temp)+rbt)+(g0/gm(temp)+f_VPD)*cube_c)/cube_e
    cube_q = (cube_d*(x1-Rd)+cube_a*cube_c+(g0/gm(temp)+f_VPD)*cube_b)/cube_e
    cube_r = -(cube_a*cube_b/cube_e)
    cube_Q = (cube_p**2.-3.*cube_q)/9.
    cube_R = (2.*cube_p**3.-9.*cube_p*cube_q+27.*cube_r)/54.
    cube_xi = arccos(max(-1.,min(1., cube_R/sqrt(cube_Q**3.))))

    Amono = -2.*sqrt(cube_Q)*cos(cube_xi/3.)-cube_p/3.

    return Amono


def compute_an_analytic(temp, vpd, x1c, x2c, x1j, x2j, x1t, x2t, Rd, psi,
                        model='misson', g0=0.019, rbt=2./3., ca=400., m0=5.278,
                        psi0=-0.1, D0=30., n=1.85):
    """
    Returns An, the net CO2 assimilation rate [umolCO2 m-2 s-1] according to the analytical solution of Evers et al. (2010) jxb.
    
    :Parameters:
    - **temp**: float, leaf temperature [degreeC]
    - **vpd**: float, vapor pressure deficit [kPa]
    - **x1c** and **x2c** are Rubisco-limiting CO2 assimilation of Farquhar model (Vcmax and Kc*(1+O/Ko), resp.)
    - **x1j** and **x2j** are electron transport-limiting CO2 assimilation of Farquhar model (J/4 and 2*gammax, resp.)
    - **x1t** and **x1t** are triose phosphate export rates of Farquhar model (3*TPU and -gamma, resp.)
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters
    
    :Return:
    - **An**: CO2 net assimilation rate [umolCO2 m-2 s-1]
    - **Cc**: CO2 chloroplast partial pressure [ubar]
    - **Ci**: intercullar CO2 partial pressure [ubar]
    - **gsw**: stomatal conductance for water [mol m-2 s-1]
    """

    f_VPD = fvpd_3(model, vpd, psi, psi_crit=psi0, c1=m0, c2=n, D0=D0)

    gammax = x2j/2.
    ci_asterisk = gammax-Rd/gm(temp)

    Ac = compute_amono_analytic(x1c,x2c,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    Aj = compute_amono_analytic(x1j,x2j,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    At = compute_amono_analytic(x1t,x2t,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    
    An = min(Ac,Aj,At)

    # chlorophyll partial pressure [ubar]
    if min(Ac,Aj,At) == Ac:
        CC=(gammax*x1c+(Ac+Rd)*x2c)/(x1c-Ac-Rd)
    elif min(Ac,Aj,At) == Aj:
        CC=(gammax*x1j+(Aj+Rd)*x2j)/(x1j-Aj-Rd)
    else:
        CC=(gammax*x1t+(At+Rd)*x2t)/(x1t-At-Rd)

    # Intercellular partial pressure [ubar]
    CI = CC + An/gm(temp);  #intercelluar co2 concentration

    # stomatal conductance for CO2 [umol m-2 s-1 ubar-1]
    GS =g0/1.6+(An + Rd)*f_VPD/(CI-ci_asterisk)
    # stomatal conductance for water [umol m-2 s-1 ubar-1]
    GSW = GS*1.6

    return An, mutils.cpa2cmol(temp,CC), mutils.cpa2cmol(temp,CI), GSW


def an_gs_ci(par_photo, meteo_leaf, psi, Tlc, model='misson', g0=0.019,rbt=2./3.,
           ca=400., m0=5.278, psi0=-0.1, D0=30., n=1.85):
    """
    Estimates simultaneously the values of net CO2 assimilation rate (An) and
    intercellular CO2 concentration (Ci), by an analytic scheme.

    :Parameters:
    - **par_photo**: a dictionary of Farquhar's model parameters
    - **meteo_leaf**: list of local (leaf-scale) meteorological data
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **Tlc**: float, leaf temperature [degreeC]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters

    :Returns:
    - **An**: net CO2 assimilation rate [umol m-2 s-1]
    - **Cc**: chloroplast CO2 concentration [umol mol-1]
    - **Cc**: intercellular CO2 concentration [umol mol-1]
    - **gs**: stomatal conductances to water [mol m-2 s-1]
    """

    Tac = meteo_leaf['Tac']
    PPFD = meteo_leaf['PPFD'] 
    hs = meteo_leaf['hs']
    
    PPFD = max(1.e-6, PPFD) # To avoid numerical instability at PPFD=0 (Eq. S9 from Evers et al., 2010 JxBot,61:2203–2216)

    VPD = mutils.VPD_leaf_air(Tac,Tlc,hs)

    x1c,x2c,x1j,x2j,x1t,x2t,Rd = compute_an_2par(par_photo, PPFD, Tlc)

    An, Cc, Ci, gs = compute_an_analytic(Tlc,VPD,x1c,x2c,x1j,x2j,x1t,x2t,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)

    return An, Cc, Ci, gs


def transpiration_rate(Tlc, ea, gs, gb, Pa = 101.3):
    """
    Returns E, leaf transpiration in [mol m-2 s-1]
    
    :Parameters:
    - **Tlc**: float, leaf temperature [degreeC]
    - **ea**: float, air vapor pressure [kPa]
    - **gs**: stomatal conductance for water [mol m-2 s-1]
    - **gb**: float, boundary layer conductance for H2O [mol m-2 s-1]
    - **Pa**: float, atmospheric pressure [kPa].
    
    :Worning:
    gb is calculated for both sides of a flat leaf.
    """

#    gva = gb*1.37 # boundary layer conductance for water vapour transport [mol m2 s-1] # R: ancienne formule : gb*1.4
    gv = 1./((2./gb)+(1./gs)) # Formulation by Pearcy 1989
    es_l = mutils.s_avpd(Tlc)  # Saturated vapor pressure at leaf surface [kPa]
    E = gv*((es_l-ea)/Pa) # Transpiration rate [mol m-2 s-1]

    return E


def VineExchange(g, par_photo, par_photo_N, par_gs, meteo, E_type2, leaf_lbl_prefix='L',
                 rbt=2./3.):
    """
    Calculates gas exchange fluxes at the leaf scale according to the analytical scheme described by Evers et al. (JxBot 2010, 2203–2216).

    :Parameters:
    - **g**: a multiscale tree graph object
    - **par_photo**: dictionary, the parameters of the Farquhar model for net CO2 assimilation (cd :func:`par_photo_default`)
    - **par_gs**: dictionary, the parameters of the stomatal conductance model (model, g0, m0, psi0, D0, n)
    - **meteo**: a `pandas.DataFrame`-like object
    - **E_type2**: string, one of two 'Ei' or 'Eabs'
    - **leaf_lbl_prefix**: string, the prefix of the label of the leaves
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    
    
    :Attaches to each leaf:
    - **An**: the net CO2 assimilation [umol m-2 s-1]
    - **Ci**: intercellular CO2 concentration [umol mol]
    - **gs**: stomatal conductance to water vapor [mol m-2 s-1]
    - **gb**: boundary layer conductance to water vapor [mol m-2 s-1]
    - **E** leaf transpiration [mol m-2 s-1].
    """
    
    model,g0max,m0,psi0,D0,n = [par_gs[ikey] for ikey in ('model','g0','m0','psi0','D0','n')]
    
    meteo_leaf = deepcopy(meteo)
    meteo_leaf = meteo_leaf.iloc[0]
    
    for vid in g:
        if vid>0:
            node=g.node(vid)
            if node.label.startswith(leaf_lbl_prefix):
                Tac = meteo_leaf.Tac
                hs = meteo_leaf.hs
                u = meteo_leaf.u
                Ca = meteo_leaf.Ca
                Pa = meteo_leaf.Pa
                
                try:
                    u *= node.u_coef
                    node.u = u
                except:
                    pass

                psi = node.properties()['psi_head'] #leaf water potential [MPa] (assumed equal to that of the petiole)
                Tlc = node.properties()['Tlc']
                PPFD_leaf = node.properties()[E_type2]
#                w = node.Length / 100. # leaf length in m                

                meteo_leaf['PPFD'] = PPFD_leaf
                meteo_leaf['Rg'] = PPFD_leaf/(0.48*4.6)
#                meteo_leaf['u'] = min(1., meteo_leaf['u']) # Hack: see Nobel p.338
                
#                if not 'par_photo' in node.properties():
                leaf_par_photo = deepcopy(par_photo)
                leaf_par_photo['Vcm25'] = par_photo_N['Vcm25_N'][0]*node.Na+par_photo_N['Vcm25_N'][1]
                leaf_par_photo['Jm25'] = par_photo_N['Jm25_N'][0]*node.Na+par_photo_N['Jm25_N'][1]
                leaf_par_photo['TPU25'] = par_photo_N['TPU25_N'][0]*node.Na+par_photo_N['TPU25_N'][1]
                leaf_par_photo['Rd'] = par_photo_N['Rd_N'][0]*node.Na+par_photo_N['Rd_N'][1]
                dHd_max = leaf_par_photo['dHd']
                dHd = dHd_sensibility(psi, Tlc, dHd_max=dHd_max, dHd_min1=195.,dHd_min2=190.,
                    psi_max=-.75, psi_min=-2., Tleaf_1=32, Tleaf_2=33)

                leaf_par_photo['dHd'] = dHd
                node.par_photo = leaf_par_photo
                
                g0 = g0max#*g0_sensibility(psi, psi_crit=-1, n=4)

                An, Cc, Ci, gs = an_gs_ci(node.par_photo, meteo_leaf, psi, Tlc,
                                           model, g0, rbt, Ca, m0, psi0, D0, n)

                # Boundary layer conductance
#                gb = mutils.computeBoundaryLayerConductance(u, w)
                l_w = node.Length/100.*0.72 # leaf length in the downwind direction [m]
                d_bl = 4.*(l_w/max(1.e-3,u))**0.5 /1000. # Boundary layer thikness in [m] (Nobel, 2009 pp.337)
                Dj0 = 2.13*1.e-5 #[m2 s-1] at P=1. atm and t=0. °C (Nobel, pp.545)
                Dj = Dj0*(101.3/Pa)*((Tac+273.15)/273.15)**1.8 #(Nobel, eq.8.8, pp.379)
                gb = Dj*(Pa * 1.e-3) / ((R*1.e-3) * (Tac+273.15) * d_bl) # [mol m-2 s-1] (Nobel, 2009 pp.337)

                # Transpiration
                es_a = mutils.s_avpd(Tac)
                ea = es_a*hs/100.
                E = transpiration_rate(Tlc, ea, gs, gb, Pa)

                node.An = An
                node.Ci = Ci
                node.gs = gs
                node.gb = gb
                node.E = max(0.,E)

    return