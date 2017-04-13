# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Hydraulic strcture module of HydroShoot.

This module computes xylem water potential value at each node of the shoot.
"""

from scipy import exp, absolute, pi, log, optimize
#from sympy.solvers.solveset import solveset_real
#from sympy import Symbol
from copy import deepcopy
import warnings

from openalea.plantgl.all import surface as surf
import openalea.mtg.traversal as traversal

# Constants
rho=997.0479    # Density of liquid water at 25 degreC temperature [kg m-3]
g_p=9.81    #Gravity acceleration [m s-2]


def k_max(diameter, a=2.8, b=0.1, min_kmax=0.):
    """
    Returns maximim stem conductivity according to te allometric relationship (Kh=a*D^b) proposed by Tyree and Zimermmann (2002, p.147)

    :Parameters:
    - **diameter**: bulk diameter of the sapwood [m]
    - **a** and **b**: respectively the slope and power parameters for the Kh(D) relationship
    - **min_kmax**: float, minimum value for the maximum conductivity [kg s-1 m MPa-1]
    
    :Notice:
    - Indicative values ranges are [2.5,2.8] for a and [2.0,5.0] for b according to Tyree and Zimermmann (2002, p.147)
    """

    return max(min_kmax, float(a*(diameter**b)))


def k_vulnerability(psi, model='tuzet', fifty_cent=-0.51,sig_slope=3):
    """
    Returns the ratio of actual to maximum stem conductance (between 1 and 0).

    :Parameters:
    - **psi**: stem water potential [MPa]
    - **model**: one of **'misson'** (logistic function with polynomial formula), **'tuzet'** (logistic function with exponential formula), or **'linear'** for linear reduction
    - **fifty_cent**: water potential at which the conductivity of the stem is reduced by 50% reltive to its maximum value [MPa]
    - **sig_slope**: a shape parameter controling the slope of the S-curve.
    """

    if model == 'misson':
        k_reduction = 1. / (1. + (psi/fifty_cent)** sig_slope)
    elif model == 'tuzet':
        k_reduction = (1. + exp(sig_slope * fifty_cent)) / (1. + exp(float(sig_slope * (fifty_cent - psi))))
    elif model == 'linear':
        k_reduction = 1 - min(0.95,float(psi/fifty_cent))
    else:
        raise ValueError ("The 'model' argument must be one of the following ('misson','tuzet', 'linear').")

    return k_reduction


def def_param_soil():
    """
    Returns a dictionary of classes of default soil hydrodynamic parameters for the
    van Genuchten-Muallem model.
    For each soil class (`dict.key`), the data are organized as follows:
    name : (theta_r, theta_s, alpha[cm-1], n, k_sat[cm d-1])
    
    :Ref:
    - Carsel and Parrish (1988) WRR 24,755–769 DOI: 10.1029/WR024i005p00755
    """

    def_dict = {'Sand': (0.045, 0.430, 0.145, 2.68, 712.8),
                'Loamy_Sand': (0.057, 0.410, 0.124, 2.28, 350.2),
                'Sandy_Loam': (0.065, 0.410, 0.075, 1.89, 106.1),
                'Loam': (0.078, 0.430, 0.036, 1.56, 24.96),
                'Silt': (0.034, 0.460, 0.016, 1.37, 6.00),
                'Silty_Loam': (0.067, 0.450, 0.020, 1.41, 10.80),
                'Sandy_Clay_Loam': (0.100, 0.390, 0.059, 1.48, 31.44),
                'Clay_Loam': (0.095, 0.410, 0.019, 1.31, 6.24),
                'Silty_Clay_Loam': (0.089, 0.430, 0.010, 1.23, 1.68),
                'Sandy_Clay': (0.100, 0.380, 0.027, 1.23, 2.88),
                'Silty_Clay': (0.070, 0.360, 0.005, 1.09, 0.48),
                'Clay': (0.068, 0.380, 0.008, 1.09, 4.80),
                'Costum': (0.065, 0.410, 0.005, 2., 106.1)}

    return def_dict


def k_soil_soil(psi, soil_class):
    """
    Returns actual soil hydraulic conductivity [cm d-1].
    
    :Parameters:
    - **psi**: float, bulk soil-matrix water potential [MPa]
    - **soil_class**: string, the name of the soil hydrodynamic class as proposed by Carsel and Parrish (1988) DOI: 10.1029/WR024i005p00755
    """

    psi *= 1.e6/(rho*g_p)*100. # conversion MPa to cm_H20
    param = def_param_soil()[soil_class]
    theta_r,theta_s,alpha,n,k_sat = [param[i] for i in range(5)]
    m = 1. - 1./n
#    p = 2.

#    k_soil = k_sat * (1./(1.+(alpha*(-psi))**n))**(p-p/n) * (1. - (1. - (1. / (1. + (alpha * (-psi))**n)))**m)**2
    
    se = 1. / (1. + abs(alpha * psi)**n)**m
    k_soil = k_sat * (se)**0.5 * (1. - (1. - se**(1./m))**m)**2

    return k_soil


def k_soil_root(k_soil, d, r):
    """
    Returns the hydraulic water conductivity at the soil-root interface according to Gardner (1960 Soil Sci 89, 63-73) [cm d-1].

    :Parameters:
    - **k_soil**: float, soil-matrix hydraulic conductivity [cm d-1]
    - **d**: float, mean distance between the neighbouring roots [m]
    - **r**: float, mean root radius [m]
    """

    return 4. * pi * k_soil / log((d/r)**2)


def soil_water_potential(psi_soil_init, flux, soil_class, soil_total_volume,
                         psi_min = -3.):
    """
    Estimates soil water potential based on Brooks and Corey (1964) van Genuchten water retention curves.

    :Parameters:
    - **psi_soil_init**: Initial soil warer potential [MPa]
    - **flux**: transpired flux [Kg T-1]
    - **soil_class**: one of the following 'Sand','Loamy Sand','Sandy Loam','Loam','Silt','Silty_Loam','Sandy_Clay_Loam','Clay_Loam','Silty_Clay_Loam','Sandy_Clay','Silty_Clay','Clay'
    - **intra_dist**,**inter_dist**,**depth**: intra-row spacing, inter-row spacing, roots depth, all in [m]

#    TODO: implement a true Richards solution
    """

    param = def_param_soil()[soil_class]
    theta_r,theta_s,alpha,n,k_sat = [param[i] for i in range(5)]
    m = 1. - 1./n

    psi = psi_soil_init*1.e6/(rho*g_p)*100. # conversion MPa to cm_H20
    theta_init = theta_r + (theta_s-theta_r) / (1. + absolute(alpha*psi)**n)**m

    F = flux * 1.e-3 # kg T-1 to m3 T-1

    # Rough estimation of total porosity volume [m3]
    porosity_volume = soil_total_volume * theta_s

    d_theta = F / porosity_volume # [m3 m-3]

    theta = max(theta_r, theta_init - d_theta)

    if theta == theta_r:
        psi_soil = psi_min
    else:
        eq = lambda psi : 1. / (1. + absolute(alpha*psi)**n)**m - (theta-theta_r)/(theta_s-theta_r)
        psi_soil = optimize.fsolve(eq, psi)[0]/(1.e6/(rho*g_p)*100)
#        psiX = Symbol('psi')
#        eq = 1. / (1. + absolute(alpha*psi)**n)**m - (theta-theta_r)/(theta_s-theta_r)
#        psi_soil = (solveset_real(eq, psiX).inf)/(1.e6/(rho*g_p)*100)
#        if not psi_soil.is_Float:
#            psi_soil.removeO()

    return float(psi_soil)


def hydraulic_prop(mtg, vtx_label='inT', MassConv=18.01528, LengthConv=1.e-2,
                   a=2.6, b=2.0, min_kmax=0.):
    """
    Returns water flux, `F` [kg s-1] and maximum stem conductivity [kg m s-1 Pa-1] of each internode.
    
    :Attention:
    1. The units of **F** and **K** are related:
     - if `F` is given as water flux [kg s-1], then `K` must be given as [kg m s-1 Pa-1] (or [kg m s-1 MPa-1])
     - else if `F` is given as water flux density [kg m-2 s-1], then `K` must be given as [kg m-1 s-1 Pa-1] (or [kg m-1 s-1 MPa-1]).
    2. Transpiration flux density per leaf surface area `E` must be in [mol m-2 s-1], otherwise, the `MassConv` value must be readapted.
    
     The resulting water portential, calculated in :func:`transient_xylem_water_potential` is then given in [MPa].
     
    :Parameters:
    - **mtg**: an MTG object
    - **vtx_label**: string, the label prefix of the basal vertex at highest scale
    - **MassConv**: molar mass of H2O [gr mol-1]
    - **LengthConv**: conversion coefficient from the length unit of the mtg to that of 1 m
    - **a** and **b**: respectively the slope and power parameters for the Kh(D) relationship
    - **min_kmax**: float, minimum value for the maximum conductivity [kg s-1 m MPa-1]
    """

    # Getting the basal vertex of the highest scale vertices
    vid_base = mtg.node(mtg.root).vid_base

    for vtx_id in traversal.post_order2(mtg,vid_base):
        n = mtg.node(vtx_id)
        if n.label.startswith('LI'):
            try:
                leaf_area = n.leaf_area*1.
            except:
                leaf_area = surf(n.geometry)*LengthConv**2    #[m2]
#                leaf_area = (0.0175*(n.Length*10.)**1.9057)*LengthConv**2 #[m2]
                n.leaf_area = leaf_area
# Note: The surface of the leaf mesh is overestimated compared to allometry results
            n.Flux = ((n.E)*MassConv*1.e-3)*leaf_area
#            n.FluxC = ((n.An)*44.0095*1.e-9)*leaf_area # [kgCO2 s-1]
            n.FluxC = (n.An)*leaf_area # [umol s-1]

        elif n.label.startswith(('in','cx','Pet')):
            n.Flux = sum([vtx.Flux for vtx in n.children()])
            Diam = 0.5*(n.TopDiameter + n.BotDiameter)*LengthConv
            n.Kmax = k_max(Diam,a,b,min_kmax)

            n.FluxC = sum([vtx.FluxC for vtx in n.children()])

        elif n.label.startswith('rhyzo'):
            n.Flux = sum([vtx.Flux for vtx in n.children()])
            n.FluxC = sum([vtx.FluxC for vtx in n.children()])
            
            n.Kmax = None
            

    return mtg


def transient_xylem_water_potential(g, model='tuzet', vtx_label='inT',
                                    LengthConv=1.e-2,psi_soil=-0.6,psi_min=-3.,
                                    fifty_cent=-0.51, sig_slope=1.,
                                    dist_roots=0.013, rad_roots=.0001,
                                    negligible_shoot_resistance = False,
                                    start_vid = None, stop_vid = None):
    """
    Returns a transient hydraulic structure of a plant shoot based on constant values of stem conductivities.
    Stems are assumed isotropic, linear, element.
    
    :Attention:
    1. the flux through the stem must have been given in [kg s-1]
    2. the hydraulic condutctivity of the stem must have been given in [kg m s-1 MPa-1]

    :Parameters:
    - **g**: an MTG object
    - **model**: either 'misson' (logistic function with polynomial formula) or 'tuzet' (logistic function with exponential formula)
    - **vtx_label**: string, the label prefix of the basal vertex at highest scale
    - **LengthConv**: conversion coefficient from the length unit used for building the mtg to 1 m
    - **psi_soil**: soil water potential [MPa]
    - **psi_min**: minimum simulated xylem water potential [MPa]
    - **fifty_cent**: water potential at which the conductivity of the stem is reduced by 50% reltive to its maximum value
    - **sig_slope**: a shape parameter controling the slope of the S-curve
    - **dist_roots**: float, the average distance between roots [m]
    - **rad_roots**: float, the average radius of individual roots [m]
    - **negligible_shoot_resistance**: logical, defalut `False`, whether shoot resistance should be considered (False) or not (True)
    - **start_vid**: `None` or integer, vertex id from which the iteration starts (if `None` it is then taken to the basal element)
    - **stop_vid**: `None` or integer, vertex id at which the iteration breaks (def `None` iteration will include all elements up until the leaves)
    """

    vid_base = g.node(g.root).vid_base
#    vid_collar = g.node(g.root).vid_collar

    if start_vid is None:
        start_vid = vid_base

    for vtx_id in traversal.pre_order2(g, start_vid):
        if vtx_id == stop_vid:
            break
        else:
            n = g.node(vtx_id)
            p = n.parent()
    
            if n.label.startswith('LI'):
                n.psi_head = p.psi_head
            else:
                F = n.properties()['Flux']
                L = n.properties()['Length']*LengthConv
                Z_head = n.properties()['TopPosition'][2]*LengthConv
                Z_base = n.properties()['BotPosition'][2]*LengthConv

    
                psi_base = psi_soil if vtx_id == vid_base else p.psi_head
    
                try:
                    psi_head = n.psi_head
                    assert (n.psi_head != None), "Water potential has a `None` value"
                except:
                    psi_head = psi_base
    
                psi = 0.5 * (psi_head + psi_base)

                if n.label.startswith('rhyzo'):
#                    # Hack
#                    k_sat = g.node(vid_collar).Kmax
#                    KL = k_sat * k_vulnerability(psi, 'misson', -.2, 5)
#                    psi_head = max(psi_min, psi_base - L*F/KL)
                    cyl_diameter = n.TopDiameter * LengthConv
                    depth = n.depth * LengthConv
                    F *= 8640./(pi *cyl_diameter* depth) # [cm d-1]

                    if n.label.startswith('rhyzo0'):
                        k_soil = k_soil_soil(psi, n.soil_class) # [cm d-1]
                        gL = k_soil_root(k_soil, dist_roots, rad_roots) #[cm d-1 m-1]
                        psi_head = max(psi_min, psi_base - (F/gL) * rho*g_p*1.e-6)
                    else:
                        KL = k_soil_soil(psi, n.soil_class) # [cm d-1]
                        psi_head = max(psi_min, psi_base - (L*F/KL) * rho*g_p*1.e-6)

                else:
                    Kmax = n.properties()['Kmax']

                    if not negligible_shoot_resistance:
                        KL = Kmax *k_vulnerability(psi, model, fifty_cent,sig_slope)
                        psi_head = max(psi_min, psi_base - L*F/KL - (rho*g_p*(Z_head-Z_base))*1.e-6)
                    else:
                        KL = None
                        psi_head = max(psi_min, psi_base - (rho*g_p*(Z_head-Z_base))*1.e-6)
#                def _psi_headX(_psi_head):
#                    psi = 0.5*(_psi_head+psi_base)
#                    _KL = Kmax*k_vulnerability(psi, model, fifty_cent,sig_slope)
#                    eq_error = _psi_head - (psi_base - L*F/_KL - (rho*g_p*(Z_head-Z_base))*1.e-6)
#                    return eq_error
#    
#                psi_head = optimize.newton_krylov(_psi_headX,psi_base).tolist()
    
    
                n.psi_head = psi_head
                n.KL = KL

    return


def xylem_water_potential(g, psi_soil=-0.8, model='tuzet', psi_min=-3.0,
                          psi_error_crit=0.001, max_iter=100, vtx_label='inT',
                          LengthConv=1.E-2, fifty_cent=-0.51,sig_slope=0.1,
                          dist_roots=0.013, rad_roots=.0001,
                          negligible_shoot_resistance = False,
                          start_vid = None, stop_vid = None,
                          psi_step = 0.5):
    """
    Returns the hydraulic structure of a plant shoot. Stems are assumed isotropic, linear, element.

    :Parameters:
    - **g**: an MTG object
    - **psi_soil**: soil water potential [MPa]
    - **model**: either 'misson' (logistic function with polynomial formula) or 'tuzet' (logistic function with exponential formula)
    - **psi_min**: minimum simulated xylem water potential [MPa]
    - **psi_error_crit**: maximum accepted error in water potential [MPa]
    - **max_iter**: maximum allowed iteration
    - **vtx_label**: string, the label prefix of the basal vertex at highest scale
    - **LengthConv**: conversion coefficient from the length unit used for building the mtg to 1 m
    - **fifty_cent** and **sig_slope**: Shape parameters of the vulnerability curve, see :func:`Kstem`
    - **start_vid**: `None` or integer, vertex id from which the iteration starts (if `None` it is then taken to the basal element)
    - **stop_vid**: `None` or integer, vertex id at which the iteration breaks (def `None` iteration will include all elements up until the leaves)
    """

    counter = 0
    
    psi_error = psi_error_crit
    psi_error_trace = []
    ipsi_step = psi_step

    while psi_error >= psi_error_crit:
        psi_error = 0.

        psi_prev = deepcopy(g.property('psi_head'))

        transient_xylem_water_potential(g, model, vtx_label, LengthConv,
                                        psi_soil,psi_min,fifty_cent,sig_slope,
                                        dist_roots,rad_roots,
                                        negligible_shoot_resistance,
                                        start_vid, stop_vid)

        psi_new = deepcopy(g.property('psi_head'))

        psi_error_trace.append(psi_error)

        if counter > max_iter:
            warnings.warn("The numerical solution of the hydraulic architecture did not converge.")
#            break
        else:
            try:
                if abs(psi_error_trace[-1] - psi_error_trace[-2]) < psi_error_crit:
                    ipsi_step = max(0.01, ipsi_step/2.)
            except:
                pass

            psi_avg_dict = {}
    
            for vtx_id in g.property('psi_head').keys():
                psi_error += abs(psi_prev[vtx_id]-psi_new[vtx_id])
                psi_avg_dict[vtx_id] = psi_prev[vtx_id] + psi_step * (psi_new[vtx_id] - psi_prev[vtx_id])
            
            g.properties()['psi_head'] = psi_avg_dict


        counter += 1

    return counter


#def soil_rhyzo_water_potential(g, psi_collar,LengthConv, dist_roots, rad_roots):
#    """
#    """
#
#    vid_collar = g.node(g.root).vid_collar
#    vid_base = g.node(g.root).vid_base
#
#    psi_soil = psi_collar
#
#    for vid in g.Ancestors(vid_collar)[1:]:
#        n = g.node(vid)
#        p = n.parent()
#
#        cyl_diameter = n.TopDiameter * LengthConv
#        depth = n.depth * LengthConv
#
#        F = n.Flux * 8640./(pi *cyl_diameter* depth) # [cm d-1]
#        L = n.Length * LengthConv
#
#        psi_head = psi_collar if n == g.node(vid_collar).parent() else n.psi_head
#        psi_base = psi_soil if vid == vid_base else p.psi_head
#
#        psi = 0.5 * (psi_head + psi_base)
#
#        if n.label.startswith('rhyzo0'):
#            k_soil = k_soil_soil(psi, n.soil_class) # [cm d-1]
#            gL = k_soil_root(k_soil, dist_roots, rad_roots) #[cm d-1 m-1]
#            psi_base = min(-0.001,max(-15, psi_head + (L*F/gL) * rho*g_p*1.e-6))
#
#        else:
#            KL = k_soil_soil(psi, n.soil_class) # [cm d-1]
#            psi_base = min(-0.001, max(-15, psi_head + (L*F/KL) * rho*g_p*1.e-6))
#
#        if vid != vid_base:
#            p.psi_head = psi_base 
#        else:
#            psi_soil = psi_base
#
#
#    print psi_soil
#
#    return

def soil_rhyzo_water_potential(flux, psi_collar_init, psi_bulk,
                               psi_error_thershold, soil_class, time_step,
                               rhyzo_rad, soil_rad, depth, psi_min=-3.):
    """
    """

    param = def_param_soil()[soil_class]
    theta_r,theta_s,alpha,n,k_sat = [param[i] for i in range(5)]

    # Rough estimation of total porosity volumes [m3]
    bulk_soil_volume = pi * (soil_rad**2 - rhyzo_rad**2) * depth * theta_s
    rhyzosphere_volume = pi * rhyzo_rad**2 * depth * theta_s

    F = flux * 1.e-3 # kg T-1 to m3 T-1

    influx = 0.
    psi_rhyzo = psi_collar_init

    for ipsi in range(100):
        psi_prev = psi_rhyzo
        F -= influx
        psi_rhyzo = soil_water_potential(psi_prev, F, soil_class,
                                         rhyzosphere_volume, psi_min)

        if abs(psi_prev - psi_rhyzo) < psi_error_thershold:

            psi_bulk = soil_water_potential(psi_bulk, influx, soil_class,
                                            bulk_soil_volume, psi_min)
            break
        else:
            k_soil = k_soil_soil(psi_bulk, soil_class) * 1.e-2 * pi * 2.*soil_rad * depth # cm d-1 to m3 d-1
            k_soil /= {'D':1., 'H':24., 'T':1440., 'S':86400}[time_step] # m3 d-1 to m3 T-1
            influx = k_soil * (psi_bulk - psi_rhyzo)* 1.e6/(rho*g_p) / (soil_rad - rhyzo_rad)

#    print soil_class, psi_rhyzo
    return psi_rhyzo, psi_bulk