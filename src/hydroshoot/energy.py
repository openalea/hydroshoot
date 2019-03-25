# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Energy balance module of HydroShoot.

This module computes leaf (and eventually other elements) tempertaure of a
given plant shoot.
"""

from scipy import optimize, mean#, spatial
from sympy.solvers import nsolve
from sympy import Symbol
from copy import deepcopy
import time


import openalea.mtg.traversal as traversal
#from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import turtle
import alinea.astk.icosphere as ico


from hydroshoot import utilities as utils
from hydroshoot import irradiance as HSCaribu
from hydroshoot.architecture import vine_orientation, vine_mtg_properties, vine_mtg_geometry, vine_transform


#def energy_params(a_PAR=0.87, a_NIR=0.35, a_glob=0.6, e_sky=1.0, e_leaf=0.96,
#                    sigma=5.670373e-8, e_soil = 0.95, lambda_=44.00e3, Cp = 29.07):
#    """
#    Returns a dictionary of spectrometric and energy balance-related properties.
#
#    Parameters:
#    - **a_PAR**: Leaf absorptance to the PAR [-]
#    - **a_NIR** : Leaf absorptance to the NIR radiation [-]
#    - **a_glob**: Leaf absorptance to the global radiation [-]
#    - **e_sky**: sky emissivity [-]
#    - **e_leaf**: leaf emissivity [-]
#    - **e_soil**: soil emissivity [-]
#    - **sigma**: Stefan-Boltzmann constant [W m-2 K-4]
#    - **lambda_**: Latent heat for evaporization [J mol-1], [W s mol-1]
#    - **Cp**: Isobaric heat capacity of the air [J mol-1 K-1]
#    """
#
#    energy_prop_dict = {
#    'a_PAR' : a_PAR,
#    'a_NIR' : a_NIR,
#    'a_glob' : a_glob,
#    'e_sky' : e_sky,
#    'e_leaf' : e_leaf,
#    'e_soil' : e_soil,
#    'sigma' : sigma,
#    'lambda_' : lambda_,
#    'Cp' : Cp
#    }
#
#    return energy_prop_dict


#def form_factors_matrix(g, pattern, LengthConv, limit=-0.01):
#    """
#    Associates form factor values to leaves and soil elements.
#
#    :Parameters:
#    - **g**: an MTG object
#    - **limit**: float (negative), a threhold for below which form factor values are ignored
#    """
##    max_height = max([coord[2] for coord in g.property('TopPosition').values()])
#    d_sphere = spatial.distance.cdist(pattern,pattern, 'euclidean').max()
#    cscene = CaribuScene(g,pattern=pattern)
#    form_factor_array=cscene.form_factors(aggregate=True, d_sphere=d_sphere)
#
#    col = form_factor_array.columns
#    for vid in g.VtxList():
#        try:
#            if g.node(vid).label.startswith(('L','soil','other')):
#
##                k_tot = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
##                        if ff < 0]))
##                k_leaves = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
##                        if ff < 0 and not g.node(col[ivid]).label.startswith('soil')]))
##                k_soil = min(1.,max(0., k_tot - k_leaves))
##                k_sky = min(1.,2.-k_tot)
## +++
##                k_sky = max(0.,min(1.,2.+ sum(form_factor_array[vid][form_factor_array[vid] < 0])))
##
##                k_leaves = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
##                        if ff < 0 and not g.node(col[ivid]).label.startswith(('soil','other'))]))
##
##                # hacks while awaiting for a complete debug of teh form factors
##                k_soil = 1. #min(1,max(2. - k_sky - k_leaves, 0.))
#                k_soil = form_factor_array[str(vid)][[vids for vids in col if g.node(int(vids)).label.startswith(('other','soil'))][0]]
#
##                g.node(vid).k_sky = k_sky
##                g.node(vid).k_leaves = k_leaves
#                g.node(vid).k_soil = k_soil
##                g.node(vid).vis_a_vis = g.node(vid).vis_a_vis = {col[ivid]:ff \
##                    for ivid, ff in enumerate(form_factor_array[vid]) if ff < limit}
#        except:
#            pass
#
#    return


def form_factors_simplified(g, pattern, leaf_lbl_prefix='L',
                            stem_lbl_prefix=('in','Pet','cx'),
                            turtle_sectors='46', icosphere_level=3,
                            unit_scene_length='cm'):
    """
    Returns sky and soil contribution factors (resp. k_sky and k_soil) to the energy budget equation.
    Both factors are calculated and attributed to each element of the scene.

    :Note 1:
    This function is a simplified approximation of the form factors matrix which is calculated by the function :func:**form_factors_matrix**.
    The canopy is turned upside-down and light is projected in each case to estimate the respective contribution of the sky ({z}>=0) and soil ({z}<=0) to energy budget calculations.

    :Note 2:
    When **icosphere_level** is defined, **turtle_sectors** is ignored.

    :Parameters:
    - **g**: an MTG object
    - **pattern**: tuple, 2D Coordinates of the domain bounding the scene for its replication.
                     (xmin, ymin, xmax, ymax) scene is not bounded along z axis.
                     Alternatively a *.8 file.
                     if `None` (default), scene is not repeated
    - **leaf_lbl_prefix**, **stem_lbl_prefix**: string, the prefices of the leaf label and stem label, resp.
    - **turtle_sectors**: string, see :func:`turtle` from `sky_tools` package
    - **icosphere_level**: integer, the level of refinement of the dual icosphere. By default 46 ^polygons are returned (level=3). See :func:`alinea.astk.icosphere.turtle_dome` for details
    - **unit_scene_length**: the unit of length used for both scene coordinates and pattern (should be one of `CaribuScene.units` default)
    """

    opt_prop_ff={'SW':{'leaf':(0.001,0.0),'stem':(0.001,),'other':(0.001,0.0)},
                 'SW':{'leaf':(0.001,0.0),'stem':(0.001,),'other':(0.001,0.0)}}

    g = HSCaribu.optical_prop(g,leaf_lbl_prefix, stem_lbl_prefix,'SW',opt_prop_ff)

    if not icosphere_level:
        energy, emission, direction, elevation, azimuth = turtle.turtle(sectors=turtle_sectors,format='uoc',energy=1.)
    else:
        vert,fac = ico.turtle_dome(icosphere_level)
        direction = ico.sample_faces(vert, fac, iter=None, spheric=False).values()
        direction = [i[0] for i in direction]
        direction = map(lambda x: tuple(list(x[:2])+[-x[2]]),direction)

    caribu_source = zip(len(direction)*[1./len(direction)],direction)        


    for s in ('pirouette', 'cacahuete'):
        print '... %s'%s
        for v in traversal.pre_order2(g,3):
            vine_orientation(g,v,180., v_axis=[1.,0.,0.], local_rotation = False)
        
        for v in traversal.iter_mtg2(g,g.root):
            vine_mtg_properties(g,v)
            vine_mtg_geometry(g,v)
            vine_transform(g,v)
        
        
        # Compute irradiance interception and absorbtion
        g, caribu_scene = HSCaribu.hsCaribu(mtg=g,meteo=None,local_date=None,
                           geo_location=None, E_type=None,
                           unit_scene_length=unit_scene_length, tzone=None,
                           wave_band='SW', source = caribu_source, direct=True,
                           infinite=True, nz=50, dz=5, ds=0.5,
                           pattern=pattern, turtle_sectors=None,
                           turtle_format=None,
                           leaf_lbl_prefix=leaf_lbl_prefix,
                           stem_lbl_prefix=stem_lbl_prefix,
                           opt_prop=opt_prop_ff, rotation_angle = 0.)
        
    #    caribu_scene.getIncidentEnergy()
        if s == 'pirouette':
            k_soil_dict = g.properties()['Ei']
            max_k_soil = float(max(k_soil_dict.values()))
            g.properties()['k_soil'] = {vid:k_soil_dict[vid]/max_k_soil for vid in k_soil_dict}
        elif s == 'cacahuete':
            k_sky_dict = g.properties()['Ei']
            max_k_sky = float(max(k_sky_dict.values()))
            g.properties()['k_sky'] = {vid:k_sky_dict[vid]/max_k_sky for vid in k_sky_dict}


    for vid in g.properties()['Ei']:
        g.node(vid).k_leaves = max(0., 2.-(g.node(vid).k_soil+g.node(vid).k_sky))

    return


#Energy_Prop = energy_params()
#nrj_Prop_tuple = ('a_PAR','a_NIR','a_glob','e_sky','e_leaf','e_soil','sigma','lambda_','Cp')
#a_PAR,a_NIR,a_glob,e_sky,e_leaf,e_soil,sigma,lambda_,Cp = [Energy_Prop[ikey] for ikey in nrj_Prop_tuple]

a_PAR=0.87
a_NIR=0.35
a_glob=0.6
e_sky=1.0
e_leaf=0.96
e_soil = 0.95
sigma=5.670373e-8
lambda_=44.0e3
Cp = 29.07


def leaf_temperature(g, macro_meteo, solo=True, simple_ff=True,
                     leaf_lbl_prefix='L', max_iter=100, t_error_crit=0.01,
                     t_step = 0.5):
    """
    Returns the "thermal structure", temperatures [degreeC] of each individual leaf and soil elements.

    :Parameters:
    - **g**: an MTG object
    - **macro_meteo**: a dictionary with keys: 'T_sky', 'T_soil', 'T_air', 'Pa', all in **absolute temperatures [K]**.
    - **solo**: logical,
        - True (default), calculates energy budget for each element, assuming the temperatures of surrounding leaves constant (from previous calculation step)
        - False, computes simultaneously all temperatures using `sympy.solvers.nsolve` (**very costly!!!**)
    - **leaf_lbl_prefix**: string, the prefix of the label of the leaves
    - **max_iter**: integer, the allowable number of itrations (for solo=True)
    - **t_error_crit**: float, the allowable error in leaf temperature (for solo=True)
    """

#   Iterative calculation of leaves temperature
    if solo:
        t_error_trace = []
        it_step = t_step
        for it in range(max_iter):
            t_prev = deepcopy(g.property('Tlc'))
#            T_leaves = mean([g.node(vid).Tlc for vid in g.property('Tlc').keys() if g.node(vid).label.startswith(leaf_lbl_prefix)]) + 273.15
            T_leaves = mean([g.node(vid).Tlc for vid in g.property('Tlc').keys()]) + 273.15
            t_dict = {}
#           +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            for vid in g.property('Tlc').keys():
                node = g.node(vid)
                E_glob = node.Ei/(0.48*4.6) # Ei not Eabs
                k_sky = node.k_sky
                k_leaves = node.k_leaves
                k_soil = node.k_soil
                u = node.u

#                gbH = node.gb*1.37*0.9184 * Cp # Boundary layer conductance for heat [mol m2 s-1. The 0.9184 see Campbell and Norman (1998) in Gutschick (2016)
#                gbH = 3.9 * (macro_meteo['u']/l_w)**0.5
                l_w = node.Length/100.*0.72 # leaf length in the downwind direction [m]
                d_bl = 4.*(l_w/max(1.e-3,u))**0.5 /1000. # Boundary layer thikness in [m] (Nobel, 2009 pp.337)
                gbH = 2.*0.026 / d_bl # Boundary layer conductance to heat [W m-2 K-1]

# TODO: Replace the constant thermal conductivity coefficient of the air (0.026 W m-1 C-1) by a model accounting for air temperature, humidity and pressure (cf. Nobel, 2009 Appendix I)
                E = node.E
                T_sky, T_air, T_soil, Pa = [macro_meteo[ikey] for ikey in ('T_sky', 'T_air', 'T_soil', 'Pa')]
                t_leaf = node.Tlc if 'Tlc' in node.properties() else T_air - 273.15
                

                if not simple_ff:
                    E_leaves = -sigma*sum([node.vis_a_vis[ivid]*(g.node(ivid).Tlc+273.15)**4 \
                            for ivid in node.vis_a_vis.keys()])
                else:
#                    E_leaves = k_leaves*sigma*(T_leaves)**4
                    E_leaves = k_leaves*sigma*(t_leaf + 273.15)**4

                def _VineEnergyX(T_leaf):
                    E_SW = a_glob*E_glob
                    delta_E_LW = e_leaf*(k_sky*e_sky*sigma*(T_sky)**4+\
                                         e_leaf*E_leaves+\
                                         k_soil*e_soil*sigma*(T_soil)**4)\
                                 - 2*e_leaf*sigma*(T_leaf)**4
                    E_Y = -lambda_*E
                    E_H = -gbH*(T_leaf-T_air)
                    E_error = E_SW + delta_E_LW + E_Y + E_H
                    return E_error

                t_leaf0 = optimize.newton_krylov(_VineEnergyX, t_leaf+273.15) - 273.15
                t_leaf0
                t_dict[vid] = t_leaf0
                
                node.gbH = gbH

            g.properties()['Tlc'] = t_dict

#           +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            t_new = deepcopy(g.property('Tlc'))

#           Evaluation of leaf temperature conversion creterion
            error_dict={vtx:abs(t_prev[vtx]-t_new[vtx]) for vtx in g.property('Tlc').keys()}
            
            t_error = max(error_dict.values())
            t_error_trace.append(t_error)

            if t_error < t_error_crit:
                break
            else:
                try:
                    if abs(t_error_trace[-1] - t_error_trace[-2]) < t_error_crit:
                        it_step = max(0.01, it_step/2.)
                except:
                    pass

                t_new_dict = {}
                for vtx_id in t_new.keys():
                    tx = t_prev[vtx_id] + it_step*(t_new[vtx_id]-t_prev[vtx_id])
                    t_new_dict[vtx_id] = tx

#               t_new_dict = {vtx_id:0.5*(t_prev[vtx_id]+t_new[vtx_id]) for vtx_id in t_new.keys()}
                g.properties()['Tlc'] = t_new_dict

                
#                g.properties()['Tlc'] = {vtx_id:0.5*(t_prev[vtx_id]+t_new[vtx_id]) for vtx_id in t_new.keys()}


#   Matrix iterative calculation of leaves temperature
    else:
        t_lst = []
        t_dict = {}

        t_dict={vid:Symbol('t%d'%vid) for vid in g.property('geometry').keys()}
    #    for vid in g.property('geometry').keys():
    ##        if g.node(vid).label.startswith(leaf_lbl_prefix):
    #        exec('t%d = %s' % (vid, None))
    #        locals()['t'+str(vid)] = Symbol('t'+str(vid))
    #        t_lst.append(locals()['t'+str(vid)])
    #        t_dict[vid] = locals()['t'+str(vid)]

        eq_lst = []
        t_leaf_lst = []
        for vid in g.property('geometry').keys():
            if g.node(vid).label.startswith(leaf_lbl_prefix):
                node = g.node(vid)
                E_glob = node.Ei/(0.48*4.6) # Ei not Eabs
                k_sky = node.k_sky
                k_leaves = node.k_leaves
                k_soil = node.k_soil
                gbH = node.gbH
                E = node.E
                T_sky, T_air, T_soil, Pa = [macro_meteo[ikey] for ikey in ('T_sky', 'T_air', 'T_soil', 'Pa')]
                t_leaf = node.Tlc if 'Tlc' in node.properties() else T_air - 273.15

                t_leaf_lst.append(t_leaf)
                t_lst.append(t_dict[vid])

        #        exec('eq%d = %s' % (vid, None))

                eq_aux = 0.
                for ivid in node.vis_a_vis.keys():
                    if not g.node(ivid).label.startswith('soil'):
                        eq_aux += -node.vis_a_vis[ivid] * ((t_dict[ivid])**4)

                eq = (a_glob * E_glob +
                    e_leaf * sigma * (k_sky * e_sky * (T_sky**4) +
                    e_leaf * eq_aux + k_soil * e_soil * (T_sky**4) -
                    2 * (t_dict[vid])**4) -
                    lambda_ * E - gbH * Cp * (t_dict[vid] - T_air))

                eq_lst.append(eq)

        tt = time.time()
        t_leaf0_lst = nsolve(eq_lst, t_lst, t_leaf_lst, verify=False) - 273.15
        print ("---%s seconds ---" % (time.time()-tt))

        ivid = 0
        for vid in g.property('geometry').keys():
    #        if g.node(vid).label.startswith(leaf_lbl_prefix):
            g.node(vid).Tlc = float(t_leaf0_lst[ivid])
            ivid += 1

        it = 1

    return it


def soil_temperature(g, meteo, T_sky, soil_lbl_prefix='other'):
    """
    Returns soil temperature based on a simplified energy budget formula.

    Parameters:
    - **t_air**: air temperature in degrees.
    """
    hs, Pa, t_air = [float(meteo[x]) for x in ('hs', 'Pa', 'Tac')]
    T_air = t_air + 273.15

    node=[g.node(vid) for vid in g.property('geometry') if g.node(vid).label.startswith('other')][0]
    T_leaf = mean(g.property('Tlc').values()) + 273.15

    E_glob = node.Ei/(0.48*4.6) # Ei not Eabs
    t_soil = node.Tsoil if 'Tsoil' in node.properties() else t_air

    def _SoilEnergyX(T_soil):
        E_SW = (1-0.25)*E_glob # 0.25 is rough estimation of albedo of a bare soil
        delta_E_LW = e_soil*sigma*(1.*e_sky*(T_sky)**4 + 1.*e_leaf*T_leaf**4- ((T_soil)**4)) # hack: 0% loss to deeper soil layers
#                             k_leaves*e_leaf*sigma*(T_leaf)**4)
        E_Y = -lambda_ * 0.06 * utils.vapor_pressure_deficit(t_air, T_soil - 273.15, hs) / Pa # 0.06 is gM from Bailey 2016 AFM 218-219:146-160
        E_H = -0.5 * Cp * (T_soil-T_air) # 0.5 is gH from Bailey 2016 AFM 218-219:146-160
        E_error = E_SW + delta_E_LW + E_Y + E_H
        return E_error

    t_soil0 = optimize.newton_krylov(_SoilEnergyX, t_soil+273.15) - 273.15
#                print t_leaf,t_leaf0
    node.Tsoil = t_soil0


    return t_soil0