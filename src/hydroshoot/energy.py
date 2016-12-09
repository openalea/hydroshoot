# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Energy balance module of HydroShoot.

This module computes leaf (and eventually other elements) tempertaure of a
given plant shoot.
"""
#import numpy as np
#import openalea.mtg.traversal as traversal

from scipy import optimize, mean
from sympy.solvers import nsolve
from sympy import Symbol
import time

from alinea.caribu.CaribuScene import CaribuScene
from hydroshoot import meteo_utils as mutils


# Global variables
##global a_PAR
##global a_NIR
##global a_glob
##global e_sky
##global e_leaf
##global e_soil
##global sigma
##global lambda_
##global Cpn
#global E_glob
#global k_sky
#global k_leaves
#global k_soil
#global gbH
#global E
#global T_sky
#global T_air
#global T_soil
#global Pa


def MTG_Energy_Prop(a_PAR=0.87, a_NIR=0.35, a_glob=0.6, e_sky=1.0, e_leaf=0.96,
                    sigma=5.670373e-8, e_soil = 0.95, lambda_=40.66e3, Cp = 29.07):
    """
    Returns a dictionary of spectrometric and energy balance-related properties.

    Parameters:
    - **a_PAR**: Leaf absorptance to the PAR [-]
    - **a_NIR** : Leaf absorptance to the NIR radiation [-]
    - **a_glob**: Leaf absorptance to the global radiation [-]
    - **e_sky**: sky emissivity [-]
    - **e_leaf**: leaf emissivity [-]
    - **e_soil**: soil emissivity [-]
    - **sigma**: Stefan-Boltzmann constant [W m-2 K-4]
    - **lambda_**: Latent heat for evaporization [J mol-1], [W s mol-1]
    - **Cp**: Isobaric heat capacity of the air [J mol-1 K-1]
    """

    energy_prop_dict = {
    'a_PAR' : a_PAR,
    'a_NIR' : a_NIR,
    'a_glob' : a_glob,
    'e_sky' : e_sky,
    'e_leaf' : e_leaf,
    'e_soil' : e_soil,
    'sigma' : sigma,
    'lambda_' : lambda_,
    'Cp' : Cp
    }

    return energy_prop_dict


def MTG_vis_a_vis(g, limit=-0.01):
    """
    Associates form factor values to leaves and soil elements.

    :Parameters:
    - **g**: an MTG object
    - **limit**: float (negative), a threhold for below which form factor values are ignored
    """

    cscene = CaribuScene(g)
    form_factor_array=cscene.form_factors(aggregate=True)

    col = form_factor_array.columns
    for vid in g.VtxList():
        try:
            if g.node(vid).label.startswith(('L','soil','other')):

#                k_tot = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
#                        if ff < 0]))
#                k_leaves = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
#                        if ff < 0 and not g.node(col[ivid]).label.startswith('soil')]))
#                k_soil = min(1.,max(0., k_tot - k_leaves))
#                k_sky = min(1.,2.-k_tot)

                k_sky = max(0.,min(1.,2.+ sum(form_factor_array[vid][form_factor_array[vid] < 0])))

                k_leaves = min(2.,sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) \
                        if ff < 0 and not g.node(col[ivid]).label.startswith(('soil','other'))]))

                # hacks while awaiting for a complete debug of teh form factors
                k_soil = min(1,max(2. - k_sky - k_leaves, 0.))

                g.node(vid).k_sky = k_sky
                g.node(vid).k_leaves = k_leaves
                g.node(vid).k_soil = k_soil
                g.node(vid).vis_a_vis = g.node(vid).vis_a_vis = {col[ivid]:ff \
                    for ivid, ff in enumerate(form_factor_array[vid]) if ff < limit}
        except:
            pass

    return


#def VineEnergyX(T_leaf):
#    """
#    """
##    E_SW = a_PAR*E_PAR + a_NIR*E_NIR
#    E_SW = a_glob*E_glob
#    delta_E_LW = e_leaf*(k_sky*e_sky*sigma*(T_sky)**4+k_leaves*e_leaf*sigma*(T_leaf)**4+k_soil*e_soil*sigma*(T_soil)**4) - 2*e_leaf*sigma*(T_leaf)**4
#    E_Y = -lambda_*E #(1./(2/(1.37*gb)+1/gs))*VPD/Pa
#    E_H = -1.37*gbH*Cp*(T_leaf-T_air)
#    E_error = E_SW + delta_E_LW + E_Y + E_H
#
#    return E_error


Energy_Prop = MTG_Energy_Prop()
nrj_Prop_tuple = ('a_PAR','a_NIR','a_glob','e_sky','e_leaf','e_soil','sigma','lambda_','Cp')
a_PAR,a_NIR,a_glob,e_sky,e_leaf,e_soil,sigma,lambda_,Cp = [Energy_Prop[ikey] for ikey in nrj_Prop_tuple]


#def leaf_temperature_solo(mtg, macro_meteo, leaf_lbl_prefix='L'):
#    """
#    """
#    global E_glob
#    global k_sky
#    global k_leaves
#    global k_soil
#    global gbH
#    global E
#    global T_sky
#    global T_air
#    global T_soil
#    global Pa
#
#    for vid in traversal.pre_order2(mtg,mtg.node(mtg.root).vid_base):
#        if mtg.node(vid).label.startswith(leaf_lbl_prefix):
#            node = mtg.node(vid)
#            E_glob = node.Eabs/(0.48*4.6)
#            k_sky, k_leaves, k_soil = [node.vis_a_vis[ikey] for ikey in ('k_sky','k_leaves','k_soil')]
#            gbH = node.gbH
#            E = node.E
#            T_sky, T_air, T_soil, Pa = [macro_meteo[ikey] for ikey in ('T_sky', 'T_air', 'T_soil', 'Pa')]
#
#            t_leaf = node.Tlc if 'Tlc' in node.properties() else T_air - 273.15
#
#            t_leaf0 = optimize.newton_krylov(VineEnergyX, t_leaf+273.15) - 273.15
#            node.Tlc = t_leaf0
#
#    return mtg


def leaf_temperature(g, macro_meteo, solo=True, leaf_lbl_prefix='L'):
    """
    Returns the "thermal structure", temperatures [degreeC] of each individual leaf and soil elements.

    :Parameters:
    - **g**: an MTG object
    - **macro_meteo**: a dictionary with keys: 'T_sky', 'T_soil', 'T_air', 'Pa', all in **absolute temperatures [K]**.
    - **solo**: logical,
        - True (default), calculates energy budget for each element, assuming the temperatures of surrounding leaves constant (from previous calculation step)
        - False, computes simultaneously all temperatures using `sympy.solvers.nsolve` (**very costly!!!**)
    - **leaf_lbl_prefix**: string, the prefix of the label of the leaves
    """
#    labels = g.property('label')
#    leaves_ids = [vid for vid, label in labels.iteritems() if label.startswith(leaf_lbl_prefix)]
#    symbols = {(vid, Symbol('t%d'%vid)) for lid in leaves_id}
    # equations = {(vid,equation(vid, symbols) for vid in leaves_id}
    # nsolve([equations[vid] for vid in leaves_id], [symbols[vid] for vid in leaves_id], verify=False)

    # Calculates leaf's temperature
    if solo:
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

                E_leaves = -sigma*sum([node.vis_a_vis[ivid]*(g.node(ivid).Tlc+273.15)**4 \
                            for ivid in node.vis_a_vis.keys()])

                def _VineEnergyX(T_leaf):
                    E_SW = a_glob*E_glob
                    delta_E_LW = e_leaf*(k_sky*e_sky*sigma*(T_sky)**4+\
                                         e_leaf*E_leaves+\
#                                         k_leaves*e_leaf*sigma*(T_leaf)**4+\
                                         k_soil*e_soil*sigma*(T_soil)**4)\
                                 - 2*e_leaf*sigma*(T_leaf)**4
                    E_Y = -lambda_*E
                    E_H = -gbH*Cp*(T_leaf-T_air)
                    E_error = E_SW + delta_E_LW + E_Y + E_H
                    return E_error

                t_leaf0 = optimize.newton_krylov(_VineEnergyX, t_leaf+273.15) - 273.15
#                print t_leaf,t_leaf0
                node.Tlc = t_leaf0
                node.SW = a_glob*E_glob
                node.LW_in = e_leaf*(k_sky*e_sky*sigma*(T_sky)**4+\
#                                         e_leaf*E_leaves+\
                                         k_leaves*e_leaf*sigma*(t_leaf0+273.15)**4+\
                                         k_soil*e_soil*sigma*(T_soil)**4)
                node.LW_out = - 2*e_leaf*sigma*(t_leaf0+273.15)**4
                node.E_Y = -lambda_*E
                node.E_H = -gbH*Cp*(t_leaf0+273.15-T_air)

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

    return

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
        E_Y = -lambda_* 0.06 * mutils.VPD_leaf_air(t_air, T_soil-273.15,hs)/Pa # 0.06 is gM from Bailey 2016 AFM 218-219:146-160
        E_H = -0.5 * Cp * (T_soil-T_air) # 0.5 is gH from Bailey 2016 AFM 218-219:146-160
        E_error = E_SW + delta_E_LW + E_Y + E_H
        return E_error

    t_soil0 = optimize.newton_krylov(_SoilEnergyX, t_soil+273.15) - 273.15
#                print t_leaf,t_leaf0
    node.Tsoil = t_soil0


    return t_soil0

#******************************************************************************
# To be deleted later (debugging stuff)
#******************************************************************************
#mtg_dict = {}
#for vid in mtg.property('geometry').keys():
#    node = mtg.node(vid)
#    E_glob = node.Eabs/(0.48*4.6)
#    k_dict = node.vis_a_vis
#    k_sky = k_dict['k_sky']
#    k_veg_soil = {ivid:k_dict[ivid] for ivid in k_dict.keys() if ivid!='k_sky'}
#    gbH = node.gbH
#    E = node.E
#    T_sky, T_air, T_soil, Pa = [macro_meteo[ikey] for ikey in ('T_sky', 'T_air', 'T_soil', 'Pa')]
#    t_leaf = node.Tlc if 'Tlc' in node.properties() else T_air - 273.15
#
#    mtg_dict[vid] = {'E_glob':E_glob, 'k_dict':k_dict, 'k_sky':k_sky,
# 'k_veg_soil':k_veg_soil, 'gbH':gbH, 'E':E, 'T_sky':T_sky, 'T_air':T_air,
# 'T_soil':T_soil, 'Pa':Pa,'t_leaf':t_leaf}
#
#from csv import writer, reader
#import ast
#
#w = writer(open("/home/albashar/Documents/Christophe_exemple/mtg2.csv", "w"))
#for vid in g2.VtxList():
#    if vid > 0:
#        if g2.node(vid).label.startswith('L'):
#            w.writerow([vid, g2.node(vid).properties()])
#
#g_dict = {}
#for key, val in csv.reader(open("/home/albashar/Documents/Christophe_exemple/mtg.csv")):
#    g_dict[key] = ast.literal_eval(val)
#
#
#    for ikey, ival in val.:
#        print ikey, ival
#        i_dict[ikey] = ival
#    g_dict[key] = i_dict
#
#import ast
#ast.literal_eval(val)
#{'muffin': 'lolz', 'foo': 'kitty'}
#
#vis_dict = {col[ivid]:ff for ivid, ff in enumerate(form_factor_array[vid]) if ff < limit}
#T_leaves = mean([g.node(col[ivid]).Tlc for ivid, ff in enumerate(form_factor_array[vid]) \
#                if ff < 0 and not g.node(col[ivid]).label.startswith('soil')]) + 273.15
#k_soil = sum([-ff for ivid, ff in enumerate(form_factor_array[vid]) if ff < 0 and g.node(col[ivid]).label.startswith('soil')])
#e_leaf*(k_sky*e_sky*sigma*(T_sky)**4+k_leaves*e_leaf*sigma*(T_leaves)**4+k_soil*e_soil*sigma*(T_soil)**4) - 2*e_leaf*sigma*(T_leaf)**4
