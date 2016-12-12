# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:17:18 2016

Application test of HydroShoot

@author: albashar
"""

from pandas import read_csv, DataFrame, date_range, DatetimeIndex #, read_excel
import datetime as dt
import matplotlib as mpl
from copy import deepcopy
import numpy as np
import time
#from scipy.stats import linregress

import openalea.mtg.traversal as traversal
from openalea.plantgl.all import Scene, surface
from openalea.mtg.plantframe import color

import hydroshoot.architecture as HSArc
import hydroshoot.irradiance as HSCaribu
import hydroshoot.exchange as HSExchange
import hydroshoot.hydraulic as HSHyd
import hydroshoot.energy as HSEnergy
import hydroshoot.display as HSVisu
from hydroshoot.data_access import get_path

mpl.style.use('ggplot')

#==============================================================================
# # Construction of the plant mock-up
#==============================================================================

# Path for plant digitalization data.
csv_file_path = get_path('grapevine_pot.csv')

g=HSArc.VineMTG(csv_file_path)

#Local Coordinates Correction
for v in traversal.iter_mtg2(g,g.root):
    n = g.node(g.Trunk(v, Scale=1)[0])
    theta = 180 if int(n.index()) < 200 else -90 if int(n.index()) < 300 else 0.
    HSArc.VineOrient(g,v,theta, local_rotation = True)

#Scene rotation
for v in traversal.iter_mtg2(g,g.root):
    HSArc.VineOrient(g,v,90., local_rotation = False)

for v in traversal.iter_mtg2(g,g.root):
    HSArc.VinePhytoModular(g,v)
    HSArc.VineMTGProp(g,v)
    HSArc.VineMTGGeom(g,v)#,theta_1=90,theta_2=180,theta_2_cv=10.)
    HSArc.VineTransform(g,v)

# Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
#MyScene = HSVisu.visu(g,def_elmnt_color_dict=True,scene=Scene())

#==============================================================================
# Divers initialisations
#==============================================================================
#sky temperature [degreeC]
t_sky = -20

# Unit length conversion (from scene unit to the standard [m]) unit)
LengthConv=1.e-2

# Topological location
latitude = 43.61
longitude = 3.87
elevation = 44.0
geo_location = (latitude, longitude, elevation)

# Attaching optical properties to MTG elements
g = HSCaribu.opticals(g, wave_band='SW')

# Default farquhar parameters
par_photo = HSExchange.par_photo_default()

# Suppression of undesired geometry for light and energy calculations
geom_prop = g.properties()['geometry']
vidkeys = []
for vid in g.properties()['geometry']:
    n = g.node(vid)
    if not n.label.startswith(('L','other','soil')):
        vidkeys.append(vid)
[geom_prop.pop(x) for x in vidkeys]
g.properties()['geometry'] = geom_prop

# Climate data
meteo_file = get_path('grapevine_pot_meteo.csv')
meteo_tab = read_csv(meteo_file, sep=';', decimal='.', header=0)
meteo_tab.time = DatetimeIndex(meteo_tab.time)
meteo_tab = meteo_tab.set_index(meteo_tab.time)

# Adding missing data
if 'Ca' not in meteo_tab.columns:
    meteo_tab['Ca'] = [400.]*len(meteo_tab) # ppm [CO2]
if 'Pa' not in meteo_tab.columns:
    meteo_tab['Pa'] = [101.3]*len(meteo_tab)# atmospheric pressure

# Soil water potential forcing
psi_pd = DataFrame([-0.17,-0.13,-0.19,-0.38,-0.61,-0.19], #1
                   index= [x.date() for x in
                   DatetimeIndex(('2012-7-20','2012-7-27','2012-8-1',
                                        '2012-8-2','2012-8-3','2012-8-9'))],
                   columns=['psi'])
psi_soil = - 0.6

# Maximum number of iterations for both temperature and hydraulic calculations
max_iter = 100

# Maximum acceptable error threshold in hydraulic calculations
psi_error_threshold = 0.05  # [MPa]

# Maximum acceptable error threshold in temperature calculations
t_error_crit = 0.02  # [MPa]

# Identifying and attaching the base node of a single MTG
vid_base = HSArc.MTGbase(g,vtx_label='inT')
g.node(g.root).vid_base = vid_base

# Initializing all xylem potential values to soil water potential
if len(g.property('psi_head'))==0:
    for vtx_id in traversal.pre_order2(g,vid_base):
        g.node(vtx_id).psi_head = psi_pd.psi[0]

# Initializing sapflow to 0
if len(g.property('Flux'))==0:
    for vtx_id in traversal.pre_order2(g,vid_base):
        g.node(vtx_id).Flux = 0.

# Computation of the form factor matrix
tt = time.time()
HSEnergy.form_factors_matrix(g, -0.000001)
print ("---%s minutes ---" % ((time.time()-tt)/60.))

#==============================================================================
# Simulations
#==============================================================================

# Determination of the simulation period
sdate = dt.datetime(2012,8,1,14,00,0,)
edate = dt.datetime(2012,8,1,14,00,0,)
deltat = edate - sdate
datet = date_range(sdate, edate, freq='H')
meteo = meteo_tab.ix[datet]

tt = time.time()

# The time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for date in meteo.time:
    print '++++++++++++++++++++++Date',date

#   Select of meteo data
    imeteo = meteo[meteo.time==date]

#   Read soil water potntial at midnight
    if date.hour == 0:
        try:
            psi_soil_init = psi_pd.ix[date.date()][0]
            psi_soil = psi_soil_init
        except KeyError:
            pass
# Estimate soil water potntial evolution due to transpiration
    else:
        psi_soil = HSHyd.soil_water_potential(psi_soil,g.node(vid_base).Flux*3600.,
                          soil_class='Sandy_Loam',
                          intra_dist=1.,inter_dist=1.8,depth=1.2)
    print 'psi_soil', round(psi_soil,3)

#   Initialize leaf [and other elements] temperature equal to air temperature
    if len(g.property('Tlc'))==0:
        g.properties()['Tlc'] = {vid:imeteo.Tac[0] for vid in g.VtxList() if vid >0 and g.node(vid).label.startswith(('L','other'))}

#   Compute irradiance distribution over the scene
    caribu_source = HSCaribu.irradiance_distribution(imeteo,geo_location,'Rg_Watt/m2',
                                              'Europe/Paris', '46', 'soc')

#   Compute irradiance interception and absorbtion
    g, caribu_scene = HSCaribu.hsCaribu(g, imeteo, source = caribu_source,
                            unit_scene_length='cm',local_date=None,
                                 geo_location=None, E_type=None, infinite=False,
                                 pattern=((-470.,-180.),(470.,180.)))

#   Hack forcing of soil temperture (model of soil temperature under development)
    dt_soil = [-2,-2,-2,-1,-1, 0, 0, 0,10,15,20,20,18,16,14,10, 6, 5, 4, 3, 2, 1, 0,-1]
    t_soil = imeteo.Tac[0] + dt_soil[date.hour]

#   Climatic data for energy balance module
    macro_meteo = {'T_sky':t_sky+273.15, 'T_soil':t_soil+273.15,
                   'T_air':imeteo.Tac[0]+273.15,'Pa':imeteo.Pa[0]}

# The t loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    t_error_list = []; t_iter_list = []
    for it in range(max_iter):
        t_prev = deepcopy(g.property('Tlc'))

# The psi loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for i in range(max_iter):
            psi_prev = deepcopy(g.property('psi_head'))

#           Computes gas-exchange fluxes. Leaf T and Psi are from prev calc loop
            HSExchange.VineExchange(g, par_photo, imeteo, psi_soil, E_type2='Ei',
                                    leaf_lbl_prefix='L', model='misson',
                                    g0=0.0, rbt=2./3., ca=400., m0=5.278,
                                    psi0=-1., D0=30., n=4.)

#           Computes sapflow and hydraulic properties
            HSHyd.hydraulic_prop(g, vtx_label='inT', MassConv=18.01528,
                              LengthConv=LengthConv,a=1.6,b=2)

#           Computes xylem water potential
            N_iter_psi=HSHyd.xylem_water_potential(g, psi_soil, psi_min=-3.0,
                             model='misson', max_iter=100, psi_error_crit=0.001,
                             vtx_label='inT',LengthConv=1.e-2,fifty_cent=-1.51,
                             sig_slope=1)

            psi_new = g.property('psi_head')

#           Evaluation of xylem conversion creterion
            psi_error_dict = {}
            for vtx_id in g.property('psi_head').keys():
                psi_error_dict[vtx_id] = abs(psi_prev[vtx_id]-psi_new[vtx_id])

            psi_error = max(psi_error_dict.values())
            print 'psi_error = ', round(psi_error,3), '::Nb_iter = ', N_iter_psi

            if psi_error < psi_error_threshold:
                break
            else:
#                try:
#                    if abs(psi_error_prev-psi_error)/max(psi_error_prev,psi_error) < 0.01:
#                        c_n = 0.75
#                        c_o = 0.25
#                    else:
#                        c_n = 0.5
#                        c_o = 0.5
#                except:
#                    c_n = 0.5
#                    c_o = 0.5
                for vtx_id in psi_new.keys():
                    g.node(vtx_id).psi_head = 0.5*(psi_prev[vtx_id]+psi_new[vtx_id])
#                    g.node(vtx_id).psi_head = (c_o*psi_prev[vtx_id]+c_n*psi_new[vtx_id])
#            print 'c_o = %f, c_n = %f'%(c_o, c_n)
            psi_error_prev = psi_error

# End psi loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#       Compute leaf temperature
        t_iter=HSEnergy.leaf_temperature(g, macro_meteo,solo=True,leaf_lbl_prefix='L',
                                  max_iter=100, t_error_crit=t_error_crit)
        t_iter_list.append(t_iter)

        t_new = deepcopy(g.property('Tlc'))

#       Evaluation of leaf temperature conversion creterion
        error_dict={vtx:abs(t_prev[vtx]-t_new[vtx]) for vtx in g.property('Tlc').keys()}

        t_error = round(max(error_dict.values()),3)
        print '********** t_error = %d'%t_error, 'max(iter) = %d'%max(t_iter_list)
#        t_error_list.append(t_error)
        if max(error_dict.values()) < t_error_crit:
            break
        else:
            for vtx_id in t_new.keys():
                g.node(vtx_id).Tlc = (t_prev[vtx_id]+t_new[vtx_id])*0.5
#            
#            if len(t_error_list) >= 10:
#                slope, intercept, r_value, p_value, std_err = linregress(x=range(10),y=t_error_list[-10:])
#                if r_value**2 < 0.9:
#                    for vtx_id in t_new.keys():
#                        g.node(vtx_id).Tlc = (t_prev[vtx_id]+t_new[vtx_id])*0.5
#                else:
#                    print 'vvvvvvvvvvvv', r_value**2
            
# End t loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# End time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#==============================================================================
# Illustrations
#==============================================================================

# Mock-ups
# --------

## Intercepted light (results in 'fig_02_light_interception.png' and 'fig_02_light_interception_colorbar.png')
#light_scene = HSVisu.visu(g,plot_prop='Ei',scene=Scene())
#
## Net photosynthesis (results in 'fig_03_An.png' and 'fig_03_An_colorbar.png')
#An_scene = HSVisu.visu(g,plot_prop='An',scene=Scene())
#
## Stomatal conductance (results in 'fig_04_gs.png' and 'fig_04_gs_colorbar.png')
#gs_scene = HSVisu.visu(g,plot_prop='gs',scene=Scene(), fmt='%0.6f')
#
## Water potential (results in 'fig_05_Tlc.png' and 'fig_05_Tlc_colorbar.png')
#Tlc_scene = HSVisu.visu(g,plot_prop='Tlc',scene=Scene(), fmt='%6.4f')


# Plots
# ----

# Hydraulic map (result in 'fig_06_psi.png')
fig_psi,ax_psi = HSVisu.property_map(g,prop='psi_head',style='r',add_head_loss=True, color='red')

# Tlc vs psi_head (result in 'fig_06_Tlc_psi.png')
fig_Tlc_psi,ax_Tlc_psi = HSVisu.prop_fun_prop(g, prop1='Tlc', prop2='psi_head',one_one=False)

# gs vs psi_head (result in 'fig_07_gs_psi.png')
fig_gs_psi,ax_gs_psi = HSVisu.prop_fun_prop(g, prop1='gs', prop2='psi_head',one_one=False)