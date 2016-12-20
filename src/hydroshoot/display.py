# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Visualization module for HydroShoot
"""

from pandas import DataFrame
import numpy as np
from scipy.spatial import distance
import matplotlib as mpl
import mpl_toolkits.mplot3d as m3d

from openalea.mtg.plantframe import color as pglcolor
import openalea.plantgl.all as pgl

mpl.style.use('ggplot')

from hydroshoot.architecture import MTGbase


def default_labels():
    """
    A dictionary of a number of default variable units often used for display of results.

    :This dictionary currently includes:
    - **'Eabs'**: [umol m-2 s-1]
    - **'Ei'**: [umol m-2 s-1]
    - **'Flux'**: [kg s-1]
    - **'KL'**: [kg s-1 m Mpa-1]
    - **'Kmax'**: [kg s-1 m Mpa-1]
    - **'Length'**: [cm]
    - **'Tlc'**: [degree C]
    - **'TopDiameter'**: [cm]
    - **'BotDiameter'**: [cm]
    - **'psi_head'**: [MPa]
    - **'An'**: [umol m-2 s-1]
    - **'gs'**: [mol m-2 s-1]
    - **'gb'**: [mol m-2 s-1]
    - **'E'**: [mol m-2 s-1]
    - **'Ci'**: [umol umol-1]
    - **'Cc'**: [umol mol-1]
    - **'k_sky'**: [-]
    - **'k_soil'**: [-]
    - **'k_leaves'**: [-]
    """

    lbl_dict = {
    'Eabs': '$\mathregular{E_{abs}\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'Ei': '$\mathregular{E_i\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'Flux': '$\mathregular{F\/[kg\/s^{-1}]}$',
    'KL': '$\mathregular{K_L\/[kg\/s^{-1}\/m\/MPa^{-1}]}$',
    'Kmax': '$\mathregular{K_{max}\/[kg\/s^{-1}\/m\/MPa^{-1}]}$',
    'Length': '$\mathregular{Diam_{bot}\/[cm]}$',
    'Tlc': '$\mathregular{T_{leaf}\/[^\circ\/C]}$',
    'TopDiameter': '$\mathregular{Diam_{bot}\/[cm]}$',
    'BotDiameter': '$\mathregular{Diam_{bot}\/[cm]}$',
    'psi_head': '$\mathregular{\Psi\/[MPa]}$',
    'An': '$\mathregular{A_n\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'gs': '$\mathregular{g_s\/[mol\/m^{-2}\/s^{-1}]}$',
    'gb': '$\mathregular{g_b\/[mol\/m^{-2}\/s^{-1}]}$',
    'E': '$\mathregular{E\/[mol\/m^{-2}\/s^{-1}]}$',
    'Ci': '$\mathregular{C_i\/[\mu mol_{CO_2}\/\mu mol]}$',
    'Cc': '$\mathregular{C_c\/[\mu mol_{CO_2}\/\mu mol]}$',
    'k_sky': '$\mathregular{k_{sky}\/[-]}$',
    'k_soil': '$\mathregular{k_{soil}\/[-]}$',
    'k_leaves': '$\mathregular{k_{leaves}\/[-]}$',
    }
    return lbl_dict


def default_color_dict(inT=(40,19,0), cx=(40,19,0),inT3y=(40,19,0),inT2y=(40,19,0),
                   in_=(40,19,0),Pet=(255,0,0),LI=(0,110,44),LII=(0,165,5),
                   other=(7,102,0),soil=(7,102,0),**kwargs):
    """
    Returns a dictionary of RGB colors for a given mtg labels.

    Default color keys are **inT**, **cx**, **inT3y**, **inT2y**, **in_**, **Pet**,
    **LI**, **LII**, **other**, **soil**.
    """

    color_dict = {
    'inT':inT,
    'cx':cx,
    'inT3y':inT3y,
    'inT2y':inT2y,
    'in':in_,
    'Pet':Pet,
    'LI':LI,
    'LII':LII,
    'other':other,
    'soil' : soil
    }

    color_dict.update(kwargs)

    return color_dict


def default_element_list(*args):
    """
    Returns a list of MTG nodes labels.
    If `*args` is not **None**, returns a list of given `*args`, otherwise it returns
    `['inT','cx','inT3y','inT2y','in','Pet','LI','LII','other', 'soil']`
    """

    if args:
        elmnt_list = list(args)
    else:
        elmnt_list = ['inT','cx','inT3y','inT2y','in','Pet','LI','LII','other', 'soil']

    return elmnt_list


def visu(g, plot_prop=None,min_value=None, max_value=None,tt=1001,cmap='jet',fmt='%6.0f',elmnt_labels=None,
         elmnt_color_dict=None,def_elmnt_color_dict=False, use_mtg_color=False,
         snap_shot_path=None, scene=None,sub_mtg_root=None):
    """
    Displays 3D moke-up using `plantgl` package.

    :Parameters:
    - **g**: an MTG object
    - **plot_prop**: string, an MTG property to plot
    - **min_value**: float, minimum property value for the color scale
    - **max_value**: float, maximum property value for the color scale
    - **cmap**: string, def('jet'), the colormap to use; only active when `plot_prop` is given or `use_mtg_color` is True
    - **fmt**: string, a legal text format, def('%6.0f'), the string format of matplotlib colorbar
    - **elmnt_labels**: a list of strings refering to the desired mtg elemnts to be displayed
    - **elmnt_colors_dict**: a dictionary of RGB color tuples for given mtg labels
    - **def_elmnt_color_dict**: logical, if `True`, uses :func:`default_color_dict()` to set a dictionary of RGB color tuples with `elmnt_labels` as keys
    - **use_mtg_color**: logical, if `True`, . If `False` uses the 'color' property from mtg nodes
    - **snap_shot_path**: string, if given, saves scene snapshot to the defined path and file name
    - **scene**: def None,if given, adds scene shapes to it
    - **sub_mtg_root**: integer, def `g.root`, an MTG node, if given, displays only its descendants
    """

    MyScene = pgl.Scene() if scene == None else scene

    if plot_prop is not None:
        prop = [g.node(vid).properties()[plot_prop] for vid in g.property(plot_prop) \
                if not g.node(vid).label.startswith('soil')]
        if min_value is None: min_value = min(prop)
        if max_value is None: max_value = max(prop)

        g, cb = pglcolor.colorbar(g, property_name=plot_prop, cmap=cmap,
                                 lognorm=False, N=6, fmt=fmt,
                                 min_value=min_value,max_value=max_value)
#        cb.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        g = pglcolor.colormap(g, property_name=plot_prop, cmap=cmap,\
                                lognorm=False,min_value=min_value,\
                                max_value=max_value)

        for vid in g.property(plot_prop).keys():
            n = g.node(vid)
            mesh = n.geometry
            label = n.label
            Scene_shape=pgl.Shape(mesh, pgl.Material(pgl.Color3(n.color)))
            MyScene.add(Scene_shape)

    else:

        if elmnt_labels is None:
            elmnt_labels = default_element_list()

        if elmnt_color_dict is None:
            if def_elmnt_color_dict:
                def_color_dict = default_color_dict()
                elmnt_color_dict = {key_: def_color_dict[key_] for key_ in elmnt_labels}
            elif use_mtg_color==False:
                raise StandardError ("Elements colors are missing. You may set elmnt_color_dict=True to use default colors provied by default_color_dict().")

        for vid in g.property('geometry'):
            n = g.node(vid)
            if n.label.startswith(tuple(elmnt_labels)):
                mesh = n.geometry
                label = n.label
                if use_mtg_color:
                    color = n.color
                else:
                    for i in tuple(elmnt_labels):
                            if i in label: color = elmnt_color_dict[i]
#                color=(50,50,50) if vid != tt else (255,255,0)
                Scene_shape=pgl.Shape(mesh, pgl.Material(pgl.Color3(color)))
                MyScene.add(Scene_shape)

    pgl.Viewer.display(MyScene)

    if snap_shot_path:
        pgl.Viewer.saveSnapshot(snap_shot_path)

    return MyScene


def hydraulic_map(g, prop='psi_head', fig=None,style=None):
    """
    Plots the hydraulic cart of an MTG.

    :Parameters:
    - **g**: an MTG object
    - **prop**: string, name of an MTG property to plot
    - **fig**: a matplotlib figure object, default None, if given adds the hydraulic cart to it
    - **style**: string, a legal matplotlib style
    """

    df = DataFrame(columns=['z',prop])
    prop_dict = g.property(prop)
    for vid in prop_dict.keys():
        if not g.node(vid).label.startswith(('soil', 'other')):
            df.loc[vid]=[g.node(vid).TopPosition[2],g.node(vid).properties()[prop]]

    if not fig:
        fig = mpl.pyplot.figure()
        ax=fig.add_subplot(111)
    else:
        ax=fig.axes[0]
    if not style: style = mpl.pyplot.rcParams['axes.color_cycle'][0]
    df.plot.scatter(x='z',y=prop,ax=ax, c=style)
    ax.set(xlabel='$\mathregular{Elevation\/[cm]}$',xlim=(min(df.z),max(df.z)),
           ylim=(min(df[prop]),max(df[prop])))#,ylabel='$\mathregular{\Psi_{stem}\/[Mpa]}$')

    return fig, ax


def property_map(g, prop='psi_head', fig=None, style=None, ylabel=None,
            add_head_loss=False, color=None):
    """
    Plots values of a given MTG property vs hight.

    :Parameters:
    - **g**: an MTG object
    - **prop**: string, name of an MTG property to plot
    - **fig**: a matplotlib figure object, default None, if given adds the hydraulic cart to it
    - **style**: string, a legal matplotlib style
    - **ylabel**: string, label of y axis, if None, tries to get a label from :func:`default_labels` based on **prop** string
    - **add_head_loss**: logical, if True, adds a dashed line representing water head loss due to elevation [-0.01 MPa m-1]
    - **color**: string, a legal matplotlib color name
    """

    if not fig:
        fig = mpl.pyplot.figure()
        ax=fig.add_subplot(111)
    else:
        ax=fig.axes[0]

    if ylabel is None:
        try:
            ylabel = default_labels()[prop]
        except:
            pass
    if color is None: color = (0,0,.5)
    vid_base = MTGbase(g)
    for vid in g.Extremities(vid_base):
        x=[g.node(ivid).TopPosition[2] for ivid in g.Ancestors(vid) if not g.node(ivid).label.startswith(('L'))]
        y=[g.node(ivid).properties()[prop] for ivid in g.Ancestors(vid) if not g.node(ivid).label.startswith(('L'))]
        ax.plot(x,y,'.-',color=color)
        ax.set(xlabel='z [cm]', ylabel=ylabel)
    #    ax.plot(x,y,'-',color=mpl.pyplot.rcParams['axes.color_cycle'][2])
    if add_head_loss:
        xlim = ax.get_xlim()
        ls = np.arange(xlim[0],xlim[1])
        ax.plot(ls,0.01*ls*(-0.01)+g.node(vid_base).psi_head, '--')

    return fig, ax


def prop_fun_prop(g, prop1='gs', prop2='psi_head', fig=None,style=None,
                  one_one=False, xlabel=None, ylabel=None):
    """
    Returns a scatter plot of two MTG properties.

    :Parameters:
    - **g**: an MTG object
    - **prop1**: string, name of the first MTG property to plot on vertical axis
    - **prop2**: string, name of the second MTG property to plot on horizontal axis
    - **fig**: a matplotlib figure object, default None, if given adds the hydraulic cart to it
    - **style**: string, a legal matplotlib style
    - **one_one**: logical, if True, adds a 1:1 line to figure
    - **xlabel**: string, label of x axis, if None, tries to get a label from :func:`default_labels` based on **prop** string
    - **ylabel**: string, label of y axis, if None, tries to get a label from :func:`default_labels` based on **prop** string
    """

    if not xlabel:
        try:
            xlabel = default_labels()[prop2]
        except:
            pass

    if not ylabel:
        try:
            ylabel = default_labels()[prop1]
        except:
            pass

    df = DataFrame(columns=[prop1, prop2])
    dict1, dict2= [g.property(prop) for prop in (prop1, prop2)]
    for vid in dict1.keys():
        try:
            df.loc[vid]=[float(dict1[vid]), float(dict2[vid])]
        except:
            pass

    if not fig:
        fig = mpl.pyplot.figure()
        ax=fig.add_subplot(111)
    else:
        ax=fig.axes[0]
    if not style: style = mpl.pyplot.rcParams['axes.color_cycle'][0]
    df.plot.scatter(prop2,prop1,ax=ax, c=style)
    ax.set(xlim=(min(df[prop2]),max(df[prop2])), xlabel = xlabel,
           ylim=(min(df[prop1]),max(df[prop1])), ylabel = ylabel)

    if one_one:
        ax.plot((min(df[prop2]),max(df[prop2])),(min(df[prop1]),max(df[prop1])),'k')

    return fig, ax
