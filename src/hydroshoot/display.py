# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Visualization module for HydroShoot
"""

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import openalea.plantgl.all as pgl
from openalea.mtg.plantframe import color as pglcolor

plt.style.use('ggplot')

from hydroshoot.hydraulic import def_param_soil

DEFAULT_LABELS = {
    'Eabs': r'$\mathregular{E_{abs}\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'Ei': r'$\mathregular{E_i\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'Ei10': r'$\mathregular{PPFD_{10}\/[mol\/m^{-2}\/d^{-1}]}$',
    'Flux': r'$\mathregular{F\/[kg\/s^{-1}]}$',
    'KL': r'$\mathregular{K_L\/[kg\/s^{-1}\/m\/MPa^{-1}]}$',
    'Kmax': r'$\mathregular{K_{max}\/[kg\/s^{-1}\/m\/MPa^{-1}]}$',
    'Length': r'$\mathregular{Diam_{bot}\/[cm]}$',
    'Tlc': r'$\mathregular{T_{leaf}\/[^\circ\/C]}$',
    'TopDiameter': r'$\mathregular{Diam_{top}\/[cm]}$',
    'BotDiameter': r'$\mathregular{Diam_{bot}\/[cm]}$',
    'psi_head': r'$\mathregular{\Psi\/[MPa]}$',
    'An': r'$\mathregular{A_n\/[\mu mol\/m^{-2}\/s^{-1}]}$',
    'gs': r'$\mathregular{g_s\/[mol\/m^{-2}\/s^{-1}]}$',
    'gb': r'$\mathregular{g_b\/[mol\/m^{-2}\/s^{-1}]}$',
    'E': r'$\mathregular{E\/[mol\/m^{-2}\/s^{-1}]}$',
    'Ci': r'$\mathregular{C_i\/[\mu mol_{CO_2}\/\mu mol]}$',
    'Cc': r'$\mathregular{C_c\/[\mu mol_{CO_2}\/\mu mol]}$',
    'k_sky': r'$\mathregular{k_{sky}\/[-]}$',
    'k_soil': r'$\mathregular{k_{soil}\/[-]}$',
    'k_leaves': r'$\mathregular{k_{leaves}\/[-]}$',
    'Na': r'$\mathregular{N_a\/[g\/m^{-2}]}$'
}

DEFAULT_COLORS = {
    'inT': (40, 19, 0),
    'cx': (40, 19, 0),
    'inT3y': (40, 19, 0),
    'inT2y': (40, 19, 0),
    'in': (40, 19, 0),
    'Pet': (255, 0, 0),
    'LI': (0, 110, 44),
    'LII': (0, 165, 5),
    'soil': (7, 102, 0),
    'other': (7, 102, 0)
}

DEFAULT_ELEMENTS = [
    'inT', 'cx', 'inT3y', 'inT2y', 'in', 'Pet', 'LI', 'LII', 'other', 'soil'
]


def visu(g, plot_prop=None, min_value=None, max_value=None, cmap='jet', fmt='%6.0f', elmnt_labels=None,
         elmnt_color_dict=None, def_elmnt_color_dict=False, use_mtg_color=False,
         snap_shot_path=None, scene=None, view_result=True):
    """Displays 3D moke-up using `plantgl` package.

    Args:
        g (openalea.mtg.mtg.MTG): an MTG object
        plot_prop (str): name of an MTG property to plot
        min_value (float): minimum value for the color scale
        max_value (float): maximum value for the color scale
        cmap (str): the colormap to use; only active when `plot_prop` is given or `use_mtg_color` is True (default='jet)
        fmt (str): format of matplotlib colorbar ticklabels (default='%6.0f')
        elmnt_labels (list): the desired MTG elemnts to be displayed
        elmnt_color_dict (dict): RGB color tuples for the given MTG labels
        def_elmnt_color_dict (bool): if True, uses DEFAULT_COLORS` to set a dictionary of RGB color tuples with
            `elmnt_labels` as keys (default=False)
        use_mtg_color (bool): if True, the 'color' property from MTG nodes is used (default=False)
        snap_shot_path (str or None): path to save scene snapshot (default=None)
        scene (Scene): if given, adds scene shapes to it (default=None)
        view_result (bool): if True, the scene is displayed (default=False)

    Returns:
        (Scene)

    """

    my_scene = pgl.Scene() if scene is None else scene

    if plot_prop is not None:
        prop = [g.node(vid).properties()[plot_prop] for vid in g.property(plot_prop)
                if not g.node(vid).label.startswith('soil')]
        if min_value is None:
            min_value = min(prop)
        if max_value is None:
            max_value = max(prop)

        g, cb = pglcolor.colorbar(g, property_name=plot_prop, cmap=cmap, lognorm=False, N=6, fmt=fmt,
                                  min_value=min_value, max_value=max_value)
        #        cb.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        g = pglcolor.colormap(g, property_name=plot_prop, cmap=cmap, lognorm=False,
                              min_value=min_value, max_value=max_value)

        for vid in g.property(plot_prop).keys():
            try:
                n = g.node(vid)
                mesh = n.geometry
                scene_shape = pgl.Shape(mesh, pgl.Material(pgl.Color3(n.color)))
                my_scene.add(scene_shape)
            except:
                pass

    else:

        if elmnt_labels is None:
            elmnt_labels = DEFAULT_ELEMENTS

        if elmnt_color_dict is None:
            if def_elmnt_color_dict:
                def_color_dict = DEFAULT_COLORS
                elmnt_color_dict = {key_: def_color_dict[key_] for key_ in elmnt_labels}
            elif not use_mtg_color:
                raise Exception("Element colors are missing. You may set elmnt_color_dict=True to use default colors"
                                "provied by DEFAULT_COLORS.")

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
                scene_shape = pgl.Shape(mesh, pgl.Material(pgl.Color3(color)))
                my_scene.add(scene_shape)

    if view_result:
        pgl.Viewer.display(my_scene)

    if snap_shot_path:
        pgl.Viewer.saveSnapshot(snap_shot_path)

    return my_scene


def hydraulic_map(g, prop='psi_head', ax=None, **kwargs):
    """Plots the hydraulic cart of an MTG.

    Args:
        g (openalea.mtg.mtg.MTG): an MTG object
        prop (str): name of an MTG property to plot
        ax (AxesSubplot or None): if given adds the hydraulic cart to it (default None)

    Returns:
        (Figure, AxesSubplot)

    """
    x = []
    y = []

    for vid in g.property(prop):
        if not g.node(vid).label.startswith(('soil', 'other')):
            x.append(g.node(vid).TopPosition[2])
            y.append(g.node(vid).properties()[prop])

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.scatter(x, y, **kwargs)
    ax.set(xlabel=r'$\mathregular{Elevation\/[cm]}$', xlim=(min(x), max(x)),
           ylabel=DEFAULT_LABELS['psi_head'], ylim=(min(y), max(y)))

    return fig, ax


def property_map(g, prop='psi_head', ax=None, label=None, add_head_loss=False, color='b', prop2=None, colormap=None,
                 add_color_bar=False):
    """Plots values of a given MTG property vs hight.

    Args:
        g (openalea.mtg.mtg.MTG): an MTG object
        prop (str): name of an MTG property to plot
        ax (AxesSubplot or None): if given adds the hydraulic cart to it (default None)
        label (str or None): label of y axis, if None, a label is sought from DEFAULT_LABELS
        add_head_loss (bool): if True, adds a dashed line representing water head loss due
            to elevation [-0.01 MPa m-1] (default False)
        color (str or None): color name
        prop2 (str or None): name of an MTG property whose values determine the color of the scattered points
        colormap (str or None): name of the colarmap
        add_color_bar (bool): if True, a colarbar is added (only considered if colormap is not None)

    Returns:
        (AxesSubplot)

    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    if label is None and prop in DEFAULT_LABELS.keys():
        label = DEFAULT_LABELS[prop]

    label2 = None
    if colormap is not None:
        assert (prop2 is not None), 'prop2 must be precised.'
        if prop2 in DEFAULT_LABELS.keys():
            label2 = DEFAULT_LABELS[prop2]
        cm = plt.cm.get_cmap(colormap)
    else:
        cm = colormap

    vid_collar = g.node(g.root).vid_collar

    if prop in g.node(vid_collar).properties() and not prop.startswith('k_'):
        for vid in g.Extremities(vid_collar):
            if prop in g.node(vid).properties():
                index = 0
            else:
                index = 1
            y = [g.node(ivid).TopPosition[2] for ivid in g.Ancestors(vid)[index:]]
            x = [g.node(ivid).properties()[prop] for ivid in g.Ancestors(vid)[index:]]
            ax.plot(x, y, '.-', label=label, zorder=0, color=color)
            ax.set(ylabel='z [cm]', xlabel=label)

        if colormap is not None:
            x = []
            y = []
            c = []
            for ivid in g.Extremities(vid_collar):
                x.append(g.node(ivid).properties()[prop])
                y.append(g.node(ivid).TopPosition[2])
                c.append(g.node(ivid).properties()[prop2])
            im = ax.scatter(x, y, c=c, vmin=min(c), vmax=max(c), cmap=cm)
            ax.im = im
            if add_color_bar:
                fig.colorbar(im, ax=ax, label=label2)

    else:
        x = [g.node(ivid).properties()[prop] for ivid in g.VtxList(Scale=3) if g.node(ivid).label.startswith('L')]
        y = [g.node(ivid).TopPosition[2] for ivid in g.VtxList(Scale=3) if g.node(ivid).label.startswith('L')]
        if colormap is None:
            ax.plot(x, y, '.', label=label, zorder=0, color=color)
        else:
            c = [g.node(ivid).properties()[prop2] for ivid in g.VtxList(Scale=3) if
                 g.node(ivid).label.startswith('L')]
            im = ax.scatter(x, y, c=c, vmin=min(c), vmax=max(c), cmap=cm)
            if add_color_bar:
                fig.colorbar(im, ax=ax, label=label2)

        ax.set(ylabel='z [cm]', xlabel=label)

    if add_head_loss:
        ylim = ax.get_ylim()
        ls = np.arange(ylim[0], ylim[1])
        ax.plot(0.01 * ls * (-0.01) + g.node(vid_collar).psi_head, ls, '--', label='Hydrostatic slope')

    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    fig.tight_layout()

    return ax


def prop_fun_prop(g, y_prop='gs', x_prop='psi_head', ax=None, color=None, one_one=False, xlabel=None, ylabel=None):
    """Returns a scatter plot of two MTG properties.

    Args:
        g (openalea.mtg.mtg.MTG): an MTG object
        y_prop (str): name of the MTG property to plot on the y-axis
        x_prop (str): name of the MTG property to plot on the x-axis
        ax (AxesSubplot or None): if given adds the hydraulic cart to it (default None)
        color (str or None): color name
        one_one (bool): if True, adds a 1:1 line to figure (default False)
        xlabel (str or None): label of x-axis, if None, a label is sought from DEFAULT_LABELS
        ylabel (str or None): label of y-axis, if None, a label is sought from DEFAULT_LABELS

    Returns:
        (Figure, AxesSubplot)

    """
    y_dict, x_dict = [g.property(prop) for prop in (y_prop, x_prop)]

    xy = []
    for vid in y_dict.keys():
        try:
            xy.append((float(x_dict[vid]), float(y_dict[vid])))
        except TypeError:
            pass
    x, y = zip(*xy)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.scatter(x, y, **{'c': color} if color is not None else {})

    kwargs = dict(
        xlabel=DEFAULT_LABELS[x_prop] if xlabel is None and x_prop in DEFAULT_LABELS.keys() else xlabel,
        ylabel=DEFAULT_LABELS[y_prop] if ylabel is None and y_prop in DEFAULT_LABELS.keys() else ylabel,
        xlim=(min(x), max(x)),
        ylim=(min(y), max(y))
    )

    ax.set(**kwargs)

    if one_one:
        ls = min(min(x, y)), max(max(x, y))
        ax.plot(ls, ls, 'k--')

    fig.tight_layout()
    return fig, ax


def retention_curve(g, ax=None):
    """Plots the retention curve of the bulk soil:

    Args:
        g (openalea.mtg.mtg.MTG): an MTG object
        ax (AxesSubplot or None): if given adds the hydraulic cart to it (default None)

    Returns:
        (AxesSubplot)

    """

    if not ax:
        fig, ax = plt.subplots()

    plt.xscale('log')
    soil_class = g.node(g.Ancestors(g.node(g.root).vid_base)[0]).soil_class
    param = def_param_soil()[soil_class]
    theta_r, theta_s, alpha, n, k_sat = [param[i] for i in range(5)]
    m = 1. - 1. / n

    psi_range = np.arange(0, 150000)
    theta_ls = []
    for psiX in psi_range:
        theta_ls.append(theta_r + (theta_s - theta_r) * 1. / ((1 + abs(alpha * psiX)) ** n) ** m)

    psi_range = np.array(psi_range) * 1.e-4
    ax.plot(psi_range, theta_ls, label=soil_class)
    ax.set(xlabel=r'$\mathregular{-\Psi_{soil}\/[MPa]}$',
           ylabel=r'$\mathregular{\Theta_{bulk\/soil}\/[-]}$')

    ax.legend()
    return ax
