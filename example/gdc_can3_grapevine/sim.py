from pathlib import Path

from openalea.mtg import traversal
from openalea.plantgl.all import Scene

from hydroshoot import architecture, display, model

if __name__ == '__main__':
    path_project = Path(__file__).parent

    # =============================================================================
    # Construct the plant mock-up
    # =============================================================================

    g = architecture.vine_mtg(path_project / 'digit.input')

    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_phyto_modular(g, v)
        architecture.vine_petiole(g, v, pet_ins=90., pet_ins_cv=0.,
                                  phyllo_angle=180.)
        architecture.vine_leaf(g, v, leaf_inc=-45., leaf_inc_cv=100.,
                               lim_max=12.5, lim_min=5., order_lim_max=2.25,
                               max_order=55, rand_rot_angle=90.,
                               cordon_vector=None)
        architecture.vine_mtg_properties(g, v)
        architecture.vine_mtg_geometry(g, v)
        architecture.vine_transform(g, v)

    scene = display.visu(g, def_elmnt_color_dict=True, scene=Scene(), view_result=True)

    # =============================================================================
    # Run HydroShoot
    # =============================================================================

    # =============================================================================
    # Run HydroShoot
    # =============================================================================

    model.run(
        g=g,
        wd=path_project,
        path_weather=path_project / 'meteo.input',
        scene=scene)
