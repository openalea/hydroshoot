"""This is an example on running HydroShoot on a potted grapevine with a
simple shoot architecture.
"""

from os import getcwd

from openalea.mtg import traversal, mtg
from openalea.plantgl.all import Scene

from hydroshoot import architecture, display, model


def build_mtg(file_name: str) -> (mtg.MTG, Scene):
    grapevine_mtg = architecture.vine_mtg(file_path=file_name)

    for v in traversal.iter_mtg2(grapevine_mtg, grapevine_mtg.root):
        architecture.vine_phyto_modular(grapevine_mtg, v)
        architecture.vine_mtg_properties(grapevine_mtg, v)
        architecture.vine_mtg_geometry(grapevine_mtg, v)
        architecture.vine_transform(grapevine_mtg, v)

    # Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
    mtg_scene = display.visu(grapevine_mtg, def_elmnt_color_dict=True, scene=Scene(), view_result=True)
    return grapevine_mtg, mtg_scene


if __name__ == '__main__':
    g, scene = build_mtg(file_name='grapevine_pot.csv')
    summary_results = model.run(g, str(getcwd()) + '/', scene, psi_soil=-0.5, gdd_since_budbreak=1000.)
