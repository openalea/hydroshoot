from pathlib import Path
from os import getcwd

from openalea.mtg import traversal
from openalea.plantgl.all import Scene

from openalea.hydroshoot import architecture, display, model

# =============================================================================
# Construct the plant mock-up
# =============================================================================


if __name__ == '__main__':
    g = architecture.vine_mtg(Path(__file__).parent / 'HSdanugue$5132.csv')

    # Local Coordinates Correction
    for v in traversal.iter_mtg2(g, g.root):
        n = g.node(g.Trunk(v, Scale=1)[0])
        theta = 180 if int(n.index()) < 200 else -90 if int(n.index()) < 300 else 0
        architecture.vine_orientation(g, v, theta, local_rotation=True)

    # Scene rotation
    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_orientation(g, v, 90., local_rotation=False)

    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_phyto_modular(g, v)
        architecture.vine_mtg_properties(g, v)
        architecture.vine_mtg_geometry(g, v)
        architecture.vine_transform(g, v)

    # Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
    scene = display.visu(g, def_elmnt_color_dict=True, scene=Scene(), view_result=True)

    # =============================================================================
    # Run HydroShoot
    # =============================================================================

    model.run(g, str(getcwd()) + '/', scene, psi_soil=-0.5, gdd_since_budbreak=1000.)




from pandas import read_csv
df = read_csv(r'C:\Users\albashar\Documents\dvp\grapesoil\grape_soil\example\data\benoit\data_meteo_with_low_WS.csv',
              sep=';', decimal='.')