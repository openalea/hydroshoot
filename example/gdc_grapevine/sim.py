"""This is an example on running HydroShoot on a potted grapevine with a
simple shoot architecture.
"""

from os import getcwd
import matplotlib as mpl

import openalea.mtg.traversal as traversal
from openalea.plantgl.all import Scene
import hydroshoot.architecture as HSArc
import hydroshoot.display as HSVisu
from hydroshoot.model import run

# =============================================================================
# Construct the plant mock-up
# =============================================================================

# Path for plant digitalization data.
g = HSArc.vine_mtg('digit.input')

# Local Coordinates Correction
for v in traversal.iter_mtg2(g, g.root):
    architecture.vine_phyto_modular(g, v)
    architecture.VinePetiole(g, v, pet_ins=90., pet_ins_cv=0., phyllo_angle=180.)
    architecture.VineLeaf(g, v, leaf_inc=-45., leaf_inc_cv=100., lim_max=12.5, lim_min=5.,
                          order_lim_max=6., max_order=55, rand_rot_angle=90.,
                          cordon_vector=None)
    architecture.VineMTGProp(g, v)
    architecture.VineMTGGeom(g, v)
    architecture.VineTransform(g, v)

# Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
# scene = HSVisu.visu(g,def_elmnt_color_dict=True,scene=Scene(),
#                    snap_shot_path='mockup.png')
scene = HSVisu.visu(g, def_elmnt_color_dict=True, scene=Scene(),
                    view_result=True)

# =============================================================================
# Run HydroShoot
# =============================================================================

run(g, str(getcwd()) + '/', scene)
