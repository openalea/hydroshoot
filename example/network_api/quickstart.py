from openalea.mtg import mtg
import json
from hydroshoot.solver import solve_interactions
from hydroshoot.initialisation import set_collar_water_potential_function
from copy import deepcopy
import openalea.mtg.traversal as traversal
from hydroshoot.params import Soil

import hydroshoot.network as network

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    def __getattr__(*args):
        val = dict.get(*args)
        return DotDict(val) if type(val) is dict else val

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __deepcopy__(self, memo=None):
        return DotDict(deepcopy(dict(self), memo=memo))


gg = network.setup_network()

g = mtg.MTG()
vid = g.add_component(g.root, label='inT')
g.node(g.root).vid_base=vid
g.node(g.root).vid_collar=vid
g.node(vid).dz = 0.1
g.node(vid).Kmax = 1
g.node(vid).length = 0.1

vid = g.add_child(vid, label='LI', edge_type='<')
g.node(vid).leaf_area = 1
g.node(vid).length = 1

# fraction of soil, leaves and sky seen by element (sum=2)
g.node(vid).ff_sky = 0.6674578719442436
g.node(vid).ff_leaves = 0.46764999380584427
g.node(vid).ff_soil = 0.8648921342499122

# micrometeo
g.node(vid).u = 0.03
g.node(vid).Tac = 20.74
g.node(vid).Tlc = 20.74
g.node(vid).vpd = 1.1335927884740755
g.node(vid).Ei = 500.0  # ppfd
g.node(vid).Eabs = 0.0  # ppfd
g.node(vid).Rg = 50.0  # W/m2

for vtx_id in traversal.pre_order2(g, g.node(g.root).vid_base):
    g.node(vtx_id).Flux = 0.
    g.node(vtx_id).FluxC = 0.

with open('params.json') as f:
    user_params = json.load(f)

user_params['simulation']['conv_to_meter']={'mm': 1.e-3, 'cm': 1.e-2, 'm': 1.}[user_params['simulation']['unit_scene_length']]
user_params['simulation']['conv_to_second']=3600.
user_params['exchange']['par_photo'].update(dict(Vcm25=83.00828070320627,
                           Jm25=180.8788133639022,
                           TPU25=13.879614097237127,
                           Rd=1.0534355642371143,
                           dHd=200.0))
for w in ('hydraulic_structure', 'negligible_shoot_resistance', 'energy_budget'):
    user_params['simulation']['is_' + w] = user_params['simulation'].pop(w)
user_params['soil']['soil_dimensions']['radius']=0.1
user_params['soil']['soil_volume'] = Soil.calc_soil_volume(soil_dimensions=user_params['soil']['soil_dimensions'])
user_params['soil']['rhyzo_volume'] = user_params['soil']['soil_volume']
params = DotDict(user_params)

meteo = dict(hs=78,
             Pa=101.325,
             Ca=600)
clearness=0.5
tsky = clearness * params.energy.t_sky + (1-clearness) * params.energy.t_cloud
calc_collar_water_potential = set_collar_water_potential_function(params=params)

solve_interactions(g, meteo, psi_soil=-0.2, t_soil=20, t_sky_eff=tsky, params=params, calc_collar_water_potential=calc_collar_water_potential)