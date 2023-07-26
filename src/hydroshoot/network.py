"""Utilities to build hydrosoot hydraulic network"""


from openalea.mtg.io import multiscale_edit
import openalea.mtg.traversal as traversal


def setup_network(n_metamer=1, metamer='I[+p<L]', i_kmax=0., i_length=0.1, i_dz=0.1, leaf_length=0.01, leaf_area=0.01,
                  leaf_photosynthesis=None):

    if leaf_photosynthesis is None:
        leaf_photosynthesis = dict(Vcm25=83.00828070320627,
                                   Jm25=180.8788133639022,
                                   TPU25=13.879614097237127,
                                   Rd=1.0534355642371143,
                                   dHd=200.0)

    s = '/' + metamer
    a_metamer = '<' + metamer
    for _ in range(n_metamer - 1):
        s += a_metamer

    g = multiscale_edit(s)

    g.node(g.root).vid_base = 1
    g.node(g.root).vid_collar = 1

    label = g.property('label')
    resistances = [vid for vid in label if label[vid] in ('I', 'p')]
    leaves = [vid for vid in label if label[vid] in ('L',)]

    g.properties()['Kmax'] = {vid: i_kmax for vid in resistances}
    g.properties()['dz'] = {vid: i_dz for vid in resistances}
    g.properties()['length'] = {vid: i_length for vid in resistances}

    g.properties()['leaf_area'] = {vid: leaf_area for vid in leaves}
    g.properties()['length'].update({vid: leaf_length for vid in leaves})
    g.properties()['photosynthesis'] = {vid: leaf_photosynthesis for vid in leaves}

    return g


def reset_fluxes(network):
    g = network
    label = g.property('label')
    g.properties()['Flux'] = {vid: 0 for vid in label}
    g.properties()['FluxC'] = {vid: 0 for vid in label}
    return g

def to_string(g):
    label = g.property('label')
    edge_type = g.property('edge_type')
    s='/'
    for vtx_id in traversal.pre_order2(g, 1):
        lab = label[vtx_id]
        if vtx_id == 1:
            et=''
        else:
            et = edge_type[vtx_id]
        if et =='+':
            s += ('[' + et + lab + ']')
        else:
            s+= (et + lab)
    return s