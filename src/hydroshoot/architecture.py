# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Plant architecture module of HydroShoot.

This module constructs a multiscale tree graph (MTG) of grapevines based on
digitalization data using the `openalea.mtg` package.

The resulting MTG incorporates geometry
"""

import scipy
from scipy.linalg import norm
#from numpy.linalg.linalg import norm
from numpy.linalg import det
from scipy.spatial import distance
from roman import toRoman
from pandas import read_csv
from re import search, findall
from itertools import product
from pickle import dump, load
from os import path, mkdir


from openalea.mtg import mtg, io
from openalea.plantgl.all import Point3Array
import openalea.plantgl.all as pgl

#==============================================================================
# Functions from TopVine package
#==============================================================================
def cart_to_pol (coordxy) :
    """
    Converts cartesian coordinates (x,y,z) to polar corrdinates (r,azi,incli).
    """

    x,y,z = coordxy[0], coordxy[1], coordxy[2]
    r = scipy.sqrt(x*x+y*y+z*z)
    if r==0 :
        incli =0
    else :
        incli = scipy.arcsin(z/r)
    if (x==0 and y==0):
        azi = 0
    elif (y>=0) :
        azi = scipy.arccos(x/scipy.sqrt(x*x+y*y))
    else :
        azi = -scipy.arccos(x/scipy.sqrt(x*x+y*y))

    return scipy.array([r,azi,incli])


def PolToXyz (coordpol) :
    """
    Converts polar corrdinates (r,azi,incli) to carthesian coordinates (x,y,z).
    """

    r,azi,incli = coordpol[0], coordpol[1], coordpol[2]
    z = r * scipy.sin(incli)
    l = scipy.sqrt (r*r-z*z)
    x = l * scipy.cos(azi)
    y = l * scipy.sin (azi)

    return scipy.array([x, y, z])


def transformation(obj, sx, sy, sz, rx, ry, rz, tx, ty, tz ): 
    """
    Returns a scaled, rotated and translated 3D object - Similar to 'transformation' in PovRay
    """

    s_obj = pgl.Scaled (pgl.Vector3(sx,sy,sz), obj)
    r_obj = pgl.EulerRotated (rx, ry, rz, s_obj)
    t_obj = pgl.Translated (pgl.Vector3(tx,ty,tz), r_obj)

    return t_obj


def leaf0(l=1.):
    points= Point3Array([ pgl.Vector3(0,9.8,0),pgl.Vector3(3.7,6.3,0),pgl.Vector3(2.3,3.9,0), 
                          pgl.Vector3(5.4,6.2,0), pgl.Vector3(5.9,3.2,0),pgl.Vector3(4.6,0.1,0),
                          pgl.Vector3(6.4,-1.6,0), pgl.Vector3(1,-3.7,0), pgl.Vector3(0,0,0),                     
                          pgl.Vector3(-1,-3.7,0),pgl.Vector3(-6.4,-1.6,0), 
                          pgl.Vector3(-4.6,0.1,0), pgl.Vector3(-5.9,3.2,0),                        
                          pgl.Vector3(-5.4,6.2,0), pgl.Vector3(-2.3,3.9,0),pgl.Vector3(-3.7,6.3,0)]) 
    
    indices= pgl.Index3Array([ pgl.Index3(0,1,8), pgl.Index3(15,0,8), pgl.Index3(8,11,14), 
                           pgl.Index3(11,12,14), pgl.Index3(12,13,14),pgl.Index3(10,11,8),                      
                           pgl.Index3(8,9,10), pgl.Index3(2,3,4),pgl.Index3(2,4,5), 
                           pgl.Index3(2,5,8), pgl.Index3(8,5,6),pgl.Index3(8,6,7)])
    
    f= pgl.TriangleSet(points, indices)
    return pgl.Scaled (pgl.Vector3(0.01*l,0.01*l,0.01*l), f)


def soil0(l=1.):
    points= Point3Array([pgl.Vector3(0.,0.,0.),pgl.Vector3(1.,0.,0.),pgl.Vector3(1.,1.,0.),pgl.Vector3(0.,1.,0.)])
    indices= pgl.Index3Array([pgl.Index3(0,1,2), pgl.Index3(0,2,3)])
    f= pgl.TriangleSet(points, indices)
    return pgl.Scaled (pgl.Vector3(l,l,l), f)

#==============================================================================
# MTG topology construction
#==============================================================================

def VineMTG(file_path):
    """
    Constructs a MultiScale Tree Graph (MTG) for digitilized grapvines.
    The data structure must be as follows:
    **Plant_ID**, **Axis_ID**, **inT_ID**, **inT3Y_ID**, **in_ID**, **in_order**, **x**, **y**, **z**[, **Diameter**]

    :Parameters:
    - **Plant_ID**: integer, the plant ID
    - **Axis_ID**: integer, the indentifier of permenant axes (0: trunk; 1, 2, etc.: cordons or arms)
    - **inT_ID**: a single internode or a "pruning complex" both on permanent axes (lignified)
        - Integer: a single internode
        - String: a pruning complex (noted as X.a, X,b, etc. or X.a.1, X.a.2, X.b.1, X.b.2, etc., where X is the number of the main axis internode)
    - **inT3Y_ID**: a single 3-year old internode (cane)
        - integer (1, 2, etc., single internode)
        - float (1.1, 1.2, etc., single of a sequence of internodes)
    - **in_ID**: Integer or string identifying a 1-year-old internode
        - integer: internode on a primary shoot
        - string: e.g. 2.5_1: a secondary shoot or order `1` attached to the fifth internode `5` of the second `2` primary internode (this holds for secondary, tertiary, etc. shoots)
    - **in_order**: integer, the order or the internodes (starts with 0 only for primary internodes, otherwise starts with 1 for secondary, tertiary, etc. internodes)
    - **x**, **y**, **z**: float, Cartesian coordinates of internode upper end
    - **Diameter**: float, internode diameter at the upper end

    :Returns:
    - An MTG object.
    """
#cep      -> plant
#tronc    -> trunk
#phyto    -> inT
#'phyto'  -> 'inT'
#'Cx'     -> 'cx'
#sarement -> inT3y
#en_s     -> inT3y
#'enS'    -> 'inT3y'
#en_t     -> inT2y
#'enT'    -> 'inT2y'
#'R'      -> 'shI'
#en_r     -> inI
#enR      -> 'inI'
#enC      -> 'inII'

#perennial_axis -> perennial_axis
#connect_s -> connect_inT3y
#connect_t -> connect_inT2y
#connect_r -> connect_inI

    table = read_csv(file_path,sep=';',decimal='.',header=0)

    g = mtg.MTG()

    plant_id_prev = -1
    trunk_id_prev = -1
    inT_id_prev = -1
    complexe_id_prev = '-1.0'
    inT3y_id_prev = -1

    trunk_lbl_prev = '-1_-1'
    inT_lbl_prev = '-1_-1_-1'
    complexe_lbl_prev = '-1_-1_-1.ab'
    inT3y_lbl_prev = '-1_-1_-1.ab_-1'

    ind_leaf_points = []

    B_A_S_E = map(''.join, product(*((c.upper(), c.lower()) for c in 'base')))

    for i in range(len(table)):
        plant_id, trunk_id, elmnt_id, inT3y_id, shoot_id, order = [table[table.columns[j]][i] for j in range(0,6)]
        vid_position = [round(table[table.columns[j]][i],2) for j in range(6,9)]
        try:
            vid_diam = table[table.columns[9]][i]
        except:
            vid_diam = None

        plant_id = int(plant_id)
        if trunk_id in B_A_S_E:
            baseXYZ = vid_position
            trunk_id = 0
        else:
            trunk_id = int(trunk_id)


        inT_cx = bool('.' in str(elmnt_id)) # TRUE if the element is a pruning complex

        if plant_id != plant_id_prev:
            plant = g.add_component(g.root, label=('plant'+str(plant_id)), edge_type='/')
            if 'baseXYZ' in locals(): g.node(plant).baseXYZ = baseXYZ
            plant_ids = (plant_id_prev, plant_id)
            plant_id_prev = plant_id

#       Perennial axes ********************************************************
        trunk_lbl = str(plant_id) + '_' + str(trunk_id)
        if trunk_lbl != trunk_lbl_prev:
            if trunk_lbl.split('_')[0] != trunk_lbl_prev.split('_')[0]:
                assert (trunk_id == 0), "Error! The first internode of the plant is not in the trunk !!!"
                perennial_axis = g.add_component(plant, label='trunk')
                trunk = perennial_axis
            else:
                perennial_axis = g.add_child(trunk, label=('arm'+str(trunk_id)), edge_type='+')

            trunk_id_prev = trunk_id
            trunk_lbl_prev = trunk_lbl

#       Internodes and pruning complexes of the permanent axes ****************
        elmnt_split = str(elmnt_id).split('.') if str(elmnt_id) not in ('0', '0.0') else []

        if len(elmnt_split) == 1: # Simple internode of perennial axes
            try: # TODO: optimize
                inT_id = float(elmnt_id)
                inT_lbl = str(plant_id) + '_' + str(trunk_id) +'_'+ str(elmnt_id)
            except ValueError:
                inT_id = elmnt_id
                inT_lbl = str(plant_id) + '_' + str(trunk_id) +'_'+ str(elmnt_id)

            if inT_lbl != inT_lbl_prev:
                if inT_lbl.split('_')[0] != inT_lbl_prev.split('_')[0]: # If not the same plant=> add first internode to trunk
                    inT = g.add_component(perennial_axis, label=('inT'+str(elmnt_id)), TopPosition=vid_position, TopDiameter=vid_diam)
                else:
                    if inT_lbl.split('_')[1] != inT_lbl_prev.split('_')[1]:
                        if inT_lbl_prev.split('_')[1] == '0' : inT_t = inT
                        inT = g.add_component(perennial_axis, label=('inT'+str(elmnt_id)), TopPosition=vid_position, TopDiameter=vid_diam)
                        if inT_lbl.split('_')[1] != '0' : inT = g.add_child(inT_t, child=inT, edge_type='+', TopPosition=vid_position, TopDiameter=vid_diam)
                    else:
                        if True in [c in str(elmnt_id) for c in ('rl','RL','Rl','rL')]:
                            connect_p = complexe
                            edge_type_e = '+'
                        else:#if any((c in ('rl','RL','Rl','rL')) for c in elmnt_id):
                            connect_p = inT
                            edge_type_e = '<'

                        inT = g.add_child(connect_p, label=('inT'+str(elmnt_id)), edge_type=edge_type_e, TopPosition=vid_position, TopDiameter=vid_diam)
                inT_id_prev = inT_id
                inT_lbl_prev = inT_lbl

#        Pruning complex ******************************************************
        elif len(elmnt_split) > 1:
            complexe_id = elmnt_id
            complexe_lbl = str(plant_id) + '_' + str(trunk_id) +'_'+ str(elmnt_id)
            if complexe_lbl != complexe_lbl_prev:
                if complexe_lbl.split('.')[0] == complexe_lbl_prev.split('.')[0] and complexe_lbl.split('.')[1] == complexe_lbl_prev.split('.')[1]:
                    connect_c = complexe
                else:
                    connect_c = inT
                complexe = g.add_child(connect_c, label=('cx'+str(elmnt_id)), edge_type='+', TopPosition=vid_position, TopDiameter=vid_diam)
                complexe_id_prev = complexe_id
                complexe_lbl_prev = complexe_lbl


#       Internodes of 3-year-old axes *****************************************
        if float(inT3y_id) != 0.0:
            inT3y_lbl = str(plant_id) + '_' + str(trunk_id) +'_'+ str(elmnt_id) +'_'+ str(inT3y_id)
            if inT3y_lbl != inT3y_lbl_prev:
                if inT3y_lbl.split('_')[:3] == inT3y_lbl_prev.split('_')[:3]:
                    if int(float(inT3y_lbl.split('_')[3])) == int(float(inT3y_lbl_prev.split('_')[3])): # In series
                        connect_inT3y = inT3y
                        edge_type_s = '<'
                else: # sarements in parallel
                    connect_inT3y = complexe if inT_cx else inT
                    edge_type_s = '+'

                inT3y = g.add_child(connect_inT3y, label=('inT3y'+str(inT3y_id)), edge_type=edge_type_s, TopPosition=vid_position, TopDiameter=vid_diam)
                inT3y_id_prev = inT3y_id
                inT3y_lbl_prev = inT3y_lbl


#       Internodes of 2- and 1-year-old axes **********************************
        if str(shoot_id) not in ['0','0.0'] and bool(search('[pPfFlL]', str(shoot_id))) == False:
            shoot_order = len(str(shoot_id).split('.'))

            if shoot_order == 1:    # Primary shoot

                connect_inT2y = -1
                label_shI = 'shI'

                if inT_cx == False and float(inT3y_id) == 0.0 and bool(search('[gG]+',str(shoot_id))):
                    connect_inI = inT       # Shoot>inT
                    label_shI = 'GI'

                elif inT_cx == False and float(inT3y_id) == 0.0 and bool(search('[gG]+',str(shoot_id))) == False:
                    connect_inT2y = inT       # Shoot>spur>inT

                elif float(inT3y_id) != 0.0 and bool(search('[gG]+',str(shoot_id))):
                    connect_inI = inT3y        # Shoot>sarement
                    label_shI = 'GI'

                elif float(inT3y_id) != 0.0 and bool(search('[gG]+',str(shoot_id))) == False:
                    connect_inT2y = inT3y        # Shoot>Spur>sarement

                elif inT_cx == True and float(inT3y_id) == 0.0:
                    connect_inI = complexe    # Shoot>complex
                    label_shI = 'GI'

                shoot_id2 = float(findall('\d+',str(shoot_id))[0])


                if order == 0:
                    if connect_inT2y != -1:
                        if str(shoot_id2) in('1','1.0'):
                            edge_type_t = '+'
                        else:
                            edge_type_t = '<'
                            connect_inT2y = inT2y
                        inT2y = g.add_child(connect_inT2y, label=('inT2y'+str(shoot_id)), edge_type=edge_type_t, TopPosition=vid_position, TopDiameter=vid_diam)
                        connect_inI = inT2y
                else:
                    if order == 1:
                        inI, shoot = g.add_child_and_complex(connect_inI, label=('inI'+str(order)), edge_type='+', TopPosition=vid_position, TopDiameter=vid_diam)
                        g.node(shoot).label = label_shI+str(shoot_id)
                        g.node(shoot).edge_type = '+'
                    else:
                        inI = g.add_child(inI, label=('inI'+str(order)), edge_type='<', TopPosition=vid_position, TopDiameter=vid_diam)


            elif shoot_order >= 1:    # Secondary, tertiary, quaternary, etc. shoots
                label_inode = 'in' + toRoman(shoot_order) + str(order)
                in_connect_orders = [int(x.split('_')[0]) for x in str(shoot_id).split('.')[1:]]

                if order == 1:    # A new ramification
                    in_iter = inI
                    for i2, in_order in enumerate(in_connect_orders):
                        connect_inode = g.Axis(in_iter)[in_order-1]
                        try :
                            in_iter = g.Sons(connect_inode, EdgeType='+')[0]
                        except:
                            pass

                    inode, shoot = g.add_child_and_complex(connect_inode, label=label_inode, edge_type='+', TopPosition=vid_position, TopDiameter=vid_diam)
                    g.node(shoot).edge_type = '+'
                    g.node(shoot).label = 'sh' + toRoman(shoot_order) + str(shoot_id).split('.')[shoot_order-1].split('_')[1]
                else:
                    inode = g.add_child(inode, label=label_inode, edge_type='<', TopPosition=vid_position, TopDiameter=vid_diam)

#       Petioles***************************************************************
        if bool(search('[pP]', str(shoot_id))):
            try: # implies that leaves must follow directly their holding internodes
                pet_connect = max(inI, inode)
            except NameError:
                pet_connect = inI
            pet_label = 'Pet' + str(g.node(pet_connect).index())
            petiole = g.add_child(pet_connect, label=pet_label, edge_type='+', TopPosition=vid_position, TopDiameter=vid_diam)

#       Leaves*****************************************************************
        if bool(search('[fFlL]', str(shoot_id))):
            ind_leaf_points.append(vid_position)
            if bool(search('5', str(shoot_id))):
                leaf_label = 'LI' + str(g.node(pet_connect).index())
                leaf = g.add_child(petiole, label=leaf_label, edge_type='+', TopPosition=ind_leaf_points[2],TopPositionPoints=scipy.array(ind_leaf_points))
                ind_leaf_points=[]

    return g


def VinePhytoModular(g,v, *args):
    """
    Identifies the type of phytometer according to Louarn et al. (2007).
    If not provided, 'P0 to P1' and 'P1 to P2' probabilites are given as:
    P0toP1 = [1,0.84,0.8,0.97,1]
    P1toP2 = [0.98,0.96,0.73,0.,1]
    """

    if 'P0toP1' not in args:
        P0toP1 = [1,0.84,0.8,0.97]
        Pzones_nb = len(P0toP1)
    if 'P1toP2' not in args:
        P1toP2 = [0.98,0.96,0.73,0.]

    n = g.node(v)
    try:
        if n.label.startswith('inI1'):
            vid_axis = g.Axis(v)
            P0_counter = 1
            P1_counter = 0
            n.PhytoType = 0
            for i in vid_axis[1:]:
                vid = g.node(i)
                rand_val = scipy.rand()
                if vid.parent().PhytoType == 0:
                    if P0_counter < Pzones_nb and rand_val > P0toP1[P0_counter]:
                            vid.PhytoType = 0
                            P0_counter += 1
                    else:
                        vid.PhytoType = 1
                        P1_counter += 1
                elif vid.parent().PhytoType == 1:
                    if P1_counter < Pzones_nb and rand_val > P1toP2[(P1_counter-1)]:
                            vid.PhytoType = 1
                            P1_counter += 1
                    else:
                        vid.PhytoType = 2
                elif vid.parent().PhytoType == 2:
                    vid.PhytoType = 0
                    P0_counter += 1
    except:
        pass
    return g


def VineNFII(in_order, pruning_type='avg_field_model',N_init=0.18,N_max=2.25,N_max_order=4,in_order_max=25,sig_slope=5.7,phyto_type='P0'):
    """
    Returns NFII, the total number of secondary phytomers per primary internode.

    :Parameters:
    - **in_order**: integer, internode order
    - **pruning_type**: string, the pruning type, one of the following: 'avg_field_model','GDC_1', 'Lyre', 'Rideau_simple', 'VSP_HL', 'Lyre_ouverte', 'GDC_2', 'Pot'
    - **N_init**: float, the average number of secondary internodes which are connected to the first primary internode
    - **N_max**: float, the average maximum number of secondary internodes
    - **N_max_order**: float, the order of the primary internode which has N_max
    - **in_order_max**: float, the order of the primary internode where the number of secondary internodes is assumed equal to zero
    - **sig_slope**: float, the slope (b) of the sigmoidal relationshipe NFII=f(in_order), a modified formula of that proposed by Louarn (2005, PhD, Eq. II.8)
    - **phyto_type**: string, the type of the primary phytomere. Can be one of 'P0', 'P1' or 'P2' (Louarn et al., 2007)

    :Notes:
    - **'GDC_1'**, **'Lyre'**, **'Rideau_simple'**, **'VSP_HL'**, are data collected in 2009 from well-irrigated field-grown Syrah vines at INRA-Montpellier.
    - **'Lyre_ouverte'**, **'GDC_2'** are data collected in 2012 from well-irrigated field-grown Syrah vines at INRA-Montpellier.
    - **'Pot'** are data collected in 2002 from well-irrigated pot-grown Syrah vines at INRA-Montpellier.
    """

    GDC_1 = {1:[0,0],2:[0.1,0.2],3:[2.9,1.35],4:[3.7,2.19],5:[4.5,2.87],6:[3.5,2.05],7:[5,3.17],8:[3.8,2.2],9:[3.1,1.61],10:[0.8,0.7],11:[0.2,0.26],12:[0.1,0.2],13:[0.1,0.2],14:[0.1,0.2]}
    Lyre = {1:[0,0],2:[0.4,0.6],3:[1.2,0.87],4:[3.4,1.94],5:[2.3,1.6],6:[3.5,1.81],7:[2.5,1.02],8:[3.4,2.27],9:[3.2,1.96],10:[2.3,1.6],11:[2.5,0.73],12:[2.8,1.82],13:[3.2,1.46],14:[2.1,1.22],15:[2.7,1.46],16:[2,1.27],17:[3,1.98],18:[1.8,1.84],19:[0.6,0.78],20:[0.6,0.67],21:[0.8,1.16],22:[0.4,0.6],23:[0.5,0.98]}
    Rideau_simple = {1:[0.05,0.1],2:[0.5,0.44],3:[1.75,0.72],4:[3.55,1.72],5:[3.2,1.06],6:[3.85,1.57],7:[3.35,0.82],8:[2.75,0.99],9:[2.47,0.99],10:[2.13,0.5],11:[2.07,0.4],12:[2,0.65],13:[1.7,0.66],14:[2,1.39],15:[1,1.13],16:[1,0],17:[1.5,0.98],18:[2,1.96]}
    VSP_HL = {1:[0.3,0.59],2:[0.3,0.59],3:[1.6,1.88],4:[5.2,3.5],5:[2.6,1.5],6:[2.4,1.68],7:[2.6,2.11],8:[3,2.79],9:[2.1,1.22],10:[6.8,3.9],11:[5,2.83],12:[4.7,3.32],13:[5,2.41],14:[6.4,4.08],15:[1.1,1.32],16:[0.3,0.59],17:[0.4,0.78],18:[0.3,0.59],19:[0.3,0.59],20:[0.2,0.39],21:[0,0],22:[0,0],23:[0,0]}
    Lyre_ouverte = {1:[0.16,0.11],2:[0.63,0.17],3:[1.27,0.26],4:[1.47,0.35],5:[1.51,0.26],6:[1.38,0.24],7:[1.67,0.37],8:[1.22,0.35],9:[1.1,0.33],10:[0.88,0.28],11:[0.73,0.28],12:[0.67,0.3],13:[0.96,0.37],14:[1.45,0.72],15:[1.61,1.29],16:[0.76,0.74],17:[1.29,1.19],18:[1.21,1.25],19:[1.1,0.9],20:[1.2,1.44],21:[6,1.96]}
    GDC_2 = {1:[0.37,0.21],2:[1.12,0.32],3:[1.95,0.36],4:[2.03,0.36],5:[1.96,0.34],6:[1.77,0.32],7:[1.7,0.36],8:[1.34,0.27],9:[1.06,0.22],10:[1.06,0.26],11:[0.75,0.21],12:[0.52,0.19],13:[0.38,0.16],14:[0.15,0.12],15:[0.1,0.11],16:[0.13,0.19]}

    NFII_dict = {'GDC_1':GDC_1, 'Lyre':Lyre, 'Rideau_simple':Rideau_simple, 'VSP_HL':VSP_HL, 'Lyre_ouverte':Lyre_ouverte, 'GDC_2':GDC_2}
    Phyto_type_dict = {'P0':1.,'P1':0.63,'P2':0.63}

    N_init,N_max,N_max_order,in_order_max = map(float,(N_init,N_max,N_max_order,in_order_max))
    a1 = (N_max - N_init)/(N_max_order - 1)
    b1 = N_init
    a2 = -N_max/(in_order_max - N_max_order)
    b2 = N_max*(1+N_max_order/(in_order_max - N_max_order))

    CI_1 = 0.42 # Average confidence interval for NFII observed in field (procedure to be improved)
    CI_2_coef = 0.007*in_order+0.193 # Average confidence interval for NFII observed in pots (procedure to be improved)

    if pruning_type in NFII_dict.keys():
        try:
            NFII = NFII_dict[pruning_type][in_order][0] + NFII_dict[pruning_type][in_order][1]*(min(1,max(-1,scipy.randn()/2.96)))
        except:
            NFII = 0
    elif pruning_type == 'Pot':
        N_reduction_slope = 0.0055 if phyto_type == 'P2' else 0.
        NFII_init = a1*in_order + b1 if in_order <= N_max_order else 2*N_max*(1-1/(1+scipy.exp(-(in_order-(N_max_order-1))/sig_slope)))
        NFII_init = NFII_init*(1+CI_2_coef*(min(1,max(-1,scipy.randn()/2.96))))
        reduction_coeff = N_reduction_slope*in_order+Phyto_type_dict[phyto_type]
        NFII=reduction_coeff*NFII_init
    else : # default_pruning_type = 'avg_field_model'
        NFII_init = a1*(in_order-1) + b1 if in_order <= N_max_order else a2*in_order + b2
        NFII_init = max(0.,NFII_init)
        NFII = NFII_init+CI_1*(min(1,max(-1,scipy.randn()/2.96)))

    NFII = int(round(NFII,0))

    return NFII


def VineLII(NFII, pruning_type='avg_field_model', a_L=43.718, b_L=-37.663, a_P=1.722, b_P=10.136, c_P=-5.435):
    """
    Returns LII [mm], the total length of secondary axis.

    :Parameters:
    - **NFII**: integer : the total number of secondary phytomers per primary internode
    - **pruning_type**: string, the pruning type, one of the following: 'avg_field_model','GDC_1', 'Lyre', 'Rideau_simple', 'VSP_HL', 'Lyre_ouverte', 'GDC_2', 'Pot'
    - **a_L**, **b_L**: float, the slope and intercept, respectively, of the linear relationship LenII=f(NFII)
    - **a_P**, **b_P**, **c_P**: floats, the coefficients of the polynomial relationship LenII=f(NFII)
    """

    if pruning_type in ('avg_field_model','GDC_1', 'Lyre', 'Rideau_simple', 'VSP_HL', 'Lyre_ouverte', 'GDC_2'):
        LII = max(0., a_L*NFII + b_L)
    elif pruning_type == 'Pot':
        LII = a_P*NFII**2 + b_P*NFII + c_P

    return LII


def VineInL(phyto_num, tot_len):
    """
    Returns the length of an individual internode.

    :Parameters:
    - **phyto_num**: integer, the total number of phytomers of an axis.
    - **tot_len**: the total length of an axis.
    TODO: calculate the length of an internode as a function of its type or order
    """

    if phyto_num != 0:
        inL = tot_len/phyto_num
    else:
        inL = 0.

    return inL


def VineAxisCurv(incli_init, length, Fifty_cent=400., sig_slope=70.,curv_type='convexe'):
    """
    Returns the total curvature [radian] of an axis based on its initial inlination and its length.

    :Parameters:
    - **incl_init**: float, the initial inclination of the axis [pi/2.,-pi/2.]
    - **length**: float, the total length of the axis [mm]
    - **Fifty_cent**: float, the abscice of inflection point of the sigmoidal curve
    - **sig_slope**: float, the slope of the sigmoidal relationship curv=f(length)
    - **curv_type**: in plant exterior direction even 'convex' or 'concave'
    """

    curv_max = scipy.pi/2. + incli_init
    curv_tot = curv_max/(1.+scipy.exp(-(length-Fifty_cent)/sig_slope)) # TODO: Rule to be confirmed
    if curv_type == 'convexe':
        curv_tot = -curv_tot

    return curv_tot


def vector_rotation(vector,axis, theta):
    """
    Returns the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.

    :Parameters:
    - **Vector**: tuple-like, the vector to be rotated
    - **axis**: tuple-like, the rotation axis vector
    - **theta**: float, the rotation angle [rad]
    """
    #cf. http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    axis = scipy.asarray(axis)
    theta = scipy.asarray(theta)
    axis = axis/scipy.sqrt(scipy.dot(axis, axis))
    a = scipy.cos(theta/2.0)
    b, c, d = -axis*scipy.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    rotation_matrix = scipy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    return scipy.dot(rotation_matrix, vector)


def VineAxeIIinsert(inI_vector, insert_angle=46.,insert_angle_CI=4.6,rot_range=180.):
    """
    Returns the azimut and inclination [rad] of the insertion unit vector of the first secondary internode

    :Parameters:
    - **inI_vector**: tuple-like, the primary node vector
    - **insert_angle**: float, the average secondary-to-primary internodes insertion angle [degrees]
    - **insert_angle_CI**: float, the confidence interval associated to insert_angle [degrees]
    - **rot_range**: float, the range [degrees] within which the axis is rotated
    """

    len_I, azi_I, incli_I = [inI_vector[i] for i in [0,1,2]]

    a_insert = insert_angle + insert_angle_CI*(min(1,max(-1,scipy.randn()/2.96)))
    a_insert = a_insert*scipy.pi/180.
    rot_range = rot_range*scipy.pi/180.
    rand_rotate = rot_range*scipy.rand()
#    azi_init = azi_I + a_insert*scipy.sin(rand_rotate)
#    incli_init = incli_I + a_insert*scipy.cos(rand_rotate)
    dx, dy, dz = PolToXyz((len_I, azi_I, incli_I))
    rot_axis_1 = scipy.cross((dx, dy, dz),(0.,0.,1.))
    vec_2 = vector_rotation((dx, dy, dz), rot_axis_1, a_insert)
    vec_2 = vector_rotation(vec_2, (dx, dy, dz), rand_rotate)
    len_init, azi_init, incli_init = cart_to_pol(vec_2)

    return azi_init, incli_init


def VineAxeII(g, vid, phyllo_angle=180., PT_init=0.5, insert_angle=46.,
              insert_angle_CI=4.6, pruning_type='avg_field_model', N_init=0.18,
              N_max=2.25, N_max_order=4, in_order_max=25, slope_nfii=5.7,
              phyto_type='P0', a_L=43.718, b_L=-37.663, a_P=1.722, b_P=10.136,
              c_P=-5.435, Fifty_cent=400., slope_curv=70.,curv_type='convexe'):

    """
    Adds secondary phytomers to an existing MTG.

    :Parameters:
    - **g**: an MTG object
    - **vid**: integer, vertex ID
    - **phyllo_angle**: float, the phyllotaxis angle [degrees]
    - **PT_init**: float, the maximum curvature point of an axis (Louarn et al., 2008)
    - **insert_angle**: float, the average secondary-to-primary internodes insertion angle [degrees]
    - **insert_angle_CI**: float, the confidence interval associated to insert_angle [degrees]
    - **pruning_type**: string, the pruning type, one of the following: 'avg_field_model','GDC_1', 'Lyre', 'Rideau_simple', 'VSP_HL', 'Lyre_ouverte', 'GDC_2', 'Pot'
    - **N_init**: float, the average number of secondary internodes which are connected to the first primary internode
    - **N_max**: float, the average maximum number of secondary internodes
    - **N_max_order**: integer, the order of the primary internode which has N_max
    - **in_order_max**: integer, the order of the primary internode where NFIIthe number of secondary internodes is assumed equal to zero
    - **slope_nfii**: float, the slope (b) of the sigmoidal relationshipe NFII=f(in_order), a modified formula of that proposed by Louarn (2005, PhD, Eq. II.8)
    - **phyto_type**: string, the type of the primary phytomere. Can be one of 'P0', 'P1' or 'P2' (Louarn et al., 2007)
    - **a_L**, **b_L**: float, the slope and intercept, respectively, of the linear relationship LenII=f(NFII)
    - **a_P**, **b_P**, **c_P**: floats, the coefficients of the polynomial relationship LenII=f(NFII)
    - **Fifty_cent**: float, the abscice of inflection point of the sigmoidal curve
    - **slope_curv**: float, the slope of the sigmoidal relationship curv=f(length)
    - **curv_type**: in plant exterior direction even 'convex' or 'concave'
    """

    phyllo_angle = phyllo_angle*scipy.pi/180.

    if not 0. <= PT_init <= 1.0 : PT_init = max(0., min(PT_init,1.))
    #insert_angle=insert_angle*pi/180.
    #insert_angle_CI=insert_angle_CI*pi/180.

# Generation of axis elements
    if vid>0:
        fatherI = g.node(vid)
        if fatherI.label.startswith('inI'):
            order = g.Class(vid).split('in')[1]
            if order == 'I':
                in_order = int(fatherI.index()) if not fatherI.label[-1].isalpha() else int(findall('\d+',str(fatherI.label))[-1])  # In case where the label ends with an alphabetical letter ('M' for Multiple internodes)
                NFII = VineNFII(in_order, pruning_type,N_init,N_max,N_max_order,in_order_max,slope_nfii,phyto_type)
                tot_len = 0.1*VineLII(NFII, pruning_type, a_L, b_L, a_P, b_P, c_P)
                length = VineInL(NFII, tot_len)

#               Generation of secondary internodes
                for vtx_id in range(NFII):
                    if vtx_id == 0:
                        en_c, axeII = g.add_child_and_complex(vid, label='inII1',edge_type='+')
                        g.node(axeII).label = 'shII'+str(fatherI.index())
                        g.node(axeII).edge_type = '+'
                    else:
                        en_c = g.add_child(en_c, label=('inII'+str(vtx_id+1)),edge_type='<')

#                   Determination of the secondary internodes TopPosition
                    if vtx_id == 0:
                        grandpa = fatherI.parent()
                        counter=0
                        for vtx in grandpa.children():
                            if vtx.label == 'inII1':
                                counter += 1

#                               Getting vector of the previous axisII insertion internode
                                sibling = vtx
                                BotPos_sg = grandpa.TopPosition
                                TopPos_sg = sibling.TopPosition
                                dx_sg, dy_sg, dz_sg = [round(float(i),2) for i in scipy.subtract(TopPos_sg,BotPos_sg)]
                                len_sg,azi_sg,incli_sg = cart_to_pol((dx_sg,dy_sg,dz_sg))

#                               Getting the vector of the primary internode holding the previous axisII
                                TopPos_gp = grandpa.TopPosition
                                BotPos_gp = grandpa.parent().TopPosition
                                dx_gp, dy_gp, dz_gp = [round(float(i),2) for i in scipy.subtract(TopPos_gp,BotPos_gp)]
                                len_gp,azi_gp,incli_gp = cart_to_pol((dx_gp,dy_gp,dz_gp))

                                # Getting the vector of the primary internode holding the actual axisII
                                TopPos_f = fatherI.TopPosition
                                BotPos_f = grandpa.TopPosition
                                dx_f, dy_f, dz_f = [round(float(i),2) for i in scipy.subtract(TopPos_f,BotPos_f)]
                                len_f,azi_f,incli_f = cart_to_pol((dx_f,dy_f,dz_f))

#                               First correction azi and incli follow the dirction of the holding primary internode.
                                dazi = azi_f - azi_gp
                                dincli = incli_f - incli_gp

                                azi_tempo1 = azi_sg + dazi
                                incli_tempo1 = incli_sg + dincli

                                if abs(azi_tempo1) >= abs(azi_f):
                                    phyllo_angle = -phyllo_angle*scipy.sign(incli_f)
                                else:
                                    phyllo_angle = phyllo_angle*scipy.sign(incli_f)

                                vector = PolToXyz((length,azi_tempo1,incli_tempo1))
                                axis_I = (dx_f, dy_f, dz_f)
                                phyllotaxis = phyllo_angle*(1.+0.1*min(1,max(-1,scipy.randn()/2.96))) # the confidence interval 0.1 is to be confirmed.

#                               Second correction: the axisII insertion internode must adhere to phyllotaxis rule.
                                dx, dy, dz = vector_rotation(vector,axis_I,phyllotaxis)
                                #len_init,azi_init,incli_init = cart_to_pol((dx,dy,dz))

#                               Third correction: Considering the uncertainty in the insertion angle
                                ins_min = insert_angle - insert_angle_CI
                                ins_max = insert_angle + insert_angle_CI
                                a_insert = scipy.arccos(min(max(-1,scipy.dot((dx, dy, dz),(axis_I))),1))
                                sup_range = max(0.,ins_max - a_insert)
                                inf_range = -max(0., a_insert - ins_min)
                                rot_range = scipy.rand()*(sup_range - inf_range) + inf_range
                                normal_vec = scipy.cross((dx, dy, dz),(axis_I)) # the normal vector to the plane defined by the primary internode and the insertion secondary internode
                                dx, dy, dz = vector_rotation((dx, dy, dz),normal_vec,rot_range)
                                len_init,azi_init,incli_init = cart_to_pol((dx,dy,dz))

                        if counter ==0:
                            TopPos_f = fatherI.TopPosition
                            BotPos_f = grandpa.TopPosition
                            dx_f, dy_f, dz_f = [round(float(i),2) for i in scipy.subtract(TopPos_f,BotPos_f)]
                            len_f,azi_f,incli_f = cart_to_pol((dx_f,dy_f,dz_f))

                            father_vec = (len_f,azi_f,incli_f)
                            azi_init, incli_init = VineAxeIIinsert(father_vec)#,insert_angle,insert_angle_CI)
                            dx, dy, dz = PolToXyz((length,azi_init, incli_init))

                    else:
                        vidII = g.node(en_c)
                        fatherII = vidII.parent()
                        grandpaII = fatherII.parent()
                        dx_prec,dy_prec,dz_prec = [round(float(i),2) for i in scipy.subtract(fatherII.TopPosition,grandpaII.TopPosition)]
                        vector_prec = (dx_prec,dy_prec,dz_prec)
                        len_prec,azi_prec,incli_prec = cart_to_pol(vector_prec)

                        curv_tot = VineAxisCurv(incli_init, tot_len*10., Fifty_cent, slope_curv,curv_type)

                        PT = PT_init #+ ((scipy.pi - curv_tot)/scipy.pi)*(1-PT_init)

                        n1 = int(round(PT*NFII,0))
                        n2 = NFII - n1

                        dcurv = curv_tot/2./n1/2. if vtx_id+1 <= n1 else curv_tot/2./n2/2.
                        incli = incli_prec + dcurv if vtx_id == 0 else incli_prec + 2*dcurv

                        TopPos_f = fatherII.TopPosition

#                       In order to avoid that all vertical ramifications goes in the direction of x axis.
                        if incli_prec == scipy.pi/2.: azi_prec = scipy.rand()*scipy.pi*2.

                        dx = length*scipy.cos(incli)*scipy.cos(azi_prec)
                        dy = length*scipy.cos(incli)*scipy.sin(azi_prec)
                        dz = length*scipy.sin(incli)

                    x_bot, y_bot,z_bot = [TopPos_f[i] for i in [0,1,2]]
                    x_n, y_n, z_n = round(x_bot + dx, 2), round(y_bot + dy,2), round(z_bot + dz,2)

                    g.node(en_c).TopPosition = [x_n, y_n, z_n]

    return g


def VinePetLen(in_order, len_max=8.,Fifty_cent=30.,sig_slope=4.2):
    """
    Returns petiole length (def in cm) of an internode or a pruning complex.

    :Parameters:
    - **in_order**: integer, internode order
    - **len_max**: float, the average maximum length of a petiole [cm]
    - **Fifty_cent**: float, the abscice of inflection point of the sigmoidal curve
    - **sig_slope**: float, the slope (b) of the sigmoidal curve Diam=f(range)
    """

    petiol_len = len_max*(1-1/(1+scipy.exp(-(in_order-Fifty_cent)/sig_slope)))
    # TODO: The confidence intervals may be added

    return petiol_len


def VinePetiole(g, vid, pet_ins=90., pet_ins_cv=10., phyllo_angle=180.,
                phyllo_angle_cv=10., len_max_I=8., len_max_II=4.,
                Fifty_cent=30., sig_slope=4.2):

    #VinePetiole(in_order, len_max=80.,Fifty_cent=30.,sig_slope=4.2):
    #VineAxeIIinsert(vec_f,insert_angle,insert_angle_CI))
    """
    Adds petiols to an existing MTG.

    :Parameters:
    - **g**: an MTG object
    - **vid**: integer, vertex ID
    - **pet_ins**: float, the average insertion angle [degree] between a petiole and its holding internode, in the absence of a secondary/tertiary axis
    - **pet_ins_cv**: float, the coefficient of variation to `pet_ins` (%)
    - **phyllo_angle**: float, the phyllotaxis angle of secondary axes
    - **phyllo_angle_cv**: float, the coefficient of variation to `phyllo_angle` (%)
    - **len_max_I**: float, the maximum length of primary petioles [cm]
    - **len_max_II**: float, the maximum length of secondary petioles [cm]
    - **Fifty_cent**: float, the abscice of inflection point of the sigmoidal curve
    - **sig_slope**: float, the slope (b) of the sigmoidal curve Diam=f(range)
    """

    pet_ins, phyllo_angle = [scipy.radians(x) for x in (pet_ins, phyllo_angle)]

    if vid > 0:
        n = g.node(vid)
        if n.label.startswith('inI'):
            father = n
            grandpa = n.parent()
            grandgrandpa = grandpa.parent() if grandpa.parent() != None else grandpa
            in_order = int(father.index()) if not father.label[-1].isalpha() else int(findall('\d+',str(father.label))[-1])  # In case where the label ends with an alphabetical letter ('M' for Multiple internodes)

            len_max = len_max_I if g.Class(vid).split('in')[1] == 'I' else len_max_II
            len_max = VineLeafLen(in_order, len_max, len_max/10.) # TODO: insert parameters
            len_petI = VinePetLen(in_order, len_max, Fifty_cent,sig_slope)

            vec_f = scipy.subtract(father.TopPosition,grandpa.TopPosition)
            len_f,azi_f,incli_f = cart_to_pol(vec_f)

            vec_gp = scipy.subtract(grandpa.TopPosition,grandgrandpa.TopPosition)
            len_gp,azi_gp,incli_gp = cart_to_pol(vec_gp)

            if in_order == 1:
                # Generating first petiole
                normal_vec = scipy.cross(vec_gp,vec_f)
                rotation_axis1 = scipy.cross(vec_f,normal_vec)
                arb = 1. #scipy.rand()
            else:
                # Determining the actual petiolar vector as a function of the previous petiolar vector (phyllotaxis)
                for sibling in grandpa.children():
                    if sibling.label.startswith('PetI'):
                        BotPos_sg = grandpa.TopPosition
                        TopPos_sg = sibling.TopPosition
                        vec_sg = scipy.subtract(TopPos_sg,BotPos_sg)
                        len_sg,azi_sg,incli_sg = cart_to_pol(vec_sg)
                rotation_axis1 = scipy.cross(vec_f, vec_sg)
                arb = 1.

#           Initiating the petiole
            pet_vec = vec_f/norm(vec_f)*len_petI

#           Insertion angle
            pet_ins = pet_ins*(1.+(pet_ins_cv/100.)*min(1,max(-1,scipy.randn()/2.96)))
            pet_vec = vector_rotation(pet_vec,rotation_axis1, pet_ins)

#           Phyllotaxis
            phyllotaxis = phyllo_angle*(1.+(phyllo_angle_cv/100.)*min(1,max(-1,scipy.randn()/2.96)))*arb
            dx_pet, dy_pet, dz_pet = vector_rotation(pet_vec, vec_f, phyllotaxis)

            # Setting the TopPosition coordinates of the primary petiole
            x_bot, y_bot, z_bot = [father.TopPosition[i] for i in [0,1,2]]
            x_pet, y_pet, z_pet = round(x_bot + dx_pet, 2), round(y_bot + dy_pet,2), round(z_bot + dz_pet,2)

            label_pet = 'Pet' + (g.Class(vid)).split('in')[1] + str(in_order)

            petiole = g.add_child(vid, label=label_pet, edge_type='+')
            g.node(petiole).TopPosition = [x_pet, y_pet, z_pet]

    return g


def VineLeafLen(in_order, lim_max=15., lim_min=5.1, order_lim_max=7, max_order=40):
    """
    Returns the length of a leaf based on the order of the holding internode.

    :Parameters:
    - **in_order**: integer, the order of the internode
    - **lim_max**, **lim_min**: float, respectively the maximum and minimum limbe lengths
    - **order_lim_max**: integer, the order at which occurs lim_max
    - **max_order**: integer, the order at which the limbe length is assumed to equal zero
    """

    if in_order <= order_lim_max:
        LimbeLen = lim_min + (in_order-1)*(lim_max-lim_min)/(order_lim_max-1)
    else:
        LimbeLen = lim_max - (in_order-order_lim_max)*lim_max/(max_order-order_lim_max)

    return LimbeLen


def VineLeaf(g, vid, leaf_inc=-45., leaf_inc_cv=10., rand_rot_angle=30.,
             lim_max=15., lim_min=5.1, order_lim_max=7, max_order=40,
             cordon_vector=None):
    """
    Determines both the length and 3D orientation of a leaf based on the order of its holding internode.

    :Parameters:
    - **g**: mtg object
    - **vid**: integer, mtg vertex ID
    - **leaf_inc**: float, leaf inclination [Degrees]
    - **leaf_inc_cv**: float, coefficient of variation to `leaf_inc` (%)
    - **rand_rot_angle**: float, maximum angle of random rotation around Z axis [Degrees]
    - **lim_max** and **lim_min**: floats, maximum and minimum alloawable leaf lengths [cm], resp.
    - **order_lim_max** and **max_order**: shape parameters
    """

    leaf_inc, rand_rot_angle = [scipy.radians(angle) for angle in leaf_inc, rand_rot_angle]

    if vid > 0:
        n = g.node(vid)
        if n.label.startswith('Pet'):
            in_order = int(n.index()) if not n.label[-1].isalpha() else int(findall('\d+',str(n.label))[-1])  # In case where the label ends with an alphabetical letter ('M' for Multiple internodes)
            order = (g.Class(vid)).split('Pet')[1]
            #if order != 'I': lim_min = lim_max*0.8 # TODO: to be written more properly

            leaf_label = 'L' + order + '_' + str(in_order)
            leaf = g.add_child(vid, label=leaf_label, edge_type='+')

            limbe_len = VineLeafLen(in_order, lim_max, lim_min, order_lim_max, max_order)
            limbe_vec = scipy.array([0.,0.,1.])*limbe_len

            if cordon_vector is None:
                vec_pet = scipy.subtract(n.TopPosition, n.parent().TopPosition)
                rotation_axis = scipy.cross([0.,0.,1.],vec_pet)
            else:
                rotation_axis = Leaf_RotationAxis(n.TopPosition,cordon_vector)
            rotation_angle = scipy.pi/2. - leaf_inc *(1.+(leaf_inc_cv/100.)*min(1,max(-1,scipy.randn()/2.96)))
            dx_limbe, dy_limbe, dz_limbe = vector_rotation(limbe_vec,rotation_axis,rotation_angle)

            rotation_angle = rand_rot_angle * min(1,max(-1,scipy.randn()/2.96))
            dx_limbe, dy_limbe, dz_limbe = vector_rotation((dx_limbe, dy_limbe, dz_limbe),(0.,0.,1.),rotation_angle)

            x_bot, y_bot, z_bot = [n.TopPosition[i] for i in [0,1,2]]
            g.node(leaf).TopPosition = [round(x_bot + dx_limbe,2), round(y_bot + dy_limbe,2), round(z_bot + dz_limbe,2)]

    return g


def VineLobesTips(pPet, lobes_tips):
    """
    Returns the cartesian coordinates of vine leaf tips.

    :Parameters:
    - **pPet**: cartesian coordinates of the ptiolar tip
    - **lobes_tips**: a list containing the cartesian coordinates of the five lobes tips `p1`, `p2`, `p3`, `p4` and `p5`, respectively in the clockwise direction of a vertically schematized grapevine leaf (petiole downwards)
    """

#    points= Point3Array([ pgl.Vector3(0,9.8,0),pgl.Vector3(3.7,6.3,0),pgl.Vector3(2.3,3.9,0),
#                          pgl.Vector3(5.4,6.2,0), pgl.Vector3(5.9,3.2,0),pgl.Vector3(4.6,0.1,0),
#                          pgl.Vector3(6.4,-1.6,0), pgl.Vector3(1,-3.7,0), pgl.Vector3(0,0,0),
#                          pgl.Vector3(-1,-3.7,0),pgl.Vector3(-6.4,-1.6,0),
#                          pgl.Vector3(-4.6,0.1,0), pgl.Vector3(-5.9,3.2,0),
#                          pgl.Vector3(-5.4,6.2,0), pgl.Vector3(-2.3,3.9,0),pgl.Vector3(-3.7,6.3,0)])

    #p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15 = [points[i] for i in range(16)]

    theta_1 = 1.06185394004232810 #arccos(scipy.dot(p9,p10)/norm(p9)/norm(p10))
    theta_2 = 0.26671436996865744 #arccos(scipy.dot(p11,p10)/norm(p11)/norm(p10))
    theta_3 = 0.74194726800591759 #arccos(scipy.dot(p12,p10)/norm(p12)/norm(p10))
    theta_4 = 0.18369778647701493 #arccos(scipy.dot(p14,p13)/norm(p14)/norm(p13))
    theta_5 = 0.18551156925222045 #arccos(scipy.dot(p15,p13)/norm(p15)/norm(p13))

    ratio_1 = 0.58098705311078391 #norm(p9)/norm(p10)
    ratio_2 = 0.69745466856698446 #norm(p11)/norm(p10)
    ratio_3 = 1.01742628725623140 #norm(p12)/norm(p10)
    ratio_4 = 0.55068542551062905 #norm(p14)/norm(p13)
    ratio_5 = 0.88861969954204889 #norm(p15)/norm(p13)

    p10, p13, p0, p3, p6 = [lobes_tips[i] for i in range(5)]

    p10_pPet = p10-pPet
    p13_pPet = p13-pPet
    p0_pPet = p0-pPet
    p3_pPet = p3-pPet
    p6_pPet = p6-pPet

    norm_Pet_10_13 = scipy.cross(p10_pPet,p13_pPet)
    norm_Pet_13_0 = scipy.cross(p13_pPet,p0_pPet)
    norm_pet_0_3 = scipy.cross(p0_pPet,p3_pPet)
    norm_pet_3_6 = scipy.cross(p3_pPet,p6_pPet)

    p9  = vector_rotation(p10_pPet * ratio_1, norm_Pet_10_13, -theta_1) + pPet
    p11 = vector_rotation(p10_pPet * ratio_2, norm_Pet_10_13, theta_2) + pPet
    p12 = vector_rotation(p10_pPet * ratio_3, norm_Pet_10_13, theta_3) + pPet
    p14 = vector_rotation(p13_pPet * ratio_4, norm_Pet_13_0, theta_4) + pPet
    p15 = vector_rotation(p13_pPet * ratio_5, norm_Pet_13_0, theta_5) + pPet
    p1  = vector_rotation(p3_pPet * ratio_5, norm_pet_0_3, -theta_5) + pPet
    p2  = vector_rotation(p3_pPet * ratio_4, norm_pet_0_3, -theta_4) + pPet
    p4  = vector_rotation(p6_pPet * ratio_3, norm_pet_3_6, -theta_3) + pPet
    p5  = vector_rotation(p6_pPet * ratio_2, norm_pet_3_6, -theta_2) + pPet
    p7  = vector_rotation(p6_pPet * ratio_1, norm_pet_3_6, theta_1) + pPet

    points = Point3Array([p0, p1, p2, p3, p4, p5, p6, p7, pPet, p9, p10, p11, p12, p13, p14, p15])

    return points


def MTGbase(g, vtx_label='inT'):
    """
    Returns the basal vertex of a given vertex user-defined type.

    :Parameters:
    - **g**: an MTG object of **only one plant**
    - **vtx_label**: string, the label prefix of the basal highest-scale stem vertex
    """

    for vid in g.VtxList(Scale=3):
        n = g.node(vid)
        try:
            if n.label.startswith((vtx_label, 'rhyzo')) and n.parent() == None:
                vid_base = vid
                break
        except:
            pass

    return vid_base


def Cordon_vector(g):
    """
    Returns the vector of a **single vine** cordon axis.
    """

    ls = []
    for vid in g.VtxList(Scale=3):
        n = g.node(vid)
        if n.label.startswith('inT') and n.complex().label != 'trunk':
            ls.append(vid)

    pos_data = {vid:g.property('TopPosition')[vid] for vid in ls}
    data = scipy.array(pos_data.values())
    datamean = data.mean(axis=0)
#   Do a singular value decomposition on the mean-centered data.
    uu, dd, vv = scipy.linalg.svd(data - datamean)

    data_dist = distance.cdist(data,data,'euclidean')
    extr_pair0 = scipy.where(data_dist == data_dist.max())
    extr_pair = extr_pair0[0] if len(extr_pair0[0])==2 else scipy.concatenate((extr_pair0[0],extr_pair0[1]), axis=0)

    extr_pos = scipy.array([data[index_] for index_ in extr_pair])
    linepts = vv[0] * extr_pos
    # shift by the mean to get the line in the right place

    return linepts, scipy.array([linepts[0][i]-linepts[1][i] for i in range(2)] + [0.])


def Leaf_RotationAxis(petiole_tip,rotation_axis):
    """
    Returns the signed rotation axis of petiole around a given axis.

    :Parameters:
    - **petiole_tip**: array-like cartesian coordinates of the upper end of a petiole
    - **rotation_axis**: rotation axis vector in a cartesian system
    """

    v1 = scipy.array(petiole_tip[:2])
    v2 = scipy.array(rotation_axis[:2])

    return scipy.sign(det((v1,v2))) * rotation_axis


def AddSoil(g, side_length=10.):
    """
    Adds a soil element to an existing MTG.

    **Needs improvement!** for soil descritization.
    """

    xls, yls =[[g.node(vid).TopPosition[i] for vid in \
    g.property('geometry').keys()] for i in (0,1)]
    x_min, x_max = min(xls), max(xls)
    y_min, y_max = min(yls), max(yls)
    nbX = int((x_max - x_min)/side_length) + 1
    nbY = int((y_max - y_min)/side_length) + 1
    x = scipy.linspace(x_min, x_min+nbX*side_length,num=nbX)
    y = scipy.linspace(y_min, y_min+nbY*side_length,num=nbY)
    xu, yu = scipy.meshgrid(x,y)
    xu, yu = [reduce(lambda ix, iy: list(ix) +list(iy), ls) for ls in (xu,yu)]
    pos = zip(xu, yu)
    Soil = g.add_component(g.root, label='Soil0', edge_type='/')
    Soil = g.add_component(g.root, label='Soil', edge_type='/')

    pos=zip(xu,yu)
    for i, ipos in enumerate(pos):
        isoil = g.add_component(Soil, label=('other'+str(i)), TopPosition=list(ipos) + [0])
        g.node(isoil).geometry = transformation(soil0(side_length),
                 1., 1., 1., 0., 0.,0.,-side_length/2.,-side_length/2.,0.)

    return

def add_soil_components(g, cylinders_number, cylinders_radii, soil_dimensions,
                        soil_class, vtx_label):
    """
    Adds concentric soil cylinders to the mtg.
    
    :Parameters:
    - **g**: a multiscale tree graphe object
    - **cylinders_number**: integer, the number of cylinders to be added
    - **cylinders_radii**: list of floats, the radii of cylinders to be used
    - **soil_dimensions**: list of floats (length, width, depth) of the soil [m]
    - **soil_class**: string, the soil class name according to Carsel and Parrish (1988) WRR 24,755–769 DOI: 10.1029/WR024i005p00755
    - **vtx_label**: string, the label prefix of the basal highest-scale stem vertex
    """
    max_radius = 0.5 * min(soil_dimensions[:2])*100 #[cm]
    assert (max(cylinders_radii) <= max_radius), 'Maximum soil radius must not exceed %d cm'%max_radius
    assert (len(cylinders_radii) == cylinders_number), 'Soil cylinders number (%d) and radii elements (%d) do not match.'%(len(cylinders_radii),cylinders_number)

    depth = soil_dimensions[2]*100. #[m]
    child = g.node(MTGbase(g,vtx_label=vtx_label))
    Length = 0.
    radius_prev = 0.

    for ivid in range(cylinders_number):
        radius = cylinders_radii[ivid]
        Length = radius - radius_prev
        label = 'rhyzo%d'%ivid
#        print radius, Length, label
        child = g.node(child._vid).insert_parent(label=label, depth=depth, Length=Length,
                       TopDiameter = radius*2., BotDiameter = radius*2.,
                       TopPosition = 3*[0], BotPosition = 3*[0],
                       soil_class = soil_class)
        radius_prev = radius

    return child._vid


#==============================================================================
#==============================================================================
# MTG geometry construction
#==============================================================================

def _distance(coord1,coord2):
    """
    Calculates the distance between two points.
    coord1 and coord2 are array-like cartesian coordinates [x,y,z].

#   TODO: replace by *scipy.spatial.distance*
    """

    try:
#        deltaZ=coord2[2]-coord1[2]
        distance = scipy.sqrt((coord2[0]-coord1[0])**2+(coord2[1]-coord1[1])**2+(coord2[2]-coord1[2])**2)
    except ValueError:
        distance = scipy.sqrt((coord2[0]-coord1[0])**2+(coord2[1]-coord1[1])**2)

    return round(distance,2)


def slim_cylinder(length, radius_base, radius_top):
    """
    Tries to construct a cylinder with a low number of triangles (hack from Adel!!)"
    """

    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2*scipy.pi/3., 4*scipy.pi/3.
    r = rb
    p1 = (r*scipy.cos(a1), r*scipy.sin(a1),0)
    p2 = (r*scipy.cos(a2), r*scipy.sin(a2),0)
    p3 = (r*scipy.cos(a3), r*scipy.sin(a3),0)
    r = rt
    q1 = (r*scipy.cos(a1+scipy.pi), r*scipy.sin(a1+scipy.pi),length)
    q2 = (r*scipy.cos(a2+scipy.pi), r*scipy.sin(a2+scipy.pi),length)
    q3 = (r*scipy.cos(a3+scipy.pi), r*scipy.sin(a3+scipy.pi),length)
    set = pgl.TriangleSet([p1, p2, p3, q1, q2, q3],
                      [(2,1,0), (3,4,5), (0,5,4), (0,4,2), (2,4,3), (3,1,2), (1,3,5), (5,0,1)])

    return set


def StemElement_mesh(length, diameter_base, diameter_top, classic = True):
    """
    Computes mesh for a stem element (from Adel).

    :Parameters:
    - length: length of the stem
    - diameter_base: lower-end diameter
    - diameter_top: upper-end diameter
    - classic: logical, if True (default) the mesh covers the entire stem volume, else (False) is uses the economic `slim_cylinder` mesh
    """

    if classic:
        solid = True
        slices = 6  # 6 is the minimal number of slices for a correct computation of star (percentage error lower than 5)
        stem = pgl.Tapered(diameter_base/2., diameter_top/2., pgl.Cylinder(1., length , solid, slices))
        tessel = pgl.Tesselator()
        stem.apply(tessel)
        mesh = tessel.triangulation
    else:
        mesh = slim_cylinder(length, diameter_base /2., diameter_top /2.)

    return mesh


def leaf_obs(points):
    """
    Builds a grapevine leaf mesh, based on :func:`leaf0` function.

    :Parameters:
    **points**: a set of 15 :func:`Point3Array` points indicating leaf blade limits

    TODO: add this function to module `primitive`
    """

    indices= pgl.Index3Array([ pgl.Index3(0,1,8), pgl.Index3(15,0,8), pgl.Index3(8,11,14),
                           pgl.Index3(11,12,14), pgl.Index3(12,13,14),pgl.Index3(10,11,8),
                           pgl.Index3(8,9,10), pgl.Index3(2,3,4), pgl.Index3(2,4,5),
                           pgl.Index3(2,5,8), pgl.Index3(8,5,6), pgl.Index3(8,6,7)])

    f= pgl.TriangleSet(points, indices)

    return f


def VineDiam(g, vid, D_trunk=5.06, D_arm=3.77, D_Cx=2.91, D_3y=1.75,
             D_spur=1.15, D_cane=0.99, Fifty_cent=5.,sig_slope=10.,D_pet=0.35):
    """
    Returns the diameter [cm] of an internode or a structural element.

    :Parameters:
    - **g**: mtg object
    - **vid**: integer, mtg vertex ID
    - **D_trunk**, **D_arm** and **D_Cx**: float, the average diameters of the trunk, arms and pruning complexes, resp.
    - **D_3y**: float, the dimater of the first internode of 3-year old canes
    - **D_spur**, **D_cane**: float, the dimater of the first internodes of canes and spurs (supposed D_max), resp.
    - **D_cane**: float, the dimater of the first internode of a cane (supposed D_max)
    - **Fifty_cent**: float, the abscice of inflection point of the sigmoidal curve
    - **sig_slope**: float, the slope (b) of the sigmoidal curve Diam=f(range)
    """

    n = g.node(vid)

    if n.label.startswith(('sh','G')):
#       Setting the axisI initial diameter randomly
        if n.label.count('I') == 1:
            init_diam = D_cane*(1+0.1*(min(1,max(-1,scipy.randn()/2.96))))

#       Setting the axisII initial diameter empirically as a function of the diameter of the holding primary internode
        elif n.label.count('I') > 1:
            init_diam = n.components()[0].parent().TopDiameter
        init_diam = init_diam*((scipy.rand()*3+1)/4.)
        n.InitDiam = init_diam
        Diam = None

    if n.label.startswith('inT'):
        if 'I' and 'V' and 'y' not in n.label:
            Diam = D_trunk if n.complex().label.startswith('trunk') else D_arm
        elif '3y' in n.label:
            Diam = D_3y
    if n.label.startswith('cx'):
        Diam = D_Cx
    if n.label.startswith(('inT3y','inT2y')):
        in_order = int(n.index()) if not n.label[-1].isalpha() else int(findall('\d+',str(n.label))[-1]) # In case where the label ends with an alphabetical letter ('M' for Multiple internodes)
        Diam = D_spur #*(1-1/(1+scipy.exp(-(in_order-Fifty_cent)/sig_slope)))
    if n.label.startswith('inI'):
        init_diam = g.node(g.Complex(vid)).InitDiam
        in_order = int(n.index()) if not n.label[-1].isalpha() else int(findall('\d+',str(n.label))[-1])  # In case where the label ends with an alphabetical letter ('M' for Multiple internodes)
        Diam = init_diam #*(1-1/(1+scipy.exp(-(in_order-Fifty_cent)/sig_slope)))
    if n.label.startswith('Pet'):
        Diam = D_pet #*max(0,(scipy.rand()+1)/2)

    return Diam


def VineMTGProp(g, vid):
    """
    Attaches geometric properties to MTG vertices.

    :Parameters:
    - **g**: mtg object
    - **vid**: integer, mtg vertex ID
    """

    n = g.node(vid)
    if vid>0:
        if n.label.startswith(('sh','G')): VineDiam(g,vid) # Sets the initial diameters of the primary and secondary axes.

#       Setting the properties of internodes, pruning complices and petioles
        if n.label.startswith(('in','cx','Pet')):
            TopPosition = n.TopPosition
            TopDiameter = n.TopDiameter if n.TopDiameter != None else max(0.05,VineDiam(g,vid))
            p = n.parent()
            try:
                p == None # First phytomere at the basis of the trunk
                if hasattr(g.node(g.Trunk(vid,Scale=1)[0]), 'baseXYZ'):
                    BotPosition = g.node(g.Trunk(vid, Scale=1)[0]).baseXYZ
                else:
                    BotPosition = list(TopPosition[:2])
                    BotPosition.append(0.)

                BotDiameter = TopDiameter
            except AttributeError:
                BotPosition = p.properties()['TopPosition']
                BotDiameter = p.properties()['TopDiameter'] if n.label.startswith(('inT','cx')) else TopDiameter
            Length = _distance(TopPosition,BotPosition)
            g.node(vid).BotPosition = BotPosition
            g.node(vid).Length = Length
            if not n.TopDiameter: g.node(vid).TopDiameter = TopDiameter
            g.node(vid).BotDiameter = BotDiameter

#       Setting the properties of leaves
        if n.label.startswith('LI'):
            TopPosition = n.TopPosition
            BotPosition = n.parent().TopPosition
            Length = _distance(TopPosition,BotPosition)
            g.node(vid).BotPosition = BotPosition
            g.node(vid).Length = Length

    return g


def VineMTGGeom(g, vid):
    """
    Adds geometry to elements of an MTG object.

    :Parameters:
    - **g**: mtg object
    - **vid**: integer, mtg vertex ID
    """
#    theta_1, theta_2 = [scipy.radians(x) for x in (theta_1, theta_2)]

    try:
        n = g.node(vid)

        if n.label.startswith(('inT','cx','inT3y','inT2y','inI','Pet')):
            length = n.properties()['Length']
            diameter_base = n.properties()['BotDiameter']
            diameter_top = n.properties()['TopDiameter']
            #mesh = pgl.Cylinder(diameter_base,length,True,6)
            mesh = StemElement_mesh(length, diameter_base, diameter_top, classic = True)
            g.node(vid).geometry = mesh

        elif n.label.startswith('LI'):
            if hasattr(n, 'TopPositionPoints'): # If leaf tips are digitized
                pPet = n.parent().TopPosition
                points = VineLobesTips(pPet, n.TopPositionPoints)
                mesh = leaf_obs(points)
            else:
                leaf_mesh = transformation(leaf0(1), 1., 1., 1., -scipy.pi/2.,0.,0.,0.,0.,0.)

                length = n.properties()['Length']
                leaf_vec = scipy.subtract(n.TopPosition, n.parent().TopPosition)
                leaf_len, leaf_azi, leaf_incli = cart_to_pol(leaf_vec)
#                theta_1 = -scipy.pi/2.
#                theta_2 = 0. #scipy.pi*(1.+(theta_2_cv/100.)*min(1,max(-1,scipy.randn()/2.96)))
                mesh = transformation(leaf_mesh, length, length, 1., 0., 0.,0.,0.,0.,0.)
                #mesh = leaf0(length)
            g.node(vid).geometry = mesh

    except AttributeError:
        pass

    return g


def VineTransform(g, vid):
    """
    Transforms all elements of an MTG to their real position.
    """

    n=g.node(vid)

    if vid>0:
#       Internodes, pruning complices and petioles
        if n.label.startswith(('inT','cx','inT3y','inT2y','inI','inII','Pet')):
            Length = n.properties()['Length']
            x,y,z = [round(float(i),2) for i in scipy.subtract(n.properties()['TopPosition'],n.properties()['BotPosition'])]

            if Length==0:
                incli =0
            else:
                incli = float(scipy.arcsin(z/Length))
            if (x==0 and y==0):
                azi = 0
            elif (y>=0) :
                azi = scipy.arccos(min(max(-1,x/scipy.sqrt(x**2+y**2)),1))
            else :
                azi = -scipy.arccos(min(max(-1,x/scipy.sqrt(x**2+y**2)),1))

            mesh = n.geometry
            mesh = pgl.EulerRotated (azi, 3.14/2.-incli, 0, mesh)
            mesh = pgl.Translated (pgl.Vector3(n.properties()['BotPosition']), mesh)
            g.node(vid).geometry = mesh

#       Leaves
        elif n.label.startswith('LI'):
            if not hasattr(n, 'TopPositionPoints'):
                x_bot,y_bot,z_bot = [n.BotPosition[i] for i in [0,1,2]]
                mesh = n.geometry
                vec_lim = scipy.subtract(n.TopPosition, n.BotPosition)
                len_lim, azi_lim, incli_lim = cart_to_pol(vec_lim)
                mesh = transformation(mesh, len_lim, len_lim, 1., azi_lim, -1*incli_lim,0.,x_bot,y_bot,z_bot)
                g.node(vid).geometry = mesh

    return g


def VineOrient(g, vid, theta, v_axis=[0.,0.,1.], local_rotation=False):
    """
    Rotates an MTG around an axis by a given angle.

    :Parameters:
    - **g**: mtg object
    - **vid**: integer, mtg vertex ID
    - **theta**: float, rotation angle [degrees]
    - **v_axis**: the rotation axis, default [0.,0.,1.]
    - **local_rotation**: logical, if False (default), the MTG is rotated about the reference [0., 0., 0.], else if True, the MTG is rotated around its own base (**works only for one plant**)
    """

    if float(theta) != 0.:
        v_axis = scipy.array(v_axis)
        n = g.node(vid)
        try:
            TopPosition = n.TopPosition
            if local_rotation:
                nBase = g.node(g.Trunk(vid, Scale=1)[0])
                TopPosition = scipy.array(TopPosition) - scipy.array(nBase.baseXYZ)
                TopPosition = vector_rotation(TopPosition,v_axis,scipy.radians(theta))
                TopPosition = TopPosition + scipy.array(nBase.baseXYZ)
            else:
                TopPosition = vector_rotation(TopPosition,v_axis,scipy.radians(theta))
            n.TopPosition = TopPosition
        except:
            pass

        if 'baseXYZ' in n.properties():
            if local_rotation == False:
                baseXYZ = n.baseXYZ
                baseXYZ = vector_rotation(baseXYZ,v_axis,scipy.radians(theta))
                n.baseXYZ = baseXYZ

        if 'TopPositionPoints' in n.properties():
            lobe_tip = []
            TPPoints = scipy.array(n.TopPositionPoints)
            if local_rotation == True:
                nBase = g.node(g.Trunk(vid, Scale=1)[0])
                TPPoints = TPPoints - scipy.array(nBase.baseXYZ)
                for tipPoint in TPPoints:
                    lobe_tip.append(vector_rotation(tipPoint,v_axis,scipy.radians(theta)))
                lobe_tip = scipy.array(lobe_tip) + scipy.array(nBase.baseXYZ)
            else:
                for i in range(len(n.TopPositionPoints)):
                    lobe_tip.append(vector_rotation(n.TopPositionPoints[i],v_axis,scipy.radians(theta)))
            n.TopPositionPoints = scipy.array(lobe_tip)

    return g


#==============================================================================
#==============================================================================
# Write output
#==============================================================================

def mtg_output(g, wd):
    """
    Writes an MTG object (**g**) to an external file given its **output_path**.
    """
    properties = [(p, 'REAL') for p in g.property_names() if p not in ['edge_type', 'index', 'label']]
    mtg_lines = io.write_mtg(g, properties)
    if not path.exists(wd+'mtg/'):
        mkdir(wd+'mtg/')
    mtg_file_path = wd + 'mtg/%s.mtg'%g.date
    f = open(mtg_file_path, 'w')
    f.write(mtg_lines)
    f.close()
    return

def mtg_save(g, scene, file_path):

    if not path.exists(file_path):
        mkdir(file_path)

    geom = {vid:g.node(vid).geometry for vid in g.property('geometry')}
    g.remove_property('geometry')

    fg = file_path + 'mtg' + g.date + '.pckl'

    f = open(fg, 'w')
    dump([g,scene], f)
    f.close()
    #restore geometry
    g.add_property('geometry')
    g.property('geometry').update(geom)
    return
    
def mtg_load(wd, index):
    
    fgeom = wd + 'geometry.bgeom'
    fg = wd + 'mtg%s.pckl'%(index)
    
    scene = pgl.Scene()
    scene.read(fgeom, 'BGEOM')
    geom = {sh.id:sh.geometry for sh in scene}

    f = open(fg)
    g2, TT = load(f)
    f.close()
    
    g2.add_property('geometry')
    g2.property('geometry').update(geom)
    
    return g2, scene

def mtg_save_geometry(scene, file_path):
    """
    Saves the geometry of a scene in an external file.

    :Parameters:
    - **scene**: scene object
    - **file_path**: path string for saving the scene object
    """

    if not path.exists(file_path):
        mkdir(file_path)

    fgeom = file_path + 'geometry.bgeom'

    scene.save(fgeom, 'BGEOM')

    return