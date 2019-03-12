# *****************************************************************************
#       Copyright (C) 2017      Tom Walther
#                     2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
import logging
import json

import hacd.analysis as ana

from sweepvolume.cell_decomposition import Cell_Decomposition
from sweepvolume.geometry import Polytope, Vertex, Hyperplane

from .geometry import PolytopeDescription

import numpy as np

from ordered_set import OrderedSet

from hacd.util.geometry import norm_halfspace, bounding_box


def load_json(filepath):
    # Open and read specified input file.
    #  Input format has to be {disj: {polyID1: [...], ..., polyIDk: [...]}}

    with open(filepath, 'r') as f:
        disj = json.load(f)
        logging.debug('Loaded disjunction from file: {}'.format(filepath))
    if len(disj) > 1:
        logging.warning("JSON input file contains multiple disjunctions -"
                        " only first will be considered!")
    return disj


def polytopes_from_json(filepath, description=PolytopeDescription.INNER_DESCRIPTION):
    """
    Method to read polytope points from JSON file.
     Input format is {disj: {polyID1: [...], ..., polyIDk: [...]}}
     If polytopes are given in outer description the elements polyID: [[a1, ..., ad, b],...,[...]]
     correspond to the halfspace a1x1 + ... + adxd + b <= 0
    """
    polys = []
    disj = load_json(filepath)
    for disjID, _ply in disj.items():
        logging.info("Reading polytopes of disjunction {}"
                     " which are given in {}".format(disjID,
                                                     description.name))

        if description == PolytopeDescription.INNER_DESCRIPTION:

            for i, (plyID, ply) in enumerate(_ply.items()):
                poly_vertices = set([Vertex.vertex_from_coordinates(np.array(coordinates))
                                     for coordinates in ply])
                polys.append(Polytope(vertices=poly_vertices))

        elif description == PolytopeDescription.OUTER_DESCRIPTION:

            for i, (plyID, ply) in enumerate(_ply.items()):
                halfspaces = [(Hyperplane(np.array(hyp_vec[:-1]), hyp_vec[-1]), -1)
                              for hyp_vec in ply]
                polys.append(Polytope(halfspaces=halfspaces))

    return polys


# Method computes cell decompositions for union of polytopes and the convex hull
def get_cell_decompositions(filepath,
                            reduce_hyperplanes=True,
                            description=PolytopeDescription.INNER_DESCRIPTION):
    polytopes = polytopes_from_json(filepath, description=description)
    hyperplanes, polytope_vectors = hyperplanes_and_polytope_vectors(polytopes)
    if reduce_hyperplanes:
        hyperplanes, polytope_vectors = ana.drop_facets(hyperplanes, polytope_vectors)
    logging.info("Computing cell decomposition for union of polytopes"
                 " from {} hyperplanes.".format(len(hyperplanes)))
    bbox = bounding_box(polytopes)
    union_cd = Cell_Decomposition(hyperplanes,
                                  polytope_vectors,
                                  bounding_box=bbox)
    conv_cd = ana.conv_hull_cell_decomposition(union_cd)
    return [union_cd, conv_cd]


def hyperplanes_and_polytope_vectors(polytopes):
    """
    Method computes hyperplanes and position vectors for given polytopes
    :param polytopes:a list of list of tuples: (sweepvolume.geometry.Hyperplane, orientation)
    :return: hyperplanes: list of hyperplanes,
     polytopes as list of tuples (hyperplane_index, orientation)
    """
    assert all([isinstance(poly, Polytope) for poly in polytopes])
    hyperplanes = OrderedSet()
    polytope_vectors = []
    for poly in polytopes:
        ply = set()
        for halfspace in poly.halfspaces:
            halfspace = norm_halfspace(halfspace)
            idx = hyperplanes.add(halfspace[0])
            ply.add((idx, halfspace[1]))
        polytope_vectors.append(ply)
    return hyperplanes, polytope_vectors
