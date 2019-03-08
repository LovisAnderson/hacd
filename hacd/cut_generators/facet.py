# *****************************************************************************
#       Copyright (C) 2017      Tom Walther
#                     2017-2018 Lovis Anderson  <lovisanderson@posteo.net>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sweepvolume.sweep import Sweep
from sweepvolume.geometry import Hyperplane
from hacd.analysis import lam_close_to_border

import numpy as np

import logging


def facet_cuts(ACDNode, nr_cuts=10):
    """
    Method returns best cuts whereby cuts are drawn from the defining Hyperplanes.
    They are scored by the derivation of the convexification error
    :param ACDNode: ACDnode which carries cell decompositions for
    :param nr_cuts: The best nr_of_cuts are returned.
    :return: best nr_cuts as list of hyperplane objects
    """

    cuts = []
    for ind, hyperplane in enumerate(ACDNode.union_cd.hyperplanes):
        union_sweep, conv_sweep = get_sweeps(ACDNode, hyperplane)

        diff_delta = diff_delta_at_facet(union_sweep, conv_sweep, ind)
        logging.debug('cut: <{}> ; difference delta: {}'.format(str(hyperplane), diff_delta))
        cuts.append((Hyperplane(hyperplane.a, hyperplane.b), diff_delta))
    cuts = sorted(cuts, key=lambda x: x[1], reverse=True)
    return list(zip(*cuts[:nr_cuts])[0])


def get_sweeps(ACDNode, hyperplane):
    """
    Method computes Sweep Objects. It pertubates the given hyperplane since Sweep only works then.
    :param ACDNode: ACDNode object
    :param hyperplane: hyperplane object as in geometry.hyperplane
    :return: union_sweep, conv_sweep Sweep objects for the pertubated hyperplane
    """
    cut_plane = hyperplane.pertubate()
    sweep_direction = cut_plane.a / np.linalg.norm(cut_plane.a)
    union_sweep = Sweep(ACDNode.union_cd.events, sweep_direction)
    conv_sweep = Sweep(ACDNode.convex_cd.events, sweep_direction)
    return union_sweep, conv_sweep


def diff_delta_at_facet(union_sweep, conv_sweep, hyperplane_index):
    """
    Method computes the quasi derivation of the sweep convexification error at the facet.
    (quasi derivation since we do not divide by epsilon).
    It returns 0 if facet is only incident to one polytope.
    If this is the case it is an outer facet (ir inner we have 2 or more components).

    :param union_sweep: sweep object for union of polytopes
    :param conv_sweep: sweep object for convex hull of union of polytopes
    :param hyperplane_index: index of hyperplane in cd object
    :return: quasi derivation of the sweep convexification error at
     hyperplane with hyperplane_index in cell decomposition object or
     0 if all events that are on hyperplane are only incident to a single polytope
    """

    epsilon = 1e-6
    scale = union_sweep.sorted_vertices[-1][1] - union_sweep.sorted_vertices[0][1]
    logging.debug('sweep direction: {}'.format(union_sweep.sweep_plane))
    logging.debug('scale: {},'
                  ' 1st lam: {:.2f},'
                  ' last lam: {:.2f}'.format(scale,
                                             union_sweep.sorted_vertices[0][1],
                                             union_sweep.sorted_vertices[-1][1]))
    for event, lam in union_sweep.sorted_vertices:
        if hyperplane_index not in event.incidences:
            continue
        logging.debug('lambda: {}'.format(lam))
        if (lam_close_to_border(lam,
                                union_sweep.sorted_vertices[0][1],
                                union_sweep.sorted_vertices[-1][1])):
            return 0
        if len(event.incident_polytopes) > 1:
            union_volumes = union_sweep.calculate_volumes((lam - epsilon,
                                                           lam,
                                                           lam + epsilon))
            conv_volumes = conv_sweep.calculate_volumes((lam - epsilon,
                                                         lam,
                                                         lam + epsilon))
            # taking the derivation of the difference into account for both directions
            derivation_backwards = (conv_volumes[1] - union_volumes[1]) - \
                                   (conv_volumes[0] - union_volumes[0])
            derivation_forward = (conv_volumes[2] - union_volumes[2]) - \
                                 (conv_volumes[1] - union_volumes[1])
            return max(derivation_backwards, derivation_forward)
    return 0
