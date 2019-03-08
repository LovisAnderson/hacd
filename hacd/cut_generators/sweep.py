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
import logging
import pandas

import numpy as np

from sweepvolume.geometry import vector_distance

import hacd.analysis as ana
from sweepvolume.geometry import Hyperplane


def sweep_cuts(node,
               n=10,
               sweeps_per_orthant=100):
    """
    Method to generate cuts via Sweep-Plane.
    :param node: Node for which cuts are to be generated
    :param n: Number of cuts to be returned.
    :param sweeps_per_orthant: Number of sweeps that are tried
    :return: a list of the best n cuts
    """

    assert isinstance(n, int) and n > 0
    nr_sweeps = sweeps_per_orthant * 2 ** node.dim
    # Generate random sweep planes.
    logging.debug("-- generating {} random sweep planes...".format(nr_sweeps))
    sweepplanes = ana.random_normed_directions(nr_sweeps, node.dim)

    # Apply sweep plane algorithm to the union of polytopes and their convex hull.
    logging.debug("-- applying sweep plane algorithm...")
    sweeps_polys, sweeps_convhull = ana.convex_hull_and_union_sweeps(node.union_cd,
                                                                     node.convex_cd,
                                                                     sweepplanes)
    cut_data = ana.cut_data(sweeps_polys, sweeps_convhull)

    # Choose cuts that are distinct (enough -> see tolerance in method)
    logging.info("-- selecting {} best cuts...".format(n))

    # tolerance dependant on dimension and nr of cuts (n)
    close_vector_tolerance = 0.005 / np.sqrt(n) * 2**node.dim
    best_cuts = get_best_distinct_cuts(cut_data, nr_of_cuts=n, tolerance=close_vector_tolerance)
    cuts_with_hyperplanes = [
        [Hyperplane(sweep['direction'], -sweep['cut_lambda']), sweep['active_hyperplanes']]
        for _, sweep in best_cuts.iterrows()
    ]
    # Move cuts to closest hyperplanes if close enough, tolerance
    cuts = move_cuts(cuts_with_hyperplanes, node.union_cd, tolerance=close_vector_tolerance * 0.5)

    # Return n best cut suggestions (sorted according to score).
    return cuts


def get_best_distinct_cuts(cut_data, nr_of_cuts=10, tolerance=0.1):
    cut_data = cut_data.sort_values(by=['diff_delta'], ascending=False)
    best_cuts = pandas.DataFrame(columns=['direction',
                                          'cut_lambda',
                                          'diff_delta',
                                          'scale',
                                          'active_hyperplanes'])
    for id, row in cut_data.iterrows():
        if best_cuts.shape[0] == nr_of_cuts:
            break
        if all([vector_distance(cut['direction'], row['direction']) > tolerance for
                _, cut in best_cuts.iterrows()]):
            best_cuts = best_cuts.append(row)
    return best_cuts


def move_cuts(cuts, cd, tolerance=0.01):
    moved_cuts = []
    for cut in cuts:
        min_dist = tolerance + 1e-3
        min_hyp = None
        if not cut[1]:
            logging.warning('no hyperplanes active at cut: <[}>'.format(str(cut[0])))
            moved_cuts.append(cut[0])
            continue
        for hyperplane_ind in cut[1]:
            hyperplane = cd.hyperplanes[hyperplane_ind]
            dist = vector_distance(hyperplane.a, cut[0].a)
            [min_dist, min_hyp] = [dist, hyperplane_ind] if dist < min_dist else [min_dist, min_hyp]
        if min_dist < tolerance:
            logging.debug('hyperplane moved: {} -> {}, dist: {:0.4f}'.format(
                str(cut[0]),
                str(cd.hyperplanes[min_hyp]),
                min_dist)
            )
            moved_cuts.append(Hyperplane(cd.hyperplanes[min_hyp].a, cd.hyperplanes[min_hyp].b))
        else:
            moved_cuts.append(cut[0])
    return moved_cuts

