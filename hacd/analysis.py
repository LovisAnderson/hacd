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
from sweepvolume.cell_decomposition import Cell_Decomposition
from sweepvolume.sweep import Sweep
from sweepvolume.geometry import Polytope, vector_distance

import numpy as np
import pandas as pd
import logging


def conv_hull_cell_decomposition(cell_decomposition, reduce_hyperplanes=True):
    """
    Method computes Cell Decomposition for convex hull of events of input cell decomposition
    :param cell_decomposition: CellDecomposition object
    :param reduce_hyperplanes: boolean if close hyperplanes should be removed from convex hull
    :return: cell decomposition object for convex hull
    """

    p = Polytope(vertices=set([e.vertex for e in cell_decomposition.events]))
    hyperplanes = p.hyperplanes
    pos_vec = [zip(range(len(hyperplanes)), [1] * len(hyperplanes))]
    if reduce_hyperplanes:
        hyperplanes, pos_vec = drop_facets(hyperplanes, pos_vec)
    hyperplanes = [h.pertubate() for h in hyperplanes]

    logging.info('create convex hull cell decomposition from {} hyperplanes'
                 .format(len(hyperplanes), len(cell_decomposition.events)))
    return Cell_Decomposition(hyperplanes,
                              pos_vec,
                              bounding_box=cell_decomposition.bbox)


def convex_hull_and_union_sweeps(union_cd, conv_hull_cd, sweep_planes):
    conv_hulls_sweeps = [Sweep(conv_hull_cd.events, sweep_plane=sweep_plane)
                         for sweep_plane in sweep_planes]
    union_sweeps = [Sweep(union_cd.events, sweep_plane=sweep_plane) for sweep_plane in sweep_planes]
    return union_sweeps, conv_hulls_sweeps


def max_diff_delta(union_sweep, conv_hull_sweep):
    """
    Method calculates the maximum of g' on the set non_conv_events.
    g := sweep_vol_conv - sweep_vol_union
    :param union_sweep: Sweep object of the union of polytopes
    :param conv_hull_sweep: Sweep object of the convex hull
    :param non_conv_events: list of events: union - convex hull
    :return: lam, diff_delta_max, active_hyperplanes : the lambda at which the maximum is attained,
                                                       the maximum,
                                                       indices of active hyperplanes at the event
                                                       corresponding to lambda
    """
    eval_times, eval_events = evaluation_points(union_sweep)
    # If there are no times for possible cuts we return 0 for max_diff_delta
    if len(eval_times) == 0:
        return None, 0, None
    t2 = union_sweep.calculate_volumes(eval_times)
    t3 = conv_hull_sweep.calculate_volumes(eval_times)
    diff_delta = diff_deltas(t3, t2)
    event, lam = eval_events[diff_delta.argmax()]
    return lam, diff_delta.max(), event.incidences


def diff_deltas(union_sweep_evaluated, conv_hull_evaluated):
    """
    Method computes quasi differential of the difference of union - conv_hull.
    Method computes differential as distance of
                [diff_1 - diff_0, diff_3 - diff_2, ..., diff_n - diff_n-1]
    :param union_sweep_evaluated: np.array of sweep evaluated shortly before and after events
    :param conv_hull_evaluated:  np.array of sweep evaluated shortly before and after events
    :return: np.array of quasi differential of the difference
    """
    diff = union_sweep_evaluated - conv_hull_evaluated
    i = 0
    diff_delta = []
    while i < len(diff):
        derivation_backwards = diff[i] - diff[i + 1]
        derivation_forwards = diff[i + 2] - diff[i + 1]
        diff_delta.append(max(derivation_backwards, derivation_forwards))
        i += 3
    return np.array(diff_delta)


def evaluation_points(union_sweep, eps=np.exp(-8)):
    """
    Method computes times lambda_i - eps, lambda_i, lambda_i + eps
     for events that are not on the boundary of the convex hull.
    [lambda_1 - eps, lambda_1, lambda_1 + eps, lambda_2 - eps, lambda_2, lambda_2 + eps, ...]
    :param union_sweep: sweep object of a union of polytopes
    :param non_conv_events: a list of events that are not part of the corresponding conv. hull sweep
    :param eps: step before/after the time at which the events occur
    :return: eval_times, eval_event ; the lambdas at which the sweep is to be evaluated
                                      and the corresponding events
    """
    eval_intervals = []
    eval_events = []

    for event, lam in union_sweep.sorted_events:
        if (
                len(event.incident_polytopes) > 1
                and not lam_close_to_border(lam,
                                            union_sweep.sorted_events[0][1],
                                            union_sweep.sorted_events[-1][1])
        ):
            eval_intervals.append([lam - eps, lam, lam + eps])
            eval_events.append((event, lam))
    eval_times = np.array(eval_intervals).flatten()
    return eval_times, eval_events


def lam_close_to_border(lam, start_lam, end_lam, tolerance=0.01):
    scale = end_lam - start_lam
    if (
            lam - start_lam < tolerance * scale or
            end_lam - lam < tolerance * scale
    ):
        return True
    else:
        return False


def union_only_events(union_sweeps):
    """
    Method computes a subset of the inner events of the union of polytopes.
    Inner events are events that are not events for the convex hull.
    All events that are met first or last by a sweep are outer events.
    We might not get all outer vertices with this procedure.
    :param union_sweeps: a list of Sweep Objects
    :return: A set of inner events
    """
    logging.debug('looking for events that only occur in the union of polytopes')
    conv_vertices = set()
    union_vertices = set([e[0] for e in union_sweeps[0].sorted_events])
    for sweep in union_sweeps:
        conv_vertices.add(sweep.sorted_events[0][0])
        conv_vertices.add(sweep.sorted_events[-1][0])
    logging.debug('union vertices: {},'.format(len(union_vertices)),
                  'conv vertices: {},'.format(len(conv_vertices)),
                  'inner_vertices:{}'.format(len(union_vertices - conv_vertices)))
    return union_vertices - conv_vertices


def detect_clusters(sweep):
    """
    Method finds clusters through the polytope incidence graph.
    :param sweep: a sweep object
    :return: a list of sets whereby each set is a list of ints,
     which describe indices of polytopes in that cluster
    """
    import networkx as nx
    g = sweep.graph()
    clusters = []
    for c in nx.connected_components(g.to_undirected()):
        cluster = set()
        for node in c:
            cluster.update(node.incident_polytopes)
        clusters.append(cluster)
    return clusters


def cut_data(union_sweeps, convex_hull_sweeps):
    """
    Method constructs data frame from non_convex and convex sweeps.
    Data frame contains information about sweep-plane
    and cut.
    :param union_sweeps: list of sweep objects
    :param convex_hull_sweeps: list of sweep objects
    :return: pandas data frame containing information about sweep and best cut
    """
    cuts = []
    for sweep_ind in range(len(union_sweeps)):
        scale = abs(union_sweeps[sweep_ind].sorted_events[-1][1] -
                    union_sweeps[sweep_ind].sorted_events[0][1])
        sweepplane = union_sweeps[sweep_ind].sweep_plane
        lam, diff, active_hyperplanes = max_diff_delta(union_sweeps[sweep_ind],
                                                       convex_hull_sweeps[sweep_ind])
        cuts.append((sweepplane, lam, diff, scale, active_hyperplanes))

    cut_dataframe = pd.DataFrame(data=cuts,
                                 columns=[
                                     'direction',
                                     'cut_lambda',
                                     'diff_delta',
                                     'scale',
                                     'active_hyperplanes'
                                 ])
    return cut_dataframe


def hyperplane_identification(hyperplanes, normal_tolerance=1e-8, offset_tolerance=1e-5):
    reduced_hyperplanes = list()
    hyperplane_identification_dict = {i: i for i in range(len(hyperplanes))}
    for i, hyperplane in enumerate(hyperplanes):
        close_hyperplane_exists = False
        for j, comparision_hyperplane in enumerate(reduced_hyperplanes):
            b_norm_diff = abs(hyperplane.b / hyperplane.a_norm -
                              comparision_hyperplane.b / comparision_hyperplane.a_norm)
            if (normal_tolerance > vector_distance(hyperplane.a,
                                                   comparision_hyperplane.a,
                                                   absolute_value=False) and
                    offset_tolerance > b_norm_diff):
                close_hyperplane_exists = True
                hyperplane_identification_dict[i] = j
                break
        if not close_hyperplane_exists:
            reduced_hyperplanes.append(hyperplane)
            hyperplane_identification_dict[i] = len(reduced_hyperplanes) - 1
    logging.debug("dropped {} out of {} hyperplanes".format(
        len(hyperplanes) - len(reduced_hyperplanes),
        len(hyperplanes))
    )
    return reduced_hyperplanes, hyperplane_identification_dict


def drop_facets(hyperplanes, polytope_vectors):
    reduced_hyperplanes, hyperplane_identification_dict = hyperplane_identification(hyperplanes)
    reduced_polytope_vectors = []
    for polytope in polytope_vectors:
        reduced_polytope_vectors.append(set([(hyperplane_identification_dict[i], orient)
                                             for (i, orient) in polytope]))
    return reduced_hyperplanes, reduced_polytope_vectors


def random_normed_directions(nr_directions, dim=2):
    np.random.seed(0)
    x = np.random.normal(size=(nr_directions, dim))
    x /= np.linalg.norm(x, axis=1)[:, np.newaxis]
    return x
