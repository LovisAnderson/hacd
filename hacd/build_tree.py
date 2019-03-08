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

from node import Node
from cut_generators.cut_generators_enum import CutGenerator


def acd_tree(union_cd,
             convex_cd,
             max_vol_error=0.05,
             max_depth=10,
             cut_generator=CutGenerator.SWEEP,
             nr_cuts=10):
    assert isinstance(max_vol_error, float) and max_vol_error >= 0.0
    assert isinstance(max_depth, int) and max_depth > 0

    # Get dimension.
    dim = union_cd.dimension

    node = Node(
        union_cd,
        convex_cd,
        id="root",
        tol_rel=0.01,
        max_depth=100,
        cut_generator=cut_generator
    )

    nodes_to_decompose = [node]
    tree = dict()

    while nodes_to_decompose:
        node_to_decompose = nodes_to_decompose.pop()
        node_to_decompose.logStatistics()
        if node_to_decompose.check_abort():
            tree[node_to_decompose.id] = node_to_decompose.as_dict()
            continue
        # Find and resolve clusters
        clusters = node_to_decompose.find_clusters()

        if len(clusters) > 1:
            node_to_decompose.children = node_to_decompose.clusters_to_nodes(clusters)
        else:
            cuts = node_to_decompose.find_cuts(nr_cuts=nr_cuts)
            best_cut = node_to_decompose.best_cut(cuts)
        nodes_to_decompose += node_to_decompose.children
        tree[node_to_decompose.id] = node_to_decompose.as_dict()
    return tree
