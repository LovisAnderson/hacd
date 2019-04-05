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
import copy
import json

from sweepvolume.sweep import Sweep
from analysis import detect_clusters

from hacd.cut_generators.cut_generators_enum import CutGenerator
from hacd.cut_generators.facet import facet_cuts
from hacd.cut_generators.sweep import sweep_cuts

from analysis import conv_hull_cell_decomposition

import random

random.seed(1)


class Node(object):
    """
    Class for a node of the Approximate Convex Decomposition tree.
    """

    def __init__(self,
                 union_cd,
                 convex_cd,
                 id=None,
                 depth=0,
                 parent_id=None,
                 tol_rel=0.01,
                 tol_abs=None,
                 max_depth=100,
                 cut_generator=None):
        """
        Node constructor.
        :param union_cd: CellDecomposition object representing a union of polytopes
        :param convex_cd:CellDecomposition object representing the convex hull of union_cd
        :param id: A string for the id
        :param depth: The depth in the ACD tree.
        :param parent_id: The id of the parent ACD node.
        :param tol_rel: Relative error tolerance.
        :param tol_abs: Absolute error tolerance.
        :param max_depth: Maximum depth of the ACD tree.
        :param cut_generator: A CutGenerator (enum) object
        """

        self.dim = union_cd.dim
        assert isinstance(id, str) or id is None
        assert isinstance(depth, int) and depth >= 0
        assert isinstance(tol_rel, float) and 0 <= tol_rel <= 1
        assert isinstance(tol_abs, float) and tol_abs >= 0 or tol_abs is None
        assert isinstance(max_depth, int) and max_depth > 0

        # Build random id hash and store parent id.
        self.id = id if id is not None else random.getrandbits(32)
        self.parent_id = parent_id
        self.depth = depth
        self.union_cd = union_cd
        self.convex_cd = convex_cd

        # Store cutGenerator
        self.cut_generator = cut_generator

        self.convex_hull_volume = self._convex_hull_volume()

        # Compute volume of union of polytopes.
        self.volume = self._union_volume()

        # Init list of child ACD nodes.
        self.children = []

        # Store tolerances.
        self.tol_rel = tol_rel
        self.tol_abs = tol_abs if tol_abs else self.volume * self.tol_rel
        self.max_depth = max_depth

        logging.info("Initialized %s" % self)

    def _convex_hull_volume(self):
        return Sweep(self.convex_cd.events).calculate_volume()

    def _union_volume(self):
        return Sweep(self.union_cd.events).calculate_volume()

    def check_abort(self):
        """
        Method to check abort criteria for ACD run.

        Returns:
            True, if some abort criterion is satisfied. False, otherwise.
        """

        if self.volume_error_small_enough():
            return True

        # Check maximum depth.
        if self.depth >= self.max_depth:
            logging.info("Maximum depth reached --> no further decomposition!")
            return True

        return False

    def volume_error_small_enough(self):
        eps = 1e-5

        # Check volume tolerance threshold.
        vol_error_abs = self.convex_hull_volume - self.volume
        if self.relative_error() <= self.tol_rel + eps or vol_error_abs <= self.tol_abs + eps:
            logging.info("Volume error is small enough --> no further decomposition!")
            return True
        return False

    def relative_error(self):
        return self.convex_hull_volume / self.volume - 1

    def find_clusters(self):
        clusters = detect_clusters(Sweep(self.union_cd.possible_events))
        logging.info("Clusters found : %s" % clusters)
        return clusters

    def cluster_to_cell_decomposition(self, cluster):
        polytope_vectors = [self.union_cd.polytope_vectors[i] for i in cluster]
        union_cd = copy.deepcopy(self.union_cd)
        for v in union_cd.possible_events:
            v.update_position_vector(union_cd.hyperplanes)
        union_cd.polytope_vectors = polytope_vectors
        union_cd.events = union_cd.find_events()
        conv_cd = conv_hull_cell_decomposition(union_cd)

        return union_cd, conv_cd

    def clusters_to_nodes(self, clusters):
        children = []
        for poly_indices in clusters:
            logging.info("processing cluster: {!s}".format(poly_indices))
            cds = self.cluster_to_cell_decomposition(poly_indices)
            # It can happen that a not full dimensional cluster is found.
            # In this case the corresponding cell decomposition should have no events.
            # We do not add that cluster.
            if len(cds[0].events) == 0:
                continue
            child_node = self.__class__(cds[0],
                                        cds[1],
                                        parent_id=self.id,
                                        depth=self.depth + 1,
                                        tol_rel=self.tol_rel,
                                        tol_abs=self.tol_abs,
                                        max_depth=self.max_depth,
                                        cut_generator=self.cut_generator)
            children.append(child_node)

        return children

    def find_cuts(self, nr_cuts):
        if self.cut_generator == CutGenerator.FACET:
            return facet_cuts(self, nr_cuts)
        elif self.cut_generator == CutGenerator.SWEEP:
            return sweep_cuts(self, nr_cuts)
        else:
            return NotImplementedError

    def as_dict(self):
        node_dict = {
            'depth': self.depth,
            'volume': self.volume,
            'cell_decomposition': self.union_cd.as_dict(),
            'convex_volume': self.convex_hull_volume,
            'parent_id': str(self.parent_id),
            'children': [str(child.id) for child in self.children],
            'total_error': self.convex_hull_volume - self.volume,
            'relative_error': self.relative_error()
        }
        return node_dict

    def to_json(self, path):
        with open(path, 'w') as outfile:
            json.dump(self.as_dict(), outfile, sort_keys=True, indent=4, separators=(',', ': '))
        return

    def best_cut(self, cuts):

        logging.info("Searching best cut out of {} cuts".format(len(cuts)))

        min_score = self.convex_hull_volume
        for i, cut in enumerate(cuts):
            cut_children = self.apply_cut(cut)
            score = sum(acdNode.convex_hull_volume for acdNode in cut_children)

            logging.info("Trying cut : <%s>" % str(cut))
            logging.info("-- volume of current ACD node      : %1.2f" % self.convex_hull_volume)
            logging.info("-- total volume of child ACD nodes : %1.2f" % score)

            # test if cut is problematic and skip this cut in that case
            if self.problematic_cut(cut, cut_children, score):
                continue
            #  default cut is first cut which should've scored best in heuristic from cut generator
            if score < min_score or i == 0:
                min_score = score
                best_cut = cut
                children = cut_children
                if all(node.volume_error_small_enough() for node in children):
                    break

        # set node children
        self.children += children
        # Return best cut.
        logging.info("Best cut : <%s>" % str(best_cut))
        return best_cut

    def problematic_cut(self, cut, cut_children, score):
        import os
        import time

        rel_tol = 1.005
        if score > rel_tol * self.convex_hull_volume:
            problem = 'cut_increased_volume'
        elif rel_tol * sum([child.convex_hull_volume for child in cut_children]) < self.volume:
            problem = 'convex_volume_too_small'
        else:
            return False
        logging.warning('Cut: <{}> -> {}'.format(str(cut), problem))
        cd_dict = self.union_cd.as_dict()
        cd_dict['problematic_cut'] = cut.a
        outfile = '{}_{}.json'.format(time.strftime("%Y%m%d-%H%M%S"), problem)
        outfile_path = os.path.join(
            os.path.abspath('/OPTI/bzfander/ACD3d/python/hyperacd/results/problem_cases'),
            outfile)
        with open(outfile_path, 'w') as out_file:
            json.dump(self.as_dict(), out_file, indent=4)
        return True

    def apply_cut(self, cut):
        """
        Apply cut to current ACD node and create child nodes.
        :param cut: A Hyperplane object.
        :return: The child ACD nodes.
        """

        logging.info("Applying cut : <%s>" % str(cut))

        # Init list of child ACD nodes.
        childACDNodes = []

        # Create child ACD nodes.
        for i in [0, 1]:

            cds = self.restrict_cds(cut, -1 if i == 0 else 1)
            if not cds:
                logging.warning("union cd is empty! cut: {},"
                                " halfspace : {}".format(str(cut), -1 if i == 0 else 1))
                continue

            childACDNodes.append(
                Node(cds[0],
                     cds[1],
                     depth=self.depth + 1,
                     parent_id=self.id,
                     tol_rel=self.tol_rel,
                     tol_abs=self.tol_abs,
                     max_depth=self.max_depth,
                     cut_generator=self.cut_generator)
            )

        return childACDNodes

    def restrict_cds(self, cut, orientation):
        union_cd = copy.deepcopy(self.union_cd)
        union_cd.restrict_to_halfspace(cut, orientation)
        if len(union_cd.events) == 0:
            return None
        conv_cd = conv_hull_cell_decomposition(union_cd)
        return [union_cd, conv_cd]

    def logStatistics(self):
        """
        Method to log some statistics of the ACD.
        """

        logging.info("ACD NODE STATISTICS")
        logging.info("Depth in ACD tree         : {}".format(self.depth))
        logging.info("Concavity tolerance (rel) : {:.2f}".format(self.tol_rel))
        logging.info("Concavity tolerance (abs) : {:.2f}".format(self.tol_abs))
        logging.info("Number of polytopes       : {}".format(len(self.union_cd.polytope_vectors)))
        logging.info("Volume (polyunion)        : {:.2f}".format(self.volume))
        logging.info("Volume (convhull)         : {:.2f}".format(self.convex_hull_volume))

    def __str__(self):
        """
        String output method.
        """

        outputStr = ("ACD[%dD] node | id : %8s | parent_id : %8s | depth : %d" %
                     (self.dim, self.id, self.parent_id, self.depth))

        return outputStr
