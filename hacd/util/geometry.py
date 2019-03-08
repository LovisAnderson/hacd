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
import numpy as np
from enum import Enum


class PolytopeDescription(Enum):
    INNER_DESCRIPTION = 'points'
    OUTER_DESCRIPTION = 'hyperplanes'


def norm_halfspace(halfspace):
    """
    Method norms halfspace in such a way that
    :param halfspace:
    :return:
    """
    hyperplane, orientation = halfspace[0], halfspace[1]
    if hyperplane.b < 0:
        hyperplane.a = -hyperplane.a
        hyperplane.b = -hyperplane.b
        hyperplane.matrix = hyperplane.as_matrix()
        return hyperplane, -orientation
    else:
        return halfspace


def bounding_box(polytopes):
    """
    Method calculates bounding box for multiple polytopes.
    :param polytopes: list of polytope objects.
    :return: (lower_bounds, upper_bounds) tuple with lower (upper) bound
     for each coordinate.
    """
    all_points = np.vstack((poly.get_vertex_coordinates() for poly in polytopes))
    # we need a slightly bigger bounding box bc dropping facets
    #  and pertubation can create vertices outside of the box
    # to do : smarter tube around bbox (taking size into consideration)
    return all_points.min(axis=0) - 10, all_points.max(axis=0) + 10