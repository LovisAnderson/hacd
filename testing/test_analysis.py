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
import pytest
import sys
print(sys.path)
from sweepvolume.cell_decomposition import Cell_Decomposition
from hacd.analysis import conv_hull_cell_decomposition, union_only_events, max_diff_delta, random_normed_directions
from sweepvolume.sweep import Sweep
from sweepvolume.geometry import Vertex, Hyperplane
import numpy as np


def test_convex_hull_cd_generation(triangles_meeting_in_one_point):
    cd_conv = conv_hull_cell_decomposition(triangles_meeting_in_one_point)
    coordinates = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
    assert all([any([np.isclose([float(e_i) for e_i in e.vertex.coordinates], coordinate, atol=1e-03).all() for
                     coordinate in coordinates]) for
                e in cd_conv.events])
    assert round(Sweep(cd_conv.events).calculate_volume(), 4) == 1


def test_union_only_events(translated_triangles):
    union_sweeps = [Sweep(translated_triangles.events, a) for
                    a in random_normed_directions(50, 2)]
    events = union_only_events(union_sweeps)
    assert len(events) == 1
    assert events.pop() == Vertex.vertex_from_coordinates(np.array([0.,0.]))


def test_max_diff_delta(translated_triangles):
    cd_conv = conv_hull_cell_decomposition(translated_triangles)
    union_sweep = Sweep(translated_triangles.events, sweep_plane=np.array([1, 0]))
    conv_sweep = Sweep(cd_conv.events, sweep_plane=np.array([1, 0]))
    lam, max, active = max_diff_delta(union_sweep, conv_sweep)
    assert set(active) == {0,2,3}
    assert lam == 0


def test_non_regular_convex_hull(cube_simplex_overlapping_3d_2):
    cd_conv = conv_hull_cell_decomposition(cube_simplex_overlapping_3d_2)


@pytest.fixture
def triangles_meeting_in_one_point():
    # geometry: |><|
    h0 = Hyperplane(np.array([1, 0]), 0)
    h1 = Hyperplane(np.array([1, 0]), -1)
    h2 = Hyperplane(np.array([1, -1]), 0)
    h3 = Hyperplane(np.array([1, 1]), -1)
    hyperplanes = [h0, h1, h2, h3]
    poly1 = {(0, 1), (2, -1), (3, -1)}
    poly2 = {(1, -1), (2, 1), (3, 1)}
    return Cell_Decomposition(hyperplanes, [poly1, poly2])


@pytest.fixture
def translated_triangles():
    # geometry: <\<\ two congruent triangles that meet in 0
    h0 = Hyperplane(np.array([0, 1]), 0)
    h1 = Hyperplane(np.array([-1, 1]), -1)
    h2 = Hyperplane(np.array([1, 1]), 0)
    h3 = Hyperplane(np.array([-1, 1]), 0)
    h4 = Hyperplane(np.array([1, 1]), -1)

    hyperplanes = [h0, h1, h2, h3, h4]
    poly1 = {(0, 1), (1, -1), (2, -1)}
    poly2 = {(0, 1), (3, -1), (4, -1)}
    return Cell_Decomposition(hyperplanes, [poly1, poly2])


@pytest.fixture
def cube_simplex_overlapping_3d_2():
    c_0 = Hyperplane(np.array([1, 0, 0]), 0)
    c_1 = Hyperplane(np.array([0, 1, 0]), 0)
    c_2 = Hyperplane(np.array([0, 0, 1]), 0)
    c_3 = Hyperplane(np.array([1, 0, 0]), -1)
    c_4 = Hyperplane(np.array([0, 1, 0]), -1)
    c_5 = Hyperplane(np.array([0, 0, 1]), -1)

    s_0 = Hyperplane(np.array([1, 1, 1]), -4)
    s_1 = Hyperplane(np.array([1, 0, 0]), -0.5)
    s_2 = Hyperplane(np.array([0, 1, 0]), -0.5)
    s_3 = Hyperplane(np.array([0, 0, 1]), -0.5)

    hyperplanes = [c_0, c_1, c_2, c_3, c_4, c_5, s_0, s_1, s_2, s_3]

    simplex = {(6, -1), (7, 1), (8, 1), (9, 1)}
    unit_cube = {(0, 1), (1, 1), (2, 1), (3, -1), (4, -1), (5, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube, simplex])