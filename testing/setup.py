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
import os
from itertools import product

from hacd.cut_generators.cut_generators_enum import CutGenerator
from hacd.util.geometry import PolytopeDescription


cut_generators = {
    'sweep': CutGenerator.SWEEP,
    'facet': CutGenerator.FACET
}


RUN_CONFIG = {
    'instance_name': None,
    'cut_generator': None,
    'maxVolError': 0.01,
    'max_depth': 3,
    'nr_cuts': 7
}


INSTANCES = {
    'test2D': {
        'path': os.path.abspath('testing/test_data/test2D.json'),
        'description': PolytopeDescription.INNER_DESCRIPTION
    }
}


def pytest_addoption(parser):
    parser.addoption("--instance", help="run tests for instance", default="all")
    parser.addoption("--cut_generator", help="run tests with cut Generator", default="sweep")


def pytest_generate_tests(metafunc):
    configurations = []
    if 'run_configuration' in metafunc.fixturenames:
        instance_arg = metafunc.config.getoption('instance')
        test_instances = INSTANCES.keys() if instance_arg == 'all' else [instance_arg]
        cut_arg = metafunc.config.getoption('cut_generator')
        if cut_arg == 'all':
            test_cut_generator = cut_generators.values()
        else:
            test_cut_generator = [cut_generators[cut_arg]]

        for (instance, cut_generator) in product(test_instances, test_cut_generator):

            configurations.append(populate_configuration(
                RUN_CONFIG,
                {'name': instance,
                 'path': INSTANCES[instance]['path'],
                 'description': INSTANCES[instance]['description'],
                 'cut_generator': cut_generator})
            )
        metafunc.parametrize("run_configuration", configurations)


def populate_configuration(configuration, conf_dict):
    from copy import deepcopy
    populated_dict = deepcopy(configuration)
    for key, item in conf_dict.items():
        populated_dict[key] = item
    return populated_dict