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
import json
from pathlib import Path

from hacd.acd_tree import build_acd
from hacd.util.data_reader import get_cell_decompositions
from hacd.util.json_encoder import NumpyEncoder
import random


def test_log_files(run_configuration):
    random.seed(1)


    # outputdir ../../test_data/output/instance_name
    filepath = Path(__file__)
    twolevelsup = Path(filepath.parents[0])
    outdir = Path(twolevelsup, 'test_data', 'output', run_configuration['name'])
    outdir.mkdir(parents=True, exist_ok=True)
    # Configure logger
    print('creating logs for run configuration {}'.format(run_configuration))
    logger = logging.getLogger()
    logpath = Path(outdir, run_configuration['cut_generator'].name.lower() + '.log')
    hdlr = logging.FileHandler(str(logpath), mode='w')
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    union_cd, convex_cd = get_cell_decompositions(
        run_configuration['path'],
        run_configuration['description']
    )
    tree = build_acd(
        union_cd,
        convex_cd,
        max_vol_error=run_configuration['maxVolError'],
        max_depth=run_configuration['max_depth'],
        cut_generator=run_configuration['cut_generator'],
        nr_cuts=run_configuration['nr_cuts']
    )

    map(logger.removeHandler, logger.handlers[:])
    map(logger.removeFilter, logger.filters[:])

    with open(str(Path(outdir, run_configuration['cut_generator'].name.lower() + '.json')), 'w') as fout:
        json.dump(tree.as_dict(), fout, indent=3, cls=NumpyEncoder)
    map(logger.removeHandler, logger.handlers[:])
    map(logger.removeFilter, logger.filters[:])