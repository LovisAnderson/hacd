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
from hacd.cut_generators.cut_generators_enum import CutGenerator
from hacd.util.geometry import PolytopeDescription
from hacd.util.data_reader import get_cell_decompositions
from hacd.util import argparse_helpers
import argparse
import logging
import os
import json

from acd_tree import build_acd, render_tree_dict


def _arguments():
    """
    Helper method to define all possible arguments.

    :return: Dictionary of all possible arguments.
    """

    # Helper class that holds an argument that will be added to the parser.
    class ArgHolder:
        def __init__(self, arg, **kwargs):
            self.arg = arg
            self.kwargs = kwargs

    result = {
        "polytopePath": ArgHolder(
            "--polytopePath",
            default="test_data/test2D.json",
            type=argparse_helpers.valid_file,
            help="File Path of json that describes polytope"
        ),
        "polytopeDescription":  ArgHolder(
            "--polytopeDescription",
            action=argparse_helpers.enum_action(PolytopeDescription),
            default="innerDescription",
            help="Description of polytope"
        ),
        "outputDir": ArgHolder(
            "--outputDir",
            default="./results",
            help="Directory for output data files"
        ),
        "cutGenerator": ArgHolder(
            "--cutGenerator",
            default='sweep',
            action=argparse_helpers.enum_action(CutGenerator),
            help="Method for generation of cuts"
        ),
        "maxVolError": ArgHolder(
            "--maxVolError",
            default=0.05,
            type=float,
            help="Maximal relative (to total) volume error tolerated"
                 " in any ACD Node"
        ),
        "reduceHyperplanes": ArgHolder(
            "--reduceHyperplanes",
            default=False,
            type=bool,
            help="Indicates if hyperplane reduction"
                 " should be used as preprocessing"),
        "maxDepth": ArgHolder(
            "--maxDepth",
            default=10,
            type=int,
            help="Maximum Depth in the ACD tree"
        ),
        "nrCuts": ArgHolder(
            "--nrCuts",
            default=10,
            type=int,
            help="Nr of cuts that are tried in each step."
        ),
    }

    # Define all possible arguments

    return result


if __name__ == "__main__":
    """
    Main function.
    """


    def convert_arg_line_to_args(arg_line):
        """
        Function to allow for argument input from file
        of format: --[argName] [value].
        """
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg


    # Create argument parser object and override read-method.
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.convert_arg_line_to_args = convert_arg_line_to_args

    arguments = _arguments()
    for k, v in arguments.iteritems():
        parser.add_argument(v.arg, **v.kwargs)
        parser.set_defaults()
    # Parse args and init parameter dictionary.
    args = parser.parse_args()
    instance_file = os.path.basename(args.polytopePath)
    instance = os.path.splitext(instance_file)[0]
    output_dir = os.path.abspath(args.outputDir + '/{}_{}_mD{}_nrC{}'.format(
        instance,
        args.cutGenerator.value,
        args.maxDepth,
        args.nrCuts,
    ))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logpath = os.path.join(output_dir, 'log.txt'.format(instance))
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    hdlr = logging.FileHandler(os.path.abspath(logpath), mode='w')
    logger.addHandler(hdlr)

    params = {k: v for (k, v) in vars(args).iteritems()}
    logger.info("Using the following parameters:\n%s" % str(params))
    union_cd, convex_cd = get_cell_decompositions(args.polytopePath,
                                                  description=args.polytopeDescription,
                                                  reduce_hyperplanes=args.reduceHyperplanes)
    tree = build_acd(union_cd,
                     convex_cd,
                     max_vol_error=args.maxVolError,
                     max_depth=args.maxDepth,
                     cut_generator=args.cutGenerator,
                     nr_cuts=args.nrCuts)

    with open(os.path.join(output_dir, 'tree.json'), 'w') as fout:
        json.dump(tree.as_dict(), fout, indent=3)
    render_tree_dict(os.path.join(output_dir, 'tree.json'))
