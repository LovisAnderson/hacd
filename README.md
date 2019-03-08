# Hierarchical Approximative Convex Decomposition
This package is still under development and far from perfect.
 Nonetheless it is somewhat stable and able to compute an hierarchical
 approximate convex decomposition for a union of polytopes. A technical description of
 our approach is in the making.
## Installation
You need to have [sweepvolume](https://gitlab.com/LovisAnderson/sweepvolume) package installed.
After installation sweepvolume and its dependencies you can install hacd through
``pip install /path/to/hacd/``

## Usage
You can call hacd through 
``
python acd.py @/path/to/some/configuration.cfg
``
An example configuration can be found under testing/test_data/test2D.cfg. In the output directory
you can then find 3 files. With help of the **log.txt** file you can get insight into the process
and the decisions that have been made. In **tree.json** the whole hierarchy is saved.
 **tree.png** is a graph visualization of the quality of the tree.

## License
sweepvolume is distributed under the terms of the GNU General Public License (GPL)
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version. See http://www.gnu.org/licenses/.
 
