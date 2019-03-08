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
from enum import Enum


class CutGenerator(Enum):
    SWEEP = 'sweep'
    FACET = 'facet'

    def __eq__(self, other):
        if self.value == other.value:
            return True
        else:
            return False