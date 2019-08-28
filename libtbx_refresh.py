from __future__ import absolute_import, division, print_function

import libtbx.pkg_utils

libtbx.pkg_utils.define_entry_points(
    {
        "dials.index.basis_vector_search": [
            "real_space_grid_smart_search = exafel_project.ADSE13_25.indexing.iota_strategies:RealSpaceGridSmartSearch",
        ],
    }
)
