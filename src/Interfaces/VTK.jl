"""
Module to export ONSAS data structures to `.vtk` files suitable [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes

@reexport import WriteVTK: vtk_grid

# Add code here.

end
