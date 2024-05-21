#
# sets input
#

input_file="demo2.jl"
#input_file="demo3.jl"

using SparseArrays
using LinearAlgebra
using Arpack
using ExtendableSparse
#using BenchmarkTools
using Pardiso
using MAT
using MATLAB
#using Plots

include("./source/EM2D_defs.jl")
include("./source/EM2D_globalprocedure.jl")
#include("./source/EM2D_globalprocedure2.jl")
#include("./source/EM2D_globalprocedure3.jl")
#include("./source/EM2D_globalprocedure4.jl")
#include("./source/EM2D_globalprocedure8.jl")
include("./source/EM2D_readgmsh.jl")
include("./source/EM2D_outgmsh.jl")
include("./source/EM2D_elemental.jl")
include("./source/EM2D_analysis.jl")
include("./source/EM2D_eig.jl")
include("./source/EM2D_param_struct.jl")
include("./source/EM2D_dpim.jl")
include("./source/EM2D_dpim_FEM.jl")
#include("./source/EM2D_prove.jl")
include("./source/EM2D_dpim_realification.jl")
include("./source/EM2D_dpim_output.jl")
include("./source/quadrature.jl")
include("./source/shape_functions.jl")
include("./source/EM2D_matcont.jl")

globalprocedure()
#globalprocedure2()
#globalprocedure3()
#globalprocedure4()
#globalprocedure8()

println("done!")




