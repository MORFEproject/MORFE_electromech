mutable struct ST6
  mat::Int64
  nodes::Vector{Int64}
end

mutable struct SB3
  nodes::Vector{Int64}
  T6::Vector{Int64}
  psidof::Int64
  psi::Float64
end

mutable struct Snode
 coor::Vector{Float64}
 dof::Vector{Int64}
 u::Vector{Float64}
end

mutable struct Sinfo
 mesh_file::String
 NN::Int64
 NE::Int64
 NL::Int64
 uneq::Int64
 uaneq::Int64
 psineq::Int64
 nK::Int64
 nA::Int64
 nMat::Int64
 nmm::Int64
 Lmm::Vector{Int64}
 nza::Int64
 nzna::Int64
 nrom:: Int64
 Ffreq::Int64
 Fmult::Float64
 alpha::Float64
 beta::Float64
 tol::Float64
 neig::Int64
 style::Char
 max_order::Int64
 max_orderNA::Int64
 Sinfo() = new()
end

