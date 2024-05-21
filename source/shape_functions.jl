
function dNda!(dNda::Array{Float64},gp::Array{Float64},::Val{:TRI6n})
  a1,a2,a3 = gp
  dNda[1,1] = 4*a1-1
  dNda[1,2] = 0.0
  dNda[2,1] = 0.0
  dNda[2,2] = 4*a2-1
  dNda[3,1] = -4*a3+1
  dNda[3,2] = -4*a3+1
  dNda[4,1] = 4*a2
  dNda[4,2] = 4*a1
  dNda[5,1] = -4*a2
  dNda[5,2] = 4*(a3-a2)
  dNda[6,1] = 4*(a3-a1)
  dNda[6,2] = -4*a1
  return nothing
end

function N!(NL::Array{Float64},gp::Array{Float64},::Val{:TRI6n})
  a1,a2,a3 = gp
  NL[:]=[a1*(2*a1-1), a2*(2*a2-1), a3*(2*a3-1), 4*a1*a2, 4*a2*a3, 4*a1*a3]  
  return nothing
end
