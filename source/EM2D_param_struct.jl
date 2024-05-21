using Combinatorics


"""
mutable struct Parametrisation
"""
mutable struct Parametrisation
  m::Int64
  mindep::Int64
  Av::Vector{Any}
  corresp::Vector{Int64}
  analysed::Vector{Int64}
  W::Matrix{ComplexF64}
  f::Matrix{ComplexF64}
  Wr::Matrix{ComplexF64}
  fr::Matrix{ComplexF64}
  R::Matrix{ComplexF64}
  Wf::Matrix{ComplexF64}
  Parametrisation() = new()
end

"""
fills parametrization structure
"""
function initParametrisation!(info::Sinfo)

  P = [Parametrisation() for i in 0:info.max_order]
  for p in 1:info.max_order
    P[p].Av,P[p].corresp,P[p].m,P[p].mindep=indexset(info,p)
    P[p].W = zeros(ComplexF64,info.nA,P[p].m)  
    P[p].Wr = zeros(ComplexF64,info.nA,P[p].m)  
    P[p].f = zeros(ComplexF64,info.nrom,P[p].m)
    P[p].Wf = zeros(ComplexF64,info.nA,P[p].m)  
    P[p].fr = zeros(ComplexF64,info.nrom,P[p].m)
    P[p].R = zeros(ComplexF64,info.nA,P[p].m)
    P[p].analysed = zeros(Int64,P[p].m)
  end
  return P
  end 


"""
creates multiexponents and looks for conjugates  
"""
function indexset(info::Sinfo,p::Int64)

nz=info.nza
nzf=info.nzna  
nzhalf=Int(nz/2)
nzfhalf=Int(nzf/2)
ndof=info.nrom

a=collect(multiexponents(ndof,p))
lena=length(a)
corresp=zeros(Int64,lena)
nc=0
#corresp=[i for i in 1:lena] 
for i in 1:lena
  orderna=sum(a[i][nz+1:ndof])
#  println(orderna)
  if orderna <= info.max_orderNA  
    b=[a[i][nzhalf+1:nz]; a[i][1:nzhalf]; a[i][nz+nzfhalf+1:ndof]; a[i][nz+1:nz+nzfhalf]]
    if corresp[i]==0
      nc+=1
      for j in i:lena
        if a[j]==b 
          corresp[j]=-i
          corresp[i]=j
          continue
        end  
      end 
    end  
  end 
end    

return a,corresp,lena,nc
end





