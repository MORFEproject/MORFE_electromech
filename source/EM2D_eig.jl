
function eigAB(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},nodes::Vector{Snode})

println("Computing eigenvalues")

nK=info.nK
nA=info.nA
uneq=info.uneq

B=spzeros(Float64,nA,nA)

B[1:uneq,1:uneq] = M
B[uneq+1:2*uneq,uneq+1:2*uneq] = M

#B[1:uneq,uneq+1:2*uneq] = M
#B[uneq+1:2*uneq,1:uneq] = M

A=spzeros(Float64,nA,nA)

A[1:uneq,1:uneq] = -C
A[uneq+1:2*uneq,1:uneq] = M
A[1:uneq,uneq+1:nA] = -K[1:uneq,:]
A[2*uneq+1:nA,uneq+1:nA] = -K[uneq+1:nK,:]

#A[uneq+1:2*uneq,1:uneq] = -C
#A[1:uneq,1:uneq] = M
#A[uneq+1:nA,uneq+1:nA] = -K

mat"""
[$VR,$DR] = eigs($A,$B,2*$info.neig,'SM');
$DR=diag($DR);
"""

#DR,VR = eigs(A,B,nev=2*info.neig,which=:SM)
println(DR)

mat"""
[$VL,$DL] = eigs(transpose($A),transpose($B),2*$info.neig,'SM');
$DL=diag($DL);
"""

#DL,VL = eigs(transpose(A),transpose(B),nev=2*info.neig,which=:SM)
println(DL)

#corresp=zeros(Int,2*info.neig)
#for i=1:2:2*info.neig
#  lam1=DR[i]
#  lam2=DL[i]
#  if abs(lam1-conj(lam2)) > 1e-6
#    temp=VL[:,i]
#    VL[:,i]=VL[:,i+1]
#    VL[:,i+1]=temp
#    println("switch: ",i,",",i+1)
#  end
#end

# normalisation
for i = 1:2*info.neig
#  for j = 1:2*info.neig
    cc = transpose(VL[:,i])*B*VR[:,i]  # compl conj
#    cd = transpose(VL[:,i])*A*VR[:,j]
#    println(i,",",j)
#    println(cc)
#    println(cd/DL[i])
#  end
  for j = 1:nA   # arbitrary: scales only VR  (seems to coincide with old approach)
    VR[j,i] /= sqrt(cc)
    VL[j,i] /= sqrt(cc)
  end
  #outgmsheig(info.NN,nodes,real(VR[uneq+1:2*uneq,i]),"_eig"*string(i))  # prints displacements
end

return DR,VR,VL

end



