
function analysis(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3}) 

uneq=0
# equation numbering 
for e in 1:info.NE           
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      for d in 1:2
        if nodes[n].dof[d]==0   
          uneq=uneq+1             
          nodes[n].dof[d]=uneq   
        end  
      end
    end
  end   
end

psineq=0    # constant Psi
for e in 1:info.NL                      
  psineq+=1
  B3[e].psidof=uneq+psineq  
end 

neq=uneq+psineq
info.uneq=uneq
info.psineq=psineq
info.nK=uneq+psineq
info.nA=2*uneq+psineq

sol=zeros(Float64,neq)
F=zeros(Float64,neq)
F0=zeros(Float64,neq)

residuum=100
iter=0
while residuum>0.0001
#while iter<3

iter+=1

global K   # to make it visible out of the while loop

# faster and safer to reallocate at every iter, as the number of nozero coeffs may vary from iter to iter
K=SparseMatrixLNK(neq,neq)
C=SparseMatrixLNK(uneq,uneq)
if iter >0 
  fill!(F,0.0)
end

###################################
# mechanics
###################################

println("Assembling mechanics")

# allocation of elements arrays
Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12)
Ue=zeros(Float64,12)
Fe=zeros(Float64,12)
Ke= zeros(Float64,(12,12))

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
      Ue[(k-1)*2+1:k*2]=nodes[n].u
    end
#    T6_Ke!(Ke,Xe,material[mat])
    T6_KeNL!(Ke,Fe,Xe,Ue,material[mat])
    for i = 1:12
      dofi=dofe[i]
      if dofi>0
        F[dofi]+=Fe[i]
        for j = 1:12
#          F[dofi]-=Ke[i,j]*Ue[j]
          dofj=dofe[j]
          if dofj>0 
            K[dofi,dofj]+=Ke[i,j]
            C[dofi,dofj]+=Ke[i,j]
          end
        end
      end
    end
  end
end  

###################################
# ES forces
###################################

# allocation of elements arrays
Xe=zeros(Float64,3,2)
dofUe=zeros(Int64,6,1)
Ue=zeros(Float64,3,2)

Fe=zeros(Float64,6,1)
Ke=zeros(Float64,6,1)

for e in 1:info.NL                     
  dofPsi=B3[e].psidof
  Psi=B3[e].psi 
  for k in 1:3                                 
    n=B3[e].nodes[k]
    Xe[k,:]=nodes[n].coor              
    dofUe[(k-1)*2+1:k*2]=nodes[n].dof
  end

  B3_F!(Fe,Ke,Xe,Psi) 

  for i in 1:6
    dofUi=dofUe[i]
    if dofUi>0
      F[dofUi]+=Fe[i]
      K[dofUi,dofPsi]+=Ke[i] 
    end
  end    

end

##############################################################
# electro mixed problem
##############################################################

# allocation of elements arrays
Xe=zeros(Float64,3,2)
dofUe=zeros(Int64,6,1)
Ue=zeros(Float64,3,2)
KPUe=zeros(Float64,6,1)
KPPe=zeros(Float64,1,1)
FPe=zeros(Float64,1,1)

for e=1:info.NL                    
  dofPsi=B3[e].psidof
  Psi=B3[e].psi
  for k=1:3
    n=B3[e].nodes[k]   
    Xe[k,:]=nodes[n].coor                    
    Ue[k,:]=nodes[n].u
    dofUe[2*k-1:k*2]=nodes[n].dof
  end

  B3_UP!(KPUe,KPPe,FPe,Xe,Ue,Psi) 
  
  F[dofPsi]+=FPe[1]
  K[dofPsi,dofPsi]+=KPPe[1]
  if iter==1
    F0[dofPsi]+=FPe[1]/VDC   # for VAC in DPIM later
  end  
  for j in 1:6
    dofU=dofUe[j]
    if dofU>0
      K[dofPsi,dofU]+=KPUe[j] 
    end  
  end  
  
end

# Solution phase
println("Solving system")
K=SparseMatrixCSC(K)  # convert to other sparse format
C=SparseMatrixCSC(C)  # convert to other sparse format
#sol = K\F
ps = MKLPardisoSolver()
solve!(ps,sol,K,F)

# fill dofs
for n in 1:info.NN 
  for d in 1:2       
    dofU=nodes[n].dof[d]
    if dofU>0                      
      nodes[n].u[d]+=sol[dofU]      
    end 
  end
end
 
for e=1:info.NL            
  dofPsi=B3[e].psidof
  B3[e].psi+=sol[dofPsi]     
end

residuum=norm(sol)
println("Residuum: ",residuum)

end # while residuum

##############################
# mass assembler
##############################

M=SparseMatrixLNK(uneq,uneq)

# allocation of elements arrays
Xe=zeros(Float64,(6,2))
dofe=zeros(Int64,12,1)
Me= zeros(Float64,(12,12))

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
    end
    T6_Me!(Me,Xe,material[mat])
    for i = 1:12
      dofi=dofe[i]
      if dofi>0
        for j = 1:12
          dofj=dofe[j]
          if dofj>0 
            M[dofi,dofj]+=Me[i,j]
          end
        end
      end
    end
  end
end

M=SparseMatrixCSC(M)
C=info.beta*C+info.alpha*M

##############################
# GMSH output of static solution
##############################

outgmsh(info.NN,nodes,T6,B3,"_u0")
#  mycommand = `../gmsh.exe post.msh`
#  run(mycommand)

return C,K,M,F0

end   


