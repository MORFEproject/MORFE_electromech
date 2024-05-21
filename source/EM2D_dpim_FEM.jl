
function assembly_quad_EM!(P::Parametrisation,pos::Int64,W1::Vector{ComplexF64},W2::Vector{ComplexF64},
                           nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3})


###################################
# ES forces
###################################

Fe=zeros(ComplexF64,6,1)
Xe=zeros(Float64,3,2)
GUe=zeros(Int64,6,1)

for e in 1:info.NL

  GPsie=B3[e].psidof
  for k in 1:3
    n=B3[e].nodes[k]   
    Xe[k,:]=nodes[n].coor                
    GUe[2*k-1:2*k]=nodes[n].dof
  end
  Psi1e=W1[GPsie]
  Psi2e=W2[GPsie]

  B3_EF_quad!(Fe,Xe,Psi1e,Psi2e) 

  for i in 1:6
    dofU=GUe[i]
    if dofU>0
      P.R[dofU,pos]+=Fe[i]
    end
  end    

end  # loop over NL


##############################################################
# electro mixed problem
##############################################################

# allocation of elements arrays
Xe=zeros(Float64,3,2)
GUe=zeros(Int64,6,1)
U1e=zeros(ComplexF64,3,2)
U2e=zeros(ComplexF64,3,2)
FPsie=zeros(ComplexF64,1,1)

for e=1:info.NL                    
  GPsie=B3[e].psidof
  for k=1:3
    n=B3[e].nodes[k]   
    Xe[k,:]=nodes[n].coor                    
    GUe[2*k-1:k*2]=nodes[n].dof
  end
  fill!(U1e,0.0)
  fill!(U2e,0.0)
  for i in 1:3
    for j in 1:2
      dofU=GUe[2*(i-1)+j]
      if dofU>0
        U1e[i,j]=W1[dofU] 
        U2e[i,j]=W2[dofU] 
      end
    end 
  end  
  Psi1e=W1[GPsie]
  Psi2e=W2[GPsie]

  B3_UP_quad!(FPsie,Xe,Psi1e,Psi2e,U1e,U2e) 

#  P.R[GPsie,pos]-=FPsie[1]
  P.R[GPsie+info.uneq,pos]-=FPsie[1]  # MODDDDDDDDDDDDDDD

end

end



# matrices for electrostatic mixed FP approach: P -> Psi F -> Phi
function B3_UP_quad!(FPsie::Matrix{ComplexF64},Xe::Matrix{Float64},
                     Psi1::ComplexF64,Psi2::ComplexF64,U1e::Matrix{ComplexF64},U2e::Matrix{ComplexF64}) 


fill!(FPsie,0.0)

qr=quadrature_points(Val(:LINE3gp))
for (w,a) in qr
  D=[a-.5  a+.5 -2*a]                # derivative of shape functions
  F=D*Xe
  J=(F[1]^2+F[2]^2)^0.5
  NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2]   # shape functions
  U1=dot(NL,U1e[:,2])
  U2=dot(NL,U2e[:,2])
  FPsie[1]-=(Psi1*U2+Psi2*U1)/2*J*w
end

end


# Electrostatic force on a B3
function B3_EF_quad!(Fe::Matrix{ComplexF64},Xe::Matrix{Float64},Psi1e::ComplexF64,Psi2e::ComplexF64)

fill!(Fe,0.0)

qr=quadrature_points(Val(:LINE3gp))
for (w,a) in qr
  D=[a-.5  a+.5 -2*a]                # derivative of shape functions
  F=D*Xe
  J=(F[1]^2+F[2]^2)^0.5
  NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2]   # shape functions
  N=[0 NL[1] 0 NL[2] 0 NL[3]]
  Fe[:]+=N'*Psi1e*Psi2e/2*eps0*J*w
end

end


#
#
# standard mechanics
#
#


function assembly_cub!(P::Parametrisation,pos::Int64,U1::Vector{ComplexF64},U2::Vector{ComplexF64},U3::Vector{ComplexF64},
                       nodes::Vector{Snode},T6::Vector{ST6})

Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12)
Fe=zeros(ComplexF64,12)
U1e=zeros(ComplexF64,12)
U2e=zeros(ComplexF64,12)
U3e=zeros(ComplexF64,12)

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
    end

    for k in 1:12
     dof=dofe[k] 
     if dof>0
        U1e[k]=U1[dof] 
        U2e[k]=U2[dof] 
        U3e[k]=U3[dof] 
      else
        U1e[k]=0.0 
        U2e[k]=0.0 
        U3e[k]=0.0 
      end
    end 

    T6_He!(Fe,Xe,U1e,U2e,U3e,material[mat])

    for i in 1:12
      dofU=dofe[i]
      if dofU>0
        P.R[dofU,pos]-=Fe[i]
      end
    end    

  end   
end  # loop over NE
end


function assembly_cub0!(P::Parametrisation,pos::Int64,U1::Vector{ComplexF64},U2::Vector{ComplexF64},
                        nodes::Vector{Snode},T6::Vector{ST6})

Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12)
Fe=zeros(ComplexF64,12)
U1e=zeros(ComplexF64,12)
U2e=zeros(ComplexF64,12)
U3e=zeros(ComplexF64,12)

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
      U3e[(k-1)*2+1:k*2]=nodes[n].u
    end

    for k in 1:12
     dof=dofe[k] 
     if dof>0
        U1e[k]=U1[dof] 
        U2e[k]=U2[dof] 
      else
        U1e[k]=0.0 
        U2e[k]=0.0 
      end
    end 

    T6_He!(Fe,Xe,U1e,U2e,U3e,material[mat])

    for i in 1:12
      dofU=dofe[i]
      if dofU>0
        P.R[dofU,pos]-=3*Fe[i]
      end
    end    

  end   
end  # loop over NE
end



function assembly_quad!(P::Parametrisation,pos::Int64,U1::Vector{ComplexF64},U2::Vector{ComplexF64},
                        nodes::Vector{Snode},T6::Vector{ST6})

Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12)
Fe=zeros(ComplexF64,12)
U1e=zeros(ComplexF64,12)
U2e=zeros(ComplexF64,12)

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
    end

    for k in 1:12
     dof=dofe[k] 
     if dof>0
        U1e[k]=U1[dof] 
        U2e[k]=U2[dof] 
      else
        U1e[k]=0.0 
        U2e[k]=0.0 
      end
    end 

    T6_Ge!(Fe,Xe,U1e,U2e,material[mat])

    for i in 1:12
      dofU=dofe[i]
      if dofU>0
        P.R[dofU,pos]-=Fe[i]
      end
    end    
  end   

end  # loop over NE
end




function T6_He!(Fe::Vector{ComplexF64},Xe::Matrix{Float64},
                U1e::Vector{ComplexF64},U2e::Vector{ComplexF64},U3e::Vector{ComplexF64},mate::Vector{Float64})
  
fill!(Fe,0.0)
  
dNda=zeros(Float64,6,2)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
G=zeros(Float64,2,2,12)

∇U1 = zeros(ComplexF64,2,2)
∇U2 = zeros(ComplexF64,2,2)
∇U3 = zeros(ComplexF64,2,2)

Bnl1=zeros(ComplexF64,3,12)
Bnl2=zeros(ComplexF64,3,12)
Bnl3=zeros(ComplexF64,3,12)

mY,nu=mate[1:2]
D=mY/((1+nu)*(1-2*nu))*[1-nu nu 0;
                       nu 1-nu 0;
                       0 0 (1-2*nu)/2] 
                                            
qr=quadrature_points(Val(:TRI6gp))
  
for (w,a) in qr

  dNda!(dNda,a,Val(:TRI6n))
  Jac[:] = Xe'*dNda
  J =  Jac[1,1]*Jac[2,2]-Jac[2,1]*Jac[1,2]
  invJac[:] =[Jac[2,2] -Jac[1,2]; 
             -Jac[2,1] Jac[1,1]]/J
  dNdx[:] = dNda*invJac   

  for i in 1:2, j in 1:2, m in 1:6
    pos=i+(m-1)*2
    G[i,j,pos]=dNdx[m,j]
  end
  for i in 1:2, j in 1:2
    ∇U1[i,j]=dot(G[i,j,:],U1e)
    ∇U2[i,j]=dot(G[i,j,:],U2e)
    ∇U3[i,j]=dot(G[i,j,:],U3e)
  end

  temp=0.5*(transpose(∇U1)*∇U2+transpose(∇U2)*∇U1)
  S12=D*[temp[1,1],temp[2,2],2*temp[1,2]]
  temp=0.5*(transpose(∇U2)*∇U3+transpose(∇U3)*∇U2)
  S23=D*[temp[1,1],temp[2,2],2*temp[1,2]]
  temp=0.5*(transpose(∇U1)*∇U3+transpose(∇U3)*∇U1)
  S13=D*[temp[1,1],temp[2,2],2*temp[1,2]]

  fill!(Bnl1,0.0)
  fill!(Bnl2,0.0)
  fill!(Bnl3,0.0)
  for m in 1:2
    Bnl1[1,:]+=G[m,1,:]*∇U1[m,1]
    Bnl1[2,:]+=G[m,2,:]*∇U1[m,2]
    Bnl1[3,:]+=G[m,1,:]*∇U1[m,2]+G[m,2,:]*∇U1[m,1]

    Bnl2[1,:]+=G[m,1,:]*∇U2[m,1]
    Bnl2[2,:]+=G[m,2,:]*∇U2[m,2]
    Bnl2[3,:]+=G[m,1,:]*∇U2[m,2]+G[m,2,:]*∇U2[m,1]

    Bnl3[1,:]+=G[m,1,:]*∇U3[m,1]
    Bnl3[2,:]+=G[m,2,:]*∇U3[m,2]
    Bnl3[3,:]+=G[m,1,:]*∇U3[m,2]+G[m,2,:]*∇U3[m,1]
  end

  Fe[:] += (1.0/6.0)*(transpose(Bnl1)*S23+transpose(Bnl2)*S13+transpose(Bnl3)*S12)*w*J
end

end




function T6_Ge!(Fe::Vector{ComplexF64},Xe::Matrix{Float64},U1e::Vector{ComplexF64},U2e::Vector{ComplexF64},mate::Vector{Float64})
  
fill!(Fe,0.0)
  
dNda=zeros(Float64,6,2)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
G=zeros(Float64,2,2,12)

∇U1 = zeros(ComplexF64,2,2)
∇U2 = zeros(ComplexF64,2,2)

B=zeros(Float64,3,12)
Bnl1=zeros(ComplexF64,3,12)
Bnl2=zeros(ComplexF64,3,12)

mY,nu=mate[1:2]
D=mY/((1+nu)*(1-2*nu))*[1-nu nu 0;
                       nu 1-nu 0;
                       0 0 (1-2*nu)/2] 
                                            
qr=quadrature_points(Val(:TRI6gp))
  
for (w,a) in qr

  dNda!(dNda,a,Val(:TRI6n))
  Jac[:]=Xe'*dNda
  J=Jac[1,1]*Jac[2,2]-Jac[2,1]*Jac[1,2]
  invJac[:]=[Jac[2,2] -Jac[1,2]; 
             -Jac[2,1] Jac[1,1]]/J
  dNdx[:,:]=dNda*invJac   

  for i in 1:2, j in 1:2, m in 1:6
    pos=i+(m-1)*2
    G[i,j,pos]=dNdx[m,j]
  end
  for i in 1:2, j in 1:2
    ∇U1[i,j]=dot(G[i,j,:],U1e)
    ∇U2[i,j]=dot(G[i,j,:],U2e)
  end

  temp=0.5*(transpose(∇U1)+∇U1)
  S1=D*[temp[1,1],temp[2,2],2*temp[1,2]]
  temp=0.5*(transpose(∇U2)+∇U2)
  S2=D*[temp[1,1],temp[2,2],2*temp[1,2]]
  temp=0.5*(transpose(∇U1)*∇U2+transpose(∇U2)*∇U1)
  S12=D*[temp[1,1],temp[2,2],2*temp[1,2]]

  B[1,:]=G[1,1,:]
  B[2,:]=G[2,2,:]
  B[3,:]=G[1,2,:]+G[2,1,:]

  fill!(Bnl1,0.0)
  fill!(Bnl2,0.0)
  for m in 1:2
    Bnl1[1,:]+=G[m,1,:]*∇U1[m,1]
    Bnl1[2,:]+=G[m,2,:]*∇U1[m,2]
    Bnl1[3,:]+=G[m,1,:]*∇U1[m,2]+G[m,2,:]*∇U1[m,1]
    Bnl2[1,:]+=G[m,1,:]*∇U2[m,1]
    Bnl2[2,:]+=G[m,2,:]*∇U2[m,2]
    Bnl2[3,:]+=G[m,1,:]*∇U2[m,2]+G[m,2,:]*∇U2[m,1]
  end
  
  Fe[:] += 0.5*(transpose(B)*S12+transpose(Bnl1)*S2+transpose(Bnl2)*S1)*w*J

end

end


