

# Electrostatic force on a B3
function B3_F!(FUe::Matrix{Float64},KUPe::Matrix{Float64},Xe::Matrix{Float64},Psi::Float64)

fill!(FUe,0.0)
fill!(KUPe,0.0)
   
qr=quadrature_points(Val(:LINE5gp))
for (w,a) in qr

  D=[a-.5  a+.5 -2*a]                # derivative of shape functions
  F=D*Xe
  J=(F[1]^2+F[2]^2)^0.5
  coeff=Psi*Psi/2
  NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2]   # shape functions
  N=[0 NL[1] 0 NL[2] 0 NL[3]]
  FUe[:]+=N'*coeff*eps0*J*w
  KUPe[:]-=N'*Psi*eps0*J*w     

end
end



# matrices for electrostatic mixed FP approach: P -> Psi F -> Phi
function B3_UP!(KPUe::Matrix{Float64},KPPe::Matrix{Float64},FPe::Matrix{Float64},
                Xe::Matrix{Float64},Ue::Matrix{Float64},psi::Float64) 

fill!(KPUe,0.0)  
KPPe[1]=0.0
FPe[1]=0.0
N=zeros(Float64,6,1)  

qr=quadrature_points(Val(:LINE5gp))
for (w,a) in qr

  D=[a-.5  a+.5 -2*a]                # derivative of shape functions
  F=D*Xe
  J=(F[1]^2+F[2]^2)^0.5
  NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2]   # shape functions
  N=[0 NL[1] 0 NL[2] 0 NL[3]]
  u=dot(NL,Ue[:,2])
  FPe[1]+=(VDC-(gap-u)*psi)*J*w
  KPPe[1]+=(gap-u)*J*w
  KPUe[:]-=N'*psi*J*w

end


end


function T6_KeNL!(Ke::Matrix{Float64},Fe::Vector{Float64},
                  Xe::Matrix{Float64},Ue::Vector{Float64},mate::Vector{Float64})

fill!(Ke,0.0)      
fill!(Fe,0.0)

#nn=6
#dim=2

dNda=zeros(Float64,6,2)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
G=zeros(Float64,2,2,12)
B=zeros(Float64,3,12)
Kp=zeros(Float64,12,12)

∇U = zeros(Float64,(2,2))
Et = zeros(Float64,(2,2))
E = zeros(Float64,3)
S = zeros(Float64,3)

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
  dNdx[:,:] = dNda*invJac   

  for i in 1:2, j in 1:2, m in 1:6
    pos=i+(m-1)*2
    G[i,j,pos]=dNdx[m,j]
  end
  for i in 1:2, j in 1:2
    ∇U[i,j]=dot(G[i,j,:],Ue)
  end

  for i in 1:2, j in i:2
    Et[i,j]=0.5*(∇U[i,j]+∇U[j,i]+dot(∇U[:,i],∇U[:,j]))
  end
  E=[Et[1,1],Et[2,2],2*Et[1,2]]

  S=D*E+S0

  B[1,:]=G[1,1,:]
  B[2,:]=G[2,2,:]
  B[3,:]=G[1,2,:]+G[2,1,:]

  fill!(Kp,0.0)

  for m in 1:2
    B[1,:]+=G[m,1,:]*∇U[m,1]
    B[2,:]+=G[m,2,:]*∇U[m,2]
    B[3,:]+=G[m,1,:]*∇U[m,2]+G[m,2,:]*∇U[m,1]
    for i in 1:12, j in 1:i
      Kp[i,j]+= S[1]*G[m,1,i]*G[m,1,j]+S[2]*G[m,2,i]*G[m,2,j]+
                S[3]*(G[m,1,i]*G[m,2,j]+G[m,2,i]*G[m,1,j])
    end        
  end

  coeff=w*J
  for i in 1:12, j in 1:i
    for m in 1:3, n in 1:3
      Kp[i,j]+=B[m,i]*D[m,n]*B[n,j]
    end
  end

  Ke[:,:]+=Kp*coeff

  for i in 1:12, m in 1:3
    Fe[i]-=B[m,i]*S[m]*coeff 
  end 

end # loopgauss

for i in 1:12, j in 1:i-1
  Ke[j,i]=Ke[i,j]
end

end



function T6_Ke!(Ke::Matrix{Float64},X::Matrix{Float64},mate::Vector{Float64})
 
fill!(Ke,0.0)

dNda=zeros(Float64,6,2)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
B=zeros(Float64,12,3)

mY,nu=mate[1:2]
D=mY/((1+nu)*(1-2*nu))*[1-nu nu 0;
                       nu 1-nu 0;
                       0 0 (1-2*nu)/2] 
                                            
qr=quadrature_points(Val(:TRI6gp))
for (w,a) in qr
   
  dNda!(dNda,a,Val(:TRI6n))
  Jac[:] = X'*dNda
  J =  Jac[1,1]*Jac[2,2]-Jac[2,1]*Jac[1,2]
  invJac[:] =[Jac[2,2] -Jac[1,2]; 
             -Jac[2,1] Jac[1,1]]/J
  dNdx[:] = dNda*invJac   
  for i = 1:6
    B[1+(i-1)*2,1] = dNdx[i,1]
    B[2+(i-1)*2,2] = dNdx[i,2]
    B[1+(i-1)*2,3] = dNdx[i,2]
    B[2+(i-1)*2,3] = dNdx[i,1]
  end
  Ke[:,:]+=(B*D*B')*J*w
  
end  # loop gauss points 

return nothing

end


# mass matrix: T6
function T6_Me!(Me::Matrix{Float64},X::Matrix{Float64},mate::Vector{Float64})

fill!(Me,0.0)
dNda=zeros(Float64,(6,2))
NL=zeros(Float64,6)

rho=mate[3]

qr=quadrature_points(Val(:TRI6gp))
for (w,a) in qr         
  dNda!(dNda,a,Val(:TRI6n))            
  F=X'*dNda                                  
  J=F[1,1]*F[2,2]-F[1,2]*F[2,1]           
  N!(NL,a,Val(:TRI6n))            
  N=[NL[1] 0 NL[2] 0 NL[3] 0 NL[4] 0 NL[5] 0 NL[6] 0;
     0 NL[1] 0 NL[2] 0 NL[3] 0 NL[4] 0 NL[5] 0 NL[6]]
  Me[:,:]+=(N'*N)*rho*J*w            
end

end
   


