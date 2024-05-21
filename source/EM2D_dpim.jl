

function dpim(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
              F0::Vector{Float64},D::Vector{ComplexF64},VR::Matrix{ComplexF64},VL::Matrix{ComplexF64},
              nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3})

  info.nmm=length(info.Lmm)   # master modes
  nmm=info.nmm

  nK=info.nK
  nA=info.nA
  uneq=info.uneq

  info.nza=2*nmm
  info.nzna=2  # imposes only two nonautonomous
  info.nrom=info.nza+info.nzna

  info.nMat=nA+info.nza  # dim of system to be solved

  println("Init Parametrisation")
  P=initParametrisation!(info)

  println("Init System")
  Rhs = Array{ComplexF64}(undef,info.nMat)
  Sol = Array{ComplexF64}(undef,info.nMat)
  Mat=spzeros(ComplexF64,info.nMat,info.nMat)

  println("Order 1")
  omega0 = zeros(Float64,nmm)
  BY = Array{ComplexF64}(undef,2*uneq,info.nza)
  XTB = Array{ComplexF64}(undef,info.nza,2*uneq)

  for i = 1:nmm
    mm=info.Lmm[i]
    omega0[i] = abs(imag(D[2*mm-1]))
    lambda1 = D[2*i-1]
    lambda2 = D[2*i]

    P[1].f[i,i] = lambda1
    P[1].f[i+nmm,i+nmm] = lambda2
    P[1].W[:,i] = VR[:,2*mm-1]
    P[1].W[:,i+nmm] = VR[:,2*mm]

    BY[1:uneq,i] = M*VR[1:uneq,2*mm-1]
    BY[uneq+1:2*uneq,i] = M*VR[uneq+1:2*uneq,2*mm-1]
    BY[1:uneq,i+nmm] = M*VR[1:uneq,2*mm]
    BY[uneq+1:2*uneq,i+nmm] = M*VR[uneq+1:2*uneq,2*mm]

    XTB[i,1:uneq] = transpose(VL[1:uneq,2*mm-1])*M
    XTB[i,uneq+1:2*uneq] = transpose(VL[uneq+1:2*uneq,2*mm-1])*M
    XTB[i+nmm,1:uneq] = transpose(VL[1:uneq,2*mm])*M
    XTB[i+nmm,uneq+1:2*uneq] = transpose(VL[uneq+1:2*uneq,2*mm])*M
   end

  Lambda=Vector{ComplexF64}(undef,info.nrom)
  for i in 1:info.nza
    Lambda[i]=P[1].f[i,i]
  end
  if info.nzna==2
   Lambda[info.nza+1]=info.Fmult*im*omega0[info.Ffreq]
   P[1].f[info.nza+1,info.nza+1] = Lambda[info.nza+1]
   Lambda[info.nza+2]=-info.Fmult*im*omega0[info.Ffreq]
   P[1].f[info.nza+2,info.nza+2] = Lambda[info.nza+2]
  end

# new first order terms due to forcing: two NA var
  P[1].R[2*uneq+1:nA,info.nza+1]=0.5*F0[uneq+1:nK]*VAC  # rhs only on lower part
  homological!(Sol,Rhs,Mat,Lambda,P[1],info.nza+1,K,M,C,BY,XTB)
  P[1].R[2*uneq+1:nA,info.nza+2]=0.5*F0[uneq+1:nK]*VAC
  homological!(Sol,Rhs,Mat,Lambda,P[1],info.nza+2,K,M,C,BY,XTB)

#  for i in 1:info.nzna
##    for j in axes(info.Fmodes)[1]
#      P[1].R[neq+1:2*neq,info.nza+i]+=0.001*M*V[:,1]
##    end
#    homological!(Sol,Rhs,Mat,Lambda,P[1],info.nza+i,M,K,BY,XTB,info)
#  end

  println("Higher orders")
  for p = 2:info.max_order
    println("Order $p")
    fillrhs_quad!(nodes,T6,B3,P,p)
    fillrhs_cub!(nodes,T6,B3,P,p)
    fillWf!(P,p)
    for i in 1:P[p].m # for every alpha vector
      corresp=P[p].corresp[i]
      if corresp>0
        fillWfnonaut!(P,p,P[p].Av[i],i)
        homological!(Sol,Rhs,Mat,Lambda,P[p],i,K,M,C,BY,XTB)
        P[p].analysed[i]=1
      elseif corresp<0
        P[p].W[:,i]=conj(P[p].W[:,-corresp])
        P[p].f[1:nmm,i]=conj(P[p].f[nmm+1:2*nmm,-corresp])
        P[p].f[nmm+1:2*nmm,i]=conj(P[p].f[1:nmm,-corresp])
        P[p].analysed[i]=1
      end
    end

  end

 return  P

end



function homological!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                      Lambda::Vector{ComplexF64},P::Parametrisation,Apos::Int64,
                      K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
                      BY::Matrix{ComplexF64},XTB::Matrix{ComplexF64})


  uneq=info.uneq
  nA=info.nA
  nK=info.nK

  σ=dot(P.Av[Apos],Lambda[1:info.nrom])
  resonant_modes = zeros(Bool,info.nza)
  fill!(resonant_modes,false)
  for i = 1:info.nza
    λ = Lambda[i]
    if (abs(σ-λ)/abs(λ)<=info.tol)
      resonant_modes[i] = true
    end
  end
  println(P.Av[Apos],"  ",resonant_modes)

  fill!(Rhs,0.0)
  Rhs[1:nA]=P.R[:,Apos]

  Rhs[1:uneq]-=M*P.Wf[1:uneq,Apos]
  Rhs[uneq+1:2*uneq]-=M*P.Wf[uneq+1:2*uneq,Apos]

  fill!(Mat.nzval,0.0)
  Mat[1:uneq,1:uneq] = σ*M+C
  Mat[uneq+1:2*uneq,uneq+1:2*uneq] = σ*M
  Mat[uneq+1:2*uneq,1:uneq] = -M
  Mat[1:uneq,uneq+1:nA] = K[1:uneq,:]
  Mat[2*uneq+1:nA,uneq+1:nA] = K[uneq+1:nK,:]
  for j = 1:info.nza
    if resonant_modes[j]  # if resonant
      Mat[1:2*uneq,nA+j]=BY[:,j]
      Mat[nA+j,1:2*uneq]=XTB[j,:]
    else
      Mat[nA+j,nA+j]=1
    end
  end

#  println("Solving system2")
  ps = MKLPardisoSolver()
  solve!(ps,Sol,Mat,Rhs)

#  file = matopen("./output/matrix.mat","w")
#  write(file, "uneq",info.uneq)
#  write(file, "psineq",info.psineq)
#  write(file, "neq",neq)
#  write(file, "nza",info.nza)
#  write(file, "nzna",info.nzna)
#  write(file, "Mat",Mat)
#  write(file, "Rhs",Rhs)
#  close(file)

#  Sol = Mat\Rhs   ## DOES IT WORKKKKKKK????

#  P.W[:,Apos]=Sol[1:2*neq]
  P.W[:,Apos]=Sol[1:nA]
#  P.f[1:info.nza,Apos]=Sol[2*neq+1:2*neq+info.nza]
  P.f[1:info.nza,Apos]=Sol[nA+1:nA+info.nza]

end


function fillrhs_quad!(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},
                       P::Vector{Parametrisation},p::Int64)

  nA=info.nA
  uneq=info.uneq
  for p1 in 1:p-1
    p2=p-p1
    for k1 in 1:P[p1].m, k2 in 1:P[p2].m
      Av=P[p1].Av[k1]+P[p2].Av[k2]
      pos=findfirst(x->x==Av,P[p].Av)
      if P[p].corresp[pos]>0
        W1 = P[p1].W[uneq+1:nA,k1] # second part is U+psi, first is V
        W2 = P[p2].W[uneq+1:nA,k2]
        assembly_quad!(P[p],pos,W1,W2,nodes,T6)
        assembly_cub0!(P[p],pos,W1,W2,nodes,T6)
        assembly_quad_EM!(P[p],pos,W1,W2,nodes,T6,B3)
      end
    end
  end

end


function fillrhs_cub!(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},
                      P::Vector{Parametrisation},p::Int64)

  nA=info.nA
  uneq=info.uneq
  for p1 in 1:p-2, p2 in 1:p-1-p1  #, p3 in 1:p-2
    p3=p-p1-p2
    for k1 in 1:P[p1].m, k2 in 1:P[p2].m, k3 in 1:P[p3].m
      Av=P[p1].Av[k1]+P[p2].Av[k2]+P[p3].Av[k3]
      pos=findfirst(x->x==Av,P[p].Av)
      if P[p].corresp[pos]>0
        W1 = P[p1].W[uneq+1:nA,k1]
        W2 = P[p2].W[uneq+1:nA,k2]
        W3 = P[p3].W[uneq+1:nA,k3]
        assembly_cub!(P[p],pos,W1,W2,W3,nodes,T6)
      end
    end
  end

end


function fillWf!(P::Vector{Parametrisation},p::Int64)

  for p1 in 2:p-1, p2 in 2:p-1
    if (p1+p2)==p+1
      for i in 1:P[p1].m
        Av1=P[p1].Av[i][:]
        for j in 1:P[p2].m
          Av2=P[p2].Av[j][:]
          for s in 1:info.nrom
            if Av1[s]>0
              Av=Av1+Av2
              Av[s]-=1
              pos=findfirst(x->x==Av,P[p].Av)
#              P[p].Wf[1:neq,pos]+=Av1[s]*P[p1].W[1:neq,i]*P[p2].f[s,j]
#              P[p].Wf[neq+1:2*neq,pos]+=Av1[s]*P[p1].W[neq+1:2*neq,i]*P[p2].f[s,j]
               P[p].Wf[:,pos]+=Av1[s]*P[p1].W[:,i]*P[p2].f[s,j]
            end
          end
        end
      end
    end
  end

end


function fillWfnonaut!(P::Vector{Parametrisation},p::Int64,Av::Vector{Int64},ind_rhs::Int64)

  for r in info.nza+1:info.nrom
    if Av[r]>0
      for s in 1:info.nza
        fs_r = P[1].f[s,r]
        if abs(fs_r)>10^(-20)    # only fills if the reduced dyn is nzero
          Av_W=Av[:]
          Av_W[s]+=1
          Av_W[r]-=1
          ind_W=findfirst(x->x==Av_W,P[p].Av)
          if P[p].analysed[ind_W]==0
            println("checkkkk")
          end
#          P[p].Wf[1:neq,ind_rhs]+=P[p].W[1:neq,ind_W]*fs_r*Av_W[s]
#          P[p].Wf[neq+1:2*neq,ind_rhs]+=P[p].W[neq+1:2*neq,ind_W]*fs_r*Av_W[s]
          P[p].Wf[:,ind_rhs]+=P[p].W[:,ind_W]*fs_r*Av_W[s]
        end
      end
    end
  end

end


