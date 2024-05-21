
function realification!(P::Vector{Parametrisation})

  println("Init realification")
  
  for p in 1:info.max_order
    println("order: ",p)
    for i in 1:P[p].m
      Av=P[p].Av[i][:]
      Ivec=zeros(Int64,p)
      counter=0
      for j in 1:info.nrom 
        Ivec[counter+1:counter+Av[j]].=j
        counter+=Av[j]
      end
      pos=1
      coeff=1.0+0.0im
      Av.=0
      recursive_C2R!(Ivec,p,pos,i,Av,coeff,P[p])
    end  
  end

end  



function recursive_C2R!(Ivec::Vector{Int64},p::Int64,pos::Int64,posinit::Int64,Av::Vector{Int64},
                        coeff::ComplexF64,P::Parametrisation)

  nzhalf=Int(info.nza/2)
  nzfhalf=Int(info.nzna/2)                           
  if pos==(p+1)
    pos1=findfirst(x->x==Av,P.Av)
    P.Wr[:,pos1]+=coeff*P.W[:,posinit]
    P.fr[:,pos1]+=coeff*P.f[:,posinit]
  else
    Av1=Av[:]   # new vectors
    Av2=Av[:]
    iz=Ivec[pos]    # z var 
   
    if iz <= nzhalf
      coeff1=0.5*coeff  
      Av1[iz]+=1
      coeff2=-0.5*im*coeff
      Av2[iz+nzhalf]+=1
    elseif iz <= info.nza
      coeff1=0.5*coeff
      Av1[iz-nzhalf]+=1
      coeff2=0.5*im*coeff
      Av2[iz]+=1
    elseif iz <= info.nza+nzfhalf  # to have cos and sin the 1/2 must avoided
      coeff1=im*coeff
      Av1[iz]+=1
      coeff2=coeff
      Av2[iz+nzfhalf]+=1
    else    
      coeff1=-im*coeff
      Av1[iz-nzfhalf]+=1
      coeff2=coeff
      Av2[iz]+=1
    end     

    pos+=1
    recursive_C2R!(Ivec,p,pos,posinit,Av1,coeff1,P)
    recursive_C2R!(Ivec,p,pos,posinit,Av2,coeff2,P)
  end

end    

