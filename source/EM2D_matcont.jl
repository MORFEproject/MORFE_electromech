
function matcont(nodes::Vector{Snode},M::SparseMatrixCSC{Float64},V::Matrix{ComplexF64},P::Vector{Parametrisation})

println("Init Matcont")
  
howmany=0
for p in 1:info.max_order   
  for i in 1:P[p].m 
    howmany+=1
  end
end

Avector=zeros(howmany,info.nrom)
fdyn=zeros(howmany,info.nza)
mappings=zeros(howmany,info.NN,3)
mappings_vel=zeros(howmany,info.NN,3)
mappings_modal=zeros(howmany,info.neig)
mappings_modal_vel=zeros(howmany,info.neig)

index=1
for p in 1:info.max_order   
  for i in 1:P[p].m 
    Avector[index,:]=P[p].Av[i]
    index+=1
  end
end

index=1
for p in 1:info.max_order  
  for i in 1:P[p].m 
    for j in 1:info.nmm
      fdyn[index,j]=2*real(P[p].fr[j,i])
    end  
    for j in 1:info.nmm
      fdyn[index,j+info.nmm]=-2*imag(P[p].fr[j,i])
    end
    index+=1
  end
end

index=1
for p in 1:info.max_order  
  for i in 1:P[p].m 
    for inode=1:info.NN 
      for idof=1:2 
        dof=nodes[inode].dof[idof] 
        if dof>0 
          mappings_vel[index,inode,idof]=real(P[p].Wr[dof,i])
          mappings[index,inode,idof]=real(P[p].Wr[info.uneq+dof,i])
        end
      end
    end
#    for imode=1:info.neig # all computed modes
#      mappings_modal_vel[index,imode]=V[:,imode]'*M*real.(P[p].Wr[1:info.neq,i])
#      mappings_modal[index,imode]=V[:,imode]'*M*real.(P[p].Wr[info.neq+1:2*info.neq,i])
#    end
    index+=1
  end      
end

file = matopen("./output/param.mat","w")
write(file, "ndof",2*info.uneq)
write(file, "nz",info.nza)
write(file, "max_order",info.max_order)
write(file, "mappings",mappings)
#write(file, "mappings_modal",mappings_modal)
write(file, "mappings_vel",mappings)
#write(file, "mappings_modal_vel",mappings_modal_vel)
write(file, "Avector",Avector)
write(file, "fdyn",fdyn)
close(file)

end

