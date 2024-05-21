function output(P::Vector{Parametrisation})

##################################################
#  OUTPUT ON FILE
##################################################

# qui lo scrivo in modo compatto
# mi segno gli alfa vector

@time howmany=count_terms_dyn(P,info)
@time mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn=store_dyn_and_map(P,info,howmany)

#export_maps(mappings,Avector,info,howmany)

#file = matopen(info.output_dir*"/param.mat", "w")
file = matopen(info.output_dir*"/param.mat", "w")
write(file, "ndof",info.nA)
write(file, "max_order",info.max_order)
write(file, "mappings",mappings)
write(file, "mappings_modal",mappings_modal)
write(file, "mappings_vel",mappings)
write(file, "mappings_modal_vel",mappings_modal_vel)
write(file, "Avector",Avector)
write(file, "fdyn",fdyn)
close(file)

println("done!")

end



function count_terms_dyn(P::Vector{Parametrisation},info::Sinfo)
  # conto quanti sono
  howmany=0
  for p in 1:info.max_order   # for every order
      for i in 1:P[p].nc # for every alpha vector
        #corresp=P[p].corresp[i] # that is resonant or under the given order NA
        #  if corresp!=0
            howmany=howmany+1
        #  end
      end
  end
  return howmany
end




function store_dyn_and_map(P::Vector{Parametrisation},info::Sinfo,howmany::Integer)

mappings=zeros(howmany,mesh.nn,3)
mappings_vel=zeros(howmany,mesh.nn,3)
mappings_modal=zeros(howmany,info.neig)
mappings_modal_vel=zeros(howmany,info.neig)
Avector=zeros(howmany,P[1].nc)
fdyn=zeros(howmany,info.nm*2)

U = Field(mesh,dim)
println("Assemblying M K")
colptr,rowval = assembler_dummy_MK(mesh,U)
val=zeros(Float64,length(rowval))
K=SparseMatrixCSC(info.nK,info.nK,colptr,rowval,val)
M = deepcopy(K)
assembler_MK!(mesh,U,K,M)

#neig = maximum([info.Φ; info.neig])   
neig = maximum([info.Lmm; info.neig])   
println("Computing eigenvalues")
λ, ϕ = eigs(K,M,nev=neig,which=:SM)
λ = real(λ)
ϕ = real(ϕ)
#set the mode directions in a predictable way
for ieig=1:neig
  if sum(ϕ[:,ieig])<0.0
     ϕ[:,ieig]=-ϕ[:,ieig]
  end
end

index=1
for p in 1:info.max_order   # for every order
for i in 1:P[p].nc # for every alpha vector
  #corresp=P[p].corresp[i] # that is resonant or under the given order
  #  if corresp!=0
      Avector[index,:]=P[p].Avector[i]
      index=index+1
  #  end    
end
end

index=1
for p in 1:info.max_order   # for every order
#write(ofile,"\nOrder "*string(p)*"\n")
for i in 1:P[p].nc # for every alpha vector
  #corresp=P[p].corresp[i] # that is resonant or under the given order
  #if corresp!=0
    for j in 1:info.nm
      fdyn[index,j]=2*real(P[p].fr[j,i])
    end  
    for j in 1:info.nm
      fdyn[index,j+info.nm]=-2*imag(P[p].fr[j,i])
    end
    index=index+1
  #end    
end
end

# mi segno le mappe
# mi sergno il numero di dofs
index=1
for p in 1:info.max_order   # for every order
  for i in 1:P[p].nc # for every alpha vector
    #corresp=P[p].corresp[i] # that is resonant or under the given order
    #if corresp!=0
      for inode=1:mesh.nn # loop nei nodi
          for idof=1:3 # lop nei dof
            dof=U.dof[(inode-1)*3+idof]# estraggo il dof
            if dof>0 # check se dof è positivo
              mappings[index,inode,idof]=real(P[p].Wr[dof,i])
              mappings_vel[index,inode,idof]=real(P[p].Wr[info.nK+dof,i])
            end
          end
      end
      for imode=1:info.neig # loop nei modi
        #println(size(ϕ[:,imode]))
        #println(size(M))
        #println(size(real.(P[p].Wr[:,i])))
        mappings_modal[index,imode]=ϕ[:,imode]'*M*real.(P[p].Wr[1:info.nK,i])
        mappings_modal_vel[index,imode]=ϕ[:,imode]'*M*real.(P[p].Wr[info.nK+1:2*info.nK,i])
      end
      index=index+1
    #end
  end      
end

return  mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn

end
