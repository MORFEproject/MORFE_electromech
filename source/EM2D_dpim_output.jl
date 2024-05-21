function output(P::Vector{Parametrisation})


# writes complex parametrization on file
ofile = open("./output/outCFULL.txt","w")

uneq=info.uneq

write(ofile,"\nParametrization f\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,string(P[p].Av[i])*"\n")
    for j in 1:info.nza
      write(ofile,string(P[p].f[j,i])*"\n")
    end  
  end
end

write(ofile,"\nParametrization W\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,"\n"*string(P[p].Av[i])*"\n")
    for j in uneq+1:2*uneq
      write(ofile,string(P[p].W[j,i])*"\n")
    end  
  end
end
close(ofile)               

# writes real parametrization on file
ofile = open("./output/outRFULL.txt","w")
write(ofile,"\nParametrization fr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,string(P[p].Av[i])*"\n")
    for j in 1:info.nmm
      write(ofile,string(2*real(P[p].fr[j,i]))*"\n")
    end  
    for j in 1:info.nmm
      write(ofile,string(2*imag(P[p].fr[j,i]))*"\n")
    end  
  end
end

write(ofile,"\nParametrization Wr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,"\n"*string(P[p].Av[i])*"\n")
    for j in uneq+1:2*uneq
      write(ofile,string(real(P[p].Wr[j,i]))*"\n")
end  
  end
end
close(ofile)               

write_rdyn(info,P)

end



function write_rdyn(info::Sinfo,P::Vector{Parametrisation})

  rdyn = ["" for i in 1:info.nza]
  for i = 1:info.nza
    rdyn[i] = "a"*string(i)*"' = "
  end
  # add unfolding parameter to first dofs identity-tangent to a modal displacement
  #nm = Int(ndofs/2)
  #for i = 1:nm
  #  rdyn[nm+i] *= "+mu*z"*string(nm+i)
  #end

  for p in 1:info.max_order
    for c = 1:P[p].m
      Av=P[p].Av[c]   
      monomial = ""
      for d = 1:info.nrom        
        if (Av[d]!=0)
          monomial *= "*a"*string(d)*"^"*string(Av[d])
        end
      end
      for j in 1:info.nmm
        rcoeff=2*real(P[p].fr[j,c])
#        
#        icoeff=2*imag(P[p].fr[j,c]) # ORIGINAL
        icoeff=-2*imag(P[p].fr[j,c])
#
        if abs(rcoeff)>1e-20
          rdyn[j] *= " + "*string(rcoeff)*monomial
        end
        if abs(icoeff)>1e-20
          rdyn[j+info.nmm] *= " + "*string(icoeff)*monomial
        end
      end
    end
  end

  ofile = open("./output/equations.txt","w")
  for i = 1:info.nza
    write(ofile,rdyn[i]*";\n")
  end
  close(ofile)  

end

