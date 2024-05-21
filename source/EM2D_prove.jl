
function  cf()

L=1
mE=154e3
h=.1
J=1/12*h^3
k=2*12*mE*J/L^3

g=.5
force=eps0/g
#force=eps0*L2*VDC^2/g^2
disp=force/k

println("displacement:  ",disp)


end



function  eigcc()

L=510
E=154e3
rho=2.33e-3
h=1.5
t=1
J=1/12*t*h^3
myeig=(4.73/L)^2*sqrt(E*J/(rho*h*t))

coeff=sqrt(rho*h*t*L^4/(E*J))

println(myeig)
println(coeff)
println(myeig*coeff)


end



function  disp()

    L=510
    E=154e3
    h=1.5
    b=100
    N=9e2
    stress=N/(h*b)
    val=stress*L/E
    println(val)
        
end
    