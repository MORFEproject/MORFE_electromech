function globalprocedure()


fname = "./input/"*input_file
include(fname)
info.mesh_file=file

T6,B3,nodes = readgmsh!()

##############################
# static solution and matrix assemblage
##############################

C,K,M,F0=analysis(nodes,T6,B3)
#cf()

##############################
# eigenvalue computation
##############################

D,VR,VL=eigAB(K,M,C,nodes)

##############################
# launch DPIM
##############################

P=dpim(K,M,C,F0,D,VR,VL,nodes,T6,B3) # computes parametrization
realification!(P)  # performs realification

##################################################
#  OUTPUT ON FILE
##################################################

output(P)  # output

##################################################
#  output for matcont
##################################################

matcont(nodes,M,VR,P)

end