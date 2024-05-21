
file = "hardening.msh";

material = [[154e3, 0.0, 2.33e-3]];

solid = [[4, 1]]

dbc = [[1, 1], [1, 2], [2, 1], [2, 2]]
dbcval = [0.0, 0.0, 0.0, 0.0]
S0 = [6.0; 0.0; 0.0]

tbc = [3]

VDC = 1.2
VAC = 0.2

gap = 1.18
eps0 = 8.854e-6

info = Sinfo()
info.alpha = 3e-3
info.beta = 0.0

info.neig = 8
info.Lmm = [1]
info.Ffreq = 1
info.Fmult = 1
info.style = 'c'
info.max_order = 7
info.max_orderNA = 1
info.tol = 1e-1




