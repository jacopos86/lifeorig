from lifeorig.gillespie_algo import chemical_kinetics_solver

CKS = chemical_kinetics_solver()
initial_state = [290, 10, 0]
#
express = []
express.append("r")
express.append("0.01 * i")
#
stoichiometry = []
stoichiometry.append([-1, 1, 0])
stoichiometry.append([ 0,-1, 1])
CKS.test(initial_state, express, stoichiometry, "simul-1.png")
#
express = []
express.append("r")
express.append("0.5 * i")
CKS.test(initial_state, express, stoichiometry, "simul-2.png")
#
express = []
express.append("s * i")
express.append("0.5 * i")
CKS.test(initial_state, express, stoichiometry, "simul-3.png")
#
express = []
express.append("s * i / 100")
express.append("0.5 * i")
CKS.test(initial_state, express, stoichiometry, "simul-4.png")