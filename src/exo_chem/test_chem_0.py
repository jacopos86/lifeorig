import easychem.easychem as ec

exo = ec.ExoAtmos()
exo.solve(1.0, 1500.0)

print("Solved OK")
print("Atoms:", exo.atoms)
print("Solved flag:", exo.solved)