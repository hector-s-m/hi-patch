import numpy as np

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3) # <-- Tune coarseness of SASA calc

num_atoms=cmd.count_atoms('all')
res=np.zeros(1)

for i in range(1, num_atoms+1):
	print(cmd.get_area('index %d' %i))

np.savetxt('areas_results.txt',res)
