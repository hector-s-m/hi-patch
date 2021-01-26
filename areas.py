import numpy as np

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 4)

num_atoms=cmd.count_atoms('all')
res=np.zeros(num_atoms)

for i in range(1, num_atoms+1):
	print(cmd.get_area('index %d' %i))
