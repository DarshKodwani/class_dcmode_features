import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/mnt/users/dkodwani/Decaying_mode/class_dcmode/outputs_decaycl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['outputs_decaycl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], curve[:, 1])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()