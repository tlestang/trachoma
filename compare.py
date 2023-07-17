import numpy as np
import matplotlib.pyplot as plt


ref = np.loadtxt("reference.csv", delimiter=",")
validate = np.loadtxt("validate.csv", delimiter=",")

fig, ax = plt.subplots()

ax.plot(validate[:,0], linestyle='-', color='b')
ax.plot(validate[:,1], linestyle='-', color='g')

ax.plot(ref[:,0], linestyle='--', color='b')
ax.plot(ref[:,1], linestyle='--', color='g')

plt.show()
