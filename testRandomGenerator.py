import RandomGenerator
import numpy as np

x = RandomGenerator.poisson(1000,0.5)
x = np.array(x)

print "Max=%i Min=%i" % (x.max(),x.min())

np.save("_gen/random_values.npy", x)