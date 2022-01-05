import sys
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.pyplot import *

import numpy as np

# Read the file.
f2 = open(sys.argv[1], 'r')

# read the whole file into a single variable, which is a list of every row of the file.
lines = f2.readlines()
f2.close()

# initialize some variable to be lists:
x1 = []
y1 = []

# scan the rows of the file stored in lines, and put the values into some variables:
for line in lines:
    p = line.split()
    x1.append(float(p[1]))
    y1.append(float(p[2]))

max = max(x1)
xv = np.array(x1)
yv = np.array(y1)

xvv = []
yvv = []
previous=0
nb=0
for i in xv:
   nb=nb+1
   sdf = int(i)
   for j in range(previous+1,sdf):
      xvv.append(j);
      yvv.append(0)
   xvv.append(sdf);
   yvv.append(yv[nb-1]);
   previous = sdf;

xvvv = np.array(xvv)
yvvv = np.array(yvv)

outfile = sys.argv[1] + ".png";

# now, plot the data:
fig = plt.gcf()
plt.plot(xvvv, yvvv)
plt.fill_between(xvvv,yvvv,0,color='0.8')
xlabel('reference(bp)')
ylabel('coverage (reads/position)')
fig.savefig(outfile)
# plt.show()
