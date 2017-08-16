"""
    Usage: python pplot_data.py input_file quantity_of_interest (0 = pH, 1 = calcite, 2 = dolomite)
    Example: python pass_fail.py time_concs.txt 0
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def exitWithError(error):
    print(error)
    sys.exit(1)

# validate input arguments
nArgsExpected = 3
if len(sys.argv) < nArgsExpected:
    exitWithError("Error: {} arguments given, expected {}\n".format(len(sys.argv), nArgsExpected) +
        "Usage: {} input_file quantity_of_interest (0 = pH, 1 = calcite, 2 = dolomite)".format(sys.argv[0]))

solutionFile1 = sys.argv[1]
channel = int(sys.argv[2])

if not os.path.isfile(solutionFile1):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile1))

f = open(solutionFile1, 'r')
content = f.read().split()
dims = np.asarray(content[0 : 3], dtype=np.int16)
data = np.asarray(content[3 : ], dtype=np.float32)
data = np.reshape(data, dims, order='F')
f.close()

x = np.linspace(0, 0.5, dims[0])
fig, ax = plt.subplots()
line, = ax.plot(x, data[:, channel, 0])

def animate(i):
    line.set_ydata(data[:, channel, i])  # update the data
    return line,

# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, dims[2]), init_func=init,
                              interval=10, blit=True)
plt.show()
