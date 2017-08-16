"""
    This program plots the output vectors (vs. time) of a given solution attribute
    for two OpenFAST solutions, with the second solution assumed to be the baseline for
    comparison. It reads two OpenFAST binary output files (.outb), and
    generates three plots of the given attribute (1) comparing the two tests' respective
    values, (2) showing the difference in values, and (3) showing relative difference,
    as compared to the baseline solution.

    Usage: python3 pass_fail.py solution1
    Example: python3 pass_fail.py output-local/Test01.outb output-baseline/Test01.outb Wind1VelX
"""
import sys, os
import numpy as np
# from numpy import linalg as LA
# from fast_io import load_output
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def exitWithError(error):
    print(error)
    sys.exit(1)

# validate input arguments
nArgsExpected = 2
if len(sys.argv) < nArgsExpected:
    exitWithError("Error: {} arguments given, expected {}\n".format(len(sys.argv), nArgsExpected) +
        "Usage: {} solution1 solution2 attribute".format(sys.argv[0]))

solutionFile1 = sys.argv[1]

if not os.path.isfile(solutionFile1):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile1))

f = open(solutionFile1, 'r')
content = f.read().split()
dims = np.asarray(content[0 : 3], dtype=np.int16)
data = np.asarray(content[3 : ], dtype=np.float32)
data = np.reshape(data, dims, order='F')
f.close()

x = np.linspace(0, 0.5, dims[0])

# plt.figure(1)
# plt.subplot(311)
# plt.title('File Comparisons for ' + testname + '\nNew: ' +
#           solutionFile1 + '\n Old: ' + solutionFile2)
# plt.grid(True)
# plt.ylabel(attribute + '\n' + '(' + info1['attribute_units'][channel] + ')')
# plt.plot(x, data[:, 0, -1])
# plt.plot(timevec, dict2[:, channel], label = 'Old')
# plt.legend()
# plt.subplot(312)
# plt.grid(True)
# plt.plot(timevec, diff)
# plt.ylabel('Difference\n' + '(Old - New)\n' + '(' + info1['attribute_units'][channel] + ')')
# plt.subplot(313)
# plt.grid(True)
# plt.plot(timevec, reldiff)
# plt.ylabel('Relative\n Difference\n' + '(%)')
# plt.xlabel('Times (s)')
# plt.show()

fig, ax = plt.subplots()

# x = np.arange(0, 2*np.pi, 0.01)
line, = ax.plot(x, data[:, 1, 0])


def animate(i):
    line.set_ydata(data[:, 1, i])  # update the data
    return line,


# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, dims[2]), init_func=init,
                              interval=10, blit=True)
plt.show()
