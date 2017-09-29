import numpy as np
from scipy.io import FortranFile
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

n = 11
time = np.zeros(n)

f = FortranFile('dolomite.con', 'r')

for i in range(0, n):
    if i == 0:
        d1 = np.array(f.read_ints(dtype='3int32, 3float32'))
        mat = np.zeros((n, d1[0][0][0]))
    d = np.array(f.read_reals(dtype='float32, int32'))
    # print (d, np.size(d))
    time[i] = d[0][0]
    loc = np.array(f.read_ints(dtype='int32'))
    # print (loc, np.size(d))
    mat[i, loc - 1] = np.array(f.read_reals(dtype='float32'))
    # print (mat[i, loc - 1], np.size(d))

f.close()
x = np.linspace(0, 0.5, d1[0][0][0])

# this is for static plotting
# plt.figure(1)
# plt.plot(x, mat[n - 1, ])
# plt.show()


# this is for animation
fig, ax = plt.subplots()
line, = ax.plot(x, mat[0, ])

def animate(i):
    line.set_ydata(mat[i, ])  # update the data
    return line,

# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    ax.set_ylim(0, np.max(mat[1, ]))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, n), init_func=init,
                              interval=200, repeat=False, blit=True)
plt.show()
