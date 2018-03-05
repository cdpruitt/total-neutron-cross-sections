import numpy as np
from matplotlib import pyplot as plot
from matplotlib import animation

fig = plot.figure()
axis = plot.axes(xlim=(0,2), ylim=(-2,2))
line, = axis.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def phaseShift(x, y):
    if(abs(y)<1 and x>1):
        return 0.1
    else:
        return 0

def animate(i):
    y = np.linspace(-2,2,1000)
    x = np.linspace(0.02*i,0.02*i,1000)

    for j, k in enumerate(y):
        x[j] = x[j] - phaseShift(x[j], k)

    line.set_data(x,y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

anim.save('sinWave.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show()
