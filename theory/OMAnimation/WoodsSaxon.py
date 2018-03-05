import numpy as np
import math
import scipy.integrate as integrate
from matplotlib import pyplot as plot
from matplotlib import mlab
from matplotlib import colors
from matplotlib import animation

X_MIN = -40.
X_MAX = 40.
X_POINTS = 500

Y_MIN = -20.
Y_MAX = 25.

Y_POINTS = 200

WELL_DEPTH = 42.8 # in MeV

NEUTRON_MASS = 938.5 # in MeV/c^2
NEUTRON_ENERGY = 25 # in MeV
PLANCK_CONSTANT = 1239.842 # in MeV*(fm/c)

NUMBER_OF_WAVEFRONTS = 10
INITIAL_DISPLACEMENT = -10

NUCLEUS_MASS = 150
R_0 = 1.4 # nuclear radius constant, in fm

def wavelength(energy):
        return(PLANCK_CONSTANT*(2*NEUTRON_MASS*energy)**(-0.5)*((NUCLEUS_MASS+1)/NUCLEUS_MASS)) # in fm

def speed(energy):
        return((2*energy/NEUTRON_MASS)**0.5) # in units of c

def WoodsSaxon(x,y):
    U_0 = WELL_DEPTH
    R = R_0*NUCLEUS_MASS**(1/3.0)
    a = 0.5

    distance = math.sqrt(x**2 + y**2)

    U = U_0/(1+math.exp((distance-R)/a))

    return U

def indexOfRefraction(x, y, energy):
        return ((energy+WoodsSaxon(x,y))/energy)**(0.5) # (unitless)

class Wavefront:

    def __init__(self,
                 init_state = -3):
        self.init_state = init_state # in fm

        self.xValues = np.linspace(init_state, init_state, Y_POINTS)
        self.yValues = np.linspace(Y_MIN+5,Y_MAX-10,Y_POINTS)

    def position(self):
        return(self.xValues, self.yValues) # in fm

    def phaseShift_dt(self, x, y):
        return (speed(NEUTRON_ENERGY)/indexOfRefraction(x,y,
            NEUTRON_ENERGY))*(indexOfRefraction(x,y, NEUTRON_ENERGY)-1)\
            /(wavelength(NEUTRON_ENERGY)/(2*math.pi)) # in units of C

    def dstate_dt(self, xValues, t):
        dstate = np.zeros_like(xValues);
        for j, k in enumerate(self.yValues):
            dstate[j] = speed(NEUTRON_ENERGY)-self.phaseShift_dt(xValues[j], self.yValues[j])
        return dstate

time_elapsed = 0
phase_difference = 0

def step(waves, dt):

    global time_elapsed, phase_difference
    for wave in waves:
        wave.xValues = integrate.odeint(wave.dstate_dt, wave.xValues, [0,dt])[1]

    time_elapsed += dt
    phase_difference = waves[0].xValues[len(waves[0].xValues)/2]\
            - waves[0].xValues[0]

fig = plot.figure(figsize=(8*(X_MAX-X_MIN)/(Y_MAX-Y_MIN),8))
axes = plot.axes(xlim=(X_MIN,X_MAX), ylim=(Y_MIN,Y_MAX),)

#plot.rc('text',usetex=True)

xRange = np.linspace(X_MIN, X_MAX, X_POINTS)
yRange = np.linspace(Y_MIN, Y_MAX, Y_POINTS)

potentialGrid = [[WoodsSaxon(x,y) for x in xRange] for y in
        yRange]

cmap = colors.LinearSegmentedColormap.from_list('custom red', ['#FFFFFF', '#FF0000'], N=256)
norm = colors.Normalize(vmax=WELL_DEPTH, vmin=0)

potential = plot.contourf(xRange, yRange, potentialGrid, 30,
        norm=norm, cmap=cmap)

waves = []
lines = []

for i in range(NUMBER_OF_WAVEFRONTS):
    waves.append(Wavefront(INITIAL_DISPLACEMENT-i*wavelength(NEUTRON_ENERGY)))
    lines.append(axes.plot([], [], "b-", lw=2)[0])

time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
phase_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)

index_text = axes.text(0.60, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.60, 0.90, '', transform=axes.transAxes)
nucleus_text = axes.text(0.60, 0.85, '', transform=axes.transAxes)

dt = 0.08/speed(NEUTRON_ENERGY)

def init():
    for line in lines:
        line.set_data([], [])

    time_text.set_text('')
    phase_text.set_text('')

    index_text.set_text('Index of refraction at core = %.2f' % (indexOfRefraction(0,0,NEUTRON_ENERGY)))
    nucleus_text.set_text('Nucleus mass = %.0f' % NUCLEUS_MASS)
    energy_text.set_text('Neutron energy = %.1f [MeV]' % NEUTRON_ENERGY)

    initObjects = tuple(lines) + (time_text, index_text, phase_text,\
            nucleus_text, energy_text)
    return initObjects

def animate(i):
    global waves, dt

    step(waves, dt)
    for index, wave in enumerate(waves):
        lines[index].set_data(*wave.position())

    time_text.set_text('Time = %.1f [fm/C]' % time_elapsed)
    phase_text.set_text('Phase difference = %.2f [radians]' % phase_difference)

    animateObjects = tuple(lines) + (time_text, phase_text)
    return animateObjects

from time import time
t0 = time()
animate(0)
t1 = time()
interval = 1 * dt - (t1-t0)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=interval, repeat=False, blit=True)

#plot.show()

anim.save('WoodsSaxon.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
