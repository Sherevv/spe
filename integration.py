#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy import integrate
from matplotlib.pylab import *


class helpers:
    def F_1(self, e):
        otv = 1 + 3 * e ** 2 + 3 * (e ** 4) / 8
        return otv

    def F_2(self, e):
        otv = 1 + (15 / 2) * e ** 2 + (45 / 8) * e ** 4 + (5 / 16) * e ** 6
        return otv

    def F_3(self, e):
        otv = 1 + (31 / 2) * e ** 2 + (255 / 8) * e ** 4 + (185 / 16) * e ** 6 + (25 / 64) * e ** 8
        return otv

    def F_4(self, e):
        otv = 9 + (135 / 4) * e ** 2 + (135 / 8) * e ** 4 + (45 / 64) * e ** 6
        return otv

    def F_5(self, e):
        otv = (11 / 2) + (33 / 4) * e ** 2 + (11 / 16) * e ** 4
        return otv


def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]


class SolveEvolution:

    def __init__(self, ecc, semi_a, incl, mp, r0, Tp, ms, time):
        self.semi_a = semi_a
        self.ecc = ecc
        self.time = 10**time

        rez = self.getParams(mp, ms, r0, semi_a, ecc, incl, Tp)
        self.p = rez[0]
        self.n0 = rez[1]
        self.I = rez[2]
        self.n = rez[3]


    def sysODE(self, t, x):
        """System ODE function"""
        # Assign some variables for convenience of notation
        n0 = x[0]
        e = x[1]
        i = x[2]
        a = x[3]

        # Algebraic equations
        koef = math.cos(i) - (self.p / (n0 ** (1.0 / 3))) * math.sqrt(1 - e ** 2)

        # Output from ODE function must be a COLUMN vector, with n rows
        dydt = np.zeros((len(x), 1))
        h = helpers()

        dn0 = -3 * (n0 ** (16.0 / 3) / ((1 - e ** 2) ** (15.0 / 2))) * (
        koef * ((1 - e ** 2) ** 1.5) * h.F_2(e) - n0 * h.F_3(e))
        de = (e * n0 ** (13.0 / 3) / ((1 - e ** 2) ** (13.0 / 2))) * (
        koef * ((1 - e ** 2) ** 1.5) * h.F_5(e) - n0 * h.F_4(e))
        di = -(n0 ** (13.0 / 3) * math.sin(i) / (2 * (1 - e ** 2) ** 5)) * h.F_1(e)
        da = -self.G0A * (2.0 / 3) * dn0 * (a ** (5.0 / 2)) / (self.f0 ** (1.0 / 2))

        # Изменяем масштаб времени
        dydt[0] = dn0 * self.time
        dydt[1] = de * self.time
        dydt[2] = di * self.time
        dydt[3] = da * self.time

        return dydt

    def getParams(self, m, mu, r, a, e, I, Tp):
        """Вычисление параметров p,n0,I для систем планета-спутник
        m - масса планеты,
        mu - масса спутника,
        r - радиус планеты,
        e - эксцентриситет,
        I - наклонение,
        Tp - период вращения планеты,
        Ts- период обращения спутника
        """


        I = I * math.pi / 180.0  # Наклонение в радианах
        w = 2.0 * math.pi / (Tp * 3600)  # угловая скорость вращения планеты
        # n=2*pi/(Ts*86400) #среднее движение спутника

        f = 6.67 * 10 ** (-11)  # Гравитационная постоянная
        f0 = f * (m + mu)
        self.f0 = f0
        n = math.sqrt(f0) / a ** (3.0 / 2)
        mr = m * (mu + 0.0) / (m + mu)
        A = (2.0 / 5) * m * r ** 2
        G0 = A * w + + mr * f0 ** (2.0 / 3) * math.sqrt(1 - e ** 2) * math.cos(I) / (n ** (1.0 / 3))
        p = ((A ** (1.0 / 3)) * (f0 ** (2.0 / 3)) * mr) / (G0 ** (4.0 / 3))
        n0 = n * A / G0
        self.G0A = G0 / A
        rez = [0]*4#np.zeros((4, 1))
        rez[0] = p
        rez[1] = n0
        rez[2] = I
        rez[3] = n
        return rez
        # w =    7.2922e-05 +
        # n =    2.6617e-06 +-
        # f0 =    4.0334e+14
        # mr =    7.2597e+22
        # A =    9.7200e+37
        # G0 =    3.5526e+34
        # G0A =    3.6550e-04
        # p =  0.15602
        # n0 =  0.0072824
        # I =  0.089884


    def integrate(self):

        # Set the time range
        t_start = 0.0
        t_final = 1.0
        delta_t = 0.001
        time = linspace(0.0, 1.0, 1000)

        # Number of time steps: 1 extra for initial condition
        num_steps = np.floor((t_final - t_start) / delta_t) + 1

        # Additional Python step: create vectors to store trajectories
        t = [0]*num_steps#np.zeros((num_steps, ))
        n = np.zeros((num_steps, ))
        n0 = np.zeros((num_steps, ))
        ecc = np.zeros((num_steps, ))
        incl = np.zeros((num_steps, ))
        semi_a = np.zeros((num_steps, ))

        # Set initial condition(s): for integrating variable and time!
        t[0] = t_start
        n[0] = self.n
        n0[0] = self.n0
        ecc[0] = self.ecc
        incl[0] = self.I
        semi_a[0] = self.semi_a


        # Start by specifying the integrator:
        # use ``vode`` with "backward differentiation formula"
        #z = integrate.odeint(self.sysODE, [n0[0], ecc[0], incl[0], semi_a[0]], time)
        r = integrate.ode(self.sysODE).set_integrator('vode', method='bdf')
        r.set_initial_value([n0[0], ecc[0], incl[0], semi_a[0]], t_start)



        #f0 = 4.033409030000001e+14
        #G0A = 3.654960197509502e-04
        #a2 = np.zeros((num_steps, 1))
        #r0 = 6.378 * 10 ** 6
        #a2[0] = (f0 ** (1 / 3)) / (r0 * (G0A * n0[0]) ** (2 / 3))

        # Integrate the ODE(s) across each delta_t timestep
        k = 1
        while r.successful() and k < num_steps:
            r.integrate(r.t + delta_t)

            # Store the results to plot later
            t[k] = r.t
            #t[k] = t[k] / 3.1556926e+7
            n0[k] = r.y[0]
            ecc[k] = r.y[1]
            incl[k] = r.y[2]
            semi_a[k] = r.y[3]
            n[k] = math.sqrt(self.f0) / (semi_a[k] ** (3.0 / 2))
            #a2[k] = (f0 ** (1 / 3)) / (r0 * (G0A * n0[k]) ** (2 / 3))
            k += 1

        ind = np.where(n0 > 0.8)
        rez = [0]*6#np.zeros(shape=(6, num_steps))
        rez[0] = t
        rez[1] = n0
        rez[2] = ecc
        rez[3] = incl
        rez[4] = semi_a
        rez[5] = n


        # if len(ind) > 0:
        #     n0 = [value for (i, value) in enumerate(n0) if i not in set(ind[0])]
        #     ecc = [value for (i, value) in enumerate(ecc) if i not in set(ind[0])]
        #     semi_a = [value for (i, value) in enumerate(semi_a) if i not in set(ind[0])]
        #     incl = [value for (i, value) in enumerate(incl) if i not in set(ind[0])]
        #     t = [value for (i, value) in enumerate(t) if i not in set(ind[0])]

        #
        # # All done!  Plot the trajectories in two separate plots:
        # #fig = figure()
        # # ax1 = subplot(211)
        # # ax1.plot(t, n0)
        # # ax1.set_xlim(t_start, t_final)
        # # ax1.set_xlabel('Time [minutes]')
        # # ax1.set_ylabel('n0')
        # # ax1.grid('on')
        #fig = figure()
        #ax2 = plt.subplot(212)


        # plt.plot(time, z[:,0], time, z[:,1]) # y[:,

        # plt.plot(extr, n0, 'r')
        #plt.plot(t, ecc, 'g')
        #plt.show()
        return rez




