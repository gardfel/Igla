import numpy as np
import math as kk
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random
import Tabl


class Target(object):

    def __init__(self):
        self.v = random.uniform(200, 400)  # сумарная скорость цели
        self.ang = 0
        self.x = random.uniform(500, 4000)
        self.y = random.uniform(300, 2000)
        self.alf = random.uniform(0, 0.3925)

    def next_coord(self, *args):

        # self.alf += random.uniform(-0.0086, +0.0086)

        self.x += self.v * kk.cos(self.alf) * dt
        self.y += self.v * kk.sin(self.alf) * dt



class Rocket(object):

    def __init__(self):
        self.m = 11.29
        self.w0 = 2.296
        self.w1 = 1.696
        self.x = 0
        self.y = 0
        self.v = 51
        self.v1 = 51 + 78.9 * dt

    def navigation(self, *args):
        a = 249.7
        a1 = 29
        target = args[0]
        t = 0
        yt, xt, xx, yy, tt = [target.y], [target.x], [0], [0], [0]
        vv = [self.v]
        while (abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5):
            alf = kk.atan((target.y - self.y) / (target.x - self.x))
            self.x += (self.v + self.v1) / 2 * kk.cos(alf) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(alf) * dt
            target.next_coord()
            t += dt
            xx.append(self.x)
            yy.append(self.y)
            xt.append(target.x)
            yt.append(target.y)
            vv.append(self.v)
            tt.append(t)
            if self.v < 650:
                self.v = self.v1
                self.v1 += a * dt
            elif self.v < 887:
                self.v = self.v1
                self.v1 += a1 * dt
        plt.plot(xx, yy)
        plt.plot(xt, yt)
        plt.show()
        plt.plot(tt, yy)
        plt.show()
        plt.plot(tt, xx)
        plt.show()
        plt.plot(tt, vv)
        plt.show()



dt = 10 ** -5
# (abs(target.y - self.y) > 25) or
vertel = Target()
igla = Rocket()

igla.navigation(vertel)

Cy = Tabl.tab_3_2(650 / 335, 12.5 / 8, 45.5 / 8)

print(Cy)
