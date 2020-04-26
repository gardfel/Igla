import numpy as np
import math as kk
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random
import Tabl
import Aerodynamic


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

        l_cil = 45.5 / 8
        l_nos = 12.5 / 8
        tan_05 = 0.307
        l_kr = 1.46
        c_k = 0.03

        b_ak_kr = 0.094
        x_otn_op_kr = 0.846 / b_ak_kr

        M = self.v / 335
        Cyy_nos = [Aerodynamic.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil),
                                  Tabl.tab_3_2(M, 1.32, 22.43), Tabl.tab_3_4(M, 0, 22.43),
                                  Tabl.tab_3_4(M, 1, 23.75))]
        Cyy_kr = [Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_22(M, x_otn_op_kr)), l_kr, c_k, tan_05)]


        #while (abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5):
        while self.v < 887:
            alf = kk.atan((target.y - self.y) / (target.x - self.x))
            self.x += (self.v + self.v1) / 2 * kk.cos(alf) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(alf) * dt
            target.next_coord()
            M = self.v / 335

            Cy_nos = Aerodynamic.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil),
                                    Tabl.tab_3_2(M, 1.32, 22.43), Tabl.tab_3_4(M, 0, 22.43),
                                    Tabl.tab_3_4(M, 1, 23.75))
            Cyy_nos.append(Cy_nos)
            Cy_kr = Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_22(M, x_otn_op_kr)), l_kr, c_k, tan_05)
            Cyy_kr.append(Cy_kr)
            t += dt
            xx.append(self.x)
            yy.append(self.y)
            xt.append(target.x)
            yt.append(target.y)
            vv.append(self.v)
            tt.append(t)
            print(t)
            if self.v < 650:
                self.v = self.v1
                self.v1 += a * dt
            elif self.v < 887:
                self.v = self.v1
                self.v1 += a1 * dt
        plt.plot(tt, Cyy_nos)
        plt.axis([-0.1, t, 0, Cy_nos + 0.01])
        plt.grid(True)
        plt.show()
        plt.plot(tt, Cyy_kr)
        plt.axis([-0.1, t, 0, 0.06])
        plt.grid(True)
        plt.show()
        '''plt.plot(xx, yy)
        plt.plot(xt, yt)
        plt.show()
        plt.plot(tt, yy)
        plt.show()
        plt.plot(tt, xx)
        plt.show()
        plt.plot(tt, vv)
        plt.show()'''



dt = 10 ** -4
# (abs(target.y - self.y) > 25) or
vertel = Target()
igla = Rocket()

igla.navigation(vertel)
