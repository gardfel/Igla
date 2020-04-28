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


        M = self.v / Tabl.tab_atm(self.y, 2)

        Cyy_nos = [Aerodynamic.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil),
                                  Tabl.tab_3_2(M, 1.32, 22.43), Tabl.tab_3_4(M, 0, 22.43),
                                  Tabl.tab_3_4(M, 1, 23.75))]
        Cyy_kr = [Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_22(M, x_otn_op_kr)), l_kr, c_k, tan_05_kr)]
        Cyy_op = [Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_21(M, l_nos)), l_op, c_op, tan_05_op)]
        Cyy1_alf = [0]

        #while (abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5):
        while self.v < 887:
            alf = kk.atan((target.y - self.y) / (target.x - self.x))
            self.x += (self.v + self.v1) / 2 * kk.cos(alf) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(alf) * dt
            target.next_coord()
            M = self.v / Tabl.tab_atm(self.y, 2)

            Cy_nos = Aerodynamic.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil),
                                    Tabl.tab_3_2(M, 1.32, 22.43), Tabl.tab_3_4(M, 0, 22.43),
                                    Tabl.tab_3_4(M, 1, 23.75))
            Cyy_nos.append(Cy_nos)
            Cy_kr = Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_22(M, x_otn_op_kr)), l_kr, c_k, tan_05_kr)
            Cyy_kr.append(Cy_kr)
            Cy_op = Tabl.tab_3_5(M * kk.sqrt(Tabl.tab_3_21(M, l_nos)), l_op, c_op, tan_05_op)
            Cyy_op.append(Cy_op)
            Cy1_alf_f = Cy_nos

            K_aa_kr = 1 + 3 * D_kr - (D_kr * (1 - D_kr)) / nu_k_kr
            eps_sr_alf = 0 # для "утки"
            Cy1_alf_kr = (Cy_kr * K_aa_kr) * (1 - eps_sr_alf)

            K_aa_op = 1 + 3 * D_op - (D_op * (1 - D_op)) / nu_k_op
            Cy1_alf_op = Cy_op * K_aa_op
            k_t_op = Tabl.tab_3_21(M, l_nos)
            k_t_kr = Tabl.tab_3_22(M, x_otn_op_kr)

            Cy1_alf = Cy1_alf_f * S_f + Cy1_alf_kr * S_kr * k_t_kr + Cy1_alf_op * S_op * k_t_op
            Cyy1_alf.append(Cy1_alf)

            ni_atm = Tabl.tab_atm(self.y, 5)
            Re_f = self.v * L_f / ni_atm

            Cx_tr = Tabl.tab_4_2(Re_f) / 2 * F_f / S_f
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
        plt.plot(tt, Cyy_op)
        plt.grid(True)
        # plt.axis([-0.1, t, -0.1, 0.1])
        plt.show()
        plt.plot(tt, Cyy1_alf)
        plt.grid(True)
        #plt.axis([-0.1, t, 0, 0.4])
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
d = 0.072 # диаметр миделевого сечения

l_cil = 45.5 / 8
l_nos = 12.5 / 8
tan_05_kr = 0.307 # тангенс среднего угла крыла
l_kr = 1.46 # относительное удлинение крыла
c_k = 0.03 # относительная толщина профиля крыла

tan_05_op = -0.012 # тангенс среднего угла оперения
l_op = 7.888 # относительное удлинение оперения l^2/S
c_op = 0.03 # относительная толщина профиля крыла

b_ak_kr = 0.094
x_otn_op_kr = 0.846 / b_ak_kr # относительное расстояние между оперением и средней хордой крыльев
S_f = 0.184 / (kk.pi * (d ** 2) / 4) # относительная площадь корпуса
S_op = 0.0078575 / (kk.pi * (d ** 2) / 4) # относительная площадь передних несущих поверхностей
S_kr = 0.015 / (kk.pi * (d ** 2) / 4) # относительная площадь задних несущих поверхностей

D_kr = 0.072 / 0.146 # отерсительный диаметр корпуса
nu_k_kr = 1.322 # относительное сужение консоли задней несущей поверхности
D_op = 0.072 / 0.25 # относительный диаметр корпуса
nu_k_op = 1.069 # относительное сужение консоли передней несущей поверхности

L_f = 1.626 # длина корпуса
Ff = 0.3619  # площадь обтекаемой потоком поверхности корпуса (без донного среза)
Sf = 0.184 #



igla.navigation(vertel)
