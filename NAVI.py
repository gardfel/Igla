import math as kk
import matplotlib.pyplot as plt
# import plotly.graph_objs as go
# import numpy as np
import random
import Tabl
import Aero


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
        yt, xt, xx, yy, tt = [], [], [], [], []
        vv = []

        cyy_nos = []
        cyy_kr = []
        cyy_op = []
        cyy1_alf = []
        cxx_tr = []
        cxx_nos = []
        cxx_korm = []

        # while (abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5):

        while self.v < 887:

            mach = self.v / Tabl.tab_atm(self.y, 2)

            cy_nos = Aero.Cy_nos(Tabl.tab_3_2(mach, l_nos, l_cil), Tabl.tab_3_4(mach, 0, l_cil),
                                 Tabl.tab_3_2(mach, 1.32, 22.43), Tabl.tab_3_4(mach, 0, 22.43),
                                 Tabl.tab_3_4(mach, 1, 23.75))
            cyy_nos.append(cy_nos)
            cy_kr = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_22(mach, x_otn_op_kr)), l_kr, c_k, tan_05_kr)
            cyy_kr.append(cy_kr)
            cy_op = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos)), l_op, c_op, tan_05_op)
            cyy_op.append(cy_op)
            cy1_alf_f = cy_nos

            k_aa_kr = 1 + 3 * D_kr - (D_kr * (1 - D_kr)) / nu_k_kr
            eps_sr_alf = 0  # для "утки"
            cy1_alf_kr = (cy_kr * k_aa_kr) * (1 - eps_sr_alf)

            k_aa_op = 1 + 3 * D_op - (D_op * (1 - D_op)) / nu_k_op
            cy1_alf_op = cy_op * k_aa_op
            k_t_op = Tabl.tab_3_21(mach, l_nos)
            k_t_kr = Tabl.tab_3_22(mach, x_otn_op_kr)

            cy1_alf = cy1_alf_f * S_f + cy1_alf_kr * S_kr * k_t_kr + cy1_alf_op * S_op * k_t_op
            cyy1_alf.append(cy1_alf)

            ni_atm = Tabl.tab_atm(self.y, 5)
            re_f = self.v * L_f / ni_atm
            # print(re_f * h / L_f)

            re_t = 15 * 10 ** 6  # test (max for current speed and form)

            # x_t = re_t * ni_atm / self.v
            # print(x_t)
            x_t = 0.08  # координата точки перехода (до точного определения)

            cx_tr = Tabl.tab_4_2(re_f, x_t) / 2 * Ff / Sf
            cxx_tr.append(cx_tr)

            cx_con = Tabl.tab_4_11(mach, l_nos_)
            cx_zat = Tabl.tab_4_13(mach, l_zat)
            cx_nos = cx_con * (1 - r_ ** 2 * kk.cos(teta) ** 2 * (3.1 - 1.4 * r_ * kk.cos(teta) - 0.7 * r_ ** 2 *
                                                                  kk.cos(teta) ** 2)) + cx_zat * r_ ** 2
            cxx_nos.append(cx_nos)

            cx_korm = Tabl.tab_4_24(mach, nu_kor, l_korm)
            cxx_korm.append(cx_korm)
            """
            re_k = self.v * l_kr / ni_atm
            re_k_t =
            x_t_kr = re_k_t / re_k
            c_f_kr = Tabl.tab_4_2(re_k, x_t_kr)
            ni_c_kr = Tabl.tab_4_28()

            re_op = self.v * l_op / ni_atm
            re_op_t = 
            x_t_op = re_op_t / re_op
            c_f_op = Tabl.tab_4_2(re_op, x_t_op)
            ni_c_op =

            cx_k_pr = 2 * c_f_kr * ni_c_kr
            cx_op_pr = c_f_op * ni_c_op
            # cx_nes ="""

            # cx_0 = 1.05 * (cx_o_f * S_f + cx_0_op * k_t_op * S_op + cx_0_kr * k_t_kr * S_kr)

            t += dt
            alf = kk.atan((target.y - self.y) / (target.x - self.x))
            self.x += (self.v + self.v1) / 2 * kk.cos(alf) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(alf) * dt
            target.next_coord()
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
        """plt.plot(tt, cyy_nos)
        plt.axis([-0.1, t, 0, 0.06])
        plt.grid(True)
        plt.show()
        plt.plot(tt, cyy_kr)
        plt.axis([-0.1, t, 0, 0.06])
        plt.grid(True)
        plt.show()
        plt.plot(tt, cyy_op)
        plt.grid(True)
        # plt.axis([-0.1, t, -0.1, 0.1])
        plt.show()
        plt.plot(tt, cyy1_alf)
        plt.grid(True)
        # plt.axis([-0.1, t, 0, 0.4])
        plt.show()"""

        plt.plot(tt, cxx_korm)
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
d = 0.072  # диаметр миделевого сечения

l_cil = 45.5 / 8
l_nos = 12.5 / 8
tan_05_kr = 0.307  # тангенс среднего угла крыла
l_kr = 1.46  # относительное удлинение крыла
c_k = 0.03  # относительная толщина профиля крыла

tan_05_op = -0.012  # тангенс среднего угла оперения
l_op = 7.888  # относительное удлинение оперения l^2/S
c_op = 0.03  # относительная толщина профиля крыла

b_ak_kr = 0.094
x_otn_op_kr = 0.846 / b_ak_kr  # относительное расстояние между оперением и средней хордой крыльев
S_f = 0.184 / (kk.pi * (d ** 2) / 4)  # относительная площадь корпуса
S_op = 0.0078575 / (kk.pi * (d ** 2) / 4)  # относительная площадь передних несущих поверхностей
S_kr = 0.015 / (kk.pi * (d ** 2) / 4)  # относительная площадь задних несущих поверхностей

D_kr = 0.072 / 0.146  # отерсительный диаметр корпуса
nu_k_kr = 1.322  # относительное сужение консоли задней несущей поверхности
D_op = 0.072 / 0.25  # относительный диаметр корпуса
nu_k_op = 1.069  # относительное сужение консоли передней несущей поверхности

L_f = 1.626  # длина корпуса
Ff = 0.3619  # площадь обтекаемой потоком поверхности корпуса (без донного среза)
Sf = 0.184  #

l_zat = 23.49 / 62.59  # относительное удлинение затупления носовой части
r_ = 2 * 0.0329 / d

l_nos_ = (l_nos - r_ / 2) / kk.sqrt(1 - r_)  # относительное удлинение носовой части без затупления
print(l_nos_, l_nos)
teta = kk.atan((1 - r_) / (l_nos - r_ / 2))  # угол наклона образующей носовой части (конуса)
print(teta * 180 / kk.pi)

h = 6.3 * 10 ** -6  # Примерная высота бугоров на поверхности корпуса (в зависимости от класса чистоты) (для 7-го)

nu_kor = 0.047 / d  # относительное сужение кормовой части
l_korm = 0.0463 / d  # относительное удлинение кормовой части

igla.navigation(vertel)
