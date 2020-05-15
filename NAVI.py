import math as kk
import matplotlib.pyplot as plt
# import plotly.graph_objs as go
# import numpy as np
import random
import Tabl
import Aero


class Target(object):

    def __init__(self):
        self.v = random.uniform(200, 320)  # сумарная скорость цели
        self.ang = 0
        self.x = random.uniform(2900, 3000)
        self.y = random.uniform(1400, 1500)
        self.alf = 0.3925  # random.uniform(0, 0.3925)

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
        self.v_ = 0
        self.v1 = 51 + 78.9 * dt
        self.alf = 0
        self.alf_potr = 0
        self.delt = 0
        self.delt_potr = 0
        self.omega = 0
        self.d_omega = 0
        self.fi_viz = 0
        self.d_fi_viz = 0
        self.p1 = 2768  # сила тяги на первом режиме
        self.p2 = 900  # сила тяги на втором режиме
        self.p3 = 0

    def navigation(self, *args):

        a = 249.7
        a1 = 29
        target = args[0]
        t = 0
        yt, xt, xx, yy, tt = [], [], [], [], []
        vv = []
        i = 0
        v_sr = 0

        cyy_nos = []
        cyy_kr = []
        cyy_op = []
        cyy1_alf = []
        cyy1_delt_op = []

        cxx_tr = []
        cxx_nos = []
        cxx_korm = []
        cxx_kr_pr = []
        cxx_op_pr = []
        cxx_kr_vol = []
        cxx_op_vol = []
        cxx_0f = []
        cxx_0_kr = []
        cxx_0_op = []
        cxx_0 = []
        f_xx = []
        f_yy = []
        cxx_zat = []
        cxx_con = []
        x_ffa_nos_cill = []
        x_ffa_f = []

        # while (abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5):

        self.omega = kk.atan((target.y - self.y) / (target.x - self.x))

        # while self.v < 887:
        while ((abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5)) and (t < 16):

            par_1 = self.fi_viz
            self.fi_viz = kk.atan((target.y - self.y) / (target.x - self.x))
            self.d_fi_viz = (self.fi_viz - par_1) / dt
            self.d_omega = a_m * self.d_fi_viz

            n_y_a = self.v * self.d_omega / g + kk.cos(self.omega)

            # alf =

            # self.delt = kk.asin()
            mach = self.v / Tabl.tab_atm(self.y, 2)

            cy_nos = Aero.Cy_nos(Tabl.tab_3_2(mach, l_nos, l_cil), Tabl.tab_3_4(mach, 0, l_cil),
                                 Tabl.tab_3_2(mach, 1.32, 22.43), Tabl.tab_3_4(mach, 0, 22.43),
                                 Tabl.tab_3_4(mach, 1, 23.75))
            cyy_nos.append(cy_nos)
            cy_korm = -0.2 * 2 / 57.3 * (1 - n_korm ** 2)
            cy_kr = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_22(mach, x_otn_op_kr)), l_kr, c_kr, tan_05_kr)
            cyy_kr.append(cy_kr)
            cy_op = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos)), l_op, c_op, tan_05_op)
            cyy_op.append(cy_op)
            cy1_alf_f = cy_nos + cy_korm

            K_aa_kr = 1 + 3 * D_kr - (D_kr * (1 - D_kr)) / nu_k_kr
            k_aa_kr = (1 + 0.41 * D_kr) ** 2 * ((1 + 3 * D_kr - 1 / nu_k_kr * D_kr * (1 - D_kr)) / (1 + D_kr) ** 2)
            eps_sr_alf = 0  # для "утки"
            cy1_alf_kr = (cy_kr * K_aa_kr) * (1 - eps_sr_alf)

            K_aa_op = 1 + 3 * D_op - (D_op * (1 - D_op)) / nu_k_op
            k_aa_op = (1 + 0.41 * D_op) ** 2 * ((1 + 3 * D_op - 1 / nu_k_op * D_op * (1 - D_op)) / (1 + D_op) ** 2)
            cy1_alf_op = cy_op * K_aa_op
            k_t_op = Tabl.tab_3_21(mach, l_nos)
            k_t_kr = Tabl.tab_3_22(mach, x_otn_op_kr)

            cy1_alf = cy1_alf_f * S_f + cy1_alf_kr * S_kr * k_t_kr + cy1_alf_op * S_op * k_t_op
            cyy1_alf.append(cy1_alf)

            K_delt_0 = k_aa_op

            if mach <= 1.4:
                k_sh = 0.85
            else:
                k_sh = 0.99
            n_ef = k_sh * kk.cos(0)
            cy1_delt_op = cy1_alf_op * K_delt_0 * n_ef
            cyy1_delt_op.append(cy1_delt_op)

            ni_atm = Tabl.tab_atm(self.y, 5)
            re_f = self.v * L_f / ni_atm
            # print(re_f * h / L_f)

            re_t = Tabl.tab_4_5(mach, re_f, 7, L_f)  # test (max for current speed and form)

            x_tt = re_t * ni_atm / self.v
            if x_tt >= (0.26 / L_f):
                x_tt = f_t / Ff
            # print(x_t)
            # x_t = 0.08  # координата точки перехода (до точного определения)

            # cx_tr = Tabl.tab_4_2(re_f, x_t) / 2 * Ff / Sf
            # cxx_tr.append(cx_tr)
            cx_tr = Tabl.tab_4_2(re_f, x_tt) / 2 * Ff / Sf
            cxx_tr.append(cx_tr)

            cx_nos = Tabl.tab_4_11(mach, l_nos_)
            cxx_con.append(cx_nos)
            """cx_zat = Tabl.tab_4_13(mach, l_zat)
            cxx_zat.append(cx_zat)
            cx_nos = cx_con * (1 - r_ ** 2 * kk.cos(teta) ** 2 * (3.1 - 1.4 * r_ * kk.cos(teta) - 0.7 * r_ ** 2 *
                                                                  kk.cos(teta) ** 2)) + cx_zat * r_ ** 2"""
            cxx_nos.append(cx_nos)

            cx_korm = Tabl.tab_4_24(mach, nu_kor, l_korm)
            cxx_korm.append(cx_korm)
            cx_0f = cx_tr + cx_nos + cx_korm
            cxx_0f.append(cx_0f)

            # профилное сопротивление несущих поверхностей
            re_k = self.v * b_kr / ni_atm
            re_k_t = Tabl.tab_4_5(mach, re_k, 5, b_kr)
            x_t_kr = re_k_t / re_k
            c_f_kr = Tabl.tab_4_2(re_k, x_t_kr)
            ni_c_kr = Tabl.tab_4_28(x_t_kr, c_kr)

            re_op = self.v * b_op / ni_atm
            re_op_t = Tabl.tab_4_5(mach, re_op, 5, b_op)
            x_t_op = re_op_t / re_op
            c_f_op = Tabl.tab_4_2(re_op, x_t_op)
            ni_c_op = Tabl.tab_4_28(x_t_op, c_op)

            cx_kr_pr = c_f_kr * ni_c_kr
            cxx_kr_pr.append(cx_kr_pr)
            cx_op_pr = c_f_op * ni_c_op
            cxx_op_pr.append(cx_op_pr)
            # волновое сопротивление несущих поверхностей
            if mach < 1.1:
                cx_op_vol = Tabl.tab_4_30(mach, nu_k_op, l_op, tan_05_op, c_op)
                cx_kr_vol = Tabl.tab_4_30(mach, nu_k_kr, l_kr, tan_05_kr, c_kr)
                cxx_kr_vol.append(cx_kr_vol)
                cxx_op_vol.append(cx_op_vol)
            elif mach >= 1.1:

                cx_kr_vol = (Tabl.tab_4_30(mach, nu_k_kr, l_kr, tan_05_kr, c_kr)) * \
                            (1 + Tabl.tab_4_32(mach, tan_05_kr) * (koef_kr - 1))
                cxx_kr_vol.append(cx_kr_vol)
                cx_op_vol = (Tabl.tab_4_30(mach, nu_k_op, l_op, tan_05_op, c_op)) * \
                            (1 + Tabl.tab_4_32(mach, tan_05_op) * (koef_op - 1))
                cxx_op_vol.append(cx_op_vol)
            # cx_nes =
            cx_0_op = cx_op_pr + cx_op_vol
            cx_0_kr = cx_kr_pr + cx_kr_vol
            cxx_0_kr.append(cx_0_kr)
            cxx_0_op.append(cx_0_op)

            cx_0 = 1.05 * (cx_0f * S_f + cx_0_op * k_t_op * S_op + cx_0_kr * k_t_kr * S_kr)
            cxx_0.append(cx_0)
            f_x = cx_0 * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_xx.append(f_x)
            f_y = cy1_alf * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_yy.append(f_y)
            # индуктивное сопротивление
            delt_cx1 = 2 * Tabl.tab_4_40(mach, l_nos, 1) * kk.sin(self.alf) ** 2
            cx_f_ind = cy1_alf_f * kk.sin(self.alf) + delt_cx1 * kk.cos(self.alf)
            # cx_op_ind = cy

            # моменты
            # фокус носовой части с цилиндром
            x_fa_nos_cill = L_nos - W_nos / S_ff + Tabl.tab_5_7(mach, l_nos, l_cil, L_nos)
            x_ffa_nos_cill.append(x_fa_nos_cill)
            x_fa_korm = L_f - 0.5 * L_korm

            x_fa_f = 1 / cy1_alf_f * (cy_nos * x_fa_nos_cill + cy_korm * x_fa_korm)
            # print(cy_korm)
            x_ffa_f.append(x_fa_f)

            # mz_alf = -cy1_alf * (x_fa_f - x_t) / L_f

            t += dt
            i += 1
            v_sr += self.v

            self.x += (self.v + self.v1) / 2 * kk.cos(self.fi_viz) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(self.fi_viz) * dt
            target.next_coord()
            xx.append(self.x)
            yy.append(self.y)
            xt.append(target.x)
            yt.append(target.y)
            vv.append(self.v)
            tt.append(t)
            print(t, self.v, self.v_)
            """if (self.v >= 650) and (self.v <= 660):
                print(t)"""
            if t <= t_st:
                self.m -= dm1 * dt
                self.v_ = 1 / self.m * (self.p1 - f_x) - g * kk.sin(self.omega)
            elif (t > t_st) and (t <= t_st + t_m):
                self.m -= dm2 * dt
                self.v_ = 1 / self.m * (self.p2 - f_x) - g * kk.sin(self.omega)
            else:
                self.v_ = 1 / self.m * (self.p3 - f_x) - g * kk.sin(self.omega)
            self.v += self.v_ * dt
            """if self.v < 650:
                self.v = self.v1
                self.v1 += a * dt
            elif self.v < 887:
                self.v = self.v1
                self.v1 += a1 * dt"""
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
        print(v_sr / i)

        plt.plot(tt, x_ffa_f)
        #plt.plot(tt, cyy_op)
        #plt.plot(tt, cyy1_delt_op)
        plt.grid(True)
        plt.show()

        plt.plot(tt, f_xx)
        plt.plot(tt, f_yy)
        plt.grid(True)
        plt.show()

        plt.plot(xx, yy)
        plt.plot(xt, yt)
        plt.show()
        plt.plot(tt, yy)
        plt.show()
        plt.plot(tt, xx)
        plt.show()
        plt.plot(tt, vv)
        plt.show()


dt = 10 ** -4
g = 9.80665
t_st = 2
t_m = 9
dm1 = (11.292 - 8.996) / t_st
dm2 = (8.996 - 7.3) / (t_m)
# (abs(target.y - self.y) > 25) or
vertel = Target()
igla = Rocket()
d = 0.072  # диаметр миделевого сечения

# параметры для метода наведения:

a_m = 2  # коэффициент быстроты реакции ракеты на маневр цели (от 1 до бесконечности)

l_cil = 45.5 / 8
l_nos = 12.5 / 8
tan_05_kr = 0.307  # тангенс среднего угла крыла
l_kr = 1.46  # относительное удлинение крыла
b_kr = 0.1061  # ширина крыла у корпуса
a_kr = 0.098  # ширина крыла у корпуса (меньшая грань трапеции)
c_kr = 0.03  # относительная толщина профиля крыла
koef_kr = 1 / (1 - a_kr / b_kr)  # коэффициент перехода от ромбовидного к шестиугольному профилю крыла

tan_05_op = -0.012  # тангенс среднего угла оперения
l_op = 7.888  # относительное удлинение оперения l^2/S
b_op = 0.032  # ширина оперения у корпуса
a_op = 0.015
c_op = 0.03  # относительная толщина профиля оперения
koef_op = 1 / (1 - a_op / b_op)
print(koef_kr, koef_op, "koef")

b_ak_kr = 0.094
x_otn_op_kr = 0.846 / b_ak_kr  # относительное расстояние между оперением и средней хордой крыльев
S_f = kk.pi * (d ** 2) / 4 / (kk.pi * (d ** 2) / 4)  # относительная площадь корпуса
S_op = 0.0078575 / (kk.pi * (d ** 2) / 4)  # относительная площадь передних несущих поверхностей
S_kr = 0.015 / (kk.pi * (d ** 2) / 4)  # относительная площадь задних несущих поверхностей

D_kr = 0.072 / 0.146  # отерсительный диаметр корпуса
nu_k_kr = 1.322  # относительное сужение консоли задней несущей поверхности
D_op = 0.072 / 0.25  # относительный диаметр корпуса
nu_k_op = 1.069  # относительное сужение консоли передней несущей поверхности

L_f = 1.626  # длина корпуса
L_nos = 0.155  # длина носовой части
W_nos = 0.00053  # Объем носовой части
L_korm = 0.160
S_dn = kk.pi * 0.05 ** 2 / 4
n_korm = 0.05 / 0.072
W_korm = 0.114 * S_dn + 0.000136
S_ff = kk.pi * d ** 2 / 4
Ff = 0.3619  # площадь обтекаемой потоком поверхности корпуса (без донного среза)
Sf = kk.pi * (d ** 2) / 4  #
f_t = 0.045216

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
