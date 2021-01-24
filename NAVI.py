import math as kk
import matplotlib.pyplot as plt
import numpy as np
import random
import Tabl
import Aero
import time
from mpl_toolkits.mplot3d import Axes3D
from tqdm import *
import cProfile


class Target(object):

    def __init__(self, *args):
        self.v = args[2]  # random.uniform(200, 200)  # сумарная скорость цели
        self.ang = 0
        self.x0 = args[0]  # random.uniform(6000, 6000)
        self.y0 = args[1]  # random.uniform(4100, 4100)
        self.alf0 = 0  # -0.0925  # random.uniform(0, 0.3925)
        self.x = self.x0
        self.y = self.y0
        self.alf = self.alf0

    def next_coord(self, *args):

        flag_t = args[0]
        # self.alf += random.uniform(-0.0086, +0.0086)
        if flag_t <= 2:
            self.x += self.v * kk.cos(self.alf) * dt
            self.y += self.v * kk.sin(self.alf) * dt
        else:
            self.v += 10 * dt
            # self.alf = -0.0925 + 0.1 * (flag_t - 2)
            self.x += self.v * kk.cos(self.alf) * dt
            self.y += self.v * kk.sin(self.alf) * dt


class Rocket(object):

    def __init__(self, *args):  # 0 - коорд у цели, 1 - коорд х цели
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
        self.omega = kk.atan((args[0] - self.y) / (args[1] - self.x))
        # print(kk.atan((args[0] - self.y) / (args[1] - self.x)), "start omega")
        self.d_omega = 0
        self.fi_viz = kk.atan((args[0] - self.y) / (args[1] - self.x))
        self.d_fi_viz = 0
        self.p = 2768
        self.p1 = 2768  # сила тяги на первом режиме
        self.p2 = 950  # 1275  # сила тяги на втором режиме
        self.p3 = 0
        self.x_ct = 0.78
        self.q = 0

    def navigation(self, *args):

        a = 249.7
        a1 = 29
        target = args[0]
        # form_nos = args[1]

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
        x_ffa_op = []
        x_ffa_kr = []
        da_ball = []
        m_zz_wz = []
        machh = []
        cyy_korm = []
        alff = []
        deltt = []
        delt_r = []
        flag_v = []
        x_ffa = []
        m_zz_wz_f = []
        m_zz_wz_op = []
        m_zz_wz_kr = []

        cyy_sum = []
        cxx_sum = []

        f_yy_exp = []
        f_xx_exp = []

        fxx_t = []
        fyy_t = []
        n_max = 0

        # ((abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5)) and
        # while (t < 16):
        ind = 0
        jam = 0
        rast_min = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
        rast_krit = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
        while (rast_min >= 5) and (t <= 16) and ((rast_krit - rast_min) <= 200):
            flag_v.append(650)
            # self.delt = kk.asin()
            mach = self.v / Tabl.tab_atm(self.y, 2)
            machh.append(mach)
            """
            значение коэффициента подъемной силы носовой части сложной формы (насадок + носовая часть корпуса)
            вычисляется через соотношения:
            """
            cy1_nos = Tabl.tab_3_2(mach, l_nos_nas, l_cil_nas)  # значение Су конической части АД насадка
            cy1_1_nos = Tabl.tab_3_4(mach, 0, l_cil_nas)  # значение Су сферической части АД насадка
            cy2_nos = 2 / 57.3 * kk.cos(Omega_con1) ** 2
            # cy2_nos = Tabl.tab_3_2(mach, l_con, l_cil)  # значение Су конической части нос. части без АД насадка
            cy2_1_nos = Tabl.tab_3_4(mach, 0, l_cil)  # значеное Су сферической части нос. части без АД насадка
            cy3_nos =  Tabl.tab_3_4(mach, 1, l_cil + l_nos)  # значение Су для цилиндра с плоским торцом (нос. часть без насадка)

            Cy1_2 = cy1_nos * (1 - r_ ** 2) + cy1_1_nos * r_ ** 2  # значение Су конуса со сферическим затуплением носовой части АД насадка

            Cy2_2 = cy2_nos * (1 - r_1 ** 2) + cy2_1_nos * r_1 ** 2  # значение Су конуса со сферическим затуплением без АД насадка
            Cy_ob = Cy2_2 * (1 - d_ ** 2) + cy3_nos * d_ ** 2  # значение Су носовой части с плоским затуплением
            cy_nos = Cy1_2 * (2 * r_n / D) ** 2 + Cy_ob  # значение Су носовой части
            cyy_nos.append(cy_nos)

            cy_korm = -0.2 * 2 / 57.3 * (1 - nu_korm ** 2)
            cyy_korm.append(cy_korm)
            cy_kr = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_22(mach, x_otn_op_kr)), l_kr, c_kr, tan_05_kr)
            cyy_kr.append(cy_kr)
            # if mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos))<= 1.05:
            cy_op = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos)), l_op, c_op, tan_05_op)

            cyy_op.append(cy_op)
            cy1_alf_f = cy_nos + cy_korm

            K_aa_kr = 1 + 3 * D_kr - (D_kr * (1 - D_kr)) / nu_k_kr
            k_aa_kr = (1 + 0.41 * D_kr) ** 2 * ((1 + 3 * D_kr - 1 / nu_k_kr * D_kr * (1 - D_kr)) / (1 + D_kr) ** 2)
            eps_sr_alf = 0  # для "утки"
            z_b = Tabl.tab_3_16(mach, l_op, nu_k_op, tan_05_op) # относительная координата вихря (по рис 3.16)
            i_ = 1 # Tabl.tab_3_17() # коэффициент интерференции вихрей
            if mach >= 1.3:
                fi_eps = 1 # отношение площади задней консоли, находящейся в конусе Маха, ко всей плозади консоли
            else:
                fi_eps = 0.3
            # eps_sr_alf = 57.3 / (2 * kk.pi) * i_ / z_b * L_k_op / L_k_kr * (cy1_alf_op / l_op) * k_aa_op / K_aa_kr * fi_eps
            cy1_alf_kr = (cy_kr * K_aa_kr) * (1 - eps_sr_alf)

            K_aa_op = 1 + 3 * D_op - (D_op * (1 - D_op)) / nu_k_op
            k_aa_op = (1 + 0.41 * D_op) ** 2 * ((1 + 3 * D_op - 1 / nu_k_op * D_op * (1 - D_op)) / (1 + D_op) ** 2)
            cy1_alf_op = cy_op * K_aa_op
            k_t_op = Tabl.tab_3_21(mach, l_nos)
            k_t_kr = Tabl.tab_3_22(mach, x_otn_op_kr)

            cy1_alf = cy1_alf_f * S__f + cy1_alf_kr * S__kr * k_t_kr + cy1_alf_op * S__op * k_t_op
            cyy1_alf.append(cy1_alf)

            K_delt_0_op = k_aa_op
            K_delt_0_kr = k_aa_kr
            k_delt_0_op = k_aa_op ** 2 / K_aa_op
            k_delt_0_kr = k_aa_kr ** 2 / K_aa_kr

            if mach <= 1:
                k_sh = 0.85
            elif mach <= 1.4 and mach > 1:
                k_sh = 0.85 + 0.14 * (mach - 1) / 0.4
            else:
                k_sh = 0.99
            n_ef = k_sh * kk.cos(0)
            cy1_delt_op = cy1_alf_op * K_delt_0_op * n_ef
            cyy1_delt_op.append(cy1_delt_op)

            ni_atm = Tabl.tab_atm(self.y, 5)
            re_f = self.v * L_f / ni_atm

            re_t = Tabl.tab_4_5(mach, re_f, 7, L_f)

            x_tt = re_t * ni_atm / self.v
            if x_tt >= (0.26 / L_f) or x_tt <= 0:
                x_tt = f_t / Ff
            #print(self.v, "V")
            #print(x_tt)
            cx_tr = Tabl.tab_4_2(re_f, x_tt) / 2 * Ff / Sf
            cxx_tr.append(cx_tr)

            cx_nos = Tabl.tab_4_11(mach, l_nos_)
            cxx_nos.append(cx_nos)

            cx_korm = Tabl.tab_4_24(mach, nu_korm, l_korm)
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
            else:

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

            cx_0 = 1.05 * (cx_0f * S__f + cx_0_op * k_t_op * S__op + cx_0_kr * k_t_kr * S__kr)
            cxx_0.append(cx_0)
            f_x = cx_0 * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_xx.append(f_x)
            f_y = cy1_alf * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_yy.append(f_y)
            # индуктивное сопротивление
            delt_cx1 = 2 * Tabl.tab_4_40(mach, l_nos, 1) * kk.sin(self.alf) ** 2
            cx_f_ind = cy1_alf_f * kk.sin(self.alf) + delt_cx1 * kk.cos(self.alf)
            D_0_op = K_aa_op - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_aa_op ** 2
            D_1_op = k_aa_op * (1 - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_delt_0_op * n_ef) + K_delt_0_op * n_ef
            D_2_op = k_delt_0_op * n_ef * (1 - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_delt_0_op * n_ef)

            D_0_kr = K_aa_kr - 57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr ** 2
            D_1_kr = -K_aa_kr + 2 * 57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr
            D_2_kr = -57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr

            if self.alf < 0.5:
                cx_op_ind = cy1_alf_op * (D_0_op + D_1_op * self.delt / 0.5 + D_2_op * (self.delt / 0.5) ** 2) * self.alf / 57.3
                cx_kr_ind = cy1_alf_kr * (D_0_kr + D_1_kr * self.delt / 0.5 + D_2_kr * (self.delt / 0.5) ** 2) * self.alf / 57.3
            else:
                cx_op_ind = cy1_alf_op * (D_0_op + D_1_op * self.delt / self.alf + D_2_op * (self.delt / self.alf) ** 2) * self.alf / 57.3
                cx_kr_ind = cy1_alf_kr * (D_0_kr + D_1_kr * self.delt / self.alf + D_2_kr * (self.delt / self.alf) ** 2) * self.alf / 57.3
            cx_ind = cx_f_ind * S__f + cx_op_ind * k_t_op * S__op + cx_kr_ind * k_t_kr * S__kr
            # cx_op_ind = cy

            # фокусы отдельных частей ЛА по углу атаки
            x_fa_nos_cill = L_nos - W_nos / S_f + Tabl.tab_5_7(mach, l_nos, l_cil, L_nos)
            x_ffa_nos_cill.append(x_fa_nos_cill)
            x_fa_korm = L_f - 0.5 * L_korm
            x_fa_f = 1 / cy1_alf_f * (cy_nos * x_fa_nos_cill + cy_korm * x_fa_korm)
            # print(cy_korm)
            x_ffa_f.append(x_fa_f)

            if mach < 1:
                x__iz_op = 0.21
                x__iz_kr = 0.15
            else:
                x__iz_op = 0.4
                x__iz_kr = 0.21

            x_f_iz_op = x_b_a_op + b_a_op * x__iz_op
            x_f_delt_op = x_f_iz_op - tan_05_op * Tabl.tab_5_11(D_op, L_k_op)
            x_f_iz_kr = x_b_a_kr + b_a_kr * x__iz_kr
            x_f_delt_kr = x_f_iz_kr - tan_05_kr * Tabl.tab_5_11(D_kr, L_k_kr)

            if mach > 1:
                b__b_op = b_op / (kk.pi / 2 * d * kk.sqrt(mach ** 2 - 1))
                b__b_kr = b_kr / (kk.pi / 2 * d * kk.sqrt(mach ** 2 - 1))
                L__hv_op = L_hv_op / (kk.pi * d * kk.sqrt(mach ** 2 - 1))
                L__hv_kr = L_hv_kr / (kk.pi * d * kk.sqrt(mach ** 2 - 1))
                F_1_op = 1 - 1 / (c_const_op * b__b_op ** 2) * (kk.e ** (-c_const_op * L__hv_op ** 2) - kk.e ** (-c_const_op * (b__b_op + L__hv_op) ** 2)) + kk.sqrt(kk.pi) / (b__b_op * kk.sqrt(c_const_op)) * Tabl.tab_int_ver(L__hv_op * kk.sqrt(2 * c_const_op))

                F_1_kr = 1 - 1 / (c_const_kr * b__b_kr ** 2) * (kk.e ** (-c_const_kr * L__hv_kr ** 2) - kk.e ** (-c_const_kr * (b__b_kr + L__hv_kr) ** 2)) + kk.sqrt(kk.pi) / (b__b_kr * kk.sqrt(c_const_kr)) * Tabl.tab_int_ver(L__hv_kr * kk.sqrt(2 * c_const_kr))

                F_op = 1 - kk.sqrt(kk.pi) / (2 * b__b_op * kk.sqrt(c_const_op)) * (Tabl.tab_int_ver((b__b_op + L__hv_op) * kk.sqrt(2 * c_const_op)) - Tabl.tab_int_ver(L__hv_op * kk.sqrt(2 * c_const_op)))
                F_kr = 1 - kk.sqrt(kk.pi) / (2 * b__b_kr * kk.sqrt(c_const_kr)) * (Tabl.tab_int_ver((b__b_kr + L__hv_kr) * kk.sqrt(2 * c_const_kr)) - Tabl.tab_int_ver(L__hv_kr * kk.sqrt(2 * c_const_kr)))

                x_f_b_op = x_f_iz_op + 0.02 * l_op * tan_05_op
                x_f_b_kr = x_f_iz_kr + 0.02 * l_kr * tan_05_kr

                x_f_ind_op = x_b_op + b_op * x_f_b_op * F_op * F_1_op
                x_f_ind_kr = x_b_kr + b_kr * x_f_b_kr * F_kr * F_1_kr
                x_fa_op = 1 / K_aa_op * (x_f_iz_op + (k_aa_op - 1) * x_f_delt_op + (K_aa_op - k_aa_op) * x_f_ind_op)
                x_fa_kr = 1 / K_aa_kr * (x_f_iz_kr + (k_aa_kr - 1) * x_f_delt_kr + (K_aa_kr - k_aa_kr) * x_f_ind_kr)
            else:
                x_f_b_op = x_f_iz_op + 0.02 * l_op * tan_05_op
                x_f_b_kr = x_f_iz_kr + 0.02 * l_kr * tan_05_kr
                x_f_ind_op = x_b_op + b_op * x_f_b_op
                x_f_ind_kr = x_b_kr + b_kr * x_f_b_kr
                x_fa_op = 1 / K_aa_op * (x_f_iz_op + (k_aa_op - 1) * x_f_delt_op + (K_aa_op - k_aa_op) * x_f_ind_op)
                x_fa_kr = 1 / K_aa_kr * (x_f_iz_kr + (k_aa_kr - 1) * x_f_delt_kr + (K_aa_kr - k_aa_kr) * x_f_ind_kr)

            x_ffa_op.append(x_fa_op)
            x_ffa_kr.append(x_fa_kr)

            x_fa = 1 / cy1_alf * ((cy1_alf_f * S__f * x_fa_f) + cy1_alf_op * S__op * x_fa_op * k_t_op + cy1_alf_kr * S__kr * x_fa_kr * k_t_kr)

            x_ffa.append(x_fa)

            # координаты фокуса рулей (передних консолей) по углам отклонения

            x_fd_op = 1 / K_delt_0_op * (k_delt_0_op * x_f_iz_op + (K_delt_0_op - k_delt_0_op) * x_f_ind_op)

            # Демпфирующие моменты АД поверхностей
            # x_c_ob - координата центра тяжести объема тела вращения
            x_c_ob = L_f * ((2 * (l_nos + l_cil)**2 - l_nos**2) / (4*(l_nos+l_cil) * (l_nos+l_cil - 2/3*l_nos)))
            m_z_wz_f = -2 * (1 - self.x_ct / L_f + (self.x_ct / L_f) ** 2 - x_c_ob / L_f)
            m_zz_wz_f.append(m_z_wz_f)

            x__ct_op = (self.x_ct - x_b_a_op) / b_a_op  # координата центра тяжести, измеренная от начала САХ рулей
            m_z_wz_op = -57.3 * (cy1_alf_op * (x__ct_op - 1 / 2) ** 2 * K_aa_op)
            m_zz_wz_op.append(m_z_wz_op)

            m_z_wz_delt_kr = -57.3 * (cy1_alf_kr * K_aa_kr) * eps_sr_alf * (self.x_ct - x_c_pl_ba) / b_a_kr * (self.x_ct - x_fa_kr) / b_a_kr
            x__ct_kr = (self.x_ct - x_b_a_kr) / b_a_kr
            m_z_wz_iz_kr = -57.3 * (cy1_alf_kr * (x__ct_kr - 1 / 2) ** 2 * K_aa_kr)
            m_z_wz_kr = m_z_wz_iz_kr * K_aa_kr + m_z_wz_delt_kr
            m_zz_wz_kr.append(m_z_wz_kr)

            m_z_wz = m_z_wz_f * S__f * L__f ** 2 + m_z_wz_op * S__op * b__a_op * kk.sqrt(k_t_op) + m_z_wz_kr * S__kr * b_a_kr * kk.sqrt(k_t_kr)
            m_zz_wz.append(m_z_wz)

            m_z_delt_op = cy1_delt_op * (self.x_ct - x_fd_op) / L_f

            m_z_alf = cy1_alf * (self.x_ct - x_fa) / L_f
            h_alf = x_fd_op - 0.333
            m_sharn_alf = -cy1_alf_op * h_alf / b_a_op
            h_delt = x_fd_op - 0.333
            m_sharn_delt = -cy1_delt_op * h_delt / b_a_op
            m_sharn = m_sharn_alf + m_sharn_delt

            da_ball.append(-m_z_alf / m_z_delt_op)

            t += dt
            i += 1
            v_sr += self.v

            # функция наведения (Метод наведения и кинематические параметры перехвата)

            par_1 = self.fi_viz
            self.fi_viz = kk.atan((target.y - self.y) / (target.x - self.x))
            self.d_fi_viz = (self.fi_viz - par_1) / dt
            self.d_omega = a_m * self.d_fi_viz
            """if (t <= 0.01) and (t >= 0):
                plt.plot([self.x, target.x], [self.y, target.y], '--')
                print(self.omega, kk.atan((target.y - self.y) / (target.x - self.x)))
            """
            if (kk.fabs(self.d_omega) <= d_omega_max):
                self.omega += self.d_omega * dt
            elif (self.d_omega >= 0):
                self.omega += d_omega_max * dt
            else:
                self.omega -= d_omega_max * dt

            if kk.fabs(self.omega - self.fi_viz) <= (15 * kk.pi / 180):
                self.delt_potr = self.fi_viz - self.omega
            elif (self.fi_viz - self.omega) < 0:
                self.delt_potr = delta_max
            else:
                self.delt_potr = - delta_max

            if (self.delt - self.delt_potr) <= 0:
                self.delt += d_delta_max * dt
            else:
                self.delt -= d_delta_max * dt
            deltt.append(self.delt_potr * 180 / kk.pi)
            delt_r.append(self.delt * 180 / kk.pi)

            n_y_a = kk.fabs(self.v * self.d_omega / g + kk.cos(self.omega))
            if kk.fabs(n_y_a) >= kk.fabs(n_max):
                n_max = n_y_a
            # print(self.omega * 180 / kk.pi,self.d_omega * 180 / kk.pi)
            # print(n_y_a)

            khi = cy1_delt_op * self.delt / cy1_alf
            # print(khi)
            q = Tabl.tab_atm(self.y, 4) * self.v ** 2 / 2
            # print(n_y_a * self.m * g)
            # print(cy1_alf * q * Sf * (1 + khi), self.p / 57.3)
            self.alf_potr = (n_y_a * self.m * g) / (cy1_alf * q * Sf * (1 + khi) + self.p / 57.3)
            #print(self.alf_potr)
            alff.append(self.alf_potr)
            if self.alf < self.alf_potr:
                self.alf += (d_delta_max) / 3 * dt
            elif self.alf > self.alf_potr:
                self.alf -= (d_delta_max) / 3 * dt
            if self.alf <= 0:
                self.alf = 0
            cy_sum = cy1_alf * self.alf * kk.pi / 180 + cy1_alf + cy1_delt_op * self.delt * kk.pi / 180
            cyy_sum.append(cy_sum)
            cx_sum = cx_0 + cx_ind
            cxx_sum.append(cx_sum)
            if jam <= cy1_alf_op * self.alf * kk.pi / 180 * self.v ** 2 * d ** 2 * kk.pi / 8 *  Tabl.tab_atm(self.y, 4):
                jam = cy1_alf_op * self.alf * kk.pi / 180 * self.v ** 2 * d ** 2 * kk.pi / 8 *  Tabl.tab_atm(self.y, 4)
            f_x_exp = cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_y_exp = cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_yy_exp.append(cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)
            f_xx_exp.append(cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)

            fxx_t.append(f_x_exp * dt)
            fyy_t.append(f_y_exp * dt)

            target.next_coord(t)
            xx.append(self.x)
            yy.append(self.y)
            xt.append(target.x)
            yt.append(target.y)
            vv.append(self.v)
            tt.append(t)
            # print(t, self.v, self.v_)
            if t <= t_st:
                self.p = self.p1
                self.x_ct += dx_cm1 * dt
                self.m -= dm1 * dt
                self.v_ = 1 / self.m * (self.p - f_x_exp) - g * kk.sin(self.omega + self.alf)
            elif (t > t_st) and (t <= t_st + t_m):
                self.p = self.p2
                self.x_ct += dx_cm2 * dt
                self.m -= dm2 * dt
                self.v_ = 1 / self.m * (self.p * kk.cos(self.alf) - f_x_exp) - g * kk.sin(self.omega + self.alf)
            else:
                self.p = self.p3
                self.v_ = 1 / self.m * (self.p * kk.cos(self.alf) - f_x_exp) - g * kk.sin(self.omega + self.alf)
            self.v1 = self.v
            self.v += self.v_ * dt
            self.x += (self.v + self.v1) / 2 * kk.cos(self.omega) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(self.omega) * dt

            if self.v >= 450:
                rast_krit = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
                if rast_min >= rast_krit:
                    rast_min = rast_krit

            ind += 1

        #print(v_sr / i, t, kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2))
        #print(jam, "Y_op")
        #print('ind =', ind, self.y, self.x)
        if fin_ch <= 0:
            return self.x, self.y, t, n_max, n_y_a, kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2), sum(fxx_t), sum(fyy_t), m_sharn
        else:
            plt.plot(tt, f_yy_exp, 'r', label = 'Подъемная сила')
            plt.plot(tt, f_xx_exp, 'g', label = 'Сила лобового сопротивления')
            plt.grid(True)
            plt.legend()
            plt.xlabel('t, [c]')
            plt.ylabel('F, [H]')
            plt.show()

            plt.plot(tt, cyy_sum, 'r', label = 'Коэффициент подъемной силы')
            plt.plot(tt, cxx_sum, 'g', label = 'Коэффициент лобового сопротивления')
            plt.grid(True)
            plt.xlabel('t, [c]')
            plt.ylabel('Cy, Cx')
            plt.legend()
            plt.show()

            plt.plot(tt, da_ball)
            plt.ylabel('delta/alf')
            plt.xlabel('t, [с]')
            plt.grid(True)
            plt.show()

            plt.plot(x_ffa_op, tt, 'r', label = 'Координата фокуса оперения')
            plt.plot(x_ffa_kr, tt, 'g', label = 'Координата фокуса крыльев')
            plt.plot(x_ffa_f, tt, 'b', label = 'Координата фокуса корпуса')
            plt.plot(x_ffa, tt, 'y', label = 'Координата фокуса Ла')
            plt.grid(True)
            plt.axis([0, 1.7, 0, 16])
            plt.ylabel('t, [с]')
            plt.xlabel('x_fa, [м]')
            plt.legend()
            plt.show()

            plt.plot(tt, machh)
            plt.grid(True)
            plt.ylabel('M')
            plt.xlabel('t, [с]')
            plt.show()


            plt.plot(xx, yy, 'r', label = 'Траектория ракеты')
            plt.plot(xt, yt, 'g', label = 'Траектория цели')
            plt.grid(True)
            plt.xlabel('x, [м]')
            plt.ylabel('y, [м]')
            plt.legend()
            plt.show()

            plt.plot(tt, alff, 'r', label = 'Угол атаки')
            plt.plot(tt, deltt, 'g', label = 'Потребный угол отклоения рулей')
            plt.plot(tt, delt_r, 'b', label = 'Реальный угол отклонения рулей')
            plt.axis([0, 16, -15, 15])
            plt.xlabel('t, [c]')
            plt.ylabel('alf, delt, [град]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, cyy_op, 'r', label = 'Cy_alf_op')
            plt.plot(tt, cyy_kr, 'g', label = 'Cy_alf_kr')
            plt.plot(tt, cyy1_alf, 'b', label = 'Cy1_alf')
            plt.plot(tt, cyy_korm, 'y', label = 'Cy_alf_korm')
            plt.plot(tt, cyy_nos, 'k', label = 'Cy_alf_nos')
            # plt.plot(tt, da_ball)
            # plt.plot(tt, m_zz_wz)
            # plt.axis([0, 16, -5, 5])
            # plt.plot(tt, cyy_op)
            # plt.plot(tt, cyy1_delt_op)
            plt.ylabel('Cy1_alf')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, cxx_nos, 'r', label = 'Сопротивление носовой части')
            plt.plot(tt, cxx_0f, 'g', label = 'Сумарное сопротивление корпуса')
            plt.plot(tt, cxx_korm, 'b', label = 'Сопротивление кормы')
            plt.plot(tt, cxx_tr, 'y', label = 'Сопротивление трению')
            plt.ylabel('Cx')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, cxx_kr_pr, 'r', label = 'Профильное сопротивление')
            plt.plot(tt, cxx_kr_vol, 'g', label = 'Волновое сопротивление')
            plt.plot(tt, cxx_0_kr, 'b', label = 'Суммарное сопротивление крыльев')
            plt.ylabel('Cx')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, cxx_op_pr, 'r', label = 'Профильное сопротивление')
            plt.plot(tt, cxx_op_vol, 'g', label = 'Волновое сопротивление')
            plt.plot(tt, cxx_0_op, 'b', label = 'Суммарное сопротивление оперения')
            plt.ylabel('Cx_0')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, cxx_0_op, 'r', label = 'Суммарное лобовое сопротивление')
            plt.plot(tt, cxx_0, 'g', label = 'Лобовое сопротивление оперения')
            plt.plot(tt, cxx_0_kr, 'b', label = 'Лобовое сопротивление крыльев')
            plt.plot(tt, cxx_0f, 'y', label = 'Лобовое сопротивление корпуса')
            plt.ylabel('Cx_0')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, m_zz_wz, 'r', label = 'Суммарный демпфирующий момент')
            plt.plot(tt, m_zz_wz_kr, 'g', label = 'Демпфирующий момент крыльев')
            plt.plot(tt, m_zz_wz_op, 'b', label = 'Демпфирующий момент оперения')
            plt.plot(tt, m_zz_wz_f, 'y', label = 'Демпфирующий момент корпуса')
            plt.ylabel('mz_wz')
            plt.xlabel('t, [с]')
            plt.legend()
            plt.grid(True)
            plt.show()

            plt.plot(tt, f_xx)
            plt.plot(tt, f_yy)
            plt.grid(True)
            plt.show()


            plt.plot(xx, yy)
            plt.plot(xt, yt)
            plt.grid(True)
            plt.show()
            plt.plot(tt, yy)
            plt.show()
            plt.plot(tt, xx)
            plt.show()
            plt.plot(tt, vv)
            plt.plot(tt, flag_v, '--')
            plt.ylabel('V, [м/с]')
            plt.xlabel('t, [с]')
            plt.grid(True)
            plt.show()


class Target2(object):

    def __init__(self, *args):
        self.v = args[2]  # random.uniform(200, 200)  # сумарная скорость цели
        self.ang = 0
        self.x0 = args[0]  # random.uniform(6000, 6000)
        self.y0 = args[1]  # random.uniform(4100, 4100)
        self.alf0 = 0.0925  # random.uniform(0, 0.3925)
        self.x = self.x0
        self.y = self.y0
        self.alf = self.alf0

    def next_coord(self, *args):

        flag_t = args[0]
        # self.alf += random.uniform(-0.0086, +0.0086)
        if flag_t <= 2:
            self.x += self.v * kk.cos(self.alf) * dt
            self.y += self.v * kk.sin(self.alf) * dt
        else:

            self.alf = 0.0925 + 0.1 * (flag_t - 2)
            self.x += self.v * kk.cos(self.alf) * dt
            self.y += self.v * kk.sin(self.alf) * dt


class Rocket2(object):

    def __init__(self, *args):  # 0 - коорд у цели, 1 - коорд х цели
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
        self.omega = kk.atan((args[0] - self.y) / (args[1] - self.x))
        # print(kk.atan((args[0] - self.y) / (args[1] - self.x)), "start omega")
        self.d_omega = 0
        self.fi_viz = kk.atan((args[0] - self.y) / (args[1] - self.x))
        self.d_fi_viz = 0
        self.p = 2768
        self.p1 = 2768  # сила тяги на первом режиме
        self.p2 = 950  # 1275  # сила тяги на втором режиме
        self.p3 = 0
        self.x_ct = 0.78
        self.q = 0

    def navigation(self, *args):

        a = 249.7
        a1 = 29
        target = args[0]
        # form_nos = args[1]

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
        x_ffa_op = []
        x_ffa_kr = []
        da_ball = []
        m_zz_wz = []
        machh = []
        cyy_korm = []
        alff = []
        deltt = []
        delt_r = []
        flag_v = []
        x_ffa = []
        m_zz_wz_f = []
        m_zz_wz_op = []
        m_zz_wz_kr = []

        cyy_sum = []
        cxx_sum = []

        f_yy_exp = []
        f_yy_kr_exp = []
        f_xx_exp = []
        f_xx_kr_exp = []

        fxx_t = []
        fyy_t = []
        n_max = 0

        # ((abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5)) and
        # while (t < 16):
        ind = 0
        jam = 0
        rast_min = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
        rast_krit = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
        while (rast_min >= 5) and (t <= 16) and ((rast_krit - rast_min) <= 200):
            flag_v.append(650)
            # self.delt = kk.asin()
            mach = self.v / Tabl.tab_atm(self.y, 2)
            machh.append(mach)
            """
            значение коэффициента подъемной силы носовой части сложной формы (насадок + носовая часть корпуса)
            вычисляется через соотношения:
            """
            cy1_nos = Tabl.tab_3_2(mach, l_nos_nas, l_cil_nas)  # значение Су конической части АД насадка
            cy1_1_nos = Tabl.tab_3_4(mach, 0, l_cil_nas)  # значение Су сферической части АД насадка
            cy2_nos = 2 / 57.3 * kk.cos(Omega_con1) ** 2
            # cy2_nos = Tabl.tab_3_2(mach, l_con, l_cil)  # значение Су конической части нос. части без АД насадка
            cy2_1_nos = Tabl.tab_3_4(mach, 0, l_cil)  # значеное Су сферической части нос. части без АД насадка
            cy3_nos =  Tabl.tab_3_4(mach, 1, l_cil + l_nos)  # значение Су для цилиндра с плоским торцом (нос. часть без насадка)

            Cy1_2 = cy1_nos * (1 - r_ ** 2) + cy1_1_nos * r_ ** 2  # значение Су конуса со сферическим затуплением носовой части АД насадка

            Cy2_2 = cy2_nos * (1 - r_1 ** 2) + cy2_1_nos * r_1 ** 2  # значение Су конуса со сферическим затуплением без АД насадка
            Cy_ob = Cy2_2 * (1 - d_ ** 2) + cy3_nos * d_ ** 2  # значение Су носовой части с плоским затуплением
            cy_nos = Cy1_2 * (2 * r_n / D) ** 2 + Cy_ob  # значение Су носовой части
            cyy_nos.append(cy_nos)

            cy_korm = -0.2 * 2 / 57.3 * (1 - nu_korm ** 2)
            cyy_korm.append(cy_korm)
            cy_kr = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_22(mach, x_otn_op_kr)), l_kr, c_kr, tan_05_kr)
            cyy_kr.append(cy_kr)
            # if mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos))<= 1.05:
            cy_op = Tabl.tab_3_5(mach * kk.sqrt(Tabl.tab_3_21(mach, l_nos)), l_op, c_op, tan_05_op)

            cyy_op.append(cy_op)
            cy1_alf_f = cy_nos + cy_korm

            K_aa_kr = 1 + 3 * D_kr - (D_kr * (1 - D_kr)) / nu_k_kr
            k_aa_kr = (1 + 0.41 * D_kr) ** 2 * ((1 + 3 * D_kr - 1 / nu_k_kr * D_kr * (1 - D_kr)) / (1 + D_kr) ** 2)
            eps_sr_alf = 0  # для "утки"
            z_b = Tabl.tab_3_16(mach, l_op, nu_k_op, tan_05_op) # относительная координата вихря (по рис 3.16)
            i_ = 1 # Tabl.tab_3_17() # коэффициент интерференции вихрей
            if mach >= 1.3:
                fi_eps = 1 # отношение площади задней консоли, находящейся в конусе Маха, ко всей плозади консоли
            else:
                fi_eps = 0.3
            # eps_sr_alf = 57.3 / (2 * kk.pi) * i_ / z_b * L_k_op / L_k_kr * (cy1_alf_op / l_op) * k_aa_op / K_aa_kr * fi_eps
            cy1_alf_kr = (cy_kr * K_aa_kr) * (1 - eps_sr_alf)

            K_aa_op = 1 + 3 * D_op - (D_op * (1 - D_op)) / nu_k_op
            k_aa_op = (1 + 0.41 * D_op) ** 2 * ((1 + 3 * D_op - 1 / nu_k_op * D_op * (1 - D_op)) / (1 + D_op) ** 2)
            cy1_alf_op = cy_op * K_aa_op
            k_t_op = Tabl.tab_3_21(mach, l_nos)
            k_t_kr = Tabl.tab_3_22(mach, x_otn_op_kr)

            cy1_alf = cy1_alf_f * S__f + cy1_alf_kr * S__kr * k_t_kr + cy1_alf_op * S__op * k_t_op
            cyy1_alf.append(cy1_alf)

            K_delt_0_op = k_aa_op
            K_delt_0_kr = k_aa_kr
            k_delt_0_op = k_aa_op ** 2 / K_aa_op
            k_delt_0_kr = k_aa_kr ** 2 / K_aa_kr

            if mach <= 1:
                k_sh = 0.85
            elif mach <= 1.4 and mach > 1:
                k_sh = 0.85 + 0.14 * (mach - 1) / 0.4
            else:
                k_sh = 0.99
            n_ef = k_sh * kk.cos(0)
            cy1_delt_op = cy1_alf_op * K_delt_0_op * n_ef
            cyy1_delt_op.append(cy1_delt_op)

            ni_atm = Tabl.tab_atm(self.y, 5)
            re_f = self.v * L_f / ni_atm

            re_t = Tabl.tab_4_5(mach, re_f, 7, L_f)

            x_tt = re_t * ni_atm / self.v
            if x_tt >= (0.26 / L_f) or x_tt <= 0:
                x_tt = f_t / Ff
            #print(self.v, "V")
            #print(x_tt)
            cx_tr = Tabl.tab_4_2(re_f, x_tt) / 2 * Ff / Sf
            cxx_tr.append(cx_tr)

            cx_nos = Tabl.tab_4_11(mach, l_nos_)
            cxx_nos.append(cx_nos)

            cx_korm = Tabl.tab_4_24(mach, nu_korm, l_korm)
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
            else:

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

            cx_0 = 1.05 * (cx_0f * S__f + cx_0_op * k_t_op * S__op + cx_0_kr * k_t_kr * S__kr)
            cxx_0.append(cx_0)
            f_x = cx_0 * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_xx.append(f_x)
            f_y = cy1_alf * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_yy.append(f_y)
            # индуктивное сопротивление
            delt_cx1 = 2 * Tabl.tab_4_40(mach, l_nos, 1) * kk.sin(self.alf) ** 2
            cx_f_ind = cy1_alf_f * kk.sin(self.alf) + delt_cx1 * kk.cos(self.alf)
            D_0_op = K_aa_op - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_aa_op ** 2
            D_1_op = k_aa_op * (1 - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_delt_0_op * n_ef) + K_delt_0_op * n_ef
            D_2_op = k_delt_0_op * n_ef * (1 - 57.3 * 0.2 * 0.3 * cy1_alf_op * k_delt_0_op * n_ef)

            D_0_kr = K_aa_kr - 57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr ** 2
            D_1_kr = -K_aa_kr + 2 * 57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr
            D_2_kr = -57.3 * 0.2 * 0.3 * cy1_alf_kr * k_aa_kr

            if self.alf < 0.5:
                cx_op_ind = cy1_alf_op * (D_0_op + D_1_op * self.delt / 0.5 + D_2_op * (self.delt / 0.5) ** 2) * self.alf / 57.3
                cx_kr_ind = cy1_alf_kr * (D_0_kr + D_1_kr * self.delt / 0.5 + D_2_kr * (self.delt / 0.5) ** 2) * self.alf / 57.3
            else:
                cx_op_ind = cy1_alf_op * (D_0_op + D_1_op * self.delt / self.alf + D_2_op * (self.delt / self.alf) ** 2) * self.alf / 57.3
                cx_kr_ind = cy1_alf_kr * (D_0_kr + D_1_kr * self.delt / self.alf + D_2_kr * (self.delt / self.alf) ** 2) * self.alf / 57.3
            cx_ind = cx_f_ind * S__f + cx_op_ind * k_t_op * S__op + cx_kr_ind * k_t_kr * S__kr
            # cx_op_ind = cy

            # фокусы отдельных частей ЛА по углу атаки
            x_fa_nos_cill = L_nos - W_nos / S_f + Tabl.tab_5_7(mach, l_nos, l_cil, L_nos)
            x_ffa_nos_cill.append(x_fa_nos_cill)
            x_fa_korm = L_f - 0.5 * L_korm
            x_fa_f = 1 / cy1_alf_f * (cy_nos * x_fa_nos_cill + cy_korm * x_fa_korm)
            # print(cy_korm)
            x_ffa_f.append(x_fa_f)

            if mach < 1:
                x__iz_op = 0.21
                x__iz_kr = 0.15
            else:
                x__iz_op = 0.4
                x__iz_kr = 0.21

            x_f_iz_op = x_b_a_op + b_a_op * x__iz_op
            x_f_delt_op = x_f_iz_op - tan_05_op * Tabl.tab_5_11(D_op, L_k_op)
            x_f_iz_kr = x_b_a_kr + b_a_kr * x__iz_kr
            x_f_delt_kr = x_f_iz_kr - tan_05_kr * Tabl.tab_5_11(D_kr, L_k_kr)

            if mach > 1:
                b__b_op = b_op / (kk.pi / 2 * d * kk.sqrt(mach ** 2 - 1))
                b__b_kr = b_kr / (kk.pi / 2 * d * kk.sqrt(mach ** 2 - 1))
                L__hv_op = L_hv_op / (kk.pi * d * kk.sqrt(mach ** 2 - 1))
                L__hv_kr = L_hv_kr / (kk.pi * d * kk.sqrt(mach ** 2 - 1))
                F_1_op = 1 - 1 / (c_const_op * b__b_op ** 2) * (kk.e ** (-c_const_op * L__hv_op ** 2) - kk.e ** (-c_const_op * (b__b_op + L__hv_op) ** 2)) + kk.sqrt(kk.pi) / (b__b_op * kk.sqrt(c_const_op)) * Tabl.tab_int_ver(L__hv_op * kk.sqrt(2 * c_const_op))

                F_1_kr = 1 - 1 / (c_const_kr * b__b_kr ** 2) * (kk.e ** (-c_const_kr * L__hv_kr ** 2) - kk.e ** (-c_const_kr * (b__b_kr + L__hv_kr) ** 2)) + kk.sqrt(kk.pi) / (b__b_kr * kk.sqrt(c_const_kr)) * Tabl.tab_int_ver(L__hv_kr * kk.sqrt(2 * c_const_kr))

                F_op = 1 - kk.sqrt(kk.pi) / (2 * b__b_op * kk.sqrt(c_const_op)) * (Tabl.tab_int_ver((b__b_op + L__hv_op) * kk.sqrt(2 * c_const_op)) - Tabl.tab_int_ver(L__hv_op * kk.sqrt(2 * c_const_op)))
                F_kr = 1 - kk.sqrt(kk.pi) / (2 * b__b_kr * kk.sqrt(c_const_kr)) * (Tabl.tab_int_ver((b__b_kr + L__hv_kr) * kk.sqrt(2 * c_const_kr)) - Tabl.tab_int_ver(L__hv_kr * kk.sqrt(2 * c_const_kr)))

                x_f_b_op = x_f_iz_op + 0.02 * l_op * tan_05_op
                x_f_b_kr = x_f_iz_kr + 0.02 * l_kr * tan_05_kr

                x_f_ind_op = x_b_op + b_op * x_f_b_op * F_op * F_1_op
                x_f_ind_kr = x_b_kr + b_kr * x_f_b_kr * F_kr * F_1_kr
                x_fa_op = 1 / K_aa_op * (x_f_iz_op + (k_aa_op - 1) * x_f_delt_op + (K_aa_op - k_aa_op) * x_f_ind_op)
                x_fa_kr = 1 / K_aa_kr * (x_f_iz_kr + (k_aa_kr - 1) * x_f_delt_kr + (K_aa_kr - k_aa_kr) * x_f_ind_kr)
            else:
                x_f_b_op = x_f_iz_op + 0.02 * l_op * tan_05_op
                x_f_b_kr = x_f_iz_kr + 0.02 * l_kr * tan_05_kr
                x_f_ind_op = x_b_op + b_op * x_f_b_op
                x_f_ind_kr = x_b_kr + b_kr * x_f_b_kr
                x_fa_op = 1 / K_aa_op * (x_f_iz_op + (k_aa_op - 1) * x_f_delt_op + (K_aa_op - k_aa_op) * x_f_ind_op)
                x_fa_kr = 1 / K_aa_kr * (x_f_iz_kr + (k_aa_kr - 1) * x_f_delt_kr + (K_aa_kr - k_aa_kr) * x_f_ind_kr)

            x_ffa_op.append(x_fa_op)
            x_ffa_kr.append(x_fa_kr)

            x_fa = 1 / cy1_alf * ((cy1_alf_f * S__f * x_fa_f) + cy1_alf_op * S__op * x_fa_op * k_t_op + cy1_alf_kr * S__kr * x_fa_kr * k_t_kr)

            x_ffa.append(x_fa)

            # координаты фокуса рулей (передних консолей) по углам отклонения

            x_fd_op = 1 / K_delt_0_op * (k_delt_0_op * x_f_iz_op + (K_delt_0_op - k_delt_0_op) * x_f_ind_op)

            # Демпфирующие моменты АД поверхностей
            # x_c_ob - координата центра тяжести объема тела вращения
            x_c_ob = L_f * ((2 * (l_nos + l_cil)**2 - l_nos**2) / (4*(l_nos+l_cil) * (l_nos+l_cil - 2/3*l_nos)))
            m_z_wz_f = -2 * (1 - self.x_ct / L_f + (self.x_ct / L_f) ** 2 - x_c_ob / L_f)
            m_zz_wz_f.append(m_z_wz_f)

            x__ct_op = (self.x_ct - x_b_a_op) / b_a_op  # координата центра тяжести, измеренная от начала САХ рулей
            m_z_wz_op = -57.3 * (cy1_alf_op * (x__ct_op - 1 / 2) ** 2 * K_aa_op)
            m_zz_wz_op.append(m_z_wz_op)

            m_z_wz_delt_kr = -57.3 * (cy1_alf_kr * K_aa_kr) * eps_sr_alf * (self.x_ct - x_c_pl_ba) / b_a_kr * (self.x_ct - x_fa_kr) / b_a_kr
            x__ct_kr = (self.x_ct - x_b_a_kr) / b_a_kr
            m_z_wz_iz_kr = -57.3 * (cy1_alf_kr * (x__ct_kr - 1 / 2) ** 2 * K_aa_kr)
            m_z_wz_kr = m_z_wz_iz_kr * K_aa_kr + m_z_wz_delt_kr
            m_zz_wz_kr.append(m_z_wz_kr)

            m_z_wz = m_z_wz_f * S__f * L__f ** 2 + m_z_wz_op * S__op * b__a_op * kk.sqrt(k_t_op) + m_z_wz_kr * S__kr * b_a_kr * kk.sqrt(k_t_kr)
            m_zz_wz.append(m_z_wz)

            m_z_delt_op = cy1_delt_op * (self.x_ct - x_fd_op) / L_f

            m_z_alf = cy1_alf * (self.x_ct - x_fa) / L_f
            h_alf = x_fd_op - 0.333
            m_sharn_alf = -cy1_alf_op * h_alf / b_a_op
            h_delt = x_fd_op - 0.333
            m_sharn_delt = -cy1_delt_op * h_delt / b_a_op
            m_sharn = m_sharn_alf + m_sharn_delt

            da_ball.append(-m_z_alf / m_z_delt_op)

            t += dt
            i += 1
            v_sr += self.v

            # функция наведения (Метод наведения и кинематические параметры перехвата)

            par_1 = self.fi_viz
            self.fi_viz = kk.atan((target.y - self.y) / (target.x - self.x))
            self.d_fi_viz = (self.fi_viz - par_1) / dt
            self.d_omega = a_m * self.d_fi_viz
            """if (t <= 0.01) and (t >= 0):
                plt.plot([self.x, target.x], [self.y, target.y], '--')
                print(self.omega, kk.atan((target.y - self.y) / (target.x - self.x)))
            """
            if (kk.fabs(self.d_omega) <= d_omega_max):
                self.omega += self.d_omega * dt
            elif (self.d_omega >= 0):
                self.omega += d_omega_max * dt
            else:
                self.omega -= d_omega_max * dt

            if kk.fabs(self.omega - self.fi_viz) <= (15 * kk.pi / 180):
                self.delt_potr = self.fi_viz - self.omega
            elif (self.fi_viz - self.omega) < 0:
                self.delt_potr = delta_max
            else:
                self.delt_potr = - delta_max

            if (self.delt - self.delt_potr) <= 0:
                self.delt += d_delta_max * dt
            else:
                self.delt -= d_delta_max * dt
            deltt.append(self.delt_potr * 180 / kk.pi)
            delt_r.append(self.delt * 180 / kk.pi)

            n_y_a = kk.fabs(self.v * self.d_omega / g + kk.cos(self.omega))
            if kk.fabs(n_y_a) >= kk.fabs(n_max):
                n_max = n_y_a
            # print(self.omega * 180 / kk.pi,self.d_omega * 180 / kk.pi)
            # print(n_y_a)

            khi = cy1_delt_op * self.delt / cy1_alf
            # print(khi)
            q = Tabl.tab_atm(self.y, 4) * self.v ** 2 / 2
            # print(n_y_a * self.m * g)
            # print(cy1_alf * q * Sf * (1 + khi), self.p / 57.3)
            self.alf_potr = (n_y_a * self.m * g) / (cy1_alf * q * Sf * (1 + khi) + self.p / 57.3)
            #print(self.alf_potr)
            alff.append(self.alf_potr)
            if self.alf < self.alf_potr:
                self.alf += (d_delta_max) / 3 * dt
            elif self.alf > self.alf_potr:
                self.alf -= (d_delta_max) / 3 * dt
            if self.alf <= 0:
                self.alf = 0
            cy_sum = cy1_alf * self.alf * kk.pi / 180 + cy1_alf + cy1_delt_op * self.delt * kk.pi / 180
            cyy_sum.append(cy_sum)
            cx_sum = cx_0 + cx_ind
            cxx_sum.append(cx_sum)
            if jam <= cy1_alf_op * self.alf * kk.pi / 180 * self.v ** 2 * d ** 2 * kk.pi / 8 *  Tabl.tab_atm(self.y, 4):
                jam = cy1_alf_op * self.alf * kk.pi / 180 * self.v ** 2 * d ** 2 * kk.pi / 8 *  Tabl.tab_atm(self.y, 4)
            f_x_exp = cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_y_exp = cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8

            f_yy_exp.append(cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)
            f_yy_kr_exp.append(cy_kr * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)
            f_xx_exp.append(cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)
            f_xx_kr_exp.append((cx_kr_ind * k_t_kr * S__kr + cx_0_kr) * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)

            fxx_t.append(f_x_exp * dt)
            fyy_t.append(f_y_exp * dt)

            target.next_coord(t)
            xx.append(self.x)
            yy.append(self.y)
            xt.append(target.x)
            yt.append(target.y)
            vv.append(self.v)
            tt.append(t)
            # print(t, self.v, self.v_)
            if t <= t_st:
                self.p = self.p1
                self.x_ct += dx_cm1 * dt
                self.m -= dm1 * dt
                self.v_ = 1 / self.m * (self.p - f_x_exp) - g * kk.sin(self.omega + self.alf)
            elif (t > t_st) and (t <= t_st + t_m):
                self.p = self.p2
                self.x_ct += dx_cm2 * dt
                self.m -= dm2 * dt
                self.v_ = 1 / self.m * (self.p * kk.cos(self.alf) - f_x_exp) - g * kk.sin(self.omega + self.alf)
            else:
                self.p = self.p3
                self.v_ = 1 / self.m * (self.p * kk.cos(self.alf) - f_x_exp) - g * kk.sin(self.omega + self.alf)
            self.v1 = self.v
            self.v += self.v_ * dt
            self.x += (self.v + self.v1) / 2 * kk.cos(self.omega) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(self.omega) * dt

            if self.v >= 450:
                rast_krit = kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2)
                if rast_min >= rast_krit:
                    rast_min = rast_krit

            ind += 1

        #print(v_sr / i, t, kk.sqrt((target.y - self.y) ** 2 + (target.x - self.x) ** 2))
        #print(jam, "Y_op")
        #print('ind =', ind, self.y, self.x)

        return tt, f_yy_exp, f_xx_exp, cyy_sum, cxx_sum,  da_ball, x_ffa_op, x_ffa_kr, x_ffa_f, x_ffa,  machh, xx, yy, \
               xt, yt, alff, deltt, delt_r, cyy_op, cyy_kr, cyy1_alf, cyy_korm, cyy_nos,  cxx_nos,  cxx_0f, cxx_korm, \
               cxx_tr,  cxx_kr_pr, cxx_kr_vol, cxx_0_kr, cxx_op_pr, cxx_op_vol, cxx_0_op, cxx_0_op, cxx_0, cxx_0_kr, \
               cxx_0f, m_zz_wz, m_zz_wz_kr, m_zz_wz_op, m_zz_wz_f, f_xx_kr_exp, f_yy_kr_exp
        plt.plot(tt, f_yy_exp, 'r', label = 'Подъемная сила')
        plt.plot(tt, f_xx_exp, 'g', label = 'Сила лобового сопротивления')
        plt.grid(True)
        plt.legend()
        plt.xlabel('t, [c]')
        plt.ylabel('F, [H]')
        plt.show()

        plt.plot(tt, cyy_sum, 'r', label = 'Коэффициент подъемной силы')
        plt.plot(tt, cxx_sum, 'g', label = 'Коэффициент лобового сопротивления')
        plt.grid(True)
        plt.xlabel('t, [c]')
        plt.ylabel('Cy, Cx')
        plt.legend()
        plt.show()

        plt.plot(tt, da_ball)
        plt.ylabel('delta/alf')
        plt.xlabel('t, [с]')
        plt.grid(True)
        plt.show()

        plt.plot(x_ffa_op, tt, 'r', label = 'Координата фокуса оперения')
        plt.plot(x_ffa_kr, tt, 'g', label = 'Координата фокуса крыльев')
        plt.plot(x_ffa_f, tt, 'b', label = 'Координата фокуса корпуса')
        plt.plot(x_ffa, tt, 'y', label = 'Координата фокуса Ла')
        plt.grid(True)
        plt.axis([0, 1.7, 0, 16])
        plt.ylabel('t, [с]')
        plt.xlabel('x_fa, [м]')
        plt.legend()
        plt.show()

        plt.plot(tt, machh)
        plt.grid(True)
        plt.ylabel('M')
        plt.xlabel('t, [с]')
        plt.show()


        plt.plot(xx, yy, 'r', label = 'Траектория ракеты')
        plt.plot(xt, yt, 'g', label = 'Траектория цели')
        plt.grid(True)
        plt.xlabel('x, [м]')
        plt.ylabel('y, [м]')
        plt.legend()
        plt.show()

        plt.plot(tt, alff, 'r', label = 'Угол атаки')
        plt.plot(tt, deltt, 'g', label = 'Потребный угол отклоения рулей')
        plt.plot(tt, delt_r, 'b', label = 'Реальный угол отклонения рулей')
        plt.axis([0, 16, -15, 15])
        plt.xlabel('t, [c]')
        plt.ylabel('alf, delt, [град]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, cyy_op, 'r', label = 'Cy_alf_op')
        plt.plot(tt, cyy_kr, 'g', label = 'Cy_alf_kr')
        plt.plot(tt, cyy1_alf, 'b', label = 'Cy1_alf')
        plt.plot(tt, cyy_korm, 'y', label = 'Cy_alf_korm')
        plt.plot(tt, cyy_nos, 'k', label = 'Cy_alf_nos')
        # plt.plot(tt, da_ball)
        # plt.plot(tt, m_zz_wz)
        # plt.axis([0, 16, -5, 5])
        # plt.plot(tt, cyy_op)
        # plt.plot(tt, cyy1_delt_op)
        plt.ylabel('Cy1_alf')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, cxx_nos, 'r', label = 'Сопротивление носовой части')
        plt.plot(tt, cxx_0f, 'g', label = 'Сумарное сопротивление корпуса')
        plt.plot(tt, cxx_korm, 'b', label = 'Сопротивление кормы')
        plt.plot(tt, cxx_tr, 'y', label = 'Сопротивление трению')
        plt.ylabel('Cx')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, cxx_kr_pr, 'r', label = 'Профильное сопротивление')
        plt.plot(tt, cxx_kr_vol, 'g', label = 'Волновое сопротивление')
        plt.plot(tt, cxx_0_kr, 'b', label = 'Суммарное сопротивление крыльев')
        plt.ylabel('Cx')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, cxx_op_pr, 'r', label = 'Профильное сопротивление')
        plt.plot(tt, cxx_op_vol, 'g', label = 'Волновое сопротивление')
        plt.plot(tt, cxx_0_op, 'b', label = 'Суммарное сопротивление оперения')
        plt.ylabel('Cx_0')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, cxx_0_op, 'r', label = 'Суммарное лобовое сопротивление')
        plt.plot(tt, cxx_0, 'g', label = 'Лобовое сопротивление оперения')
        plt.plot(tt, cxx_0_kr, 'b', label = 'Лобовое сопротивление крыльев')
        plt.plot(tt, cxx_0f, 'y', label = 'Лобовое сопротивление корпуса')
        plt.ylabel('Cx_0')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, m_zz_wz, 'r', label = 'Суммарный демпфирующий момент')
        plt.plot(tt, m_zz_wz_kr, 'g', label = 'Демпфирующий момент крыльев')
        plt.plot(tt, m_zz_wz_op, 'b', label = 'Демпфирующий момент оперения')
        plt.plot(tt, m_zz_wz_f, 'y', label = 'Демпфирующий момент корпуса')
        plt.ylabel('mz_wz')
        plt.xlabel('t, [с]')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.plot(tt, f_xx)
        plt.plot(tt, f_yy)
        plt.grid(True)
        plt.show()


        plt.plot(xx, yy)
        plt.plot(xt, yt)
        plt.grid(True)
        plt.show()
        plt.plot(tt, yy)
        plt.show()
        plt.plot(tt, xx)
        plt.show()
        plt.plot(tt, vv)
        plt.plot(tt, flag_v, '--')
        plt.ylabel('V, [м/с]')
        plt.xlabel('t, [с]')
        plt.grid(True)
        plt.show()


def krit_maxi(m_0_w, Sfy, Sfx):
    return (m_0_w / m_0_w_0) ** (0) * (Sfy / pp2[k_paramz]) ** (2) * (Sfx / pp1[k_paramz]) ** (-4)

fin_ch = 0
g = 9.80665

t_st = 2.7  # время стартового участка
t_m = 8.3  # время маршевого участка

dm1 = (11.292 - 8.996) / t_st  # производная изменения массы на стартовом участке
dm2 = (8.996 - 7.3) / t_m  # производная изменения массы на маршевом участке

dx_cm1 = (0.671 - 0.78) / t_st  # скорость изменения цм на стартовом участке
dx_cm2 = (0.637 - 0.671) / t_m  # скорость изменеия цт на маршевом учестке

d = 0.072  # диаметр миделевого сечения

# параметры для метода наведения:

a_m = 2.8  # коэффициент быстроты реакции ракеты на маневр цели (от 1 до бесконечности)
eps_krit = 38 * kk.pi / 180  # предельное отклонение координатора головки самонаведения
d_omega_max = 40 * kk.pi / 180
delta_max = 15 * kk.pi / 180  # предельный угол отклонения рулей
d_delta_max = 35 * kk.pi / 180 # * 20 * 4  # скорость отклонения рулей

# геометрические параметры ракеты:

# геометрические параметры корпуса
L_f = 1.626  # длина корпуса
L__f = L_f / L_f  # относительная длина корпуса
L_nos = 0.0974  # длина носовой части
W_nos = 0.00053  # Объем носовой части
L_korm = 0.160  # длина кормовой части
S_dn = kk.pi * 0.05 ** 2 / 4  # площадь донного среза
nu_korm = 0.05 / 0.072  # относительное сужение кормовой части
W_korm = 0.114 * S_dn + 0.000136  # объем кормовой части
S_f = kk.pi * d ** 2 / 4  # площадь миделевого сечения корпуса
Ff = 0.3619  # площадь обтекаемой потоком поверхности корпуса (без донного среза)
Sf = kk.pi * (d ** 2) / 4  # площадь миделя
f_t = 0.045216

l_cil_nas = 45.5 / 8  # относительное удлинение циллиндрической части АД насадка
l_nos_nas = 12.5 / 8  # относительное удлинение носовой части АД насадка
l_korp = 1650 / 72 # относительное удлинение корпуса
l_nos = 97.4 / 72  # относительное удлинение носовой части корпуса
l_cil = l_korp - 160 / 72 - l_nos  # относительное удлинение циллиндрической части корпуса

Omega_con1 = 3.64  # угол полураствора конуса носовой части (градусы)
Omega_con2 = 1.58  # предельный угол полураствора конуса носовой части (градусы)

l_con = 565 / 72  # относительное удлинение конической части НЧ без сферического затупления
l_con_max = 1305 / 72  # максимальное относ. удлинение конческой части НЧ без сферического затупления

"""L_kr = 1.46  # размах крыла
L_k_kr = L_kr - d"""

# геометрические параметры АД насадка
r_n = 4  # радиус сферического наконечника АД насадка
r1_n = 32.86  # ражиус сферической части носовой части
d_n = 14  # диаметр
d1_n = 8  # диаметр циллиндрической части АД насадка
D = 72  # миделевый диаметр в мм

r_ = 2 * r_n / d_n
r_1 = 2 * r1_n / D
d_ = d1_n / D


tan_05_kr = 0.307  # тангенс среднего угла крыла
l_kr = 1.46  # относительное удлинение крыла
L_k_kr = 0.146 - 0.072  # размах консолей крыла
L_kr = 0.146
b_kr = 0.1061  # ширина крыла у корпуса
b_a_kr = 0.101  # САХ консоли крыльев
b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
x_b_a_kr = 1.183  # координата начала САХ крыла
x_b_kr = 1.178  # координата начала бортовой хорды крыла
a_kr = 0.098  # ширина крыла у корпуса (меньшая грань трапеции)
c_kr = 0.03  # относительная толщина профиля крыла
L_hv_kr = 0.015  # Расстояние от конца бортовой хорды оперения до кормового среза корпуса
L1_kr = 1.2311  # Расстояние от носика корпуса до серидины бортовой хорды оперения
koef_kr = 1 / (1 - a_kr / b_kr)  # коэффициент перехода от ромбовидного к шестиугольному профилю крыла

tan_05_op = -0.012  # тангенс среднего угла оперения
l_op = 7.888  # относительное удлинение оперения l^2/S
L_k_op = 0.25 - 0.72  # размах консолей оперения
b_op = 0.032  # ширина оперения у корпуса
b_a_op = 0.029  # САХ консоли оперения
b__a_op = b_a_op / L_f  # относительное значение САХ консоли оперения
x_c_pl_ba = 0.349  # координата цт площади передних консолей (середина САХ консолей)
x_b_a_op = 0.334  # коордигата начала САХ консоли опрения
x_b_op = 0.334  # координата начала бортовой хорды оперения
a_op = 0.015  # ширина оперения у корпуса
c_op = 0.03  # относительная толщина профиля оперения
L_hv_op = 0.933  # Расстояние от конца бортовой хорды оперения до кормового среза корпуса
L1_op = 0.35  # Расстояние от носика корпуса до серидины бортовой хорды оперения
koef_op = 1 / (1 - a_op / b_op)  # коэффициент перехода от ромбовидного к шестиугольному профилю оперения
print(koef_kr, koef_op, "koef")

b_ak_kr = 0.094  # средняя ад хорда консоли крыла
x_otn_op_kr = 0.846 / b_ak_kr  # относительное расстояние между оперением и средней хордой крыльев
S__f = kk.pi * (d ** 2) / 4 / (kk.pi * (d ** 2) / 4)  # относительная площадь корпуса
S__op = 0.0078575 / (kk.pi * (d ** 2) / 4)  # относительная площадь передних несущих поверхностей
S__kr = 0.0146 / (kk.pi * (d ** 2) / 4)  # относительная площадь задних несущих поверхностей
S_kr = 0.0146
S_op = 0.0078575

D_kr = 0.072 / 0.146  # отерсительный диаметр корпуса
nu_k_kr = 1.322  # относительное сужение консоли задней несущей поверхности
c_const_kr = (4 + 1 / nu_k_kr) * (1 + 8 * D_kr ** 2)  # коэффициент формы эпюры погонной нагрузки для крыльев
D_op = 0.072 / 0.25  # относительный диаметр корпуса
nu_k_op = 1.069  # относительное сужение консоли передней несущей поверхности
c_const_op = (4 + 1 / nu_k_op) * (1 + 8 * D_op ** 2)  # коэффициент формы эпюры погонной нагрузки для оперения

l_zat = 23.49 / 62.59  # относительное удлинение затупления носовой части
r_ = 2 * 0.0329 / d

l_nos_ = (l_nos - r_ / 2) / kk.sqrt(1 - r_)  # относительное удлинение носовой части без затупления
print(l_nos_, l_nos)
teta = kk.atan((1 - r_) / (l_nos - r_ / 2))  # угол наклона образующей носовой части (конуса)
print(teta * 180 / kk.pi)

h = 6.3 * 10 ** -6  # Примерная высота бугоров на поверхности корпуса (в зависимости от класса чистоты) (для 7-го)

nu_korm = 0.047 / d  # относительное сужение кормовой части
l_korm = 0.0463 / d  # относительное удлинение кормовой части

dt = 10 ** -3
vertel = Target2(6000, 4100, 220)
igla = Rocket2(vertel.y, vertel.x)
tt_1, f_yy_exp_1, f_xx_exp_1, cyy_sum_1, cxx_sum_1, da_ball_1, x_ffa_op_1, x_ffa_kr_1, x_ffa_f_1, x_ffa_1, machh_1, xx_1, yy_1, xt_1, yt_1, alff_1, deltt_1, delt_r_1, cyy_op_1, cyy_kr_1, cyy1_alf_1, cyy_korm_1, cyy_nos_1, cxx_nos_1, cxx_0f_1, cxx_korm_1, cxx_tr_1, cxx_kr_pr_1, cxx_kr_vol_1, cxx_0_kr_1, cxx_op_pr_1, cxx_op_vol_1, cxx_0_op_1, cxx_0_op_1, cxx_0_1, cxx_0_kr_1, cxx_0f_1, m_zz_wz_1, m_zz_wz_kr_1, m_zz_wz_op_1, m_zz_wz_f_1, f_xx_kr_exp_1, f_yy_kr_exp_1 = igla.navigation(vertel)

v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
param_x = 6000
param_y = 3500
param_v = v_target_max
num_x = 2
num_y = 2
num_v = 2
d_x = (param_x - 500) / (num_x - 1)
d_y = (param_y - 10) / (num_y - 1)
d_v = (v_target_max - v_target_min) / (num_v - 1)
# vertel = Target(param_x, param_y)
param_y_0 = param_y
param_x_0 = param_x
start_time = time.time()

koord_iniz_x = np.zeros((num_x, num_y))
koord_iniz_y = np.zeros((num_x, num_y))
koord_nach_x = np.zeros((num_x, num_y))
koord_nach_y = np.zeros((num_x, num_y))

dt = 10 ** -2
koord_v = [] * 10



"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')"""
t_in = np.zeros((num_x, num_y))
dict_igla = []
# for k_v in tqdm(range(0, num_v)):
pp1 = []
pp2 = []
pp3 = []
for k_v in tqdm(range(0, num_v)):
    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))
    t_in = np.zeros((num_x, num_y))
    n_max_ = np.zeros((num_x, num_y))
    n_kon = np.zeros((num_x, num_y))
    param_x = param_x_0
    koord_v = []
    for i in range(0, num_x):
        param_y = param_y_0
        for j in range(0, num_y):
            koord_v.append(param_v)
            vertel = Target(param_x, param_y, param_v)
            start_time_var = time.time()
            igla = Rocket(vertel.y, vertel.x)
            # print(param_x, param_y)
            koord_nach_x[i][j] = param_x
            koord_nach_y[i][j] = param_y
            koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][j], rast, fx_t, fy_t, m_shar = igla.navigation(vertel)
            #print("Time", time.time() - start_time_var)
            pp1.append(fx_t)
            pp2.append(fy_t)
            pp3.append(m_shar)
            param_y -= d_y
        param_x -= d_x
    """x1 = [koord_nach_x[num_x - 1][0], koord_nach_x[0][0], koord_nach_x[0][num_y - 1],
          koord_nach_x[num_x - 1][num_y - 1], koord_nach_x[num_x - 1][0]]
    y1 = [koord_nach_y[num_x - 1][0], koord_nach_y[0][0], koord_nach_y[0][num_y - 1],
          koord_nach_y[num_x - 1][num_y - 1], koord_nach_y[num_x - 1][0]]
    z1 = [koord_v[0], koord_v[1], koord_v[2], koord_v[3], koord_v[0]]
    x2 = [koord_iniz_x[num_x - 1][0], koord_iniz_x[0][0], koord_iniz_x[0][num_y - 1],
          koord_iniz_x[num_x - 1][num_y - 1], koord_iniz_x[num_x - 1][0]]
    y2 = [koord_iniz_y[num_x - 1][0], koord_iniz_y[0][0], koord_iniz_y[0][num_y - 1],
          koord_iniz_y[num_x - 1][num_y - 1], koord_iniz_y[num_x - 1][0]]
    z2 = [koord_v[0], koord_v[1], koord_v[2], koord_v[3], koord_v[0]]

    ax.plot(z2, x2, y2, 'g')
    ax.plot(z1, x1, y1, 'r')"""
    param_v -= d_v
    # print(t_in, n_max_, n_kon)
    paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
    dict_igla.append(paramz)
"""
v_target_min = -120
v_target_max = -400
param_x = 9000
param_y = 3500
param_v = v_target_max
num_x = 2
num_y = 2
num_v = 9
d_x = (param_x - 2000) / (num_x - 1)
d_y = (param_y - 10) / (num_y - 1)
d_v = (v_target_max - v_target_min) / (num_v - 1)
# vertel = Target(param_x, param_y)
param_y_0 = param_y
param_x_0 = param_x
start_time = time.time()

koord_v = [] * 10
"""
"""
for k_v in range(0, num_v):
    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))
    t_in = np.zeros((num_x, num_y))
    n_max_ = np.zeros((num_x, num_y))
    n_kon = np.zeros((num_x, num_y))
    param_x = param_x_0
    koord_v = []
    for i in range(0, num_x):
        param_y = param_y_0
        for j in range(0, num_y):
            koord_v.append(param_v)
            vertel = Target(param_x, param_y, param_v)
            start_time_var = time.time()
            igla = Rocket(vertel.y, vertel.x)
            # print(param_x, param_y)
            koord_nach_x[i][j] = param_x
            koord_nach_y[i][j] = param_y
            koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][j], rast = igla.navigation(vertel)
            # print("Time", time.time() - start_time_var)
            param_y -= d_y
            print(param_y, param_x, param_v)
        param_x -= d_x
    x1 = [-koord_nach_x[num_x - 1][0], -koord_nach_x[0][0], -koord_nach_x[0][num_y - 1],
          -koord_nach_x[num_x - 1][num_y - 1], -koord_nach_x[num_x - 1][0]]
    y1 = [koord_nach_y[num_x - 1][0], koord_nach_y[0][0], koord_nach_y[0][num_y - 1],
          koord_nach_y[num_x - 1][num_y - 1], koord_nach_y[num_x - 1][0]]
    # z1 = [koord_v[0], koord_v[1], koord_v[2], koord_v[3], koord_v[0]]
    x2 = [-koord_iniz_x[num_x - 1][0], -koord_iniz_x[0][0], -koord_iniz_x[0][num_y - 1],
          -koord_iniz_x[num_x - 1][num_y - 1], -koord_iniz_x[num_x - 1][0]]
    y2 = [koord_iniz_y[num_x - 1][0], koord_iniz_y[0][0], koord_iniz_y[0][num_y - 1],
          koord_iniz_y[num_x - 1][num_y - 1], koord_iniz_y[num_x - 1][0]]
    z2 = [-koord_v[0], -koord_v[1], -koord_v[2], -koord_v[3], -koord_v[0]]

    ax.plot(z2, x2, y2, 'g')
    ax.plot(z2, x1, y1, 'r')
    param_v -= d_v
print("Time", time.time() - start_time)
plt.show()
print(dict_igla[0]['t_kon'])"""

krit_max = 0

m_0_w = 7.3
m_nos = 1.3

m_0_w_0 = m_0_w

kritt = []
l_cil += l_nos
l_nos = 1
l_cil -= l_nos
lambdd = []
m_0_w = m_0_w - m_nos
m_nos = m_nos / W_nos
W_nos = W_nos / L_nos

L_nos = l_nos * 0.072
W_nos = W_nos * L_nos
m_nos = m_nos * W_nos
m_0_w = m_0_w + m_nos

krit_l = []
print("l_nos")
for k_l_nos in tqdm(range(0, 19)):
    l_cil += l_nos # относительное удлинение циллиндрической части корпуса
    if i >= 1:
        l_nos += 0.1
    else:
        l_nos = l_nos
    l_cil -= l_nos
    m_0_w = m_0_w - m_nos
    m_nos = m_nos / W_nos
    W_nos = W_nos / L_nos

    L_nos = l_nos * 0.072
    W_nos = W_nos * L_nos
    m_nos = m_nos * W_nos
    m_0_w = m_0_w + m_nos

    l_cil += l_korm

    l_korm = 0.278
    krit_korm = []
    lambd_k = []
    print("l_korm:")
    for k_l_korm in tqdm(range(0, 19)):
        l_cil += l_korm
        if i >= 1:
            l_korm += 0.05
        else:
            l_korm = l_korm
        l_cil -= l_korm

        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))

        t_in = np.zeros((num_x, num_y))
        dict_igla = []

        v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
        v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
        param_x = 6000
        param_y = 3500
        param_v = v_target_max
        num_x = 2
        num_y = 2
        num_v = 2
        d_x = (param_x - 500) / (num_x - 1)
        d_y = (param_y - 10) / (num_y - 1)
        d_v = (v_target_max - v_target_min) / (num_v - 1)
        param_y_0 = param_y
        param_x_0 = param_x
        start_time = time.time()

        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))

        koord_v = [] * 10
        krit_ind = 0
        k_paramz = 0
        for k_v in range(0, num_v):
            koord_iniz_x = np.zeros((num_x, num_y))
            koord_iniz_y = np.zeros((num_x, num_y))
            koord_nach_x = np.zeros((num_x, num_y))
            koord_nach_y = np.zeros((num_x, num_y))
            t_in = np.zeros((num_x, num_y))
            n_max_ = np.zeros((num_x, num_y))
            n_kon = np.zeros((num_x, num_y))
            param_x = param_x_0
            koord_v = []
            for i in range(0, num_x):
                param_y = param_y_0
                for j in range(0, num_y):
                    koord_v.append(param_v)
                    vertel = Target(param_x, param_y, param_v)
                    start_time_var = time.time()
                    igla = Rocket(vertel.y, vertel.x)
                    koord_nach_x[i][j] = param_x
                    koord_nach_y[i][j] = param_y
                    koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][j], rast, Sfx, Sfy, m_shar = igla.navigation(
                        vertel)
                    krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                    k_paramz += 1
                    param_y -= d_y
                param_x -= d_x
            param_v -= d_v
        krit_korm.append(krit_ind / k_paramz)
        lambd_k.append(l_korm)
    for i in range(len(krit_korm)):
        if krit_korm[i] >= krit_max:
            l_korm_par = lambd_k[i]
            krit_max = krit_korm[i]
    l_cil += l_korm
    l_cil -= l_korm_par

    l_korm = l_korm_par
    print("l_korm", l_korm, "l_cil", l_cil)
    krit_l.append(krit_ind / k_paramz)
    lambdd.append(l_nos)


    for i in range(len(krit_l)):
        if krit_l[i] >= krit_max:
            l_nos_par = 1 + 0.1 * i
            krit_max = krit_l[i]
"""
plt.plot(lambdd, krit_l)
plt.grid(True)
plt.show()"""

krit_max = 0
for i in range(len(krit_l)):
    if krit_l[i] >= krit_max:
        l_nos_par = 1 + 0.1 * i
        krit_max = krit_l[i]

l_cil += l_nos
l_cil -= l_nos_par
m_0_w = m_0_w - m_nos
m_nos = m_nos / W_nos
W_nos = W_nos / L_nos

L_nos = l_nos_par * 0.072
W_nos = W_nos * L_nos
m_nos = m_nos * W_nos
m_0_w = m_0_w + m_nos
l_nos = l_nos_par

print("l_cil", l_cil, "L_nos", L_nos, "m_0_w", m_0_w, "m_nos", m_nos, "W_nos", W_nos)

krit_korm = []
lambd_k = []
l_korm = 0.278

print("l_korm:")
for k_l_korm in tqdm(range(0, 19)):
    l_cil += l_korm
    if i >= 1:
        l_korm += 0.05
    else:
        l_korm = l_korm
    l_cil -= l_korm

    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 9
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10

    t_in = np.zeros((num_x, num_y))
    dict_igla = []
    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 2
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10
    krit_ind = 0
    k_paramz = 0
    for k_v in range(0, num_v):
        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))
        t_in = np.zeros((num_x, num_y))
        n_max_ = np.zeros((num_x, num_y))
        n_kon = np.zeros((num_x, num_y))
        param_x = param_x_0
        koord_v = []
        for i in range(0, num_x):
            param_y = param_y_0
            for j in range(0, num_y):
                koord_v.append(param_v)
                vertel = Target(param_x, param_y, param_v)
                start_time_var = time.time()
                igla = Rocket(vertel.y, vertel.x)
                # print(param_x, param_y)
                koord_nach_x[i][j] = param_x
                koord_nach_y[i][j] = param_y
                koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][
                    j], rast, Sfx, Sfy, m_shar = igla.navigation(
                    vertel)
                #print("Time", time.time() - start_time_var)
                krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                k_paramz += 1
                param_y -= d_y
            param_x -= d_x

        param_v -= d_v
        #print(t_in, n_max_, n_kon)
        paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
        dict_igla.append(paramz)
    krit_korm.append(krit_ind / k_paramz)
    lambd_k.append(l_korm)


krit_max = 0
for i in range(len(krit_korm)):
    if krit_korm[i] >= krit_max:
        l_korm_par = lambd_k[i]
        krit_max = krit_korm[i]
l_cil += l_korm
l_cil -= l_korm_par

l_korm = l_korm_par
print("l_korm", l_korm, "l_cil", l_cil)


"""plt.plot(lambdd, krit_l)
plt.grid(True)
plt.show()
plt.plot(lambd_k, krit_korm)
plt.grid(True)
plt.show()"""

L_kr = kk.sqrt(S_kr * l_kr)
tan_pk = 0.577
krit_kr = []
lambdd_kr = []

l_kr_kons = l_kr - d
b_0_kr = 0.12
nu_kr = (S_kr / (L_kr * b_0_kr) * 2 - 1) ** (-1)
S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d/L_kr) # S консоли крыла (const)

l_kr = 1
L_kr = kk.sqrt(S_kr * l_kr)
S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d/L_kr) # S консоли крыла (const)
L_k_kr = L_kr - d
b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
x_b_kr += b_kr
L1_kr += b_kr / 2
b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
x_b_a_kr += b_a_kr
b_a_kr = 4/3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1)**2) # САХ консоли крыльев
b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
x_b_a_kr -= b_a_kr   # координата начала САХ крыла
x_b_kr -= b_kr  # координата начала бортовой хорды крыла
L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

krit_l_kr = []
lambdd_kr = []
print("l_kr")
for k_l_kr in tqdm(range(0, 10)):
    L_kr = kk.sqrt(S_kr * l_kr)
    L_k_kr = L_kr - d
    S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d / L_kr)  # S консоли крыла (const)

    b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
    nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
    l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
    tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
    x_b_kr += b_kr
    L1_kr += b_kr / 2
    b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
    x_b_a_kr += b_a_kr
    b_a_kr = 4 / 3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1) ** 2)  # САХ консоли крыльев
    b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
    x_b_a_kr -= b_a_kr  # координата начала САХ крыла
    x_b_kr -= b_kr  # координата начала бортовой хорды крыла
    L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 9
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10

    t_in = np.zeros((num_x, num_y))
    dict_igla = []
    # for k_v in tqdm(range(0, num_v)):
    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 2
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10
    krit_ind = 0
    k_paramz = 0
    for k_v in range(0, num_v):
        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))
        t_in = np.zeros((num_x, num_y))
        n_max_ = np.zeros((num_x, num_y))
        n_kon = np.zeros((num_x, num_y))
        param_x = param_x_0
        koord_v = []
        for i in range(0, num_x):
            param_y = param_y_0
            for j in range(0, num_y):
                koord_v.append(param_v)
                vertel = Target(param_x, param_y, param_v)
                start_time_var = time.time()
                igla = Rocket(vertel.y, vertel.x)

                koord_nach_x[i][j] = param_x
                koord_nach_y[i][j] = param_y
                koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][
                    j], rast, Sfx, Sfy, m_shar = igla.navigation(
                    vertel)

                krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                k_paramz += 1
                param_y -= d_y
            param_x -= d_x

        param_v -= d_v

        paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
        dict_igla.append(paramz)
    krit_l_kr.append(krit_ind / k_paramz)
    lambdd_kr.append(l_kr)
    l_kr += 0.05


for i in range(len(krit_l_kr)):
    if krit_l_kr[i] >= krit_max:
        l_kr = lambdd_kr[i]
        krit_max = krit_l_kr[i]

L_kr = kk.sqrt(S_kr * l_kr)
S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d/L_kr) # S консоли крыла (const)
L_k_kr = L_kr - d
b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
x_b_kr += b_kr
L1_kr += b_kr / 2
b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
x_b_a_kr += b_a_kr
b_a_kr = 4/3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1)**2) # САХ консоли крыльев
b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
x_b_a_kr -= b_a_kr   # координата начала САХ крыла
x_b_kr -= b_kr  # координата начала бортовой хорды крыла
L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

print("L_kr", L_kr, "S_k_kr", S_k_kr, "L_k_kr", L_k_kr, "b_0_kr", b_0_kr, "nu_k_kr", nu_k_kr, "l_k_kr", l_k_kr)
print("tan_05_kr", tan_05_kr, "b_kr", b_kr, "b_a_kr", b_a_kr,  "b__a_kr", b__a_kr, "x_b_a_kr", x_b_a_kr, "x_b_kr", x_b_kr, "L1_kr", L1_kr)

nu_kr = 1

S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d/L_kr) # S консоли крыла (const)
b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
x_b_kr += b_kr
L1_kr += b_kr / 2
b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
x_b_a_kr += b_a_kr
b_a_kr = 4/3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1)**2) # САХ консоли крыльев
b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
x_b_a_kr -= b_a_kr   # координата начала САХ крыла
x_b_kr -= b_kr  # координата начала бортовой хорды крыла
L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения


krit_nu_kr = []
nuu_kr = []
print("nu_kr:")
for k_nu_kr in tqdm(range(0, 16)):

    S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d / L_kr)  # S консоли крыла (const)

    b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
    nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
    l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
    tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
    x_b_kr += b_kr
    L1_kr += b_kr / 2
    b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
    x_b_a_kr += b_a_kr
    b_a_kr = 4 / 3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1) ** 2)  # САХ консоли крыльев
    b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
    x_b_a_kr -= b_a_kr  # координата начала САХ крыла
    x_b_kr -= b_kr  # координата начала бортовой хорды крыла
    L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 9
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))


    koord_v = [] * 10

    t_in = np.zeros((num_x, num_y))
    dict_igla = []
    # for k_v in tqdm(range(0, num_v)):
    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 2
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))


    koord_v = [] * 10
    krit_ind = 0
    k_paramz = 0
    for k_v in range(0, num_v):
        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))
        t_in = np.zeros((num_x, num_y))
        n_max_ = np.zeros((num_x, num_y))
        n_kon = np.zeros((num_x, num_y))
        param_x = param_x_0
        koord_v = []
        for i in range(0, num_x):
            param_y = param_y_0
            for j in range(0, num_y):
                koord_v.append(param_v)
                vertel = Target(param_x, param_y, param_v)
                start_time_var = time.time()
                igla = Rocket(vertel.y, vertel.x)
                # print(param_x, param_y)
                koord_nach_x[i][j] = param_x
                koord_nach_y[i][j] = param_y
                koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][
                    j], rast, Sfx, Sfy, m_shar = igla.navigation(
                    vertel)
                # print("Time", time.time() - start_time_var)
                krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                k_paramz += 1
                param_y -= d_y
            param_x -= d_x

        param_v -= d_v
        # print(t_in, n_max_, n_kon)
        paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
        dict_igla.append(paramz)
    krit_nu_kr.append(krit_ind / k_paramz)
    nuu_kr.append(nu_kr)
    nu_kr += 0.1

"""plt.plot(nuu_kr, krit_nu_kr)
plt.grid(True)
plt.show()"""

for i in range(len(krit_nu_kr)):
    if krit_nu_kr[i] >= krit_max:
        nu_kr = nuu_kr[i]
        krit_max = krit_nu_kr[i]
S_k_kr = S_kr * (1 - (nu_kr - 1) / (nu_kr + 1) * d / L_kr) * (1 - d/L_kr) # S консоли крыла (const)
b_0_kr = S_kr / L_kr * (2 * nu_kr) / (nu_kr + 1)
nu_k_kr = nu_kr - d / L_kr * (nu_kr - 1)
l_k_kr = l_kr * (1 - d / L_kr) / (1 + (nu_kr - 1) / (nu_kr + 1) * d / L_kr)
tan_05_kr = 0 + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)
x_b_kr += b_kr
L1_kr += b_kr / 2
b_kr = S_k_kr / L_k_kr * (2 * nu_k_kr) / (nu_k_kr + 1)  # ширина крыла у корпуса
x_b_a_kr += b_a_kr
b_a_kr = 4/3 * S_kr / L_kr * (1 - nu_kr / (nu_kr + 1)**2) # САХ консоли крыльев
b__a_kr = b_a_kr / L_f  # относительное значение САХ консоли крыльев
x_b_a_kr -= b_a_kr   # координата начала САХ крыла
x_b_kr -= b_kr  # координата начала бортовой хорды крыла
L1_kr -= b_kr / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

print("S_k_kr", S_k_kr, "L_k_kr", L_k_kr, "b_0_kr", b_0_kr, "nu_k_kr", nu_k_kr, "l_k_kr", l_k_kr)
print("tan_05_kr", tan_05_kr, "b_kr", b_kr, "b_a_kr", b_a_kr, "b__a_kr", b__a_kr, "x_b_a_kr", x_b_a_kr, "x_b_kr", x_b_kr, "L1_kr", L1_kr)

# print("Time", time.time() - start_time)

L_op = kk.sqrt(S_op * l_op)
tan_pk = 0.577

l_op_kons = l_op - d
b_0_op = 0.12
nu_op = (S_op / (L_op * b_0_op) * 2 - 1) ** (-1)
S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d/L_op) # S консоли крыла (const)

l_op = 2
L_op = kk.sqrt(S_op * l_op)
S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d/L_op) # S консоли крыла (const)
L_k_op = L_op - d
b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
nu_k_op = nu_op - d / L_op * (nu_op - 1)
l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)

tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

b_a_op = 4/3 * S_op / L_op * (1 - nu_op / (nu_op + 1)**2) # САХ консоли крыльев
b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев


krit_l_op = []
lambdd_op = []

print("l_op:")
for k_l_op in tqdm(range(0, 12)):
    L_op = kk.sqrt(S_op * l_op)
    S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d / L_op)  # S консоли крыла (const)
    L_k_op = L_op - d
    b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
    nu_k_op = nu_op - d / L_op * (nu_op - 1)
    l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)

    tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

    b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

    b_a_op = 4 / 3 * S_op / L_op * (1 - nu_op / (nu_op + 1) ** 2)  # САХ консоли крыльев
    b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев


    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 9
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))
    koord_v = [] * 10
    t_in = np.zeros((num_x, num_y))
    dict_igla = []
    # for k_v in tqdm(range(0, num_v)):
    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 2
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))


    koord_v = [] * 10
    krit_ind = 0
    k_paramz = 0
    for k_v in range(0, num_v):
        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))
        t_in = np.zeros((num_x, num_y))
        n_max_ = np.zeros((num_x, num_y))
        n_kon = np.zeros((num_x, num_y))
        param_x = param_x_0
        koord_v = []
        for i in range(0, num_x):
            param_y = param_y_0
            for j in range(0, num_y):
                koord_v.append(param_v)
                vertel = Target(param_x, param_y, param_v)
                start_time_var = time.time()
                igla = Rocket(vertel.y, vertel.x)
                # print(param_x, param_y)
                koord_nach_x[i][j] = param_x
                koord_nach_y[i][j] = param_y
                koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][
                    j], rast, Sfx, Sfy, m_shar = igla.navigation(
                    vertel)
                # print("Time", time.time() - start_time_var)
                krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                k_paramz += 1
                param_y -= d_y
            param_x -= d_x

        param_v -= d_v
        # print(t_in, n_max_, n_kon)
        paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
        dict_igla.append(paramz)
    krit_l_op.append(krit_ind / k_paramz)
    lambdd_op.append(l_op)
    l_op += 0.5

"""plt.plot(lambdd_kr, krit_l_kr)
plt.grid(True)
plt.show()"""


for i in range(len(krit_l_op)):
    if krit_l_op[i] >= krit_max:
        l_op = lambdd_op[i]
        krit_max = krit_l_op[i]
l_op = 5
L_op = kk.sqrt(S_op * l_op)

S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d/L_op) # S консоли крыла (const)
L_k_op = L_op - d
b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
nu_k_op = nu_op - d / L_op * (nu_op - 1)
l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)
tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

b_a_op = 4/3 * S_op / L_op * (1 - nu_op / (nu_op + 1)**2) # САХ консоли крыльев
b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев
x_b_a_op -= b_a_op   # координата начала САХ крыла
x_b_op -= b_op  # координата начала бортовой хорды крыла
L1_op -= b_op / 2  # Расстояние от носика корпуса до серидины бортовой хорды оперения

print("L_op", L_op, "S_k_op", S_k_op, "L_k_op", L_k_op, "b_0_op", b_0_op, "nu_k_op", nu_k_op, "l_k_op", l_k_op)
print("tan_05_op", tan_05_op, "b_op", b_op, "b_a_op", b_a_op,  "b__a_op", b__a_op, "x_b_a_op", x_b_a_op, "x_b_op", x_b_op, "L1_op", L1_op)

nu_op = 1

S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d/L_op) # S консоли крыла (const)
L_k_op = L_op - d
b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
nu_k_op = nu_op - d / L_op * (nu_op - 1)
l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)

tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

b_a_op = 4/3 * S_op / L_op * (1 - nu_op / (nu_op + 1)**2) # САХ консоли крыльев
b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев


krit_nu_op = []
nuu_op = []
print("nu_op:")
for k_nu_op in tqdm(range(0, 8)):

    S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d / L_op)  # S консоли крыла (const)
    L_k_op = L_op - d
    b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
    nu_k_op = nu_op - d / L_op * (nu_op - 1)
    l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)

    tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

    b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

    b_a_op = 4 / 3 * S_op / L_op * (1 - nu_op / (nu_op + 1) ** 2)  # САХ консоли крыльев
    b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев


    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 9
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10

    t_in = np.zeros((num_x, num_y))
    dict_igla = []
    # for k_v in tqdm(range(0, num_v)):
    v_target_max = 320  # скорость и направление (+ догонный курс, - навстречу) цели
    v_target_min = 120  # скорость и направление (+ догонный курс, - навстречу) цели
    param_x = 6000
    param_y = 3500
    param_v = v_target_max
    num_x = 2
    num_y = 2
    num_v = 2
    d_x = (param_x - 500) / (num_x - 1)
    d_y = (param_y - 10) / (num_y - 1)
    d_v = (v_target_max - v_target_min) / (num_v - 1)
    # vertel = Target(param_x, param_y)
    param_y_0 = param_y
    param_x_0 = param_x
    start_time = time.time()

    koord_iniz_x = np.zeros((num_x, num_y))
    koord_iniz_y = np.zeros((num_x, num_y))
    koord_nach_x = np.zeros((num_x, num_y))
    koord_nach_y = np.zeros((num_x, num_y))

    koord_v = [] * 10
    krit_ind = 0
    k_paramz = 0
    for k_v in range(0, num_v):
        koord_iniz_x = np.zeros((num_x, num_y))
        koord_iniz_y = np.zeros((num_x, num_y))
        koord_nach_x = np.zeros((num_x, num_y))
        koord_nach_y = np.zeros((num_x, num_y))
        t_in = np.zeros((num_x, num_y))
        n_max_ = np.zeros((num_x, num_y))
        n_kon = np.zeros((num_x, num_y))
        param_x = param_x_0
        koord_v = []
        for i in range(0, num_x):
            param_y = param_y_0
            for j in range(0, num_y):
                koord_v.append(param_v)
                vertel = Target(param_x, param_y, param_v)
                start_time_var = time.time()
                igla = Rocket(vertel.y, vertel.x)

                koord_nach_x[i][j] = param_x
                koord_nach_y[i][j] = param_y
                koord_iniz_x[i][j], koord_iniz_y[i][j], t_in[i][j], n_max_[i][j], n_kon[i][
                    j], rast, Sfx, Sfy, m_shar = igla.navigation(vertel)

                krit_ind += krit_maxi(m_0_w, Sfy, Sfx)
                k_paramz += 1
                param_y -= d_y
            param_x -= d_x

        param_v -= d_v

        paramz = {'t_kon': t_in, 'n_max': n_max_, 'n_kon': n_kon}
        dict_igla.append(paramz)
    krit_nu_op.append(krit_ind / k_paramz)
    nuu_op.append(nu_op)
    nu_op += 0.05


for i in range(len(krit_nu_op)):
    if krit_nu_op[i] >= krit_max:
        nu_kr = nuu_op[i]
        krit_max = krit_nu_op[i]
"""print(krit_max)
print(l_nos, l_korm, l_kr, nu_kr, l_op, nu_op)
print(krit_nu_op)"""
S_k_op = S_op * (1 - (nu_op - 1) / (nu_op + 1) * d / L_op) * (1 - d / L_op)  # S консоли крыла (const)
L_k_op = L_op - d
b_0_op = S_op / L_op * (2 * nu_op) / (nu_op + 1)
nu_k_op = nu_op - d / L_op * (nu_op - 1)
l_k_op = l_op * (1 - d / L_op) / (1 + (nu_op - 1) / (nu_op + 1) * d / L_op)

tan_05_op = 0 - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)

b_op = S_k_op / L_k_op * (2 * nu_k_op) / (nu_k_op + 1)  # ширина крыла у корпуса

b_a_op = 4 / 3 * S_op / L_op * (1 - nu_op / (nu_op + 1) ** 2)  # САХ консоли крыльев
b__a_op = b_a_op / L_f  # относительное значение САХ консоли крыльев

print("S_k_op", S_k_op, "L_k_op", L_k_op, "b_0_op", b_0_op, "nu_k_op", nu_k_op, "l_k_op", l_k_op)
print("tan_05_op", tan_05_op, "b_op", b_op, "b_a_op", b_a_op, "b__a_op", b__a_op)

fin_ch += 10
dt = 10 ** -3
vertel = Target2(6000, 4100, 220)
igla = Rocket2(vertel.y, vertel.x)
tt_2, f_yy_exp_2, f_xx_exp_2, cyy_sum_2, cxx_sum_2, da_ball_2, x_ffa_op_2, x_ffa_kr_2, x_ffa_f_2, x_ffa_2, machh_2, xx_2, yy_2, xt_2, yt_2, alff_2, deltt_2, delt_r_2, cyy_op_2, cyy_kr_2, cyy1_alf_2, cyy_korm_2, cyy_nos_2, cxx_nos_2, cxx_0f_2, cxx_korm_2, cxx_tr_2, cxx_kr_pr_2, cxx_kr_vol_2, cxx_0_kr_2, cxx_op_pr_2, cxx_op_vol_2, cxx_0_op_2, cxx_0_op_2, cxx_0_2, cxx_0_kr_2, cxx_0f_2, m_zz_wz_2, m_zz_wz_kr_2, m_zz_wz_op_2, m_zz_wz_f_2, f_xx_kr_exp_2, f_yy_kr_exp_2 = igla.navigation(vertel)

plt.plot(tt_2, f_xx_kr_exp_2, linestyle = "-", color = 'g', label = 'Сила лобового сопротивления')
plt.plot(tt_2, f_yy_kr_exp_2, linestyle = "-", color = 'r', label = 'Подъемная сила')
plt.legend()
plt.grid(True)
plt.xlabel('t, [c]')
plt.ylabel('F, [H]')
plt.show()

#plt.plot(tt_1, f_yy_exp_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, f_xx_exp_1, linestyle = "--", color = 'g')
plt.plot(tt_2, f_yy_exp_2, linestyle = "-", color = 'r', label = 'Подъемная сила')
plt.plot(tt_2, f_xx_exp_2, linestyle = "-", color = 'g', label = 'Сила лобового сопротивления')
plt.grid(True)
plt.legend()
plt.xlabel('t, [c]')
plt.ylabel('F, [H]')
plt.show()

#plt.plot(tt_1, cyy_sum_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cxx_sum_1, linestyle = "--", color = 'g')
plt.plot(tt_2, cyy_sum_2, linestyle = "-", color = 'r', label = 'Коэффициент подъемной силы')
plt.plot(tt_2, cxx_sum_2, linestyle = "-", color = 'g', label = 'Коэффициент лобового сопротивления')
plt.grid(True)
plt.xlabel('t, [c]')
plt.ylabel('Cy, Cx')
plt.legend()
plt.show()

#plt.plot(tt_1, da_ball_1, "--")
plt.plot(tt_2, da_ball_2)
plt.ylabel('delta/alf')
plt.xlabel('t, [с]')
plt.grid(True)
plt.show()

#plt.plot(x_ffa_op_1, tt_1, linestyle = "--", color = 'r')
#plt.plot(x_ffa_kr_1, tt_1, linestyle = "--", color = 'g')
#plt.plot(x_ffa_f_1, tt_1, linestyle = "--", color = 'b')
#plt.plot(x_ffa_1, tt_1,linestyle = "--", color = 'y')
plt.plot(x_ffa_op_2, tt_2, linestyle = "-", color = 'r', label = 'Координата фокуса оперения')
plt.plot(x_ffa_kr_2, tt_2,linestyle = "-", color = 'g', label = 'Координата фокуса крыльев')
plt.plot(x_ffa_f_2, tt_2, linestyle = "-", color = 'b', label = 'Координата фокуса корпуса')
plt.plot(x_ffa_2, tt_2, linestyle = "-", color ='y', label = 'Координата фокуса Ла')
plt.grid(True)
plt.axis([-0.3, 1.7, 0, 16])
plt.ylabel('t, [с]')
plt.xlabel('x_fa, [м]')
plt.legend()
plt.show()

#plt.plot(tt_1, machh_1, "--", label = 'первоначального варианта')
plt.plot(tt_2, machh_2, label = 'финального варианта')
plt.grid(True)
plt.ylabel('M')
plt.xlabel('t, [с]')
plt.show()


#plt.plot(xx_1, yy_1, linestyle = "--", color = 'r')
#plt.plot(xt_1, yt_1, linestyle = "--", color = 'g')
plt.plot(xx_2, yy_2, linestyle = "-", color = 'r', label = 'Траектория ракеты')
plt.plot(xt_2, yt_2, linestyle = "-", color = 'g', label = 'Траектория цели')
plt.grid(True)
plt.xlabel('x, [м]')
plt.ylabel('y, [м]')
plt.legend()
plt.show()

#plt.plot(tt_1, alff_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, deltt_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, delt_r_1, linestyle = "--", color = 'b')
plt.plot(tt_2, alff_2, linestyle = "-", color = 'r', label = 'Угол атаки')
plt.plot(tt_2, deltt_2, linestyle = "-", color = 'g', label = 'Потребный угол отклоения рулей')
plt.plot(tt_2, delt_r_2, linestyle = "-", color = 'b', label = 'Реальный угол отклонения рулей')
plt.axis([0, 16, -15, 15])
plt.xlabel('t, [c]')
plt.ylabel('alf, delt, [град]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, cyy_op_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cyy_kr_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, cyy1_alf_1, linestyle = "--", color = 'b')
#plt.plot(tt_1, cyy_korm_1, linestyle = "--", color = 'y')
#plt.plot(tt_1, cyy_nos_1, linestyle = "--", color = 'k')
plt.plot(tt_2, cyy_op_2, linestyle = "-", color = 'r', label = 'Cy_alf_op')
plt.plot(tt_2, cyy_kr_2, linestyle = "-", color = 'g', label = 'Cy_alf_kr')
plt.plot(tt_2, cyy1_alf_2, linestyle = "-", color = 'b', label = 'Cy1_alf')
plt.plot(tt_2, cyy_korm_2, linestyle = "-", color = 'y', label = 'Cy_alf_korm')
plt.plot(tt_2, cyy_nos_2, linestyle = "-", color = 'k', label = 'Cy_alf_nos')
plt.ylabel('Cy1_alf')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, cxx_nos_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cxx_0f_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, cxx_korm_1, linestyle = "--", color = 'b')
#plt.plot(tt_1, cxx_tr_1, linestyle = "--", color = 'y')

plt.plot(tt_2, cxx_nos_2, linestyle = "-", color = 'r', label = 'Сопротивление носовой части')
plt.plot(tt_2, cxx_0f_2, linestyle = "-", color = 'g', label = 'Сумарное сопротивление корпуса')
plt.plot(tt_2, cxx_korm_2, linestyle = "-", color = 'b', label = 'Сопротивление кормы')
plt.plot(tt_2, cxx_tr_2, linestyle = "-", color = 'y', label = 'Сопротивление трению')

plt.ylabel('Cx')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, cxx_kr_pr_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cxx_kr_vol_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, cxx_0_kr_1, linestyle = "--", color = 'b')

plt.plot(tt_2, cxx_kr_pr_2, linestyle = "-", color = 'r', label = 'Профильное сопротивление')
plt.plot(tt_2, cxx_kr_vol_2, linestyle = "-", color = 'g', label = 'Волновое сопротивление')
plt.plot(tt_2, cxx_0_kr_2, linestyle = "-", color = 'b', label = 'Суммарное сопротивление крыльев')
plt.ylabel('Cx')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, cxx_op_pr_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cxx_op_vol_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, cxx_0_op_1, linestyle = "--", color = 'b')

plt.plot(tt_2, cxx_op_pr_2, linestyle = "-", color = 'r', label = 'Профильное сопротивление')
plt.plot(tt_2, cxx_op_vol_2, linestyle = "-", color = 'g', label = 'Волновое сопротивление')
plt.plot(tt_2, cxx_0_op_2,  linestyle = "-", color ='b', label = 'Суммарное сопротивление оперения')
plt.ylabel('Cx_0')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, cxx_0_op_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, cxx_0_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, cxx_0_kr_1, linestyle = "--", color = 'b')
#plt.plot(tt_1, cxx_0f_1, linestyle = "--", color = 'y')

plt.plot(tt_2, cxx_0_op_2, linestyle = "-", color = 'r', label = 'Суммарное лобовое сопротивление')
plt.plot(tt_2, cxx_0_2, linestyle = "-", color = 'g', label = 'Лобовое сопротивление оперения')
plt.plot(tt_2, cxx_0_kr_2, linestyle = "-", color = 'b', label = 'Лобовое сопротивление крыльев')
plt.plot(tt_2, cxx_0f_2, linestyle = "-", color = 'y', label = 'Лобовое сопротивление корпуса')
plt.ylabel('Cx_0')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

#plt.plot(tt_1, m_zz_wz_1, linestyle = "--", color = 'r')
#plt.plot(tt_1, m_zz_wz_kr_1, linestyle = "--", color = 'g')
#plt.plot(tt_1, m_zz_wz_op_1, linestyle = "--", color = 'b')
#plt.plot(tt_1, m_zz_wz_f_1, linestyle = "--", color = 'y')

plt.plot(tt_2, m_zz_wz_2, linestyle = "-", color = 'r', label = 'Суммарный демпфирующий момент')
plt.plot(tt_2, m_zz_wz_kr_2, linestyle = "-", color = 'g', label = 'Демпфирующий момент крыльев')
plt.plot(tt_2, m_zz_wz_op_2, linestyle = "-", color = 'b', label = 'Демпфирующий момент оперения')
plt.plot(tt_2, m_zz_wz_f_2, linestyle = "-", color = 'y', label = 'Демпфирующий момент корпуса')
plt.ylabel('mz_wz')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

# _______________________________


plt.plot(tt_2, f_xx_kr_exp_2, linestyle = "-", color = 'g', label = 'Сила лобового сопротивления')
plt.plot(tt_2, f_yy_kr_exp_2, linestyle = "-", color = 'r', label = 'Подъемная сила')
plt.legend()
plt.grid(True)
plt.xlabel('t, [c]')
plt.ylabel('F, [H]')
plt.show()

plt.plot(tt_1, f_yy_exp_1, linestyle = "--", color = 'r')
plt.plot(tt_1, f_xx_exp_1, linestyle = "--", color = 'g')
plt.plot(tt_2, f_yy_exp_2, linestyle = "-", color = 'r', label = 'Подъемная сила')
plt.plot(tt_2, f_xx_exp_2, linestyle = "-", color = 'g', label = 'Сила лобового сопротивления')
plt.grid(True)
plt.legend()
plt.xlabel('t, [c]')
plt.ylabel('F, [H]')
plt.show()

plt.plot(tt_1, cyy_sum_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cxx_sum_1, linestyle = "--", color = 'g')
plt.plot(tt_2, cyy_sum_2, linestyle = "-", color = 'r', label = 'Коэффициент подъемной силы')
plt.plot(tt_2, cxx_sum_2, linestyle = "-", color = 'g', label = 'Коэффициент лобового сопротивления')
plt.grid(True)
plt.xlabel('t, [c]')
plt.ylabel('Cy, Cx')
plt.legend()
plt.show()

plt.plot(tt_1, da_ball_1, "--")
plt.plot(tt_2, da_ball_2)
plt.ylabel('delta/alf')
plt.xlabel('t, [с]')
plt.grid(True)
plt.show()

plt.plot(x_ffa_op_1, tt_1, linestyle = "--", color = 'r')
plt.plot(x_ffa_kr_1, tt_1, linestyle = "--", color = 'g')
plt.plot(x_ffa_f_1, tt_1, linestyle = "--", color = 'b')
plt.plot(x_ffa_1, tt_1,linestyle = "--", color = 'y')
plt.plot(x_ffa_op_2, tt_2, linestyle = "-", color = 'r', label = 'Координата фокуса оперения')
plt.plot(x_ffa_kr_2, tt_2,linestyle = "-", color = 'g', label = 'Координата фокуса крыльев')
plt.plot(x_ffa_f_2, tt_2, linestyle = "-", color = 'b', label = 'Координата фокуса корпуса')
plt.plot(x_ffa_2, tt_2, linestyle = "-", color ='y', label = 'Координата фокуса Ла')
plt.grid(True)
plt.axis([-0.3, 1.7, 0, 16])
plt.ylabel('t, [с]')
plt.xlabel('x_fa, [м]')
plt.legend()
plt.show()

plt.plot(tt_1, machh_1, "--", label = 'первоначального варианта')
plt.plot(tt_2, machh_2, label = 'финального варианта')
plt.grid(True)
plt.ylabel('M')
plt.xlabel('t, [с]')
plt.show()


plt.plot(xx_1, yy_1, linestyle = "--", color = 'r')
plt.plot(xt_1, yt_1, linestyle = "--", color = 'g')
plt.plot(xx_2, yy_2, linestyle = "-", color = 'r', label = 'Траектория ракеты')
plt.plot(xt_2, yt_2, linestyle = "-", color = 'g', label = 'Траектория цели')
plt.grid(True)
plt.xlabel('x, [м]')
plt.ylabel('y, [м]')
plt.legend()
plt.show()

plt.plot(tt_1, alff_1, linestyle = "--", color = 'r')
plt.plot(tt_1, deltt_1, linestyle = "--", color = 'g')
plt.plot(tt_1, delt_r_1, linestyle = "--", color = 'b')
plt.plot(tt_2, alff_2, linestyle = "-", color = 'r', label = 'Угол атаки')
plt.plot(tt_2, deltt_2, linestyle = "-", color = 'g', label = 'Потребный угол отклоения рулей')
plt.plot(tt_2, delt_r_2, linestyle = "-", color = 'b', label = 'Реальный угол отклонения рулей')
plt.axis([0, 16, -15, 15])
plt.xlabel('t, [c]')
plt.ylabel('alf, delt, [град]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, cyy_op_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cyy_kr_1, linestyle = "--", color = 'g')
plt.plot(tt_1, cyy1_alf_1, linestyle = "--", color = 'b')
plt.plot(tt_1, cyy_korm_1, linestyle = "--", color = 'y')
plt.plot(tt_1, cyy_nos_1, linestyle = "--", color = 'k')
plt.plot(tt_2, cyy_op_2, linestyle = "-", color = 'r', label = 'Cy_alf_op')
plt.plot(tt_2, cyy_kr_2, linestyle = "-", color = 'g', label = 'Cy_alf_kr')
plt.plot(tt_2, cyy1_alf_2, linestyle = "-", color = 'b', label = 'Cy1_alf')
plt.plot(tt_2, cyy_korm_2, linestyle = "-", color = 'y', label = 'Cy_alf_korm')
plt.plot(tt_2, cyy_nos_2, linestyle = "-", color = 'k', label = 'Cy_alf_nos')
plt.ylabel('Cy1_alf')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, cxx_nos_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cxx_0f_1, linestyle = "--", color = 'g')
plt.plot(tt_1, cxx_korm_1, linestyle = "--", color = 'b')
plt.plot(tt_1, cxx_tr_1, linestyle = "--", color = 'y')

plt.plot(tt_2, cxx_nos_2, linestyle = "-", color = 'r', label = 'Сопротивление носовой части')
plt.plot(tt_2, cxx_0f_2, linestyle = "-", color = 'g', label = 'Сумарное сопротивление корпуса')
plt.plot(tt_2, cxx_korm_2, linestyle = "-", color = 'b', label = 'Сопротивление кормы')
plt.plot(tt_2, cxx_tr_2, linestyle = "-", color = 'y', label = 'Сопротивление трению')

plt.ylabel('Cx')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, cxx_kr_pr_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cxx_kr_vol_1, linestyle = "--", color = 'g')
plt.plot(tt_1, cxx_0_kr_1, linestyle = "--", color = 'b')

plt.plot(tt_2, cxx_kr_pr_2, linestyle = "-", color = 'r', label = 'Профильное сопротивление')
plt.plot(tt_2, cxx_kr_vol_2, linestyle = "-", color = 'g', label = 'Волновое сопротивление')
plt.plot(tt_2, cxx_0_kr_2, linestyle = "-", color = 'b', label = 'Суммарное сопротивление крыльев')
plt.ylabel('Cx')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, cxx_op_pr_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cxx_op_vol_1, linestyle = "--", color = 'g')
plt.plot(tt_1, cxx_0_op_1, linestyle = "--", color = 'b')

plt.plot(tt_2, cxx_op_pr_2, linestyle = "-", color = 'r', label = 'Профильное сопротивление')
plt.plot(tt_2, cxx_op_vol_2, linestyle = "-", color = 'g', label = 'Волновое сопротивление')
plt.plot(tt_2, cxx_0_op_2,  linestyle = "-", color ='b', label = 'Суммарное сопротивление оперения')
plt.ylabel('Cx_0')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, cxx_0_op_1, linestyle = "--", color = 'r')
plt.plot(tt_1, cxx_0_1, linestyle = "--", color = 'g')
plt.plot(tt_1, cxx_0_kr_1, linestyle = "--", color = 'b')
plt.plot(tt_1, cxx_0f_1, linestyle = "--", color = 'y')

plt.plot(tt_2, cxx_0_op_2, linestyle = "-", color = 'r', label = 'Суммарное лобовое сопротивление')
plt.plot(tt_2, cxx_0_2, linestyle = "-", color = 'g', label = 'Лобовое сопротивление оперения')
plt.plot(tt_2, cxx_0_kr_2, linestyle = "-", color = 'b', label = 'Лобовое сопротивление крыльев')
plt.plot(tt_2, cxx_0f_2, linestyle = "-", color = 'y', label = 'Лобовое сопротивление корпуса')
plt.ylabel('Cx_0')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(tt_1, m_zz_wz_1, linestyle = "--", color = 'r')
plt.plot(tt_1, m_zz_wz_kr_1, linestyle = "--", color = 'g')
plt.plot(tt_1, m_zz_wz_op_1, linestyle = "--", color = 'b')
plt.plot(tt_1, m_zz_wz_f_1, linestyle = "--", color = 'y')

plt.plot(tt_2, m_zz_wz_2, linestyle = "-", color = 'r', label = 'Суммарный демпфирующий момент')
plt.plot(tt_2, m_zz_wz_kr_2, linestyle = "-", color = 'g', label = 'Демпфирующий момент крыльев')
plt.plot(tt_2, m_zz_wz_op_2, linestyle = "-", color = 'b', label = 'Демпфирующий момент оперения')
plt.plot(tt_2, m_zz_wz_f_2, linestyle = "-", color = 'y', label = 'Демпфирующий момент корпуса')
plt.ylabel('mz_wz')
plt.xlabel('t, [с]')
plt.legend()
plt.grid(True)
plt.show()



plt.plot(lambdd, krit_l)
plt.grid(True)
plt.show()

plt.plot(lambd_k, krit_korm)
plt.grid(True)
plt.show()

plt.plot(lambdd_kr, krit_l_kr)
plt.grid(True)
plt.show()

plt.plot(nuu_kr, krit_nu_kr)
plt.grid(True)
plt.show()

plt.plot(lambdd_op, krit_l_op)
plt.grid(True)
plt.show()

plt.plot(nuu_op, krit_nu_op)
plt.grid(True)
plt.show()



