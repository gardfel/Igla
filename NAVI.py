import math as kk
import matplotlib.pyplot as plt
import random
import Tabl
import Aero


class Target(object):

    def __init__(self):
        self.v = random.uniform(200, 200)  # сумарная скорость цели
        self.ang = 0
        self.x = random.uniform(6000, 6000)
        self.y = random.uniform(4100, 4100)
        self.alf = -0.0925  # random.uniform(0, 0.3925)

    def next_coord(self, *args):

        flag_t = args[0]
        # self.alf += random.uniform(-0.0086, +0.0086)
        if flag_t <= 2:
            self.x += self.v * kk.cos(self.alf) * dt
            self.y += self.v * kk.sin(self.alf) * dt
        else:
            self.v += 10 * dt
            self.alf = -0.0925 + 0.1 * (flag_t - 2)
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
        self.p = 2768
        self.p1 = 2768  # сила тяги на первом режиме
        self.p2 = 900  # сила тяги на втором режиме
        self.p3 = 0
        self.x_ct = 0.78
        self.q = 0

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



        self.omega = kk.atan((target.y - self.y) / (target.x - self.x))
        # ((abs(target.y - self.y) > 5) or (abs(target.x - self.x) > 5)) and
        while (t < 16):
            flag_v.append(650)
            # self.delt = kk.asin()
            mach = self.v / Tabl.tab_atm(self.y, 2)
            machh.append(mach)

            cy_nos = Aero.Cy_nos(Tabl.tab_3_2(mach, l_nos, l_cil), Tabl.tab_3_4(mach, 0, l_cil),
                                 Tabl.tab_3_2(mach, 1.32, 22.43), Tabl.tab_3_4(mach, 0, 22.43),
                                 Tabl.tab_3_4(mach, 1, 23.75))
            cyy_nos.append(cy_nos)
            cy_korm = -0.2 * 2 / 57.3 * (1 - n_korm ** 2)
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
            eps_sr_alf = 57.3 / (2 * kk.pi) * i_ / z_b * L_k_op / L_k_kr * (cy1_alf_op / l_op) * k_aa_op / K_aa_kr * fi_eps
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
            if x_tt >= (0.26 / L_f):
                x_tt = f_t / Ff

            cx_tr = Tabl.tab_4_2(re_f, x_tt) / 2 * Ff / Sf
            cxx_tr.append(cx_tr)

            cx_nos = Tabl.tab_4_11(mach, l_nos_)
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

            da_ball.append(-m_z_alf / m_z_delt_op)

            t += dt
            i += 1
            v_sr += self.v

            # функция наведения (Метод наведения и кинематические параметры перехвата)

            par_1 = self.fi_viz
            self.fi_viz = kk.atan((target.y - self.y) / (target.x - self.x))
            self.d_fi_viz = (self.fi_viz - par_1) / dt
            self.d_omega = a_m * self.d_fi_viz

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
            # print(self.omega * 180 / kk.pi,self.d_omega * 180 / kk.pi)
            # print(n_y_a)

            khi = cy1_delt_op * self.delt / cy1_alf
            # print(khi)
            q = Tabl.tab_atm(self.y, 4) * self.v ** 2 / 2
            # print(n_y_a * self.m * g)
            # print(cy1_alf * q * Sf * (1 + khi), self.p / 57.3)
            self.alf_potr = (n_y_a * self.m * g) / (cy1_alf * q * Sf * (1 + khi) + self.p / 57.3)
            alff.append(self.alf_potr)
            if self.alf < self.alf_potr:
                self.alf += (d_delta_max) / 3 * dt
            elif self.alf > self.alf_potr:
                self.alf -= (d_delta_max) / 3 * dt
            cy_sum = cy1_alf * self.alf * kk.pi / 180 + cy1_alf + cy1_delt_op * self.delt * kk.pi / 180
            cyy_sum.append(cy_sum)
            cx_sum = cx_0 + cx_ind
            cxx_sum.append(cx_sum)
            f_x_exp = cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_y_exp = cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8
            f_yy_exp.append(cy_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)
            f_xx_exp.append(cx_sum * Tabl.tab_atm(self.y, 4) * self.v ** 2 * d ** 2 * kk.pi / 8)

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
                self.v_ = 1 / self.m * (self.p - f_x_exp) - g * kk.sin(self.omega)
            elif (t > t_st) and (t <= t_st + t_m):
                self.p = self.p2
                self.x_ct += dx_cm2 * dt
                self.m -= dm2 * dt
                self.v_ = 1 / self.m * (self.p - f_x_exp) - g * kk.sin(self.omega)
            else:
                self.p = self.p3
                self.v_ = 1 / self.m * (self.p - f_x_exp) - g * kk.sin(self.omega)
            self.v1 = self.v
            self.v += self.v_ * dt
            self.x += (self.v + self.v1) / 2 * kk.cos(self.omega) * dt
            self.y += (self.v + self.v1) / 2 * kk.sin(self.omega) * dt
            print(self.v)

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

        """plt.plot(tt, da_ball)
        plt.ylabel('delta/alf')
        plt.xlabel('t, [с]')
        plt.grid(True)
        plt.show()"""

        """plt.plot(x_ffa_op, tt, 'r', label = 'Координата фокуса оперения')
        plt.plot(x_ffa_kr, tt, 'g', label = 'Координата фокуса крыльев')
        plt.plot(x_ffa_f, tt, 'b', label = 'Координата фокуса корпуса')
        plt.plot(x_ffa, tt, 'y', label = 'Координата фокуса Ла')
        plt.grid(True)
        plt.axis([0, 1.7, 0, 16])
        plt.ylabel('t, [с]')
        plt.xlabel('x_fa, [м]')
        plt.legend()
        plt.show()"""

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


dt = 10 ** -4
g = 9.80665

t_st = 2.7
t_m = 8.3

dm1 = (11.292 - 8.996) / t_st
dm2 = (8.996 - 7.3) / t_m

dx_cm1 = (0.671 - 0.78) / t_st
dx_cm2 = (0.637 - 0.671) / t_m

vertel = Target()
igla = Rocket()
d = 0.072  # диаметр миделевого сечения

# параметры для метода наведения:

a_m = 2  # коэффициент быстроты реакции ракеты на маневр цели (от 1 до бесконечности)
eps_krit = 38 * kk.pi / 180  # предельное отклонение координатора головки самонаведения
d_omega_max = 40 * kk.pi / 180
delta_max = 15 * kk.pi / 180
d_delta_max = 35 * kk.pi / 180

L_f = 1.626  # длина корпуса
L__f = L_f / L_f
L_nos = 0.155  # длина носовой части
W_nos = 0.00053  # Объем носовой части
L_korm = 0.160
S_dn = kk.pi * 0.05 ** 2 / 4
n_korm = 0.05 / 0.072
W_korm = 0.114 * S_dn + 0.000136
S_f = kk.pi * d ** 2 / 4
Ff = 0.3619  # площадь обтекаемой потоком поверхности корпуса (без донного среза)
Sf = kk.pi * (d ** 2) / 4  # площадь миделя
f_t = 0.045216

l_cil = 45.5 / 8
l_nos = 12.5 / 8

tan_05_kr = 0.307  # тангенс среднего угла крыла
l_kr = 1.46  # относительное удлинение крыла
L_k_kr = 0.146 - 0.072  # размах консолей крыла
b_kr = 0.1061  # ширина крыла у корпуса
b_a_kr = 0.101  # САХ консоли крыльев
b__a_kr = b_a_kr / L_f
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
b__a_op = b_a_op / L_f
x_c_pl_ba = 0.349  # координата цт площади передних консолей (середина САХ консолей)
x_b_a_op = 0.334  # коордигата начала САХ консоли опрения
x_b_op = 0.334  # координата начала бортовой хорды оперения
a_op = 0.015
c_op = 0.03  # относительная толщина профиля оперения
L_hv_op = 0.933  # Расстояние от конца бортовой хорды оперения до кормового среза корпуса
L1_op = 0.35  # Расстояние от носика корпуса до серидины бортовой хорды оперения
koef_op = 1 / (1 - a_op / b_op)
print(koef_kr, koef_op, "koef")

b_ak_kr = 0.094
x_otn_op_kr = 0.846 / b_ak_kr  # относительное расстояние между оперением и средней хордой крыльев
S__f = kk.pi * (d ** 2) / 4 / (kk.pi * (d ** 2) / 4)  # относительная площадь корпуса
S__op = 0.0078575 / (kk.pi * (d ** 2) / 4)  # относительная площадь передних несущих поверхностей
S__kr = 0.015 / (kk.pi * (d ** 2) / 4)  # относительная площадь задних несущих поверхностей

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

nu_kor = 0.047 / d  # относительное сужение кормовой части
l_korm = 0.0463 / d  # относительное удлинение кормовой части

igla.navigation(vertel)
