import matplotlib.pyplot as plt
import math as kk


def procent(*args):

    var1 = args[0]
    tab1 = args[1]
    tab2 = args[2]

    return (var1 - tab1) / (tab2 - tab1)


def interpol(*args):

    tab1 = args[0]
    tab2 = args[1]
    proc = args[2]

    return tab2 + (tab1 - tab2) * proc


def tab_3_2(*args):
    """
    Функция для вывода Су для конической ГЧ
    :param args: число Маха, относительное удлинение носка и цилиндрической части
    :return: Значение Су ГЧ
    """

    mah = args[0]
    lambd_nos = args[1]
    lambd_cil = args[2]

    cy1iz_alf_0 = [0.0350, 0.0350, 0.0350, 0.0348, 0.0338, 0.0335, 0.0333, 0.0327, 0.0323, 0.0321, 0.0321, 0.0321,
                   0.0320, 0.0320, 0.0320, 0.0320, 0.0320, 0.0320]
    cy1iz_alf_05 = [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0426, 0.0455, 0.0465, 0.0460, 0.0448, 0.0435, 0.0425,
                    0.0417, 0.0410, 0.0402, 0.0400, 0.0395, 0.0390]
    cy1iz_alf_1 = [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0505, 0.0512, 0.0507, 0.0500, 0.0490,
                   0.0480, 0.0473, 0.0465, 0.0460, 0.0450, 0.0448]
    cy1iz_alf_2 = [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0540, 0.0560, 0.0565, 0.0560,
                   0.0555, 0.0550, 0.0545, 0.0537, 0.0530, 0.0525]
    cy1iz_alf_3 = [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0548, 0.0570, 0.0580, 0.0585,
                   0.0587, 0.0587, 0.0585, 0.0583, 0.0580, 0.0575]
    cy1iz_alf_4 = [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0548, 0.0570, 0.0584, 0.0594,
                   0.0599, 0.0600, 0.0600, 0.0600, 0.0600, 0.0600]
    razm = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6]

    if (mah ** 2 - 1) >= 0:
        razmm = kk.sqrt(mah ** 2 - 1) / lambd_nos
    else:
        razmm = -kk.sqrt(1 - mah ** 2) / lambd_nos

    otnos = lambd_cil / lambd_nos

    k = int(razmm // 0.2 + 5)

    if k >= 17:
        k = 17

    if otnos == 0:
        cy1 = interpol(cy1iz_alf_0[k], cy1iz_alf_0[k - 1], procent(razmm, razm[k - 1], razm[k]))
    elif (otnos <= 0.5) and (otnos >= 0):
        cy1 = interpol(interpol(cy1iz_alf_05[k], cy1iz_alf_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       interpol(cy1iz_alf_0[k], cy1iz_alf_0[k - 1], procent(razmm, razm[k - 1], razm[k])),(otnos / 0.5))
    elif (otnos <= 1) and (otnos >= 0.5):
        cy1 = interpol(interpol(cy1iz_alf_1[k], cy1iz_alf_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       interpol(cy1iz_alf_05[k], cy1iz_alf_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       (otnos - 0.5) / 0.5)
    elif (otnos <= 2) and (otnos >= 1):
        cy1 = interpol(interpol(cy1iz_alf_2[k], cy1iz_alf_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       interpol(cy1iz_alf_1[k], cy1iz_alf_1[k - 1], procent(razmm, razm[k - 1], razm[k])), otnos - 1)
    elif (otnos <= 3) and (otnos >= 2):
        cy1 = interpol(interpol(cy1iz_alf_3[k], cy1iz_alf_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       interpol(cy1iz_alf_2[k], cy1iz_alf_2[k - 1], procent(razmm, razm[k - 1], razm[k])), otnos - 2)
    elif (otnos <= 4) and (otnos >= 3):
        cy1 = interpol(interpol(cy1iz_alf_4[k], cy1iz_alf_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                       interpol(cy1iz_alf_3[k], cy1iz_alf_3[k - 1], procent(razmm, razm[k - 1], razm[k])), otnos - 3)
    else:
        cy1 = interpol(cy1iz_alf_4[k], cy1iz_alf_4[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return cy1


def tab_3_4(*args):
    """
    Функция для вывода Су плоской и сферической ГЧ
    :param args: число маха, 0 - для сферической, 1 - для плоской, относительное удлинение цилиндрической части
    :return: Значение Су ГЧ
    """

    mah = args[0]
    krit = args[1]
    lambd_cil = args[2]

    cy1iz_sph = [0.0345, 0.0345, 0.0350, 0.0355, 0.0383, 0.0433, 0.0440, 0.0427, 0.0408, 0.0390, 0.0370,
                 0.0355, 0.0335, 0.0318, 0.0300]
    cy1iz_cil = [0.0345, 0.0345, 0.0350, 0.0350, 0.0365, 0.0375, 0.0370, 0.0355, 0.0333, 0.0315, 0.0300,
                 0.0283, 0.0270, 0.0255, 0.0240]
    razm = [-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

    if (mah ** 2 - 1) >= 0:
        razmm = kk.sqrt(mah ** 2 - 1) / lambd_cil
    else:
        razmm = -kk.sqrt(1 - mah ** 2) / lambd_cil

    k = int(razmm // 0.1 + 5)

    if krit == 0:
        cy1 = interpol(cy1iz_sph[k], cy1iz_sph[k - 1], procent(razmm, razm[k - 1], razm[k]))
    elif krit == 1:
        cy1 = interpol(cy1iz_cil[k], cy1iz_cil[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return cy1


def tab_3_5(*args):

    """
    Функиця для вывода значения Су крыльев
    :param args: Число Маха, удлинение крыльев, относительная толщина профиля крыла, тангенс угла средней стреловидности
    :return: Значение Су для изолированных крыльев
    """

    mah = args[0]
    lambd_k = args[1]
    c_ = args[2]
    tg_khi_05 = args[3]

    # lambd*tan_0.5 = 0
    cy1_iz_kr_a_000 = [0.0176, 0.0187, 0.0200, 0.0216, 0.0235, 0.0258, 0.0288, 0.0325, 0.0365, 0.0350, 0.0306, 0.0265,
                       0.0228, 0.0200, 0.0178, 0.0158, 0.0145]
    cy1_iz_kr_a_025 = [0.0176, 0.0187, 0.0200, 0.0216, 0.0235, 0.0258, 0.0288, 0.0321, 0.0344, 0.0330, 0.0298, 0.0260,
                       0.0227, 0.0200, 0.0178, 0.0158, 0.0145]
    cy1_iz_kr_a_050 = [0.0176, 0.0187, 0.0200, 0.0216, 0.0235, 0.0258, 0.0288, 0.0321, 0.0320, 0.0304, 0.0278, 0.0249,
                       0.0218, 0.0196, 0.0175, 0.0155, 0.0140]
    cy1_iz_kr_a_100 = [0.0176, 0.0187, 0.0200, 0.0216, 0.0238, 0.0265, 0.0280, 0.0285, 0.0277, 0.0264, 0.0248, 0.0229,
                       0.0208, 0.0188, 0.0168, 0.0152, 0.0138]
    cy1_iz_kr_a_150 = [0.0170, 0.0183, 0.0198, 0.0222, 0.0243, 0.0253, 0.0250, 0.0245, 0.0236, 0.0225, 0.0214, 0.0200,
                       0.0188, 0.0173, 0.0160, 0.0146, 0.0134]

    # lambd*tan_0.5 = 1
    cy1_iz_kr_b_000 = [0.0165, 0.0175, 0.0190, 0.0205, 0.0220, 0.0240, 0.0265, 0.0295, 0.0330, 0.0330, 0.0300, 0.0260,
                       0.0225, 0.0200, 0.0175, 0.0158, 0.0145]
    cy1_iz_kr_b_025 = [0.0165, 0.0175, 0.0190, 0.0205, 0.0222, 0.0240, 0.0264, 0.0290, 0.0315, 0.0310, 0.0282, 0.0248,
                       0.0220, 0.0194, 0.0172, 0.0155, 0.0141]
    cy1_iz_kr_b_050 = [0.0168, 0.0179, 0.0194, 0.0207, 0.0224, 0.0244, 0.0270, 0.0296, 0.0302, 0.0288, 0.0263, 0.0235,
                       0.0211, 0.0190, 0.0170, 0.0155, 0.0141]
    cy1_iz_kr_b_100 = [0.0165, 0.0175, 0.0190, 0.0205, 0.0225, 0.0251, 0.0266, 0.0270, 0.0262, 0.0248, 0.0233, 0.0215,
                       0.0200, 0.0184, 0.0170, 0.0156, 0.0142]
    cy1_iz_kr_b_150 = [0.0168, 0.0179, 0.0193, 0.0211, 0.0232, 0.0244, 0.0240, 0.0235, 0.0225, 0.0216, 0.0205, 0.0194,
                       0.0180, 0.0170, 0.0158, 0.0146, 0.0134]

    # lambd*tan_0.5 = 2
    cy1_iz_kr_v_000 = [0.0155, 0.0165, 0.0177, 0.0190, 0.0202, 0.0220, 0.0240, 0.0276, 0.0296, 0.0298, 0.0281, 0.0254,
                       0.0220, 0.0194, 0.0175, 0.0157, 0.0143]
    cy1_iz_kr_v_025 = [0.0155, 0.0165, 0.0177, 0.0190, 0.0202, 0.0220, 0.0244, 0.0274, 0.0285, 0.0279, 0.0256, 0.0233,
                       0.0208, 0.0190, 0.0172, 0.0155, 0.0142]
    cy1_iz_kr_v_050 = [0.0156, 0.0167, 0.0178, 0.0190, 0.0207, 0.0225, 0.0248, 0.0271, 0.0272, 0.0259, 0.0242, 0.0221,
                       0.0201, 0.0184, 0.0169, 0.0153, 0.0141]
    cy1_iz_kr_v_100 = [0.0155, 0.0164, 0.0176, 0.0190, 0.0205, 0.0227, 0.0243, 0.0247, 0.0242, 0.0232, 0.0218, 0.0205,
                       0.0189, 0.0175, 0.0162, 0.0150, 0.0137]
    cy1_iz_kr_v_150 = [0.0152, 0.0162, 0.0173, 0.0190, 0.0210, 0.0220, 0.0222, 0.0219, 0.0212, 0.0202, 0.0192, 0.0183,
                       0.0172, 0.0163, 0.0153, 0.0142, 0.0133]

    # lambd*tan_0.5 = 3
    cy1_iz_kr_g_000 = [0.0135, 0.0144, 0.0158, 0.0170, 0.0184, 0.0200, 0.0220, 0.0245, 0.0265, 0.0274, 0.0265, 0.0245,
                       0.0220, 0.0197, 0.0178, 0.0160, 0.0146]
    cy1_iz_kr_g_025 = [0.0135, 0.0144, 0.0158, 0.0170, 0.0184, 0.0200, 0.0215, 0.0236, 0.0250, 0.0245, 0.0226, 0.0208,
                       0.0190, 0.0175, 0.0162, 0.0150, 0.0139]
    cy1_iz_kr_g_050 = [0.0135, 0.0144, 0.0158, 0.0170, 0.0184, 0.0200, 0.0216, 0.0233, 0.0239, 0.0231, 0.0218, 0.0202,
                       0.0185, 0.0172, 0.0158, 0.0147, 0.0138]
    cy1_iz_kr_g_100 = [0.0135, 0.0144, 0.0158, 0.0170, 0.0184, 0.0204, 0.0215, 0.0221, 0.0217, 0.0209, 0.0199, 0.0187,
                       0.0175, 0.0163, 0.0153, 0.0143, 0.0134]
    cy1_iz_kr_g_150 = [0.0136, 0.0145, 0.0158, 0.0174, 0.0187, 0.0197, 0.0200, 0.0198, 0.0194, 0.0188, 0.0181, 0.0172,
                       0.0162, 0.0153, 0.0144, 0.0136, 0.0127]

    razm = [-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]

    if (mah ** 2 - 1) <= 0:
        razmm = -lambd_k * (1 - mah ** 2) ** 0.5
    elif (lambd_k * kk.sqrt(mah ** 2 - 1)) >= 5.1:
        return 4 / (57.3 * (kk.sqrt(mah ** 2 - 1)))
    else:
        razmm = lambd_k * (mah ** 2 - 1) ** 0.5

    if razmm >= 5:
        razmm = 5
    if razmm <= -3.5:
        razmm = -3.5

    k = int(razmm // 0.5 + 8)
    if k >= 16:
        k = 16
    if k <= 0:
        k = 0

    otnos = lambd_k * tg_khi_05

    krit = lambd_k * c_ ** (1 / 3)

    if otnos <= 0:
        if krit == 0:
            cy1 = interpol(cy1_iz_kr_a_000[k], cy1_iz_kr_a_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
        elif (krit > 0) and (krit <= 0.25):
            cy1 = interpol(interpol(cy1_iz_kr_a_025[k], cy1_iz_kr_a_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_a_000[k], cy1_iz_kr_a_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0, 0.25))
        elif (krit > 0.25) and (krit <= 0.5):
            cy1 = interpol(interpol(cy1_iz_kr_a_050[k], cy1_iz_kr_a_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_a_025[k], cy1_iz_kr_a_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0.25, 0.5))
        elif (krit > 0.5) and (krit <= 1):
            cy1 = interpol(interpol(cy1_iz_kr_a_100[k], cy1_iz_kr_a_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_a_050[k], cy1_iz_kr_a_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0.5, 1))
        elif (krit > 1) and (krit <= 1.5):
            cy1 = interpol(interpol(cy1_iz_kr_a_150[k], cy1_iz_kr_a_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_a_100[k], cy1_iz_kr_a_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 1, 1.5))
        else:
            cy1 = interpol(cy1_iz_kr_a_150[k], cy1_iz_kr_a_150[k - 1], procent(razmm, razm[k - 1], razm[k]))

    elif (otnos > 0) and (otnos <= 1):
        param = procent(otnos, 0, 1)
        if krit == 0:
            cy1a = interpol(cy1_iz_kr_a_000[k], cy1_iz_kr_a_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_b_000[k], cy1_iz_kr_b_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0) and (krit <= 0.25):
            cy1a = interpol(interpol(cy1_iz_kr_a_025[k], cy1_iz_kr_a_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_a_000[k], cy1_iz_kr_a_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1b = interpol(interpol(cy1_iz_kr_b_025[k], cy1_iz_kr_b_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_000[k], cy1_iz_kr_b_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.25) and (krit <= 0.5):
            cy1a = interpol(interpol(cy1_iz_kr_a_050[k], cy1_iz_kr_a_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_a_025[k], cy1_iz_kr_a_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1b = interpol(interpol(cy1_iz_kr_b_050[k], cy1_iz_kr_b_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_025[k], cy1_iz_kr_b_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.5) and (krit <= 1):
            cy1a = interpol(interpol(cy1_iz_kr_a_100[k], cy1_iz_kr_a_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_a_050[k], cy1_iz_kr_a_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1b = interpol(interpol(cy1_iz_kr_b_100[k], cy1_iz_kr_b_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_050[k], cy1_iz_kr_b_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 1) and (krit <= 1.5):
            cy1a = interpol(interpol(cy1_iz_kr_a_150[k], cy1_iz_kr_a_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_a_100[k], cy1_iz_kr_a_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1b = interpol(interpol(cy1_iz_kr_b_150[k], cy1_iz_kr_b_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_100[k], cy1_iz_kr_b_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1 = interpol(cy1b, cy1a, param)
        else:
            cy1a = interpol(cy1_iz_kr_a_150[k], cy1_iz_kr_a_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_b_150[k], cy1_iz_kr_b_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)

    elif (otnos > 1) and (otnos <= 2):
        param = procent(otnos, 1, 2)
        if krit == 0:
            cy1a = interpol(cy1_iz_kr_b_000[k], cy1_iz_kr_b_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_v_000[k], cy1_iz_kr_v_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0) and (krit <= 0.25):
            cy1a = interpol(interpol(cy1_iz_kr_b_025[k], cy1_iz_kr_b_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_000[k], cy1_iz_kr_b_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1b = interpol(interpol(cy1_iz_kr_v_025[k], cy1_iz_kr_v_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_000[k], cy1_iz_kr_v_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.25) and (krit <= 0.5):
            cy1a = interpol(interpol(cy1_iz_kr_b_050[k], cy1_iz_kr_b_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_025[k], cy1_iz_kr_b_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1b = interpol(interpol(cy1_iz_kr_v_050[k], cy1_iz_kr_v_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_025[k], cy1_iz_kr_v_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.5) and (krit <= 1):
            cy1a = interpol(interpol(cy1_iz_kr_b_100[k], cy1_iz_kr_b_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_050[k], cy1_iz_kr_b_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1b = interpol(interpol(cy1_iz_kr_v_100[k], cy1_iz_kr_v_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_050[k], cy1_iz_kr_v_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 1) and (krit <= 1.5):
            cy1a = interpol(interpol(cy1_iz_kr_b_150[k], cy1_iz_kr_b_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_b_100[k], cy1_iz_kr_b_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1b = interpol(interpol(cy1_iz_kr_v_150[k], cy1_iz_kr_v_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_100[k], cy1_iz_kr_v_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1 = interpol(cy1b, cy1a, param)
        else:
            cy1a = interpol(cy1_iz_kr_b_150[k], cy1_iz_kr_b_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_v_150[k], cy1_iz_kr_v_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)

    elif (otnos > 2) and (otnos <= 3):
        param = procent(otnos, 2, 3)
        if krit == 0:
            cy1a = interpol(cy1_iz_kr_v_000[k], cy1_iz_kr_v_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_g_000[k], cy1_iz_kr_g_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0) and (krit <= 0.25):
            cy1a = interpol(interpol(cy1_iz_kr_v_025[k], cy1_iz_kr_v_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_000[k], cy1_iz_kr_v_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1b = interpol(interpol(cy1_iz_kr_g_025[k], cy1_iz_kr_g_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_g_000[k], cy1_iz_kr_g_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0, 0.25))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.25) and (krit <= 0.5):
            cy1a = interpol(interpol(cy1_iz_kr_v_050[k], cy1_iz_kr_v_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_025[k], cy1_iz_kr_v_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1b = interpol(interpol(cy1_iz_kr_g_050[k], cy1_iz_kr_g_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_g_025[k], cy1_iz_kr_g_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.25, 0.5))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 0.5) and (krit <= 1):
            cy1a = interpol(interpol(cy1_iz_kr_v_100[k], cy1_iz_kr_v_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_050[k], cy1_iz_kr_v_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1b = interpol(interpol(cy1_iz_kr_g_100[k], cy1_iz_kr_g_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_g_050[k], cy1_iz_kr_g_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 0.5, 1))
            cy1 = interpol(cy1b, cy1a, param)
        elif (krit > 1) and (krit <= 1.5):
            cy1a = interpol(interpol(cy1_iz_kr_v_150[k], cy1_iz_kr_v_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_v_100[k], cy1_iz_kr_v_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1b = interpol(interpol(cy1_iz_kr_g_150[k], cy1_iz_kr_g_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cy1_iz_kr_g_100[k], cy1_iz_kr_g_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(krit, 1, 1.5))
            cy1 = interpol(cy1b, cy1a, param)
        else:
            cy1a = interpol(cy1_iz_kr_v_150[k], cy1_iz_kr_v_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1b = interpol(cy1_iz_kr_g_150[k], cy1_iz_kr_g_150[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cy1 = interpol(cy1b, cy1a, param)

    else:
        if krit == 0:
            cy1 = interpol(cy1_iz_kr_g_000[k], cy1_iz_kr_g_000[k - 1], procent(razmm, razm[k - 1], razm[k]))
        elif (krit > 0) and (krit <= 0.25):
            cy1 = interpol(interpol(cy1_iz_kr_g_025[k], cy1_iz_kr_g_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_g_000[k], cy1_iz_kr_g_000[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0, 0.25))
        elif (krit > 0.25) and (krit <= 0.5):
            cy1 = interpol(interpol(cy1_iz_kr_g_050[k], cy1_iz_kr_g_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_g_025[k], cy1_iz_kr_g_025[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0.25, 0.5))
        elif (krit > 0.5) and (krit <= 1):
            cy1 = interpol(interpol(cy1_iz_kr_g_100[k], cy1_iz_kr_g_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_g_050[k], cy1_iz_kr_g_050[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 0.5, 1))
        elif (krit > 1) and (krit <= 1.5):
            cy1 = interpol(interpol(cy1_iz_kr_g_150[k], cy1_iz_kr_g_150[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(cy1_iz_kr_g_100[k], cy1_iz_kr_g_100[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(krit, 1, 1.5))
        elif krit > 1.5:
            cy1 = interpol(cy1_iz_kr_g_150[k], cy1_iz_kr_g_150[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return cy1 * lambd_k


def tab_3_16(*args):
    """
    Определение относительной координаты вихря
    :param args: число Маха, относительное удлинение, сужение и тангенс угла средней стреловидности передних консолей
    :return:
    """

    z_v_1_0 = []
    z_v_1_2 = []
    z_v_1_4 = []

    z_v_2_0 = []
    z_v_2_066 = []
    z_v_2_2 = []

    z_v_inf_0 = []
    z_v_inf_2 = []

    razm = []


def tab_3_17(*args):
    """
    Функция вывода коэффициента интерференции передних и задних несущих поверхностей
    :param args: относительное сужение задних поверхностей,
    :return:
    """


def tab_3_21(*args):

    """
    Функция вывода значения коэффициента торможения потока при носовой части
    :param args: Число Маха, относительное удлинение носовой части
    :return:
    """

    mah = args[0]
    lambd = args[1]

    x_t_2 = [0.920, 0.950, 0.974, 0.983, 0.988, 0.992, 0.994, 0.996, 1.000, 1.000, 1.000, 1.000, 1.000]
    x_t_3 = [0.780, 0.865, 0.925, 0.955, 0.970, 0.980, 0.985, 0.990, 0.992, 0.996, 1.000, 1.000, 1.000]
    x_t_4 = [0.640, 0.768, 0.845, 0.900, 0.933, 0.950, 0.970, 0.980, 0.986, 0.990, 0.995, 1.000, 1.000]
    x_t_5 = [0.540, 0.650, 0.750, 0.835, 0.885, 0.920, 0.940, 0.960, 0.970, 0.980, 0.987, 0.992, 1.000]
    lambd_tab = [0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]

    k = int(lambd // 0.2 - 2)

    if mah <= 2:
        k_t = interpol(interpol(x_t_2[k], x_t_2[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])), 1,
                       procent(mah, 0, 2))
    elif (mah <= 3) and (mah >= 2):
        k_t = interpol(interpol(x_t_3[k], x_t_3[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       interpol(x_t_2[k], x_t_2[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       procent(mah, 2, 3))
    elif (mah <= 4) and (mah >= 3):
        k_t = interpol(interpol(x_t_4[k], x_t_4[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       interpol(x_t_3[k], x_t_3[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       procent(mah, 3, 4))
    elif (mah <= 5) and (mah >= 4):
        k_t = interpol(interpol(x_t_5[k], x_t_5[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       interpol(x_t_4[k], x_t_4[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k])),
                       procent(mah, 4, 5))
    else:
        k_t = interpol(x_t_5[k], x_t_5[k - 1], procent(lambd, lambd_tab[k - 1], lambd_tab[k]))

    return k_t


def tab_3_22(*args):

    """
    Функция вывода коэффициента торможения потока, вызванного обтеканием передних несущих поверхностей
    :param args: Число Маха, относительное расстояние между оперением и крыльями
    :return:
    """

    mah = args[0]
    x_ = args[1]

    k_t_00 = [0.967, 0.968, 0.967, 0.960, 0.890, 0.775, 0.730, 0.705, 0.690, 0.680, 0.672, 0.670, 0.670]
    k_t_02 = [0.972, 0.974, 0.973, 0.965, 0.940, 0.855, 0.823, 0.810, 0.801, 0.798, 0.792, 0.789, 0.788]
    k_t_04 = [0.977, 0.979, 0.978, 0.975, 0.960, 0.905, 0.880, 0.866, 0.860, 0.858, 0.855, 0.850, 0.850]
    k_t_06 = [0.980, 0.982, 0.981, 0.978, 0.965, 0.925, 0.900, 0.890, 0.885, 0.882, 0.880, 0.880, 0.880]
    k_t_08 = [0.984, 0.986, 0.985, 0.981, 0.970, 0.933, 0.910, 0.905, 0.900, 0.897, 0.895, 0.895, 0.895]
    k_t_10 = [0.987, 0.989, 0.988, 0.985, 0.975, 0.940, 0.920, 0.910, 0.906, 0.905, 0.900, 0.900, 0.900]
    mah_tab = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

    if mah <= 1:
        k = int(mah // 0.25 + 1)
    else:
        k = int(mah // 0.5 + 3)

    if x_ == 0:
        k_t = interpol(k_t_00[k], k_t_00[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))
    elif (x_ <= 0.2) and (x_ >= 0):
        k_t = interpol(interpol(k_t_02[k], k_t_02[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       interpol(k_t_00[k], k_t_00[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       procent(x_, 0, 0.2))
    elif (x_ <= 0.4) and (x_ >= 0.2):
        k_t = interpol(interpol(k_t_04[k], k_t_04[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       interpol(k_t_02[k], k_t_02[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       procent(x_, 0.2, 0.4))
    elif (x_ <= 0.6) and (x_ >= 0.4):
        k_t = interpol(interpol(k_t_06[k], k_t_06[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       interpol(k_t_04[k], k_t_04[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       procent(x_, 0.4, 0.6))
    elif (x_ <= 0.8) and (x_ >= 0.6):
        k_t = interpol(interpol(k_t_08[k], k_t_08[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       interpol(k_t_06[k], k_t_06[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       procent(x_, 0.6, 0.8))
    elif (x_ <= 1) and (x_ >= 0.8):
        k_t = interpol(interpol(k_t_10[k], k_t_10[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       interpol(k_t_08[k], k_t_08[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                       procent(x_, 0.8, 1))
    else:
        k_t = interpol(k_t_10[k], k_t_10[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))

    return k_t


def tab_4_2(*args):
    """
    Вывод коэффициента трения плоской пластины в зависимости от числа Re и относительной координаты перехода ламинарного
    пограничного слоя в турбулентный
    :param args: число Re, относительная координата перехода ламинарного пограничного слоя в турбулентный
    :return:
    """

    re_ = args[0]
    x_t = args[1]

    d_cf_0 = [0.00900, 0.00788, 0.00730, 0.00695, 0.00670, 0.00650, 0.00635, 0.00620, 0.00610, 0.00600, 0.00535,
              0.00505, 0.00485, 0.00470, 0.00457, 0.00445, 0.00435, 0.00430, 0.00422, 0.00380, 0.00365, 0.00350,
              0.00340]
    d_cf_01 = [0.00855, 0.00750, 0.00690, 0.00660, 0.00635, 0.00615, 0.00595, 0.00580, 0.00570, 0.00560, 0.00500,
               0.00465, 0.00448, 0.00433, 0.00420, 0.00410, 0.00402, 0.00398, 0.00390, 0.00350, 0.00332, 0.00320,
               0.00310]
    d_cf_02 = [0.00800, 0.00695, 0.00642, 0.00606, 0.00585, 0.00565, 0.00550, 0.00535, 0.00525, 0.00515, 0.00455,
               0.00430, 0.00410, 0.00395, 0.00385, 0.00375, 0.00368, 0.00360, 0.00355, 0.00325, 0.00305, 0.00295,
               0.00285]
    d_cf_03 = [0.00750, 0.00645, 0.00590, 0.00555, 0.00530, 0.00510, 0.00495, 0.00485, 0.00475, 0.00465, 0.00410,
               0.00383, 0.00365, 0.00355, 0.00345, 0.00335, 0.00330, 0.00325, 0.00320, 0.00290, 0.00274, 0.00265,
               0.00255]
    d_cf_04 = [0.00690, 0.00585, 0.00535, 0.00500, 0.00478, 0.00460, 0.00445, 0.00435, 0.00425, 0.00415, 0.00365,
               0.00340, 0.00325, 0.00310, 0.00300, 0.00295, 0.00290, 0.00285, 0.00280, 0.00250, 0.00240, 0.00230,
               0.00225]
    d_cf_05 = [0.00635, 0.00530, 0.00480, 0.00450, 0.00425, 0.00410, 0.00395, 0.00380, 0.00375, 0.00365, 0.00315,
               0.00295, 0.00280, 0.00265, 0.00260, 0.00250, 0.00248, 0.00245, 0.00240, 0.00220, 0.00210, 0.00200,
               0.00195]
    d_cf_06 = [0.00570, 0.00470, 0.00425, 0.00390, 0.00370, 0.00350, 0.00340, 0.00330, 0.00320, 0.00310, 0.00270,
               0.00245, 0.00235, 0.00225, 0.00220, 0.00215, 0.00210, 0.00205, 0.00200, 0.00185, 0.00174, 0.00168,
               0.00162]
    d_cf_1 = [0.00275, 0.00200, 0.00160, 0.00140, 0.00125, 0.00115, 0.00105, 0.00100, 0.00095, 0.00090, 0.00065,
              0.00052, 0.00049, 0.00040, 0.00038, 0.00035, 0.00032, 0.00030, 0.00028, 0.00020, 0.00015, 0.00010,
              0.00008]
    re_mas = [1.0E+06, 2.0E+06, 3.0E+06, 4.0E+06, 5.0E+06, 6.0E+06, 7.0E+06, 8.0E+06, 9.0E+06, 1.0E+07, 2.0E+07,
              3.0E+07, 4.0E+07, 5.0E+07, 6.0E+07, 7.0E+07, 8.0E+07, 9.0E+07, 1.0E+08, 2.0E+08, 3.0E+08, 4.0E+08,
              5.0E+08]

    k = 0

    for i in range(len(re_mas)):
        if re_ <= re_mas[i]:
            if re_ >= re_mas[i - 1]:
                k = i

    if x_t == 0:
        d_cf = interpol(d_cf_0[k], d_cf_0[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                          kk.log(re_mas[k], 10)))
    elif (x_t <= 0.1) and (x_t >= 0):
        d_cf = interpol(interpol(d_cf_01[k], d_cf_01[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_0[k], d_cf_0[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                   kk.log(re_mas[k], 10))), procent(x_t, 0, 0.1))
    elif (x_t <= 0.2) and (x_t >= 0.1):
        d_cf = interpol(interpol(d_cf_02[k], d_cf_02[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_01[k], d_cf_01[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.1, 0.2))
    elif (x_t <= 0.3) and (x_t >= 0.2):
        d_cf = interpol(interpol(d_cf_03[k], d_cf_03[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_02[k], d_cf_02[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.2, 0.3))
    elif (x_t <= 0.4) and (x_t >= 0.3):
        d_cf = interpol(interpol(d_cf_04[k], d_cf_04[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_03[k], d_cf_03[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.3, 0.4))
    elif (x_t <= 0.5) and (x_t >= 0.4):
        d_cf = interpol(interpol(d_cf_05[k], d_cf_05[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_04[k], d_cf_04[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.4, 0.5))
    elif (x_t <= 0.6) and (x_t >= 0.5):
        d_cf = interpol(interpol(d_cf_06[k], d_cf_06[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))),
                        interpol(d_cf_05[k], d_cf_05[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.5, 0.6))
    elif (x_t <= 1) and (x_t >= 0.6):
        d_cf = interpol(interpol(d_cf_1[k], d_cf_1[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                   kk.log(re_mas[k], 10))),
                        interpol(d_cf_06[k], d_cf_06[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                                     kk.log(re_mas[k], 10))), procent(x_t, 0.6, 1))
    else:
        d_cf = interpol(d_cf_1[k], d_cf_1[k - 1], procent(kk.log(re_, 10), kk.log(re_mas[k - 1], 10),
                                                          kk.log(re_mas[k], 10)))

    return d_cf


def tab_4_3(*args):

    """
    Функция для вывода коэффициента трения плоской пластины в зависимости от числа Маха
    :param args: Число Маха, координата перехода ламинарного в турбулентный слой
    :return: значение коэффициента трения плоской пластины
    """

    mah = args[0]
    x_t = args[1]

    nu_m_0 = [1.00, 0.93, 0.77, 0.63, 0.51, 0.42, 0.35]
    nu_m_02 = [1.00, 0.94, 0.78, 0.64, 0.52, 0.43, 0.36]
    nu_m_05 = [1.00, 0.95, 0.80, 0.66, 0.55, 0.47, 0.41]
    nu_m_06 = [1.00, 0.955, 0.82, 0.68, 0.57, 0.50, 0.45]
    nu_m_08 = [1.00, 0.96, 0.85, 0.74, 0.64, 0.58, 0.52]
    nu_m_1 = [1.00, 0.99, 0.97, 0.93, 0.89, 0.84, 0.80]
    mah_tab = [0, 1, 2, 3, 4, 5, 6]

    k = int(mah // 1 + 1)

    if x_t == 0:
        nu_m = interpol(nu_m_0[k], nu_m_0[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))
    elif (x_t <= 0.2) and (x_t >= 0):
        nu_m = interpol(interpol(nu_m_02[k], nu_m_02[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        interpol(nu_m_0[k], nu_m_0[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        procent(x_t, 0, 0.2))
    elif (x_t <= 0.5) and (x_t >= 0.2):
        nu_m = interpol(interpol(nu_m_05[k], nu_m_05[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        interpol(nu_m_02[k], nu_m_02[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        procent(x_t, 0.2, 0.5))
    elif (x_t <= 0.6) and (x_t >= 0.5):
        nu_m = interpol(interpol(nu_m_06[k], nu_m_06[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        interpol(nu_m_05[k], nu_m_05[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        procent(x_t, 0.5, 0.6))
    elif (x_t <= 0.8) and (x_t >= 0.6):
        nu_m = interpol(interpol(nu_m_08[k], nu_m_08[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        interpol(nu_m_06[k], nu_m_06[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        procent(x_t, 0.6, 0.8))
    elif (x_t <= 1) and (x_t >= 0.8):
        nu_m = interpol(interpol(nu_m_1[k], nu_m_1[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        interpol(nu_m_08[k], nu_m_08[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                        procent(x_t, 0.8, 1))
    else:
        nu_m = interpol(nu_m_1[k], nu_m_1[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))

    return nu_m



def tab_4_5(*args):
    """
    Вывод критического числа Рейнольдса при шероховатой поверхности тела
    :param args: 1 - число Маха, 2 - число Рейнольдса, 3 - класс чистоты (от 4 до 9), 4 - длина тела
    :return: критическое число Рейнольдса
    """

    mah = args[0]
    re_ = args[1]
    cla = args[2]
    l_t = args[3]

    h = [40e-06, 20e-06, 10e-06, 6.3e-06, 3.2e-06, 1.6e-06]

    re_0 = [3e+06, 1.85e+06, 1e+06, 7e+05, 5.25e+05, 4e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05]
    re_1 = [3.8e+06, 3e+06, 2e+06, 1.6e+06, 1.3e+06, 9.4e+06, 2e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05]
    re_2 = [4e+06, 3.8e+06, 3.2e+06, 3e+06, 2.7e+06, 2.5e+06, 1.35e+06, 1.5e+05, 1e+05, 1e+05, 1e+05, 1e+05, 1e+05]
    re_3 = [3.3e+06, 3.2e+06, 3.1e+06, 3e+06, 2.8e+06, 2.7e+06, 2.1e+06, 1.5e+06, 9e+05, 4.8e+05, 1e+05, 1e+05, 1e+05]
    re_4 = [3.2e+06, 3.2e+06, 3.2e+06, 3.2e+06, 3.2e+06, 3.15e+06, 2.9e+06, 2.6e+06, 2e+06, 1.7e+06, 1.4e+06, 1e+05, 1e+05]
    re_5 = [3.8e+06, 3.8e+06, 3.8e+06, 3.8e+06, 3.8e+06, 3.8e+06, 3.7e+06, 3.4e+06, 3e+06, 2.8e+06, 2.6e+06, 1.3e+06, 1e+05]
    razm = [10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000, 2000, 4000]

    razmm = re_ * h[cla - 4] / l_t
    k = 0

    for i in range(len(razm)):
        if (razmm >= razm[i - 1]) and (razmm <= razm[i]):
            k = i
        elif (razmm >= razm[i]):
            k = i

    if mah == 0:
        re_t = interpol(re_0[k], re_0[k - 1], procent(razmm, razm[k - 1], razm[k]))
    elif (mah > 0) and (mah <= 1):
        re_t = interpol(interpol(re_1[k], re_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                        interpol(re_0[k], re_0[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(mah, 0, 1))
    elif (mah > 1) and (mah <= 2):
        re_t = interpol(interpol(re_2[k], re_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                        interpol(re_1[k], re_1[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(mah, 1, 2))
    elif (mah > 2) and (mah <= 3):
        re_t = interpol(interpol(re_3[k], re_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                        interpol(re_2[k], re_2[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(mah, 2, 3))
    elif (mah > 3) and (mah <= 4):
        re_t = interpol(interpol(re_4[k], re_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                        interpol(re_3[k], re_3[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(mah, 3, 4))
    elif (mah > 4) and (mah <= 5):
        re_t = interpol(interpol(re_5[k], re_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                        interpol(re_4[k], re_4[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(mah, 4, 5))
    else:
        re_t = interpol(re_5[k], re_5[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return re_t


def tab_4_11(*args):

    """
    Вывод Сх для конической ГЧ
    :param args: на вход число Маха и относительное удлинение носовой части
    :return: Значение Сх
    """

    mah = args[0]
    lambd = args[1]

    cx_nos_15 = [0.044, 0.061, 0.110, 0.193, 0.320, 0.373, 0.400, 0.408, 0.412, 0.408, 0.400, 0.390, 0.378, 0.355,
                 0.335, 0.319, 0.306, 0.296, 0.290, 0.262, 0.250, 0.242, 0.237, 0.235, 0.233, 0.232]
    cx_nos_20 = [0.017, 0.032, 0.055, 0.109, 0.230, 0.260, 0.273, 0.276, 0.265, 0.252, 0.238, 0.228, 0.220, 0.206,
                 0.197, 0.189, 0.184, 0.180, 0.179, 0.164, 0.153, 0.146, 0.141, 0.138, 0.136, 0.134]
    cx_nos_25 = [0.014, 0.018, 0.032, 0.064, 0.162, 0.210, 0.216, 0.211, 0.197, 0.186, 0.175, 0.167, 0.160, 0.150,
                 0.142, 0.137, 0.133, 0.130, 0.129, 0.119, 0.111, 0.105, 0.101, 0.098, 0.096, 0.095]
    cx_nos_30 = [0.004, 0.008, 0.017, 0.035, 0.121, 0.167, 0.168, 0.158, 0.147, 0.139, 0.133, 0.127, 0.122, 0.115,
                 0.110, 0.106, 0.103, 0.100, 0.099, 0.091, 0.084, 0.079, 0.075, 0.073, 0.072, 0.071]
    cx_nos_40 = [0.000, 0.004, 0.009, 0.021, 0.075, 0.091, 0.091, 0.088, 0.085, 0.082, 0.080, 0.078, 0.076, 0.073,
                 0.070, 0.068, 0.066, 0.065, 0.064, 0.059, 0.055, 0.053, 0.051, 0.050, 0.049, 0.047]
    cx_nos_50 = [0.000, 0.002, 0.006, 0.010, 0.035, 0.057, 0.059, 0.058, 0.057, 0.056, 0.056, 0.055, 0.054, 0.053,
                 0.051, 0.050, 0.049, 0.049, 0.048, 0.045, 0.042, 0.040, 0.038, 0.037, 0.036, 0.035]
    mah_tab = [0.6, 0.7, 0.8, 0.9, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5,
               3, 3.5, 4, 4.5, 5, 5.5]
    k = 0

    if mah <= 0.6:
        k = 1
    elif mah <= 1:
        k = int(mah // 0.1 - 5)
    elif mah <= 1.4:
        k = int(mah // 0.05 - 15)
    elif mah <= 2:
        k = int(mah // 0.1 - 1)
    elif mah <= 5.5:
        k = int(mah // 0.5 + 15)
    elif mah > 5.5:
        k = len(mah_tab) - 1

    if lambd <= 1.5:
        cx_nos = interpol(cx_nos_15[k], cx_nos_15[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))
    elif (lambd <= 2) and (lambd > 1.5):
        cx_nos = interpol(interpol(cx_nos_20[k], cx_nos_20[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_15[k], cx_nos_15[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 1.5, 2))
    elif (lambd <= 2.5) and (lambd > 2):
        cx_nos = interpol(interpol(cx_nos_25[k], cx_nos_25[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_20[k], cx_nos_20[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 2, 2.5))
    elif (lambd <= 3) and (lambd > 2.5):
        cx_nos = interpol(interpol(cx_nos_30[k], cx_nos_30[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_25[k], cx_nos_25[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 2.5, 3))
    elif (lambd <= 4) and (lambd > 3):
        cx_nos = interpol(interpol(cx_nos_40[k], cx_nos_40[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_30[k], cx_nos_30[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 3, 4))
    elif (lambd <= 5) and (lambd > 4):
        cx_nos = interpol(interpol(cx_nos_50[k], cx_nos_50[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_40[k], cx_nos_40[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 4, 5))
    else:
        cx_nos = interpol(cx_nos_50[k], cx_nos_50[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))

    return cx_nos


def tab_4_13(*args):
    """
    Вывод значения Сх носовой части с эллиптической образующей
    :param args: Число Маха, относительное удлинение носовой части
    :return:
    """

    mah = args[0]
    lambd = args[1]

    cx_nos_000 = [0.79, 0.84, 0.92, 1.05, 1.19, 1.29, 1.37, 1.43, 1.47, 1.51, 1.54, 1.56, 1.59, 1.60, 1.61, 1.62, 1.62]
    cx_nos_025 = [0.110, 0.210, 0.400, 0.600, 0.740, 0.850, 0.930, 0.980, 1.010, 1.040, 1.060, 1.080, 1.095, 1.105,
                  1.110, 1.115, 1.115]
    cx_nos_050 = [0.06, 0.04, 0.12, 0.34, 0.51, 0.62, 0.69, 0.74, 0.77, 0.80, 0.81, 0.82, 0.83, 0.84, 0.84, 0.84, 0.84]
    cx_nos_100 = [0.000, 0.001, 0.050, 0.180, 0.300, 0.390, 0.450, 0.480, 0.510, 0.520, 0.525, 0.525, 0.530, 0.530,
                  0.530, 0.525, 0.520]
    cx_nos_200 = [-0.02, -0.02, -0.02, 0.04, 0.16, 0.21, 0.23, 0.24, 0.24, 0.24, 0.24, 0.24, 0.23, 0.23, 0.23, 0.23,
                  0.22]
    mah_tab = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6]

    k = 0

    if mah <= 0.6:
        k = 0
    else:
        for i in range(len(mah_tab)):
            if (mah <= mah_tab[i]) and (mah >= mah_tab[i - 1]):
                k = i

    if lambd == 0:
        cx_nos = interpol(cx_nos_000[k], cx_nos_000[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))
    elif (lambd <= 0.25) and (lambd >= 0):
        cx_nos = interpol(interpol(cx_nos_025[k], cx_nos_025[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_000[k], cx_nos_000[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 0, 0.25))
    elif (lambd <= 0.5) and (lambd >= 0.25):
        cx_nos = interpol(interpol(cx_nos_050[k], cx_nos_050[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_025[k], cx_nos_025[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 0.25, 0.5))
    elif (lambd <= 1) and (lambd >= 0.5):
        cx_nos = interpol(interpol(cx_nos_100[k], cx_nos_100[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_050[k], cx_nos_050[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 0.5, 1))
    elif (lambd <= 2) and (lambd >= 1):
        cx_nos = interpol(interpol(cx_nos_200[k], cx_nos_200[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          interpol(cx_nos_100[k], cx_nos_100[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k])),
                          procent(lambd, 1, 2))
    else:
        cx_nos = interpol(cx_nos_200[k], cx_nos_200[k - 1], procent(mah, mah_tab[k - 1], mah_tab[k]))

    return cx_nos


def tab_4_24(*args):
    """
    Вывод значения Сх для кормовой части
    :param args: число Маха, относительное сужение кормовой части, относительное удлинение кормовой части
    :return: значение коэффициента Сх кормовой части
    """

    mah = args[0]
    nu_korm = args[1]
    lambd = args[2]

    cx_korm_nu_000_20 = [0.068, 0.075, 0.085, 0.100, 0.123, 0.140, 0.143, 0.143, 0.137, 0.130, 0.122, 0.116, 0.111,
                         0.107, 0.103, 0.100, 0.085, 0.074, 0.066, 0.060, 0.055, 0.051, 0.049]
    cx_korm_nu_000_30 = [0.061, 0.067, 0.075, 0.086, 0.103, 0.107, 0.108, 0.106, 0.102, 0.097, 0.092, 0.089, 0.086,
                         0.083, 0.080, 0.078, 0.069, 0.063, 0.058, 0.054, 0.051, 0.049, 0.048]

    cx_korm_nu_050_15 = [0.089, 0.089, 0.091, 0.093, 0.097, 0.099, 0.102, 0.105, 0.104, 0.102, 0.099, 0.095, 0.092,
                         0.088, 0.083, 0.080, 0.064, 0.052, 0.043, 0.036, 0.030, 0.026, 0.022]
    cx_korm_nu_050_20 = [0.060, 0.060, 0.062, 0.064, 0.068, 0.070, 0.073, 0.076, 0.075, 0.073, 0.070, 0.066, 0.063,
                         0.060, 0.056, 0.052, 0.044, 0.038, 0.033, 0.028, 0.024, 0.020, 0.017]
    cx_korm_nu_050_25 = [0.053, 0.053, 0.053, 0.053, 0.053, 0.050, 0.054, 0.054, 0.054, 0.052, 0.050, 0.048, 0.045,
                         0.043, 0.041, 0.040, 0.033, 0.028, 0.024, 0.021, 0.018, 0.016, 0.014]

    cx_korm_nu_075_10 = [0.050, 0.051, 0.052, 0.053, 0.055, 0.056, 0.057, 0.058, 0.059, 0.059, 0.057, 0.055, 0.051,
                         0.048, 0.044, 0.040, 0.030, 0.024, 0.020, 0.017, 0.014, 0.012, 0.011]
    cx_korm_nu_075_15 = [0.043, 0.042, 0.039, 0.036, 0.034, 0.030, 0.034, 0.034, 0.034, 0.033, 0.032, 0.031, 0.030,
                         0.029, 0.027, 0.025, 0.020, 0.016, 0.013, 0.011, 0.009, 0.008, 0.007]
    cx_korm_nu_075_20 = [0.038, 0.036, 0.032, 0.029, 0.026, 0.025, 0.024, 0.023, 0.022, 0.022, 0.021, 0.021, 0.020,
                         0.019, 0.019, 0.018, 0.013, 0.010, 0.008, 0.007, 0.006, 0.006, 0.006]

    mah_tab = [0.6, 0.7, 0.8, 0.9, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]

    k = 0
    cx_korm = 0
    if mah <= 0.6:
        k = 0
    else:
        for i in range(len(mah_tab)):
            if (mah <= mah_tab[i]) and (mah >= mah_tab[i - 1]):
                k = i

    if nu_korm == 0:
        cx_korm = interpol(interpol(cx_korm_nu_000_30[k], cx_korm_nu_000_30[k - 1], procent(mah, mah_tab[k - 1],
                                                                                            mah_tab[k])),
                           interpol(cx_korm_nu_000_20[k], cx_korm_nu_000_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                            mah_tab[k])),
                           procent(lambd, 2, 3))

    elif (nu_korm <= 0.5) and (nu_korm >= 0):
        proc = procent(nu_korm, 0, 0.5)
        if lambd <= 1.5:
            cx_korm = interpol(interpol(cx_korm_nu_050_15[k], cx_korm_nu_050_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                mah_tab[k])),
                               interpol(cx_korm_nu_000_20[k], cx_korm_nu_000_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                mah_tab[k])), proc)
        elif (lambd <= 2) and (lambd >= 1.5):
            proc1 = procent(lambd, 1.5, 2)
            cx_nu_05 = interpol(interpol(cx_korm_nu_050_20[k], cx_korm_nu_050_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_050_15[k], cx_korm_nu_050_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc1)
            cx_nu_00 = interpol(cx_korm_nu_000_20[k], cx_korm_nu_000_20[k - 1],
                                procent(mah, mah_tab[k - 1], mah_tab[k]))
            cx_korm = interpol(cx_nu_05, cx_nu_00, proc)
        elif (lambd <= 2.5) and (lambd >= 2):
            proc1 = procent(lambd, 2, 2.5)
            proc2 = procent(lambd, 2, 3)
            cx_nu_05 = interpol(interpol(cx_korm_nu_050_25[k], cx_korm_nu_050_25[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_050_20[k], cx_korm_nu_050_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc1)
            cx_nu_00 = interpol(interpol(cx_korm_nu_000_30[k], cx_korm_nu_000_30[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_000_20[k], cx_korm_nu_000_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc2)
            cx_korm = interpol(cx_nu_05, cx_nu_00, proc)
        elif (lambd >= 2.5) and (lambd <= 3):
            proc1 = procent(lambd, 2, 3)
            cx_nu_05 = interpol(cx_korm_nu_050_25[k], cx_korm_nu_050_25[k - 1],
                                procent(mah, mah_tab[k - 1], mah_tab[k]))
            cx_nu_00 = interpol(interpol(cx_korm_nu_000_30[k], cx_korm_nu_000_30[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_000_20[k], cx_korm_nu_000_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc1)
            cx_korm = interpol(cx_nu_05, cx_nu_00, proc)
        else:
            cx_nu_05 = interpol(cx_korm_nu_050_25[k], cx_korm_nu_050_25[k - 1],
                                procent(mah, mah_tab[k - 1], mah_tab[k]))
            cx_nu_00 = interpol(cx_korm_nu_000_30[k], cx_korm_nu_000_30[k - 1],
                                procent(mah, mah_tab[k - 1], mah_tab[k]))
            cx_korm = interpol(cx_nu_05, cx_nu_00, proc)
    elif (nu_korm <= 0.75) and (nu_korm >= 0.5):
        proc = procent(nu_korm, 0.5, 0.75)
        if lambd <= 1:
            cx_korm = interpol(interpol(cx_korm_nu_075_10[k], cx_korm_nu_075_10[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                mah_tab[k])),
                               interpol(cx_korm_nu_050_15[k], cx_korm_nu_050_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                mah_tab[k])), proc)
        elif (lambd >= 1) and (lambd <= 1.5):
            proc1 = procent(lambd, 1, 1.5)
            cx_nu_075 = interpol(interpol(cx_korm_nu_075_15[k], cx_korm_nu_075_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                  mah_tab[k])),
                                 interpol(cx_korm_nu_075_10[k], cx_korm_nu_075_10[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                  mah_tab[k])), proc1)
            cx_nu_05 = interpol(cx_korm_nu_050_15[k], cx_korm_nu_050_15[k - 1], mah, mah_tab[k - 1], mah_tab[k])
            cx_korm = interpol(cx_nu_075, cx_nu_05, proc)
        elif (lambd >= 2) and (lambd <= 1.5):
            proc1 = procent(lambd, 1.5, 2)
            cx_nu_075 = interpol(interpol(cx_korm_nu_075_20[k], cx_korm_nu_075_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                  mah_tab[k])),
                                 interpol(cx_korm_nu_075_15[k], cx_korm_nu_075_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                  mah_tab[k])), proc1)
            cx_nu_05 = interpol(interpol(cx_korm_nu_050_20[k], cx_korm_nu_050_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_050_15[k], cx_korm_nu_050_15[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc1)
            cx_korm = interpol(cx_nu_075, cx_nu_05, proc)
        elif (lambd >= 2.5) and (lambd <= 2):
            proc1 = procent(lambd, 2, 2.5)
            cx_nu_075 = interpol(cx_korm_nu_075_20[k], cx_korm_nu_075_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                         mah_tab[k]))
            cx_nu_05 = interpol(interpol(cx_korm_nu_050_25[k], cx_korm_nu_050_25[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])),
                                interpol(cx_korm_nu_050_20[k], cx_korm_nu_050_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                                 mah_tab[k])), proc1)
            cx_korm = interpol(cx_nu_075, cx_nu_05, proc)
        else:
            cx_nu_075 = interpol(cx_korm_nu_075_20[k], cx_korm_nu_075_20[k - 1], procent(mah, mah_tab[k - 1],
                                                                                         mah_tab[k]))
            cx_nu_05 = interpol(cx_korm_nu_050_25[k], cx_korm_nu_050_25[k - 1], procent(mah, mah_tab[k - 1],
                                                                                        mah_tab[k]))
            cx_korm = interpol(cx_nu_075, cx_nu_05, proc)

    return cx_korm


def tab_4_28(*args):
    """
    Вывод поправочного коэффициента, учитывающего влияние толщины профиля
    :param args: относительное положение точки перехода ламинарного пограничного слоя в турбулентный (Х_т_),
    относительная толщина профиля
    :return: Значение поправочного коэффициента"""

    x_t = args[0]
    c_ = args[1]

    nu_t_00 = [1.00, 1.03, 1.05, 1.08, 1.11, 1.13, 1.16, 1.19, 1.22, 1.25, 1.29, 1.33, 1.37]
    nu_t_02 = [1.000, 1.020, 1.040, 1.060, 1.080, 1.104, 1.127, 1.155, 1.180, 1.205, 1.235, 1.260, 1.295]
    nu_t_04 = [1.00, 1.01, 1.03, 1.04, 1.05, 1.07, 1.09, 1.10, 1.12, 1.14, 1.16, 1.17, 1.20]

    c_mas = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]

    k = int(c_ // 0.01 + 1)

    if x_t == 0:
        nu_t = interpol(nu_t_00[k], nu_t_00[k - 1], procent(c_, c_mas[k - 1], c_mas[k]))
    elif (x_t >= 0) and (x_t <= 0.2):
        nu_t = interpol(interpol(nu_t_02[k], nu_t_02[k - 1], procent(c_, c_mas[k - 1], c_mas[k])),
                        interpol(nu_t_00[k], nu_t_00[k - 1], procent(c_, c_mas[k - 1], c_mas[k])),
                        procent(x_t, 0, 0.2))
    elif (x_t >= 0.2) and (x_t <= 0.4):
        nu_t = interpol(interpol(nu_t_04[k], nu_t_04[k - 1], procent(c_, c_mas[k - 1], c_mas[k])),
                        interpol(nu_t_02[k], nu_t_02[k - 1], procent(c_, c_mas[k - 1], c_mas[k])),
                        procent(x_t, 0.2, 0.4))
    else:
        nu_t = interpol(nu_t_04[k], nu_t_04[k - 1], procent(c_, c_mas[k - 1], c_mas[k]))

    return nu_t


def tab_4_30(*args):
    """
    Функция вывода волнового сопротивления крыльев с ромбовидным профилем
    :param args: число Маха, относительное сужение крыла, относительное удлинение крыла, тангенс угла передней
     стреловидности, относительная толщина профиля крыла,
    :return: Знаение волнового сопротивления для профиля крыла
    """

    mah = args[0]
    nu_k = args[1]
    lambd_k = args[2]
    tan_05x = args[3]
    c_ = args[4]

# а - для относительного сужения крыла = 1
    cx_v_a_0_05 = [2.81, 3.20, 3.10, 2.50, 1.90, 1.60, 1.33, 1.13, 0.97, 0.83, 0.73, 0.67, 0.62, 0.58, 0.55, 0.52]
    cx_v_a_0_10 = [2.73, 2.71, 2.65, 2.37, 1.90, 1.60, 1.33, 1.13, 0.97, 0.83, 0.73, 0.67, 0.62, 0.58, 0.55, 0.52]
    cx_v_a_0_15 = [2.10, 2.09, 2.07, 1.95, 1.75, 1.58, 1.33, 1.13, 0.97, 0.83, 0.73, 0.67, 0.62, 0.58, 0.55, 0.52]
    cx_v_a_1 = [1.520, 1.680, 1.880, 2.070, 2.000, 1.700, 1.430, 1.220, 1.050, 0.910, 0.802, 0.750, 0.680, 0.650,
                0.600, 0.570]
    cx_v_a_2 = [0.85, 0.90, 0.98, 1.13, 1.37, 1.50, 1.50, 1.40, 1.19, 1.02, 0.90, 0.80, 0.72, 0.66, 0.60, 0.57]
    cx_v_a_3 = [0.440, 0.470, 0.510, 0.600, 0.720, 0.850, 1.040, 1.180, 1.190, 1.130, 1.000, 0.880, 0.790, 0.700,
                0.650, 0.600]
    cx_v_a_4 = [0.22, 0.24, 0.26, 0.30, 0.38, 0.45, 0.57, 0.70, 0.82, 0.92, 0.99, 0.97, 0.89, 0.79, 0.70, 0.64]
# б - для относительного сужения крыла = 5
    cx_v_b_0_05 = [2.80, 3.12, 2.95, 2.35, 1.90, 1.58, 1.34, 1.15, 1.00, 0.88, 0.80, 0.70, 0.64, 0.60, 0.55, 0.52]
    cx_v_b_0_10 = [2.68, 2.63, 2.56, 2.25, 1.90, 1.58, 1.34, 1.15, 1.00, 0.88, 0.80, 0.70, 0.64, 0.60, 0.55, 0.52]
    cx_v_b_0_15 = [2.12, 2.09, 2.05, 1.95, 1.76, 1.52, 1.34, 1.15, 1.00, 0.88, 0.80, 0.70, 0.64, 0.60, 0.55, 0.52]
    cx_v_b_1 = [1.830, 1.930, 1.990, 1.990, 1.850, 1.600, 1.400, 1.220, 1.070, 0.930, 0.830, 0.740, 0.680, 0.640,
                0.580, 0.550]
    cx_v_b_2 = [1.08, 1.15, 1.24, 1.37, 1.46, 1.46, 1.40, 1.28, 1.17, 1.02, 0.91, 0.81, 0.74, 0.66, 0.60, 0.56]
    cx_v_b_3 = [0.600, 0.650, 0.700, 0.800, 0.880, 1.000, 1.080, 1.100, 1.090, 1.040, 0.970, 0.880, 0.800, 0.720,
                0.650, 0.600]
    cx_v_b_4 = [0.330, 0.360, 0.385, 0.450, 0.490, 0.560, 0.650, 0.720, 0.800, 0.850, 0.890, 0.890, 0.860, 0.820,
                0.750, 0.650]
# в - для относительного сужения крыла = бесконечности
    cx_v_v_0_05 = [2.85, 2.98, 2.70, 2.25, 1.85, 1.59, 1.39, 1.19, 1.04, 0.91, 0.81, 0.74, 0.67, 0.61, 0.57, 0.52]
    cx_v_v_0_10 = [2.60, 2.54, 2.44, 2.13, 1.84, 1.59, 1.39, 1.19, 1.04, 0.91, 0.81, 0.74, 0.67, 0.61, 0.57, 0.52]
    cx_v_v_0_15 = [2.10, 2.05, 2.00, 1.85, 1.70, 1.51, 1.32, 1.15, 1.00, 0.90, 0.81, 0.74, 0.67, 0.61, 0.57, 0.52]
    cx_v_v_1 = [1.930, 1.925, 1.910, 1.850, 1.730, 1.570, 1.375, 1.215, 1.090, 0.970, 0.870, 0.780, 0.710, 0.640,
                0.580, 0.550]
    cx_v_v_2 = [1.31, 1.35, 1.37, 1.40, 1.38, 1.33, 1.26, 1.18, 1.09, 1.00, 0.90, 0.80, 0.73, 0.66, 0.60, 0.57]
    cx_v_v_3 = [0.800, 0.815, 0.840, 0.890, 0.930, 0.970, 0.990, 0.995, 0.980, 0.940, 0.890, 0.840, 0.780, 0.710,
                0.650, 0.600]
    cx_v_v_4 = [0.46, 0.47, 0.48, 0.50, 0.53, 0.57, 0.61, 0.65, 0.70, 0.73, 0.76, 0.77, 0.76, 0.75, 0.71, 0.64]

    razm = [0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]

    if (mah ** 2 - 1) >= 0:
        razmm = lambd_k * (mah ** 2 - 1) ** 0.5
    else:
        razmm = -lambd_k * (1 - mah ** 2) ** 0.5
        return 0

    if razmm >= 7:
        razmm = 7

    otnos = lambd_k * tan_05x
    param = lambd_k * c_ ** (1 / 3)

    if razmm <= 0.5:
        k = int(razmm // 0.25 + 1)
    else:
        k = int(razmm // 0.5 + 2)

    if k >= 15:
        k = 15

    if nu_k == 0:
        if otnos == 0:
            if param <= 0.5:
                cx_v = interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
            elif (param <= 1) and (param >= 0.5):
                cx_v = interpol(interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                procent(param, 0.5, 1))
            elif (param <= 1.5) and (param >= 1):
                cx_v = interpol(interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                procent(param, 1, 1.5))
            else:
                cx_v = interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))

        elif (otnos <= 1) and (otnos >= 0):
            proc2 = procent(otnos, 0, 1)
            if param <= 0.5:
                cx_v = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                proc2)
            elif (param <= 1) and (param >= 0.5):
                cx_v1 = interpol(interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_v = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1,
                                proc2)
            elif (param <= 1.5) and (param >= 1):
                cx_v1 = interpol(interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_v = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1,
                                proc2)
            else:
                cx_v = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                proc2)

        elif (otnos <= 2) and (otnos >= 1):
            cx_v = interpol(interpol(cx_v_a_2[k], cx_v_a_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 1, 2))
        elif (otnos <= 3) and (otnos >= 2):
            cx_v = interpol(interpol(cx_v_a_3[k], cx_v_a_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_a_2[k], cx_v_a_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 2, 3))
        elif (otnos <= 4) and (otnos >= 3):
            cx_v = interpol(interpol(cx_v_a_4[k], cx_v_a_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_a_3[k], cx_v_a_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 3, 4))
        else:
            cx_v = interpol(cx_v_a_4[k], cx_v_a_4[k - 1], procent(razmm, razm[k - 1], razm[k]))

    elif (nu_k <= 5) and (nu_k >= 0):
        proc1 = procent(nu_k, 0, 5)
        if otnos == 0:
            if param <= 0.5:
                cx_va = interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_vb = interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_v = interpol(cx_vb, cx_va, proc1)
            elif (param <= 1) and (param >= 0.5):
                cx_va = interpol(interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_vb = interpol(interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_v = interpol(cx_vb, cx_va, proc1)
            elif (param <= 1.5) and (param >= 1):
                cx_va = interpol(interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_vb = interpol(interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_v = interpol(cx_vb, cx_va, proc1)
            else:
                cx_va = interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_vb = interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_v = interpol(cx_vb, cx_va, proc1)

        elif (otnos <= 1) and (otnos >= 0):
            proc2 = procent(otnos, 0, 1)
            if param <= 0.5:
                cx_va = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_v = interpol(cx_vb, cx_va, proc1)
            elif (param <= 1) and (param >= 0.5):
                cx_v1a = interpol(interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_a_0_05[k], cx_v_a_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 0.5, 1))
                cx_va = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1a,
                                 proc2)
                cx_v1b = interpol(interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 0.5, 1))
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1b,
                                 proc2)
                cx_v = interpol(cx_vb, cx_va, proc1)
            elif (param <= 1.5) and (param >= 1):
                cx_v1a = interpol(interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_a_0_10[k], cx_v_a_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 1, 1.5))
                cx_va = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1a,
                                 proc2)
                cx_v1b = interpol(interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 1, 1.5))
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1b,
                                 proc2)
                cx_v = interpol(cx_vb, cx_va, proc1)
            else:
                cx_va = interpol(interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_a_0_15[k], cx_v_a_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_v = interpol(cx_vb, cx_va, proc1)

        elif (otnos <= 2) and (otnos >= 1):
            cx_va = interpol(interpol(cx_v_a_2[k], cx_v_a_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_a_1[k], cx_v_a_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 1, 2))
            cx_vb = interpol(interpol(cx_v_b_2[k], cx_v_b_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 1, 2))
            cx_v = interpol(cx_vb, cx_va, proc1)
        elif (otnos <= 3) and (otnos >= 2):
            cx_va = interpol(interpol(cx_v_a_3[k], cx_v_a_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_a_2[k], cx_v_a_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_vb = interpol(interpol(cx_v_b_3[k], cx_v_b_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_2[k], cx_v_b_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_v = interpol(cx_vb, cx_va, proc1)
        elif (otnos <= 4) and (otnos >= 3):
            cx_va = interpol(interpol(cx_v_a_4[k], cx_v_a_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_a_3[k], cx_v_a_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 3, 4))
            cx_vb = interpol(interpol(cx_v_b_4[k], cx_v_b_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_3[k], cx_v_b_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_v = interpol(cx_vb, cx_va, proc1)
        else:
            cx_va = interpol(cx_v_a_4[k], cx_v_a_4[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cx_vb = interpol(cx_v_b_4[k], cx_v_b_4[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cx_v = interpol(cx_vb, cx_va, proc1)

    elif (nu_k <= 1000) and (nu_k >= 5):
        proc1 = procent(nu_k, 5, 1000)
        if otnos == 0:
            if param <= 0.5:
                cx_vb = interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_vv = interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_v = interpol(cx_vv, cx_vb, proc1)
            elif (param <= 1) and (param >= 0.5):
                cx_vb = interpol(interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_vv = interpol(interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_v = interpol(cx_vv, cx_vb, proc1)
            elif (param <= 1.5) and (param >= 1):
                cx_vb = interpol(interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_vv = interpol(interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_v = interpol(cx_vv, cx_vb, proc1)
            else:
                cx_vb = interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_vv = interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))
                cx_v = interpol(cx_vv, cx_vb, proc1)

        elif (otnos <= 1) and (otnos >= 0):
            proc2 = procent(otnos, 0, 1)
            if param <= 0.5:
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_vv = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_v = interpol(cx_vv, cx_vb, proc1)
            elif (param <= 1) and (param >= 0.5):
                cx_v1b = interpol(interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_b_0_05[k], cx_v_b_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 0.5, 1))
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1b,
                                 proc2)
                cx_v1v = interpol(interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 0.5, 1))
                cx_vv = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1v,
                                 proc2)
                cx_v = interpol(cx_vv, cx_vb, proc1)
            elif (param <= 1.5) and (param >= 1):
                cx_v1b = interpol(interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_b_0_10[k], cx_v_b_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 1, 1.5))
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1b,
                                 proc2)
                cx_v1v = interpol(interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                  procent(param, 1, 1.5))
                cx_vv = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1v,
                                 proc2)
                cx_v = interpol(cx_vv, cx_vb, proc1)
            else:
                cx_vb = interpol(interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_b_0_15[k], cx_v_b_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_vv = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 proc2)
                cx_v = interpol(cx_vv, cx_vb, proc1)

        elif (otnos <= 2) and (otnos >= 1):
            cx_vb = interpol(interpol(cx_v_b_2[k], cx_v_b_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_1[k], cx_v_b_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 1, 2))
            cx_vv = interpol(interpol(cx_v_v_2[k], cx_v_v_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 1, 2))
            cx_v = interpol(cx_vv, cx_vb, proc1)
        elif (otnos <= 3) and (otnos >= 2):
            cx_vb = interpol(interpol(cx_v_b_3[k], cx_v_b_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_2[k], cx_v_b_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_vv = interpol(interpol(cx_v_v_3[k], cx_v_v_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_v_2[k], cx_v_v_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_v = interpol(cx_vv, cx_vb, proc1)
        elif (otnos <= 4) and (otnos >= 3):
            cx_vb = interpol(interpol(cx_v_b_4[k], cx_v_b_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_b_3[k], cx_v_b_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 3, 4))
            cx_vv = interpol(interpol(cx_v_v_4[k], cx_v_v_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             interpol(cx_v_v_3[k], cx_v_v_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                             procent(otnos, 2, 3))
            cx_v = interpol(cx_vv, cx_vb, proc1)
        else:
            cx_vb = interpol(cx_v_b_4[k], cx_v_b_4[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cx_vv = interpol(cx_v_v_4[k], cx_v_v_4[k - 1], procent(razmm, razm[k - 1], razm[k]))
            cx_v = interpol(cx_vv, cx_vb, proc1)

    else:
        if otnos == 0:
            if param <= 0.5:
                cx_v = interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
            elif (param <= 1) and (param >= 0.5):
                cx_v = interpol(interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                procent(param, 0.5, 1))
            elif (param <= 1.5) and (param >= 1):
                cx_v = interpol(interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                procent(param, 1, 1.5))
            else:
                cx_v = interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k]))

        elif (otnos <= 1) and (otnos >= 0):
            proc2 = procent(otnos, 0, 1)
            if param <= 0.5:
                cx_v = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                proc2)
            elif (param <= 1) and (param >= 0.5):
                cx_v1 = interpol(interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_05[k], cx_v_v_0_05[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 0.5, 1))
                cx_v = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1,
                                proc2)
            elif (param <= 1.5) and (param >= 1):
                cx_v1 = interpol(interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 interpol(cx_v_v_0_10[k], cx_v_v_0_10[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                 procent(param, 1, 1.5))
                cx_v = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])), cx_v1,
                                proc2)
            else:
                cx_v = interpol(interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                interpol(cx_v_v_0_15[k], cx_v_v_0_15[k - 1], procent(razmm, razm[k - 1], razm[k])),
                                proc2)

        elif (otnos <= 2) and (otnos >= 1):
            cx_v = interpol(interpol(cx_v_v_2[k], cx_v_v_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_v_1[k], cx_v_v_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 1, 2))
        elif (otnos <= 3) and (otnos >= 2):
            cx_v = interpol(interpol(cx_v_v_3[k], cx_v_v_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_v_2[k], cx_v_v_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 2, 3))
        elif (otnos <= 4) and (otnos >= 3):
            cx_v = interpol(interpol(cx_v_v_4[k], cx_v_v_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(cx_v_v_3[k], cx_v_v_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(otnos, 3, 4))
        else:
            cx_v = interpol(cx_v_v_4[k], cx_v_v_4[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return cx_v * lambd_k * c_ ** 2


def tab_4_32(*args):

    """
    Функция вывода коэффициента фи
    :param args: число маха,
    :return:
    """

    mah = args[0]
    tan_x = args[1]

    fi_tab = [0.00, 0.07, 0.20, 0.38, 0.54, 0.67, 0.77, 0.85, 0.91, 0.96, 0.99]
    razm = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]

    if (mah ** 2 - 1) <= 0:
        return 0
    else:
        razmm = (mah ** 2 - 1) ** 0.5 - tan_x
        k = int(razmm // 0.25 + 1)
        if k >= 10:
            k = 10
        fi = interpol(fi_tab[k], fi_tab[k - 1], procent(razmm, razm[k - 1], razm[k]))
        return fi


def tab_4_40(*args):

    """
    Функция вывода коэффициента кси (для вычисления индуктивного сопротивления корпуса)
    :param args: число Маха, относительное удлинение носовой части, 0 - оживальная носовая часть / 1 - коническая
    :return: Значение коэффициента для данных условий
    """

    mah = args[0]
    lambd_nos = args[1]
    param = args[2]

    ksi_ozhiv = [-0.32, -0.28, 0.00, 0.32, 0.51, 0.62, 0.69, 0.72, 0.75, 0.76]
    ksi_konic = [-0.200, -0.200, -0.150, -0.080, 0.000, 0.110, 0.210, 0.300, 0.370, 0.410]
    razm = [-0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8]

    if (mah ** 2 - 1) <= 0:
        razmm = -kk.sqrt(1 - mah ** 2) / lambd_nos
    else:
        razmm = kk.sqrt(mah ** 2 - 1) / lambd_nos

    k = int(razmm // 0.4 + 3)
    ksi = 0

    if param == 0:
        ksi = interpol(ksi_ozhiv[k], ksi_ozhiv[k - 1], procent(razmm, razm[k - 1], razm[k]))
    if param == 1:
        ksi = interpol(ksi_konic[k], ksi_konic[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return ksi


def tab_4_43(*args):
    """

    :param args:
    :return:
    """


def tab_5_7(*args):

    """
    Определение фокуса корпуса
    :param args: число Маха, относительное удлинение носовой части и корпуса, длина носовой части
    :return:
    """

    mah = args[0]
    lambd_nos = args[1]
    lambd_korp = args[2]
    l_nos = args[3]

    dx_05 = [0.040, 0.045, 0.052, 0.070, 0.090, 0.110, 0.125, 0.134, 0.138, 0.139, 0.140, 0.140, 0.138, 0.135, 0.130,
             0.120, 0.110]
    dx_1 = [0.040, 0.045, 0.052, 0.070, 0.090, 0.125, 0.170, 0.220, 0.250, 0.265, 0.270, 0.270, 0.270, 0.270, 0.270,
            0.270, 0.270]
    dx_2 = [0.040, 0.045, 0.052, 0.070, 0.090, 0.125, 0.170, 0.240, 0.320, 0.400, 0.450, 0.485, 0.507, 0.525, 0.540,
            0.550, 0.557]
    dx_3 = [0.040, 0.045, 0.052, 0.070, 0.090, 0.125, 0.170, 0.240, 0.320, 0.410, 0.480, 0.540, 0.576, 0.605, 0.630,
            0.650, 0.667]
    dx_4 = [0.040, 0.045, 0.052, 0.070, 0.090, 0.125, 0.170, 0.240, 0.320, 0.410, 0.490, 0.560, 0.615, 0.660, 0.690,
            0.720, 0.745]
    razm = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4]

    if (mah ** 2 - 1) <= 0:
        razmm = -kk.sqrt(1 - mah ** 2) / lambd_nos
    else:
        razmm = kk.sqrt(mah ** 2 - 1) / lambd_nos

    k = int(razmm // 0.2 + 5)

    param = lambd_korp / lambd_nos

    if param <= 0.5:
        dx = interpol(dx_05[k], dx_05[k - 1], procent(razmm, razm[k - 1], razm[k]))
    elif (param <= 1) and (param >= 0.5):
        dx = interpol(interpol(dx_1[k], dx_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                      interpol(dx_05[k], dx_05[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(param, 0.5, 1))
    elif (param <= 2) and (param >= 1):
        dx = interpol(interpol(dx_2[k], dx_2[k - 1], procent(razmm, razm[k - 1], razm[k])),
                      interpol(dx_1[k], dx_1[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(param, 1, 2))
    elif (param <= 3) and (param >= 2):
        dx = interpol(interpol(dx_3[k], dx_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                      interpol(dx_2[k], dx_2[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(param, 2, 3))
    elif (param <= 4) and (param >= 3):
        dx = interpol(interpol(dx_4[k], dx_4[k - 1], procent(razmm, razm[k - 1], razm[k])),
                      interpol(dx_3[k], dx_3[k - 1], procent(razmm, razm[k - 1], razm[k])), procent(param, 3, 4))
    else:
        dx = interpol(dx_4[k], dx_4[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return dx * l_nos


def tab_5_8(*args):
    """
    Определение координаты фокуса изолированных крыльев
    :param args: Число маха, относительный размах крыльев, тангенс угла средней стреловидности,
    относительное сужение крыльев
    :return: координата фокуса изолированных крыльев
    """

    mah = args[0]
    lambd_k = args[1]
    tan_05 = args[2]
    nu_k = args[3]

    xf_0_1 = [0.230, 0.230, 0.228, 0.221, 0.211, 0.193, 0.180, 0.163, 0.143, 0.140, 0.160, 0.200, 0.239, 0.274, 0.346,
              0.395, 0.430, 0.441, 0.448, 0.448, 0.449]
    xf_0_3 = [0.245, 0.245, 0.244, 0.238, 0.230, 0.220, 0.212, 0.200, 0.196, 0.204, 0.235, 0.270, 0.305, 0.332, 0.380,
              0.412, 0.440, 0.448, 0.451, 0.453, 0.454]
    xf_0_5 = [0.250, 0.252, 0.254, 0.252, 0.250, 0.241, 0.237, 0.233, 0.234, 0.250, 0.275, 0.310, 0.335, 0.355, 0.395,
              0.425, 0.447, 0.451, 0.454, 0.457, 0.458]
    xf_0_inf = [0.259, 0.262, 0.264, 0.266, 0.264, 0.261, 0.258, 0.252, 0.260, 0.276, 0.305, 0.335, 0.355, 0.378, 0.410,
                0.435, 0.453, 0.459, 0.464, 0.468, 0.470]

    xf_1_1 = [0.220, 0.218, 0.215, 0.208, 0.200, 0.186, 0.177, 0.170, 0.165, 0.170, 0.184, 0.215, 0.250, 0.290, 0.346,
              0.390, 0.418, 0.435, 0.450, 0.455, 0.458]
    xf_1_3 = [0.250, 0.254, 0.256, 0.256, 0.255, 0.252, 0.250, 0.252, 0.260, 0.272, 0.295, 0.320, 0.340, 0.360, 0.392,
              0.415, 0.433, 0.446, 0.455, 0.458, 0.460]
    xf_1_5 = [0.275, 0.272, 0.274, 0.276, 0.285, 0.298, 0.306, 0.315, 0.327, 0.340, 0.359, 0.372, 0.390, 0.402, 0.425,
              0.435, 0.445, 0.449, 0.456, 0.459, 0.461]
    xf_1_inf = [0.282, 0.284, 0.288, 0.295, 0.305, 0.322, 0.333, 0.345, 0.360, 0.375, 0.390, 0.400, 0.411, 0.420, 0.435,
                0.446, 0.450, 0.452, 0.458, 0.460, 0.462]

    xf_2_1 = [0.215, 0.220, 0.218, 0.214, 0.209, 0.200, 0.195, 0.196, 0.202, 0.213, 0.242, 0.275, 0.306, 0.335, 0.388,
              0.420, 0.450, 0.460, 0.465, 0.466, 0.469]
    xf_2_3 = [0.280, 0.280, 0.282, 0.284, 0.290, 0.300, 0.306, 0.314, 0.330, 0.348, 0.370, 0.390, 0.410, 0.425, 0.454,
              0.470, 0.482, 0.484, 0.484, 0.484, 0.484]
    xf_2_5 = [0.310, 0.315, 0.318, 0.327, 0.340, 0.355, 0.365, 0.378, 0.390, 0.405, 0.420, 0.437, 0.450, 0.460, 0.477,
              0.485, 0.487, 0.488, 0.488, 0.488, 0.488]
    xf_2_inf = [0.326, 0.333, 0.342, 0.355, 0.367, 0.387, 0.399, 0.410, 0.422, 0.440, 0.450, 0.460, 0.470, 0.480, 0.489,
                0.491, 0.491, 0.492, 0.493, 0.492, 0.491]

    xf_3_1 = [0.207, 0.209, 0.207, 0.200, 0.190, 0.192, 0.197, 0.200, 0.210, 0.225, 0.250, 0.290, 0.325, 0.357, 0.410,
              0.450, 0.475, 0.490, 0.500, 0.504, 0.505]
    xf_3_3 = [0.300, 0.300, 0.300, 0.305, 0.312, 0.327, 0.335, 0.348, 0.363, 0.380, 0.397, 0.415, 0.430, 0.446, 0.470,
              0.485, 0.500, 0.507, 0.510, 0.510, 0.510]
    xf_3_5 = [0.349, 0.350, 0.357, 0.368, 0.380, 0.397, 0.408, 0.418, 0.433, 0.445, 0.460, 0.473, 0.485, 0.490, 0.505,
              0.514, 0.520, 0.519, 0.517, 0.516, 0.515]
    xf_3_inf = [0.368, 0.375, 0.382, 0.396, 0.410, 0.430, 0.440, 0.450, 0.462, 0.472, 0.485, 0.492, 0.500, 0.507, 0.517,
                0.525, 0.528, 0.527, 0.525, 0.524, 0.523]

    razm = [-4, -3.5, -3, -2.5, -2, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

    if (mah ** 2 - 1) <= 0:
        razmm = -kk.sqrt(1 - mah ** 2) * lambd_k
    else:
        razmm = kk.sqrt(mah ** 2 - 1) * lambd_k

    param = lambd_k * tan_05
    k = 0
    for i in range(len(razm)):
        if (razmm <= razm[i]) and (razmm >= razm[i - 1]):
            k = i

    if param == 0:
        if nu_k == 1:
            x_f = interpol(xf_0_1[k], xf_0_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
        elif (nu_k <= 3) and (nu_k >= 1):
            x_f = interpol(interpol(xf_0_3[k], xf_0_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_0_1[k], xf_0_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 1, 3))
        elif (nu_k <= 5) and(nu_k >= 3):
            x_f = interpol(interpol(xf_0_5[k], xf_0_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_0_3[k], xf_0_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 3, 5))
        elif (nu_k <= 1000) and (nu_k >= 5):
            x_f = interpol(interpol(xf_0_inf[k], xf_0_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_0_5[k], xf_0_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 5, 1000))
        else:
            x_f = interpol(xf_0_inf[k], xf_0_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))

    elif (param <= 1) and (param >= 0):
        proc1 = procent(param, 0, 1)
        if nu_k == 1:
            x_f1 = interpol(xf_0_1[k], xf_0_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_1_1[k], xf_1_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 3) and (nu_k >= 1):
            x_f1 = interpol(interpol(xf_0_3[k], xf_0_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_0_1[k], xf_0_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f2 = interpol(interpol(xf_1_3[k], xf_1_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_1[k], xf_1_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 5) and (nu_k >= 3):
            x_f1 = interpol(interpol(xf_0_5[k], xf_0_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_0_3[k], xf_0_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f2 = interpol(interpol(xf_1_5[k], xf_1_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_3[k], xf_1_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 1000) and (nu_k >= 5):
            x_f1 = interpol(interpol(xf_0_inf[k], xf_0_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_0_5[k], xf_0_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f2 = interpol(interpol(xf_1_inf[k], xf_1_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_5[k], xf_1_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f = interpol(x_f2, x_f1, proc1)
        else:
            x_f1 = interpol(xf_0_inf[k], xf_0_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_1_inf[k], xf_1_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)

    elif (param <= 2) and (param >= 1):
        proc1 = procent(param, 1, 2)
        if nu_k == 1:
            x_f1 = interpol(xf_1_1[k], xf_1_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_2_1[k], xf_2_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 3) and (nu_k >= 1):
            x_f1 = interpol(interpol(xf_1_3[k], xf_1_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_1[k], xf_1_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f2 = interpol(interpol(xf_2_3[k], xf_2_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_1[k], xf_2_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 5) and (nu_k >= 3):
            x_f1 = interpol(interpol(xf_1_5[k], xf_1_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_3[k], xf_1_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f2 = interpol(interpol(xf_2_5[k], xf_2_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_3[k], xf_2_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 1000) and (nu_k >= 5):
            x_f1 = interpol(interpol(xf_1_inf[k], xf_1_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_1_5[k], xf_1_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f2 = interpol(interpol(xf_2_inf[k], xf_2_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_5[k], xf_2_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f = interpol(x_f2, x_f1, proc1)
        else:
            x_f1 = interpol(xf_1_inf[k], xf_1_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_2_inf[k], xf_2_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)

    elif (param <= 3) and (param >= 2):
        proc1 = procent(param, 2, 3)
        if nu_k == 1:
            x_f1 = interpol(xf_2_1[k], xf_2_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_3_1[k], xf_3_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 3) and (nu_k >= 1):
            x_f1 = interpol(interpol(xf_2_3[k], xf_2_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_1[k], xf_2_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f2 = interpol(interpol(xf_3_3[k], xf_3_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_3_1[k], xf_3_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 1, 3))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 5) and (nu_k >= 3):
            x_f1 = interpol(interpol(xf_2_5[k], xf_2_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_3[k], xf_2_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f2 = interpol(interpol(xf_3_5[k], xf_3_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_3_3[k], xf_3_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 3, 5))
            x_f = interpol(x_f2, x_f1, proc1)
        elif (nu_k <= 1000) and (nu_k >= 5):
            x_f1 = interpol(interpol(xf_2_inf[k], xf_2_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_2_5[k], xf_2_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f2 = interpol(interpol(xf_3_inf[k], xf_3_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            interpol(xf_3_5[k], xf_3_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                            procent(nu_k, 5, 1000))
            x_f = interpol(x_f2, x_f1, proc1)
        else:
            x_f1 = interpol(xf_2_inf[k], xf_2_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f2 = interpol(xf_3_inf[k], xf_3_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))
            x_f = interpol(x_f2, x_f1, proc1)

    else:
        if nu_k == 1:
            x_f = interpol(xf_3_1[k], xf_3_1[k - 1], procent(razmm, razm[k - 1], razm[k]))
        elif (nu_k <= 3) and (nu_k >= 1):
            x_f = interpol(interpol(xf_3_3[k], xf_3_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_3_1[k], xf_3_1[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 1, 3))
        elif (nu_k <= 5) and(nu_k >= 3):
            x_f = interpol(interpol(xf_3_5[k], xf_3_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_3_3[k], xf_3_3[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 3, 5))
        elif (nu_k <= 1000) and (nu_k >= 5):
            x_f = interpol(interpol(xf_3_inf[k], xf_3_inf[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           interpol(xf_3_5[k], xf_3_5[k - 1], procent(razmm, razm[k - 1], razm[k])),
                           procent(nu_k, 5, 1000))
        else:
            x_f = interpol(xf_3_inf[k], xf_3_inf[k - 1], procent(razmm, razm[k - 1], razm[k]))

    return x_f


def tab_5_9(*args):

    mah = args[0]
    lambd_k = args[1]
    c_ = args[2]

    x_f_025 = []
    x_f_050 = []
    x_f_075 = []
    x_f_100 = []
    x_f_125 = []
    x_f_150 = [0.21, 0.225, 0.21, 0.225, 0.35]
    x_f_185 = [0.21, 0.21, 0.21, 0.25, 0.2, 0.325, 0.37, 0.38, 0.39, 0.4]

    razm = [-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1]



def tab_5_11(*args):
    """
    Определение расстояния между фокусом изолированного крыла и точкой приложения дополнительной нормальной силы консоли
    :param args: , размах консоли
    :return:
    """

    d_ = args[0]
    l_k = args[1]

    f1 = [0.0150, 0.0279, 0.0314, 0.0328, 0.0332, 0.0325, 0.0315, 0.0300, 0.0285, 0.0240, 0.0190, 0.0140, 0.0095,
          0.0050, 0.0005]
    d_tab = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

    if d_ <= 0.4:
        k = int(d_ // 0.05 + 1)
    else:
        k = int(d_ // 0.1 + 5)

    return interpol(f1[k], f1[k - 1], procent(d_, d_tab[k - 1], d_tab[k])) * l_k / 2


def tab_atm(*args):
    """
    вывод параметров стандартной атмосферы при определенной высоте
    :param args: высота [м], выбор параметра (1 - Температура [к], 2 - местная скорость звука [м / с],
    3 - давление [Па], 4 - плотность [кг / м^3], 5 - кинематическая вязкость [м ^ 2 / с])
    :return:
    """

    h = args[0]
    param = args[1]

    h_tab = [-2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]

    t_atm = [301.2, 297.9, 294.7, 291.4, 288.2, 284.9, 281.7, 278.4, 275.2, 271.9, 268.7, 262.2, 255.7, 249.2, 242.7,
             236.2, 229.7, 223.3]
    a_atm = [347.9, 346.0, 344.1, 342.2, 340.3, 338.4, 336.4, 334.5, 332.5, 330.6, 328.6, 324.6, 320.6, 316.5, 312.3,
             308.1, 303.9, 299.6]
    p_atm = [127783, 120696, 113931, 107478, 101330, 95464, 89877, 84559, 79499, 74690, 70123, 61661, 54052, 47217,
             41106, 35653, 30801, 26500]
    ro_atm = [1.48, 1.41, 1.35, 1.28, 1.23, 1.17, 1.11, 1.06, 1.01, 0.96, 0.91, 0.82, 0.74, 0.66, 0.59, 0.53, 0.47,
              0.41]
    ni_atm = [1.253 * 10 ** -5, 1.301 * 10 ** -5, 1.352 * 10 ** -5, 1.405 * 10 ** -5, 1.46 * 10 ** -5, 1.52 * 10 ** -5,
              1.58 * 10 ** -5, 1.65 * 10 ** -5, 1.71 * 10 ** -5, 1.79 * 10 ** -5, 1.86 * 10 ** -5, 2.03 * 10 ** -5,
              2.21 * 10 ** -5, 2.42 * 10 ** -5, 2.65 * 10 ** -5, 2.9 * 10 ** -5, 3.2 * 10 ** -5, 3.53 * 10 ** -5]

    if h <= 3000:
        k = int(h // 500 + 5)
    else:
        k = int(h // 1000 + 8)

    if param == 1:
        return interpol(t_atm[k], t_atm[k - 1], procent(h, h_tab[k - 1], h_tab[k]))
    elif param == 2:
        return interpol(a_atm[k], a_atm[k - 1], procent(h, h_tab[k - 1], h_tab[k]))
    elif param == 3:
        return interpol(p_atm[k], p_atm[k - 1], procent(h, h_tab[k - 1], h_tab[k]))
    elif param == 4:
        return interpol(ro_atm[k], ro_atm[k - 1], procent(h, h_tab[k - 1], h_tab[k]))
    elif param == 5:
        return interpol(ni_atm[k], ni_atm[k - 1], procent(h, h_tab[k - 1], h_tab[k]))
    else:
        print("Ошибка: неверное значение при выборе параметра")


def tab_int_ver(*args):

    """
    вывод значения функции Лапласа-Гаусса по "Бронштейн И.Н.,Семендяев К.А.Справочник по математике. М., "Наука", с.81"
    :param args:
    :return:
    """

    x = args[0]

    F_x = [0.0000, 0.0399, 0.0797, 0.1192, 0.1585, 0.1974, 0.2358, 0.2737, 0.3108, 0.3473, 0.3829, 0.4177,
           0.4515, 0.4843, 0.5161, 0.5467, 0.5763, 0.6047, 0.6319, 0.6579, 0.6827, 0.7063, 0.7287, 0.7499,
           0.7699, 0.7887, 0.8064, 0.8230, 0.8385, 0.8529, 0.8664, 0.8789, 0.8904, 0.9011, 0.9109, 0.9199,
           0.9281, 0.9357, 0.9426, 0.9488, 0.9545, 0.9756, 0.9876, 0.9940, 0.9973, 0.99953, 0.99994]
    xx = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
          0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75,
          1.80, 1.85, 1.90, 1.95, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00]
    if x <= 2:
        k = int(x // 0.05 + 1)
    elif (x > 2) and (x <= 3):
        k = int(x // 0.25 + 33)
    elif (x > 3) and (x < 4):
        k = int(x // 0.5 + 39)
    else:
        return 0.99999

    return interpol(F_x[k], F_x[k - 1], procent(x, xx[k - 1], xx[k]))
