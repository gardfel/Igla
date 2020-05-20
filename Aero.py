import math as kk


def Cy_nos(*args):

    """
    Расчет Су головной части
    :param args:
    :return:
    """

    cy1 = args[0]
    cy1_1 = args[1]
    cy2 = args[2]
    cy2_1 = args[3]
    cy3 = args[4]

    r = 4
    r1 = 32.86
    d = 14
    d1 = 8
    D = 72

    r_ = 2 * r / d
    r_1 = 2 * r1 / D
    d_ = d1 / D


    Cy1_2 = cy1 * (1 - r_ ** 2) + cy1_1 * r_ ** 2
    Cy2_2 = cy2 * (1 - r_1 ** 2) + cy2_1 * r_1 ** 2
    Cy_ob = Cy2_2 * (1 - d_ ** 2) + cy3 * d_ ** 2
    Cy_iz = Cy1_2 * (2 * r / D) ** 2 + Cy_ob

    return Cy_iz


"""def Cy_(*args):"""



