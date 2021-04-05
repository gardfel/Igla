import random
import Tabl
import Aero
import math as kk
import matplotlib.pyplot as plt
import time
import numpy
import scipy


"""x = time.time()

dt = 10 ** -4
c = 335
v = 51
a = 249.7
a1 = 29
M = v / c

l_cil = 45.5 / 8
l_nos = 12.5 / 8
tan_05 = 0.307
l_kr = 1.46
c_k = 0.03
t = 0
tt = [0]
v = 1
Mah = v / 334
cy = [Tabl.tab_3_5(Mah, l_kr, c_k, tan_05)]
"""

'''if (Mah ** 2 - 1) <= 0:
    razmm = [-l_kr * (-Mah ** 2 + 1) ** 0.5]
else:
    razmm = [l_kr * (Mah ** 2 - 1) ** 0.5]

while v < 800:
    Mah = v / 334
    cyy = Tabl.tab_3_5(Mah, l_kr, c_k, tan_05)
    cy.append(cyy)
    if (Mah ** 2 - 1) <= 0:
        razm = -l_kr * (-Mah ** 2 + 1) ** 0.5
    else:
        razm = l_kr * (Mah ** 2 - 1) ** 0.5
    razmm.append(razm)
    v += 1
plt.plot(razmm, cy)
plt.show()
'''
'''
Cyy_nos = [Aero.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil), Tabl.tab_3_2(M, 1.32, 22.43),
                              Tabl.tab_3_4(M, 0, 22.43), Tabl.tab_3_4(M, 1, 23.75))]
Cyy_kr = [Tabl.tab_3_5(M, l_kr, c_k, tan_05)]

while v < 887:
    M = v / c
    Cy_nos = Aero.Cy_nos(Tabl.tab_3_2(M, l_nos, l_cil), Tabl.tab_3_4(M, 0, l_cil), Tabl.tab_3_2(M, 1.32, 22.43),
                                Tabl.tab_3_4(M, 0, 22.43), Tabl.tab_3_4(M, 1, 23.75))
    Cyy_nos.append(Cy_nos)
    Cy_kr = Tabl.tab_3_5(M, l_kr, c_k, tan_05)
    Cyy_kr.append(Cy_kr)
    t += dt
    tt.append(t)
    if v < 650:
        v += a * dt
    elif v < 887:
        v += a1 * dt
    print(t)

print(time.time() - x)
plt.plot(tt, Cyy_nos)
plt.axis([-0.1, t, 0, Cy_nos + 0.01])
plt.grid(True)
plt.show()
plt.plot(tt, Cyy_kr)
plt.axis([-0.1, t, 0, 0.04])
plt.grid(True)
plt.show()'''
"""
print(180 * kk.pi / 180)
print(Tabl.tab_3_5(2.648, 7.888, 0.03, -0.012))
M = 2
print()
h = [40e-06, 20 * 10 ** -6, 10 * 10 ** -6, 6.3 * 10 ** -6, 3.2 * 10 ** -6, 1.6 * 10 ** -6]
print(int(2.01 // 0.05 + 1))
print(int(3.01 // 0.25 + 33))
print(int(3.01 // 0.5 + 39))
y1 = [26.85, 23.90, 21.16, 18.70, 16.49, 14.33, 12.11, 10, 8.1, 6.50]
for i in range(len(y1)):
    y1[i] = y1[i] * 0.002 / 8.6918 + 0.006
print(y1)
mach_per = kk.sqrt((10 / 7.888) ** 2 + 1)
print(mach_per, "mach")
cy_per = 4 / (57.3 * kk.sqrt(mach_per ** 2 - 1)) / 7.888
print(cy_per, "cy_per")
"""
'''
Cy_iz = [[-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6],
         [0.0350, 0.0350, 0.0350, 0.0348, 0.0338, 0.0335, 0.0333, 0.0327, 0.0323, 0.0321, 0.0321, 0.0321, 0.0320, 0.0320, 0.0320, 0.0320, 0.0320, 0.0320],
         [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0426, 0.0455, 0.0465, 0.0460, 0.0448, 0.0435, 0.0425, 0.0417, 0.0410, 0.0402, 0.0400, 0.0395, 0.0390],
         [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0505, 0.0512, 0.0507, 0.0500, 0.0490, 0.0480, 0.0473, 0.0465, 0.0460, 0.0450, 0.0448],
         [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0540, 0.0560, 0.0565, 0.0560, 0.0555, 0.0550, 0.0545, 0.0537, 0.0530, 0.0525],
         [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0548, 0.0570, 0.0580, 0.0585, 0.0587, 0.0587, 0.0585, 0.0583, 0.0580, 0.0575],
         [0.0350, 0.0350, 0.0350, 0.0358, 0.0380, 0.0430, 0.0477, 0.0517, 0.0548, 0.0570, 0.0584, 0.0594, 0.0599, 0.0600, 0.0600, 0.0600, 0.0600, 0.0600]]

otnos = [0, 0.5, 1, 2, 3, 4]
print(Cy_iz[3][8])
if (M ** 2 - 1) >= 0:
    razmm = kk.sqrt(M ** 2 - 1) / l_nos
else:
    razmm = -kk.sqrt(1 - M ** 2) / l_nos
k = int(razmm // 0.2 + 5)

otn = l_cil / l_nos

if otn <= 1:
    jj = int(otn // 0.5 + 2)
else:
    jj = int(otn // 1 + 3)

Cy = Tabl.interpol(Tabl.interpol(Cy_iz[jj][k], Cy_iz[jj][k - 1], Tabl.procent(razmm, Cy_iz[0][k - 1], Cy_iz[0][k])), Tabl.interpol(Cy_iz[jj - 1][k], Cy_iz[jj - 1][k - 1], Tabl.procent(razmm, Cy_iz[0][k - 1], Cy_iz[0][k])), Tabl.procent(otn, otnos[jj - 2], otnos[jj - 1]))
print("new", Cy)
print("old ", Tabl.tab_3_2(M, l_nos, l_cil))
'''
# Корпус:
# const:

d = 0.072
d_korm = 0.05

# Крылья:
# const:

tan_1_kr = 0

# var:

L_kr = 0.146  # удлинение крыльев
S_kr = 0.0146  # площадь
b_0_kr = 0.12  # корневая хорда

# вычисления формы:

l_kr = L_kr ** 2 / S_kr  # относительное удлинение крыльев
nu_kr = ((S_kr / (L_kr * b_0_kr) - 0.5) ** (-1) / 2)  # относительное сужение крыльев
b_1_kr = b_0_kr / nu_kr  # концевая хорда крыльев
print((b_0_kr - b_1_kr) / (L_kr / 2) - 2 / l_kr * ((nu_kr - 1)/(nu_kr + 1)))
b_a_kr = 4 / 3 * S_kr / L_kr * (1 - (nu_kr / (nu_kr + 1) ** 2))  # САХ
print(b_a_kr)
z_a_kr = L_kr / 6 * ((nu_kr + 2) / (nu_kr + 1))  # расстояние от САХ до оси ЛА
# консоли:
b_kr = b_0_kr * (1 - (nu_kr - 1) / nu_kr * d_korm / L_kr)  # бортовая хорда крыльев
S_k_kr = S_kr * (1 - ((nu_kr - 1) / (nu_kr + 1)) * d_korm / L_kr) * (1 - d_korm / L_kr)  # площадь консолей крыльев
l_k_kr = l_kr * ((1 - d_korm / L_kr) / (1 - ((nu_kr - 1) / (nu_kr + 1) * d_korm / L_kr)))  # относительное удлинение консолей
L_k_kr = L_kr - d_korm  # размах консолей крыльев
nu_k_kr = nu_kr - d_korm / L_kr * (nu_kr - 1)
tan_05_kr = tan_1_kr + 2 / l_k_kr * (nu_k_kr - 1) / (nu_k_kr + 1)  # средняя стреловидность
b_a_k_kr = 4 / 3 * S_k_kr / L_k_kr * (1 - nu_k_kr / (nu_k_kr + 1) ** 2)  # САХ консоли
z_a_k_kr = L_k_kr / 6 * (nu_k_kr + 2) / (nu_k_kr + 1)


# tan_0_kr = tan_1_kr + 4 / l_kr *

print(tan_05_kr)

print(L_k_kr)
print()
print(S_k_kr, L_k_kr)
print(b_kr)
print(b_1_kr)
print(nu_kr)
print(l_kr)


# Рули:
# const:

tan_0_op = 0

# var:

L_op = 0.25  # удлинение рулей
S_op = 0.0078575  # площадь рулей
b_0_op = 0.033  # корневая хорда рулей

# вычисления формы:

l_op = L_op ** 2 / S_op  # относительное удлинение крыльев
nu_op = ((S_op / (L_op * b_0_op) - 0.5) ** (-1) / 2)  # относительное сужение крыльев
b_1_op = b_0_op / nu_op  # концевая хорда крыльев
print((b_0_op - b_1_op) / (L_op / 2) - 2 / l_op * ((nu_op - 1)/(nu_op + 1)))
b_a_op = 4 / 3 * S_op / L_op * (1 - (nu_op / (nu_op + 1) ** 2))  # САХ
print(b_a_op)
z_a_op = L_op / 6 * ((nu_op + 2) / (nu_op + 1))  # расстояние от САХ до оси ЛА
# консоли:
b_op = b_0_op * (1 - (nu_op - 1) / nu_op * d / L_op)  # бортовая хорда крыльев
S_k_op = S_op * (1 - ((nu_op - 1) / (nu_op + 1)) * d / L_op) * (1 - d / L_op)  # площадь консолей крыльев
l_k_op = l_op * ((1 - d / L_op) / (1 - ((nu_op - 1) / (nu_op + 1) * d / L_op)))  # относительное удлинение консолей
L_k_op = L_op - d  # размах консолей крыльев
nu_k_op = nu_op - d / L_op * (nu_op - 1)
tan_05_op = tan_0_op - 2 / l_k_op * (nu_k_op - 1) / (nu_k_op + 1)  # средняя стреловидность
b_a_k_op = 4 / 3 * S_k_op / L_k_op * (1 - nu_k_op / (nu_k_op + 1) ** 2)  # САХ консоли
z_a_k_op = L_k_op / 6 * (nu_k_op + 2) / (nu_k_op + 1)
print(tan_05_op)

print(L_k_op)
print()
print(S_k_op, L_k_op)
print(b_op)
print(b_1_op)
print(nu_op)
print(l_op)
