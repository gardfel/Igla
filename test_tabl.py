from pytest import approx
from Tabl import *
from math import *
import random



def test_tab_3_5():
    def convert_foo_3_5(a1, a2, a3):
        lambd_k = 1.1
        while a1 >= lambd_k:
            lambd_k *= 1.5 
        mach = sqrt(1 - (a1 / lambd_k)**2) 
        c_ = (a2 / lambd_k)**3 
        tg_khi_05 = a3 / lambd_k 
        return tab_3_5(mach, lambd_k, c_, tg_khi_05) / lambd_k

    db = [(0.019,0,1.6,1.5),
        (0.029,0,1.1,0.4),
        (0.014,0,1.7,4),
        (0.028,0,0.3,1.7),
        (0.021,1,0.9,2.2),
        (0.029,1,0.2,1.5),
        (0.024,1,1.5,-1.2),
        (0.021,1,1.7,0.6),
        (0.014,2,1.6,3.5),
        (0.025,2,0.7,0.8),
        (0.022,2,1.4,-1),
        (0.017,2,1.6,-2.5),
        (0.018,3,1,-1.7),
        (0.011,3,1.7,5.2),
        (0.017,3,0.9,2.4),
        (0.026,3,0.2,1)]
    for answ, a1, a2, a3 in db:
        assert convert_foo_3_5(a1, a2, a3) == approx(answ, abs=0.005)


# def test_tab_3_2():
#     def convert_foo_3_2(a1, a2):
#         lambd_nos = random.uniform(1.2, 2.5)
#         mach = sqrt(1+(lambd_nos*a1)**2)
#         lambd_cil = a2 * lambd_nos
#         return tab_3_2(mach, lambd_nos, lambd_cil)
    
#     assert convert_foo_3_2() == approx(,abs=0.05)


