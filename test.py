import Tabl as tb
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt


def graph_test_3_2():
    a2 = 1
    a3 = [0, 0.5, 1, 2, 3, 4]
    x1 = np.linspace(-0.8, 2.6, 49)
    y1 = [[], [], [], [], [], []]
    for i in range(len(a3)):
        a1 = 0.36
        k = 0
        while a1 < 2.8:
            y1[i].append(tb.tab_3_2(a1, a2, a3[i]))
            a1 += 0.05
            k += 1
        # print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([-0.8, 3, 0, 0.07])
    plt.show()


def graph_test_3_4():
    x1 = np.linspace(-0.8, 2.6, 52)
    a2 = [0, 1]
    y1 = [[], []]
    a3 = 1
    for i in range(len(a2)):
        a1 = 0.9
        k = 0
        while a1 <= 1.42:
            y1[i].append(tb.tab_3_4(a1, a2[i], a3))
            a1 += 0.01
            k += 1
        # print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([-0.4, 1.2, 0, 0.05])
    plt.show()


def graph_test_3_5():
    from math import sqrt
    import numpy as np

    def convert_foo_3_5(a1, a2, a3):
        lambd_k = 1.2
        if a1 < 0:
            while (a1 / lambd_k)**2 > 1:
                lambd_k *= 1.5
            mach = sqrt(1 - (a1 / lambd_k)**2)
        else:
            mach = sqrt((a1 / lambd_k)**2 + 1)
        c_ = (a2 / lambd_k)**3
        tg_khi_05 = a3 / lambd_k
        # print(lambd_k)
        return tb.tab_3_5(mach, lambd_k, c_, tg_khi_05) / lambd_k

    def plot_line(ax, a2, a3):
        xs = np.linspace(-3, 9, 1000)
        ys = [convert_foo_3_5(a1, a2, a3) for a1 in xs]
        ax.plot(xs, ys, label=f'a2 = {a2}')
        ax.grid(True)
        ax.legend()

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    a3 = 0
    ax = ax1
    plot_line(ax, 0.25, a3)
    plot_line(ax, 0.5, a3)
    plot_line(ax, 1, a3)
    plot_line(ax, 1.5, a3)
    plot_line(ax, 1.7, a3)

    a3 = 1
    ax = ax2
    plot_line(ax, 0.25, a3)
    plot_line(ax, 0.5, a3)
    plot_line(ax, 1, a3)
    plot_line(ax, 1.5, a3)
    plot_line(ax, 1.7, a3)

    a3 = 2
    ax = ax3
    plot_line(ax, 0.25, a3)
    plot_line(ax, 0.5, a3)
    plot_line(ax, 1, a3)
    plot_line(ax, 1.5, a3)
    plot_line(ax, 1.7, a3)

    a3 = 3
    ax = ax4
    plot_line(ax, 0.25, a3)
    plot_line(ax, 0.5, a3)
    plot_line(ax, 1, a3)
    plot_line(ax, 1.5, a3)
    plot_line(ax, 1.7, a3)

    plt.show()


def graph_test_3_21():
    a1 = [2, 3, 4, 5]
    y1 = [[], [], [], []]
    x1 = np.linspace(0.6, 3, 49)

    for i in range(len(a1)):
        a2 = 0.6
        k = 0

        while a2 <= 3:
            y1[i].append(tb.tab_3_21(a1[i], a2))
            a2 += 0.05
            k += 1
        plt.plot(x1, y1[i])
        print(k)
    plt.grid(True)
    plt.axis([0.6, 3, 0.6, 1])
    plt.show()


def graph_test_3_22():
    a2 = [1, 0.8, 0.6, 0.4, 0.2, 0]
    y1 = [[], [], [], [], [], []]
    x1 = np.linspace(0, 5, 51)

    for i in range(len(a2)):
        a1 = 0
        k = 0
        while a1 < 5:
            y1[i].append(tb.tab_3_22(a1, a2[i]))
            a1 += 0.1
            k += 1
        print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0, 5, 0.65, 1])
    plt.show()


def graph_test_4_2():
    a2 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1]
    x_n = 10 ** 6
    x_k = 3 * 10 ** 8
    x1 = np.linspace(x_n, x_k, 60)
    y1 = [[], [], [], [], [], [], [], []]
    for i in range(len(a2)):
        a1 = 10 ** 6
        k = 0
        while a1 <= 3 * 10 ** 8:
            y1[i].append(tb.tab_4_2(a1, a2[i]) * 2)
            a1 *= 1.1
            k += 1
        print(k)
        print(len(y1[i]))
        print(len(x1))
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([10 ** 6, 3 * 10 ** 8, 0, 0.009])
    plt.show()


def graph_test_4_3():

    a2 = [0, 0.2, 0.5, 0.6, 0.8, 1]
    y1 = [[], [], [], [], [], []]
    x1 = np.linspace(0, 6, 61)
    for i in range(len(a2)):
        a1 = 0
        k = 0
        while a1 < 6:
            y1[i].append(tb.tab_4_3(a1, a2[i]))
            a1 += 0.1
            k += 1
        print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0, 6, 0, 1])
    plt.show()


def graph_test_4_11():
    a2 = [1.5, 2, 2.5, 3, 4, 5]
    x1 = np.linspace(0.6, 5.5, 99)
    y1 = [[], [], [], [], [], []]

    for i in range(len(a2)):
        a1 = 0.6
        k = 0
        while a1 <= 5.5:
            y1[i].append(tb.tab_4_11(a1, a2[i]))
            a1 += 0.05
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0.6, 5.5, 0, 0.44])
    plt.show()


def graph_test_4_13():
    a2 = [0, 0.25, 0.5, 1, 2]
    y1 = [[], [], [], [], []]
    x1 = np.linspace(0.4, 3.6, 321)

    for i in range(len(a2)):
        a1 = 0.4
        k = 0
        while a1 <= 3.6:
            y1[i].append(tb.tab_4_13(a1, a2[i]))
            a1 += 0.01
            k += 1
        print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0.4, 3.6, 0, 1.8])
    plt.show()


def graph_test_4_24():
    a2 = [[0, 2], [0, 3], [0.5, 1.5], [0.5, 2], [0.5, 2.5], [0.75, 1], [0.75, 1.5], [0.75, 2]]
    x1 = np.linspace(0.6, 5.5, 50)
    y1 = [[], [], [], [], [], [], [], []]
    for i in range(len(a2)):
        a1 = 0.6
        k = 0
        while a1 <= 5.5:
            y1[i].append(tb.tab_4_24(a1, a2[i][0], a2[i][1]))
            a1 += 0.1
            k += 1
        print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0.6, 5.5, 0, 0.17])
    plt.show()


def graph_test_4_28():
    a1 = [0, 0.2, 0.4]
    y1 = [[], [], []]
    x1 = np.linspace(0, 0.12, 60)
    for i in range(len(a1)):
        a2 = 0
        k = 0
        while a2 <= 0.12:
            y1[i].append(tb.tab_4_28(a1[i], a2))
            a2 += 0.002
            k += 1
        print(k)
        plt.plot(x1, y1[i])
    plt.grid(True)
    plt.axis([0, 0.12, 1, 1.4])
    plt.show()


def graph_test_4_30():
    x1 = np.linspace(0, 7, 61)
    a2 = [1, 5, 1000]
    a3 = [1, 1, 1, 1, 1, 1, 1]
    a4 = [0, 0, 0, 1, 2, 3, 4]
    a5 = [0.125, 1, 3.375, 0, 0, 0, 0]
    y1 = [[], [], [], [], [], [], []]

    for i in range(len(a2)):
        for j in range(len(a3)):
            a1 = 1
            k = 0
            y1 = [[], [], [], [], [], [], []]
            while a1 <= 7.07:
                y1[j].append(tb.tab_4_30(a1, a2[i], a3[j], a4[j], a5[j]))
                a1 += 0.1
                k += 1
            print(k)
            plt.plot(x1, y1[j])
        plt.grid(True)
        plt.axis([0, 7, 0, 4.1])
        plt.show()


def graph_test_4_32():
    a1 = 1.42
    a2 = 1
    y1 = []
    x1 = np.linspace(0, 2.5, 45)
    k = 0
    while a1 < 3.64:
        y1.append(tb.tab_4_32(a1, a2))
        a1 += 0.05
        k += 1
    plt.plot(x1, y1)
    plt.grid(True)
    plt.axis([0, 2.5, 0, 1])
    print(k)
    plt.show()


def graph_test_4_40():
    a2 = 1
    a3 = [0, 1]
    x1 = np.linspace(-0.8, 2.8, 31)
    y1 = [[], []]
    for i in range(len(a3)):
        k = 0
        a1 = 0.44
        print(sqrt(1 - a1 ** 2) / 1)
        while a1 < 1.95:
            y1[i].append(tb.tab_4_40(a1, a2, a3[i]))
            a1 += 0.05
            k += 1
        plt.plot(x1, y1[i])
        print(k)
    plt.grid(True)
    plt.axis([-0.8, 2.8, -0.4, 0.8])
    plt.show()


def graph_test_5_7():
    a2 = 1
    a3 = [0.5, 1, 2, 3, 4]
    a4 = 1
    x1 = np.linspace(-0.8, 2.4, 41)
    y1 = [[], [], [], [], []]
    for i in range(len(a3)):
        a1 = 0.6
        k = 0
        while a1 < 2.6:
            y1[i].append(tb.tab_5_7(a1, a2, a3[i], a4))
            a1 += 0.05
            k += 1
        plt.plot(x1, y1[i])
        print(k)
    plt.grid(True)
    plt.axis([-0.8, 2.4, 0, 0.8])
    plt.show()


def graph_test_5_8():
    a2 = 4
    a3 = [0, 0.25, 0.5, 0.75]
    a4 = [1, 3, 5, 1000]
    x1 = np.linspace(-4, 4, 29)
    y1 = [[], [], [], []]

    for i in range(len(a3)):
        y1 = [[], [], [], []]
        for j in range(len(a4)):
            k = 0
            a1 = 0
            while a1 < 1.42:
                y1[j].append(tb.tab_5_8(a1, a2, a3[i], a4[j]))
                a1 += 0.05
                k += 1
            print(k)
            plt.plot(x1, y1[j])
        plt.grid(True)
        plt.axis([-4, 4, 0.1, 0.6])
        plt.show()


def graph_test_5_11():
    a1 = 0
    a2 = 1
    y1 = []
    x1 = np.linspace(0, 1, 50)
    k = 0
    while a1 <= 1:
        y1.append(tb.tab_5_11(a1, a2) / a2 * 2)
        k += 1
        a1 += 0.02
    plt.plot(x1, y1)
    print(k)
    plt.grid(True)
    plt.axis([0, 1, 0, 0.04])
    plt.show()


# graph_test_3_2()
# graph_test_3_4()
# graph_test_3_5()
# graph_test_3_21()
# graph_test_3_22()
graph_test_4_2()
# graph_test_4_3()
# graph_test_4_11()
# graph_test_4_13()
# graph_test_4_24()
# graph_test_4_28()
# graph_test_4_30()
# graph_test_4_32()
# graph_test_4_40()
# graph_test_5_7()
# graph_test_5_8()
# graph_test_5_11()

