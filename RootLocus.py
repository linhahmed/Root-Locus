import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import *


def center_of_asymptotes(poles):
    n = len(poles)
    return (np.sum(poles)) / (n)


def angle_of_asymptote(n, h):
    angle = wrapangle((180 * (2 * h + 1)) / (n))
    return angle


def angles_of_all_asymptotes(poles):
    n = len(poles)
    num_asymptotes = n
    angles = []
    for h in range(0, num_asymptotes):
        angles.append(angle_of_asymptote(n, h))
    return np.array(angles)


def equation(poles, s):
    fun = s - poles[0]
    for i in range(1, len(poles)):
        fun *= (s - poles[i])
    eqn = simplify(fun)
    return eqn


def breakAwaypts(poles, s):
    eqn = equation(poles, s)
    derivative = sympy.diff(eqn, s)
    diff = simplify(derivative)
    sol = solve(diff)
    angle = 0
    angles = []
    thepoint = []
    for j in range(0, len(sol)):
        for i in range(0, len(poles)):
            k = -1 * calculateArg(sol[j] - poles[i])
            angle += k
        angles.append(angle)
    for h in range(0, len(angles)):
        if round(angles[h], 1) == -180 or round(angles[h], 1) == 180:
            thepoint.append(sol[h])
    return thepoint[0]


def calculateArg(number):
    number = complex(number)
    arg = 0
    x = number.real
    y = number.imag
    if x > 0 and round(y, 1) == 0:
        arg = 0
    elif x < 0 and round(y, 1) == 0:
        arg = 180
    elif y != 0:
        arg = 2 * np.degrees(np.arctan2(y, (((x ** 2 + y ** 2) ** (1 / 2)) + x)))
    return wrapangle(arg)


def routh(poles):
    s, k = symbols('s k')
    eqn = equation(poles, s) + k
    eqn = expand(eqn)
    arr = []
    a = Poly(eqn, s)
    ar = a.coeffs()
    for i in range(0, len(ar)):
        arr.append([])
        if i % 2 == 0:
            arr[0].append(ar[i])
        else:
            arr[1].append(ar[i])
    arr[1].append(0)
    o = 2
    for l in range(1, len(arr) - 1):
        for j in range(0, len(arr[l]) - 1):
            if arr[l][j] == 0:
                arr[o].append(0)
            else:
                arr[o].append((arr[l][j] * arr[l - 1][j + 1] - arr[l - 1][j] * arr[l][j + 1]) / arr[l][j])
        arr[o].append(0)
        o += 1
    sol = solve(arr[len(arr) - 2][0])
    fun = arr[2][0] * s ** 2 + sol[0]
    solfun = solve(fun)
    return solfun


def wrapangle(angle):
    if angle >= 180:
        return angle - 360
    else:
        return angle % 360


def get_angle(point1, point2):
    delta_real = point2.real - point1.real
    delta_j = point2.imag - point1.imag
    angle = np.arctan2(delta_j, delta_real)
    return wrapangle(np.degrees(angle))


def angleofDeparture(poles, index):
    p = poles[index]
    pole_angles = [get_angle(pole, p) for pole in poles if pole != p]
    return wrapangle(180 - np.sum(pole_angles))


def getallanglesofDeparture(poles):
    angles = []
    for index in range(len(poles)):
        angles.append(angleofDeparture(poles, index))
    return angles


def drawPoles(poles):
    fig, ax = plt.subplots()
    ax.set_xlabel('Re')
    ax.set_ylabel('Im')
    ax.axvline(x=0, color='g', lw=0.5)
    ax.grid(True, which='both')
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    drawPolespts(poles)
    drawAsymptoteline(poles)
    drawbreakAwaypts(poles)
    drawlinereal(poles)
    drawRouthpts(poles)
    drawAngleofDepatureline(poles)
    drawcurvereal(poles)
    drawcurveimag(poles)
    plt.show()


def drawPolespts(poles):
    x = []
    y = []
    for i in range(0, len(poles)):
        ele = complex(poles[i])
        x.append(ele.real)
        y.append(ele.imag)
    plt.scatter(x, y, marker='x', color='r', lw=0.5)


def drawAsymptoteline(poles):
    angles = angles_of_all_asymptotes(poles)
    elo = complex(center_of_asymptotes(poles))
    x, y = symbols('x y')
    for i in angles:
        a = np.tan(i * np.pi / 180)
        b = elo.imag - a * elo.real
        xx = np.linspace(-125, 75, 100)
        yy = a * xx + b
        plt.plot(xx, yy, linestyle=':', lw=0.5, color='b')


def drawbreakAwaypts(poles):
    s = symbols('s')
    pt = complex(breakAwaypts(poles, s))
    plt.scatter(pt.real, pt.imag, marker='o', color='r', lw=0.5)


def drawlinereal(poles):
    point1 = [poles[0], 0]
    point2 = [poles[1], 0]
    x_values = [point1[0], point2[0]]
    y_values = [point1[1], point2[1]]
    plt.plot(x_values, y_values, lw=2, color='b')


def drawRouthpts(poles):
    routhpt = routh(poles)
    x = []
    y = []
    for i in routhpt:
        i = complex(i)
        x.append(i.real)
        y.append(i.imag)
    plt.scatter(x, y, marker='o', color='r', lw=0.1)


def drawAngleofDepatureline(poles):
    # returns the inretsected point with x axis
    angles = getallanglesofDeparture(poles)
    for i in range(2, len(angles)):
        a = np.tan(angles[i] * np.pi / 180)
        elo = complex(poles[i])
        b = elo.imag - a * elo.real
        xx = np.linspace(-100, -50, 100)
        yy = a * xx + b
        plt.plot(xx, yy, linestyle='-', lw=0.5, color='g')


def drawcurvereal(poles):
    h = complex(center_of_asymptotes(poles)).real
    s = symbols('s')
    v = complex(breakAwaypts(poles, s)).real
    a = v - h
    b = h - v
    y = np.arange(-100, 100, 0.01)
    temp = 1 + (y ** 2 / b ** 2)
    temp = temp * a ** 2
    x = temp ** (1 / 2) + h
    plt.plot(x, y, lw=1)


def drawcurveimag(poles):
    h = complex(center_of_asymptotes(poles)).real
    s = symbols('s')
    pole1 = complex(poles[2])
    a = ((pole1.real - h) ** 2 - pole1.imag ** 2) ** (1 / 2)
    b = ((pole1.real - h) ** 2 - pole1.imag ** 2) ** (1 / 2)
    y = np.arange(10, 75, 0.01)
    temp = 1 + ((y) ** 2 / b ** 2)
    temp = temp * a ** 2
    x = -temp ** (1 / 2) + h
    plt.plot(x, y, lw=1, color='y')
    plt.plot(x, -1 * y, lw=1, color='y')


poles = [0, -25, -50 + 10j, -50 - 10j]
print("Center of Asymtotes is:")
print(center_of_asymptotes(poles))
print("Angles of Asymtotes are:")
print(angles_of_all_asymptotes(poles))
s = symbols('s')
print("The eqution is:")
print(equation(poles, s))
print("break Away points are:")
print(breakAwaypts(poles, s))
print("Imaginary crossings points are:")
print(routh(poles))
print("Angles of Depature are:")
print(getallanglesofDeparture(poles))
drawPoles(poles)
