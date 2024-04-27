import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi


plt.rcParams["text.usetex"] = True
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}

plt.rcParams.update(tex_fonts)
plt.rc('text', usetex=True)







def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





epsilon = 60

#3.2
# erster punkt, 2 messungen für rechts und links
delta_min1 = 33.5 #grad
aichung1 = -5.4
delta_min2 = 291.6
aichung2 = 330.4

delta_min1 = delta_min1 - aichung1
delta_min2 = abs(delta_min2 - aichung2)

n1 = sin(2*pi/360*((delta_min1 + epsilon)/2))/sin(2*pi/360*(epsilon/2))
sigma_n1 = gauss("sin(2*pi/360*((delta_min1 + epsilon)/2))/sin(2*pi/360*(epsilon/2))")


#2. punkt
Nullpunkt = 339
gelb = 19.4

gelb = gelb + 360 - Nullpunkt

#3. punkt
#Nullpunkt wie bei punkt 2

#Prisma I

W1_Hg = matrix("""
2 19.3;
3 19.4;
5 19.5;
7 19.7;
8 20.2;
9 20.4;
10 20.45
""")# Nr. aus der tabelle, Wellenlänge,gemessener winkel

#Prisma II
W2_Hg = matrix("""
2 623.4 27.8;
3 578.0 28.2;
5 546.1 28.4;
6 499.2 29;
7 491.6 29.1;
8 435.8 30.1;
9 407.8 30.7;
10 404.7 30.8

""")

#Prisma III
W3_Hg = matrix("""
2 623.8 38.2;
3 579.1 38.88;
4 577 38.9;
5 546.1 39.4;
6 499.2 40.5;
7 491.6 40.6;
8 435.8 42.7;
9 407.8 44.3;
10 404.7 44.5
""")

#4,Punkt Helium Lampe, Prisma III
#rot, gelb, zyan, zyan, violett, violett
W_He = matrix("""
37.8;
38.7;
40.4;
40.6;
42.25;
42.5
""")