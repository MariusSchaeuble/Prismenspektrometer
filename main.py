import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
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


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)




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
sigma_epsilon = 0

#3.2
# erster punkt, 2 messungen für rechts und links
delta_min1 = 33.5 #grad
sigma_delta_min1 = 0.05
aichung1 = -5.4
delta_min2 = 291.6
sigma_delta_min2 = 0.05
aichung2 = 330.4

delta_min1 = delta_min1 - aichung1
delta_min2 = abs(delta_min2 - aichung2)

n1 = sin(2*pi/360*((delta_min1 + epsilon)/2))/sin(2*pi/360*(epsilon/2))
sigma_n1 = gauss("sin(2*(pi/360)*((delta_min1 + epsilon)/2))/sin(2*(pi/360)*(epsilon/2))")

n2 = sin(2*pi/360*((delta_min2 + epsilon)/2))/sin(2*pi/360*(epsilon/2))
sigma_n2 = gauss("sin(2*(pi/360)*((delta_min2 + epsilon)/2))/sin(2*(pi/360)*(epsilon/2))")

RCP(n1, sigma_n1)
RCP(n2, sigma_n2)

#2. punkt
Nullpunkt = 339
gelb = 19.4

gelb = gelb + 360 - Nullpunkt
gelb = 2*pi/360*gelb
sigma_gelb = 0.1
sigma_gelb = gauss("2*pi/360*gelb")
theta = 60
sigma_theta = 1
theta = 2*pi/360*theta
sigma_theta = gauss("2*pi/360*theta")

epsilon = 2*pi*epsilon/360
n_gelb = sqrt(2*cos(epsilon)*sin(epsilon + gelb - theta)*sin(theta) + sin(theta)**2 + sin(epsilon + gelb - theta)**2)/sin(epsilon)
sigma_n_gelb = gauss("sqrt(2*cos(epsilon)*sin(epsilon + gelb - theta)*sin(theta) + sin(theta)**2 + sin(epsilon + gelb - theta)**2)/sin(epsilon)")




#3. punkt
#Nullpunkt wie bei punkt 2

#Prisma I

def parabel(x, a, b, c):
    return a*x**2 + b*x + c


W1_Hg = matrix("""
2 623.4 19.3;
3 578 19.4;
5 546.1 19.5;
7 491.6 19.7;
8 435.8 20.2;
9 407.8 20.4;
10 404.7 20.45
""")# Nr. aus der tabelle, Wellenlänge,gemessener winkel

delta1 = 360 + W1_Hg[:, 2] - Nullpunkt
delta1 = delta1*2*pi/360
n1_Hg = ones(len(delta1))
for i in range(len(delta1)):
    n1_Hg[i] = (sqrt(2*cos(epsilon)*sin(epsilon + delta1[i, 0] - theta)*sin(theta) + sin(theta)**2 + sin(epsilon + delta1[i, 0] - theta)**2)/sin(epsilon))

lamda1_Hg = toArray(W1_Hg[:, 1])

errorbar(lamda1_Hg, n1_Hg, sigma_n_gelb*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, n1_Hg,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Brechungsindex', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('n1_Hg.png')
show()

al1 = optimizedParameters[0]
bl1 = optimizedParameters[1]
opn1 = optimizedParameters

delta1 = toArray(360 + W1_Hg[:, 2] - Nullpunkt)
errorbar(lamda1_Hg, delta1, 1*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Ablenkwinkel in °', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('delta1.png')
show()

ad1 = optimizedParameters[0]*2*pi/360
bd1 = optimizedParameters[1]*2*pi/360

def ableitungsQuozient(x):
    return (2*ad1*x + bd1)/(2*al1*x + bl1)
    #return bd1/bl1

Term11 = ones(len(n1_Hg))
sigma_Term11 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term11[i] = 2*sin(epsilon/2)/sqrt(1 - n1_Hg[i]**2*sin(epsilon/2)**2)
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term11[i] = gauss("2*sin(epsilon/2)/sqrt(1 - temp1**2*sin(epsilon/2)**2)")


Term12 = ones(len(n1_Hg))
sigma_Term12 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term12[i] = sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term12[i] = gauss("sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))".replace("n1_Hg[i]", "temp1"))


errorbar(lamda1_Hg, Term11, sigma_Term11, None,'x', label='Genäherter Term')
errorbar(lamda1_Hg, Term12, sigma_Term12, None, 'x', label='Allgemeiner Term')
plot(lamda1_Hg, ableitungsQuozient(lamda1_Hg), label='AbleitungsQuozient')
#optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
#plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Wert der Quozienten', fontsize=20)
legend(fontsize=15)
plt.tight_layout()
savefig('Term1X.png')
show()


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

delta1 = 360 + W2_Hg[:, 2] - Nullpunkt
delta1 = delta1*2*pi/360
n1_Hg = ones(len(delta1))
for i in range(len(delta1)):
    n1_Hg[i] = (sqrt(2*cos(epsilon)*sin(epsilon + delta1[i, 0] - theta)*sin(theta) + sin(theta)**2 + sin(epsilon + delta1[i, 0] - theta)**2)/sin(epsilon))

lamda1_Hg = toArray(W2_Hg[:, 1])

errorbar(lamda1_Hg, n1_Hg, sigma_n_gelb*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, n1_Hg,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Brechungsindex', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('n2_Hg.png')
show()

al1 = optimizedParameters[0]
bl1 = optimizedParameters[1]
opn2 = optimizedParameters

delta1 = toArray(360 + W2_Hg[:, 2] - Nullpunkt)
errorbar(lamda1_Hg, delta1, 1*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Ablenkwinkel in °', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('delta2.png')
show()

ad1 = optimizedParameters[0]*2*pi/360
bd1 = optimizedParameters[1]*2*pi/360

def ableitungsQuozient(x):
    return (2*ad1*x + bd1)/(2*al1*x + bl1)
    #return bd1/bl1

Term11 = ones(len(n1_Hg))
sigma_Term11 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term11[i] = 2*sin(epsilon/2)/sqrt(1 - n1_Hg[i]**2*sin(epsilon/2)**2)
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term11[i] = gauss("2*sin(epsilon/2)/sqrt(1 - temp1**2*sin(epsilon/2)**2)")


Term12 = ones(len(n1_Hg))
sigma_Term12 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term12[i] = sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term12[i] = gauss("sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))".replace("n1_Hg[i]", "temp1"))


errorbar(lamda1_Hg, Term11, sigma_Term11, None,'x', label='Genäherter Term')
errorbar(lamda1_Hg, Term12, sigma_Term12, None, 'x', label='Allgemeiner Term')
plot(lamda1_Hg, ableitungsQuozient(lamda1_Hg), label='AbleitungsQuozient')
#optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
#plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Wert der Quozienten', fontsize=20)
legend(fontsize=15)
plt.tight_layout()
savefig('Term2X.png')
show()


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

delta1 = 360 + W3_Hg[:, 2] - Nullpunkt
delta1 = delta1*2*pi/360
n1_Hg = ones(len(delta1))
for i in range(len(delta1)):
    n1_Hg[i] = (sqrt(2*cos(epsilon)*sin(epsilon + delta1[i, 0] - theta)*sin(theta) + sin(theta)**2 + sin(epsilon + delta1[i, 0] - theta)**2)/sin(epsilon))

lamda1_Hg = toArray(W3_Hg[:, 1])

errorbar(lamda1_Hg, n1_Hg, sigma_n_gelb*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, n1_Hg,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Brechungsindex', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('n3_Hg.png')
show()

al1 = optimizedParameters[0]
bl1 = optimizedParameters[1]
opn3 = optimizedParameters

delta1 = toArray(360 + W3_Hg[:, 2] - Nullpunkt)
errorbar(lamda1_Hg, delta1, 1*ones(len(lamda1_Hg)), None,'x', label='Messwerte')
optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Ablenkwinkel in °', fontsize=20)
legend(fontsize=20)
plt.tight_layout()
savefig('delta3.png')
show()

ad1 = optimizedParameters[0]*2*pi/360
bd1 = optimizedParameters[1]*2*pi/360
opDelta3 = optimizedParameters


def ableitungsQuozient(x):
    return (2*ad1*x + bd1)/(2*al1*x + bl1)
    #return bd1/bl1

Term11 = ones(len(n1_Hg))
sigma_Term11 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term11[i] = 2*sin(epsilon/2)/sqrt(1 - n1_Hg[i]**2*sin(epsilon/2)**2)
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term11[i] = gauss("2*sin(epsilon/2)/sqrt(1 - temp1**2*sin(epsilon/2)**2)")


Term12 = ones(len(n1_Hg))
sigma_Term12 = ones(len(n1_Hg))
for i in range(len(n1_Hg)):
    Term12[i] = sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))
    temp1 = n1_Hg[i]
    sigma_temp1 = sigma_n_gelb
    sigma_Term12[i] = gauss("sin(epsilon)*n1_Hg[i]/sqrt((n1_Hg[i]**2 - sin(theta)**2)*(1 - (sin(epsilon)*sqrt(n1_Hg[i]**2 - sin(theta)**2) - cos(epsilon)*sin(theta))**2))".replace("n1_Hg[i]", "temp1"))


errorbar(lamda1_Hg, Term11, sigma_Term11, None,'x', label='Genäherter Term')
errorbar(lamda1_Hg, Term12, sigma_Term12, None, 'x', label='Allgemeiner Term')
plot(lamda1_Hg, ableitungsQuozient(lamda1_Hg), label='AbleitungsQuozient')
#optimizedParameters, s = opt.curve_fit(parabel,lamda1_Hg, delta1,  maxfev=10000)
#plt.plot(lamda1_Hg, parabel(lamda1_Hg, *optimizedParameters), label="fit")
xlabel('Wellenlänge in nm', fontsize=20)
ylabel('Wert der Quozienten', fontsize=20)
legend(fontsize=15)
plt.tight_layout()
savefig('Term3X.png')
show()



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

delta_He = 360 + W_He - Nullpunkt
delta_He = toArray(delta_He)

def disp3(lamda):
    return lamda**2*opDelta3[0] + lamda*opDelta3[1] + opDelta3[2]

lamda_He = ones(len(delta_He))
sigma_lamda_He = ones(len(delta_He))
for i in range(len(delta_He)):
    def disp3(lamda):
        return lamda ** 2 * opDelta3[0] + lamda * opDelta3[1] + opDelta3[2] - delta_He[i]
    lamda_He[i] = opt.fsolve(disp3, 500)
    def disp3UP(lamda):
        return lamda ** 2 * opDelta3[0] + lamda * opDelta3[1] + opDelta3[2] - delta_He[i] + 0.1
    def disp3DOWN(lamda):
        return lamda ** 2 * opDelta3[0] + lamda * opDelta3[1] + opDelta3[2] - delta_He[i] - 0.1
    up = opt.fsolve(disp3UP, 500)
    down = opt.fsolve(disp3DOWN, 500)
    sigma_lamda_He[i] = max([abs(lamda_He[i] - up), abs(lamda_He[i] - down)])


#Auflösungen
b = 0.05
sigma_b = 0.03
A1 = abs(b*(2*opn1[0]*578 + opn1[1]))
sigma_A1 = -sigma_b*(2*opn1[0]*578 + opn1[1])

A2 = -b*(2*opn2[0]*578 + opn2[1])
sigma_A2 = -sigma_b*(2*opn2[0]*578 + opn2[1])

A3 = -b*(2*opn3[0]*578 + opn3[1])
sigma_A3 = -sigma_b*(2*opn3[0]*578 + opn3[1])


