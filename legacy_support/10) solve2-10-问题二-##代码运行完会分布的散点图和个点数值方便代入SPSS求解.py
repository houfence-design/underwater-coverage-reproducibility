# Historical source extracted from supplied notebook; not used for revised results.

##代码运行完会出现命中概率分布的散点图和个点数值方便代入SPSS求解
from pycurl import INFILE
from scipy import integrate
import math

from scipy import integrate
from scipy.stats import multivariate_normal, sigmaclip
import numpy as np
from statsmodels.sandbox.distributions.gof_new import v_st70_upp
from statsmodels.sandbox.stats.multicomp import maxzero
from sympy.abc import sigma
from sympy.solvers.diophantine.diophantine import transformation_to_normal


sigma_x = 120
sigma_y = 120
sigma_z = 40
L = 100
B = 20
x0 = 0
y0 = 0
v0 = 0
mu_x = 0
mu_y = 0
h0 = 150
l = 120
H = 25



def bounds_y(*args):
    return [(y0 - (B / 2 + 20)), (y0 + B / 2 + 20)]


def bounds_x(*args):
    return [(x0 - (L / 2 + 20)), (x0 + L / 2 + 20)]

def bounds_v(*args):
    return [(v0 - (H / 2 + 20)), (v0 + H/ 2 + 20)]

def bounds_tv(*args):
    return [float('inf'),float('inf')]

def normal_fv(v):
    return 1/np.sqrt(2*np.pi)*np.exp((-1/2)*((v-h0)/sigma_z)*((v-h0)/sigma_z))

def truncated_normal_f(v):
    return normal_fv(v)/((1-phi_integrate(v)) * sigma_z)
def bivariate_normal_f(x, y):
    #return multivariate_normal.pdf([x, y], mean=[mu_x, mu_y], cov=covariance_matrix)
    coeff = 1 / (2 * np.pi * sigma_x * sigma_y)
    exponent = -0.5 * (((x - mu_x) / sigma_x) ** 2 + ((y - mu_y) / sigma_y) ** 2)
    return coeff * np.exp(exponent)
def normal_f(v):
    #return multivariate_normal.pdf([x, y], mean=[mu_x, mu_y], cov=covariance_matrix)
    coeff = 1 / (np.sqrt(2 * np.pi))
    exponent = -0.5 * (v ** 2)
    return coeff * np.exp(exponent)
def f(x, y, v):
    return bivariate_normal_f(x, y) * truncated_normal_f(v)


def phi_integrate(v):
    #integrate.nquad(normal_f, [[-float('inf'), (l-h0)/sigma_z]])
    result, _ = integrate.nquad(normal_f, [[-float('inf'), (l-h0)/sigma_z]])
    return result


def solve(x, y, v):
    #sigma_x = 120
    #sigma_y = 120
    #L = 100
    #B = 20
    global x0
    global y0
    global v0
    y0 = y
    x0 = x
    v0 = v

    # 均值
    ##mu_y = 0
    
    
    result, error = integrate.nquad(f, [bounds_x, bounds_y , bounds_v])
    #print(f"积分结果: {result}, 估计误差: {error}")
    return result

def solve2(v):
    global v0
    v0 = v
    result, error = integrate.nquad(truncated_normal_f, [bounds_v,])
    return result
#solution = solve(0, 0, 150)



max = 0.00
maxx = 0.0
maxy = 0.0
maxv = 0.0
resultx = []
resulty = []
resultz = []
resultv = []
i = 0
for x0 in range(-500, 501, 100):
    for y0 in range(-500, 501, 100):
        for v0 in range(120, 300, 30):
            z0 = solve(x0, y0, v0)
            #i = i + 1
        #print(i)
            resultx.append(x0)
            resulty.append(y0)
            resultz.append(z0)
            resultv.append(v0)
        if max <= z0:
            max = z0
            maxx = x0
            maxy = y0
            maxv = v0

import math

import matplotlib.pyplot as plt
import numpy as np
from sympy.stats.sampling.sample_numpy import numpy


def display(arrayx, arrayy, arrayz, arrayv):
    if type(arrayx) != list or type(arrayy) != list or type(arrayz) != list or type(arrayv) != list:
        print("参数必须是list")
        return 1

    x = np.array(arrayx)
    y = np.array(arrayy)
    z = np.array(arrayz)
    v = np.array(arrayv)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, v, c=z)
    #plt.contour([x, y,z], cmap=plt.get_cmap('jet'))
    #ax.scatter(x, y, v, z)
    plt.show()


#display(resultx, resulty, resultz)

display(resultx, resulty, resultz, resultv)

for x in range(len(resultx)):
    print(resultx[x], resulty[x], resultv[x], resultz[x])


