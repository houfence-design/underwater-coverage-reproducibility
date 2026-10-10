# Historical source extracted from supplied notebook; not used for revised results.

##代码运行完会出现命中概率分布的散点图和个点数值方便代入SPSS求解
from scipy import integrate
import math

from scipy import integrate
from scipy.stats import multivariate_normal
import numpy as np

sigma_x = 120
sigma_y = 120
L = 100
B = 20
x0 = 0
y0 = 0
mu_x = 0
mu_y = 0


def bounds_y():
    return [(y0 - (B / 2 + 20)), (y0 + B / 2 + 20)]


def bounds_x(*args):
    return [(x0 - (L / 2 + 20)), (x0 + L / 2 + 20)]


def bivariate_normal_f(x, y):
    #return multivariate_normal.pdf([x, y], mean=[mu_x, mu_y], cov=covariance_matrix)
    coeff = 1 / (2 * np.pi * sigma_x * sigma_y)
    exponent = -0.5 * (((x - mu_x) / sigma_x) ** 2 + ((y - mu_y) / sigma_y) ** 2)
    return coeff * np.exp(exponent)


def solve(x, y):
    #sigma_x = 120
    #sigma_y = 120
    #L = 100
    #B = 20
    global x0
    global y0
    y0 = y
    x0 = x

    # 均值
    ##mu_y = 0
    
    
    result, error = integrate.nquad(bivariate_normal_f, [bounds_x, bounds_y])
    #print(f"积分结果: {result}, 估计误差: {error}")
    return result


max = 0.00
maxx = 0.0
maxy = 0.0
resultx = []
resulty = []
resultz = []
i = 0
for x0 in range(-500, 501, 10):
    for y0 in range(-500, 501, 10):
        z0 = solve(x0, y0)
        i = i + 1
        #print(i)
        resultx.append(x0)
        resulty.append(y0)
        resultz.append(z0)
        if max <= z0:
            max = z0
            maxx = x0
            maxy = y0

import math

import matplotlib.pyplot as plt
import numpy as np
from sympy.stats.sampling.sample_numpy import numpy


def display(arrayx, arrayy, arrayz):
    if type(arrayx) != list or type(arrayy) != list or type(arrayz) != list:
        print("参数必须是list")
        return 1

    x = np.array(arrayx)
    y = np.array(arrayy)
    z = np.array(arrayz)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plt.contour([x, y,z], cmap=plt.get_cmap('jet'))
    ax.scatter(x, y, z)
    plt.show()


#display(resultx, resulty, resultz)

display(resultx, resulty, resultz)

for x in range(len(resultx)):
    print(resultx[x], resulty[x], resultz[x])

