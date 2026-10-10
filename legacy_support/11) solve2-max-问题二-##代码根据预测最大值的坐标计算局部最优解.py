# Historical source extracted from supplied notebook; not used for revised results.

##代码根据预测最大值的坐标计算局部最优解
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
from scipy.optimize import minimize


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


def solve(x):
    #sigma_x = 120
    #sigma_y = 120
    #L = 100
    #B = 20
    global x0
    global y0
    global v0
    y0 = x[1]
    x0 = x[0]
    v0 = x[2]

    # 均值
    ##mu_y = 0
    
    
    result, error = integrate.nquad(f, [bounds_x, bounds_y , bounds_v])
    #print(f"积分结果: {result}, 估计误差: {error}")
    return -result

def solve2(v):
    global v0
    v0 = v
    result, error = integrate.nquad(truncated_normal_f, [bounds_v,])
    return result
#solution = solve(0, 0, 150)


x0 = np.asarray((0, 0, 150))
minimize(fun=solve, x0=x0, method='Powell')






