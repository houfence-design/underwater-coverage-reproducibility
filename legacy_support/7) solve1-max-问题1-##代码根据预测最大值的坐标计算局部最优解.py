# Historical source extracted from supplied notebook; not used for revised results.

##代码根据预测最大值的坐标计算局部最优解
from scipy import integrate
import math

from scipy import integrate
from scipy.stats import multivariate_normal
import numpy as np
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
from scipy import integrate
import math

from scipy import integrate
from scipy.stats import multivariate_normal, sigmaclip
import numpy as np
from statsmodels.sandbox.distributions.gof_new import v_st70_upp
from sympy.abc import sigma
from sympy.solvers.diophantine.diophantine import transformation_to_normal
from scipy.optimize import minimize

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


def solve(x):
    #sigma_x = 120
    #sigma_y = 120
    #L = 100
    #B = 20
    global x0
    global y0
    y0 = x[0]
    x0 = x[1]

    # 均值
    ##mu_y = 0
    
    
    result, error = integrate.nquad(bivariate_normal_f, [bounds_x, bounds_y])
    #print(f"积分结果: {result}, 估计误差: {error}")
    return -result
x0 = np.asarray((0, 0))
minimize(fun=solve, x0=x0, method='Powell')





