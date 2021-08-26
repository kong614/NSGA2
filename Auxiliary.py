from importing import *

def area_under_curve(x, curve):
    area = 0
    for i in range(1, len(x)):
        area += (x[i] - x[i - 1])*(curve[i] + curve[i - 1])/2
    return area

def generate_random_parameter(limits):
    random.seed()
    parameter = random.uniform(limits[0], limits[1])
    return parameter

def newt_rhap(function, derivative, ans = [], tolerance = 0.0001, init = 0.01):
    xi = init
    if len(ans) == 0:
        for _ in range(infinite):
            xn = xi - function(xi)/derivative(xi)
            if abs(xn - xi) < tolerance:
                break
            xi = xn
    else:
        for _ in range(infinite):
            xn = xi - function(ans, xi)/derivative(ans, xi)
            if abs(xn - xi) < tolerance:
                break
            xi = xn
    return xn

def bisection(ans, function, tolerance = 0.0001, interval = [0, 1]):
    a, b = interval[0], interval[1]
    val_a, val_b = function(ans, a), function(ans, b)
    sign_a, sign_b = nup.sign(val_a), nup.sign(val_b)
    max_iterations = math.log((b - a)/tolerance, 2)
    for _ in range(max_iterations + 2):
        c = (a + b)/2
        val_c = function(ans, c)
        sign_c = nup.sign(val_c)
        if (sign_a != sign_c):
            b = c
            val_b = val_c
            sign_b = sign_c
        elif (sign_b != sign_c):
            a = c
            val_a = val_c
            sign_a = sign_c
    return (a + b)/2

def return_true(_ = 0):
    return True

def linear_interpolate(x0, x1, y0, y1, x):
    if x1 - x0 == 0:
        return (y0 + y1)/2
    return y0 + ((x - x0)*(y1 - y0)/(x1 - x0))