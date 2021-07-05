import math
from sympy import *
import collections

triple = collections.namedtuple('triple', 'middle f_middle simpson')
x = Symbol("x")
global counter
counter = 1


def simpson_calc(f, a, fa, b, fb):
    """
        Calculation of function based on Simpson Method.
    :param f: given function
    :param a: starting point
    :param fa: function value at point a
    :param b: end point
    :param fb: function value at point b
    :return: middle point, function value at middle point, simpson value.
    """
    middle = (a + b) / 2
    f_middle = f.subs(x, middle)
    simpson = abs(b - a) / 6 * (fa + 4 * f_middle + fb)
    return triple(middle, f_middle, simpson)


def approx_check(f, a, fa, b, fb, eps, current_value, middle, f_middle, side = ""):
    """
    Recursive function that will be called recursively until the approximation reaches the
    accepted epsilion error.
    :param f: given function
    :param a: starting point
    :param fa: function value at point a
    :param b: end point
    :param fb: function value at point b
    :param eps: accepted error
    :param current_value: current value of the calculation
    :param middle: point between start and end points
    :param f_middle: function value at middle point
    :param side: STRING - indicated witch side of the integral we are calculating - left or right.
    :return: integral value.
    """
    global counter
    left_res = simpson_calc(f, a, fa, middle, f_middle)
    right_res = simpson_calc(f, middle, f_middle, b, fb)
    delta = (left_res.simpson + right_res.simpson - current_value) / 15

    print("Iteration number", counter, side + ": ", left_res.simpson + right_res.simpson + delta / 15)
    if abs(delta) <= eps:
        counter += 1
        return left_res.simpson + right_res.simpson + delta / 15
    else:
        counter += 1
        return approx_check(f, a, fa, middle, f_middle, eps / 2, left_res.simpson, left_res.middle, left_res.f_middle, "left side") +\
            approx_check(f, middle, f_middle, b, fb, eps / 2, right_res.simpson, right_res.middle, right_res.f_middle, "right side")


def adaptive_quadrature(f, a, b, eps):
    """
    The main function starts all the calculation
    :param f: given function
    :param a: starting point in the function
    :param b: end point in the function
    :param eps: accepted error
    :return: integral value
    """
    fa = f.subs(x, a)
    fb = f.subs(x, b)
    res = simpson_calc(f, a, fa, b, fb)
    final_result = approx_check(f, a, fa, b, fb, eps, res.simpson, res.middle, res.f_middle, "(based on simpson method)")
    return final_result


# (a, b) = (0.0, math.pi / 4)
# f = math.e ** (3 * x) * sin(2 * x)
# res = adaptive_quadrature(f, a, b, 10 ** -4)
# print("Adaptive Quadrature - integration with", counter, "iterations {} to {} = {}\n".format(a, b, res))
