from sympy import *

# used symbols, (ei, ej) -> multiindex of the micro element, n -> number of cell per direction
x, y, ei, ej, n = symbols('x y ei ej n')
e = Matrix([[ei], [ej]])

# basis functions
phi_1d = [1 - x, x]
phi = [px * py for py in [p.subs(x, y) for p in phi_1d] for px in phi_1d]
phi_grad = [Matrix([[diff(f, x)],[diff(f, y)]]) for f in phi]

# quadrature point
q = [0.25, 0.25]

# coordinates of macro element corners
c = [Symbol('c_{}'.format(i)) for i in range(4)]

# corners of reference element
w = [Matrix([[i], [j]]) for j in range(2) for i in range(2)]

# corner of micro element in reference element
F = [(w_i + e)/n for w_i in w]

# coordinates of micro element corners
v = [sum(ci * phi_i.subs(([x, f[0]], [y, f[1]])) for ci, phi_i in zip(c, phi)) for f in F]

# first column of DX(q)
J_x = sum((vi * grad[0].subs(([x, q[0]], [y, q[1]])) for vi, grad in zip(v, phi_grad)))
J_x = J_x.simplify().nsimplify()  # without nsimplify there are terms like 1.0 * c_i
print(J_x, 'number of ops = {}'.format(count_ops(J_x)))

# computing DX by transforming q

# transformed q in reference element of macro element
T = Matrix([[(q[0] + e[0]) / n], [(q[1] + e[1]) / n]])

# first column of DX(q)
J2_x = sum((ci * grad[0].subs(([x, T[0]], [y, T[1]])) / n for ci, grad in zip(c, phi_grad)))
J2_x = J2_x.simplify().nsimplify()
print(J2_x, 'number of ops = {}'.format(count_ops(J2_x)))

# this is false but by looking at both expressions this should be true
assert J2_x == J_x


