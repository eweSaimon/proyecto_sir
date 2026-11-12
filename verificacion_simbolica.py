#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificacion_simbolica.py
Comprobación SIMBÓLICA (SymPy) del análisis matemático del proyecto:
conservación, R0, jacobiano, valores propios, pico, integral primera y umbral u*.
Complementa la verificación NUMÉRICA de programa_malware_sir.py.
"""
import sympy as sp

S, I, R, N, beta, gamma, u = sp.symbols("S I R N beta gamma u", positive=True)
b = (1 - u) * beta                                  # beta efectiva

# Sistema con intervención
f = sp.Matrix([-b * S * I / N,
               b * S * I / N - gamma * I,
               gamma * I])
print("1) Conservación: d(S+I+R)/dt =", sp.simplify(sum(f)))

# Jacobiano
J = f.jacobian([S, I, R])
print("\n2) Jacobiano J(S,I,R):")
sp.pprint(J)
print("   Vector propio izquierdo (1,1,1):  (1,1,1)·J =", sp.simplify(sp.Matrix([[1, 1, 1]]) * J))

# Jacobiano en E0 = (N, 0, 0) y valores propios
J0 = sp.simplify(J.subs({S: N, I: 0, R: 0}))
print("\n3) J(E0), E0 = (N,0,0):")
sp.pprint(J0)
print("   Valores propios:", sp.roots(sp.Poly(J0.charpoly().as_expr(), sp.Symbol("lambda")), multiple=False))
lam = sp.Symbol("lambda")
print("   Polinomio característico:", sp.factor(J0.charpoly(lam).as_expr()))
print("   Espacio nulo de J(E0):", [list(v) for v in J0.nullspace()])
print("   -> dimensión del núcleo =", len(J0.nullspace()),
      "(multiplicidad algebraica de 0 = 2): una dirección tangente a la familia de equilibrios "
      "y otra que sale del plano S+I+R=N")

# Jacobiano en un equilibrio genérico (S*, 0, R*) con S*+R*=N
Ss = sp.Symbol("S_e", positive=True)
Jg = sp.simplify(J.subs({S: Ss, I: 0}))
print("\n4) J en equilibrio genérico (S_e, 0, R_e):")
sp.pprint(Jg)
print("   Polinomio característico:", sp.factor(Jg.charpoly(lam).as_expr()))
print("   Valores propios: 0, 0, (1-u)*beta*S_e/N - gamma  (estabilidad transversal si S_e/N < gamma/((1-u)beta))")

# Reducción por conservación: R = N - S - I  => sistema 2D (S, I)
f2 = sp.Matrix([-b * S * I / N, b * S * I / N - gamma * I])
J2 = f2.jacobian([S, I])
J2_E0 = J2.subs({S: N, I: 0})
print("\n5) Sistema reducido 2D (S,I) en E0:")
sp.pprint(J2_E0)
print("   Valores propios:", J2_E0.eigenvals())

# Pico: dI/dt = 0  =>  S = gamma*N/beta_ef
print("\n6) Pico: dI/dt=0 con I>0 =>", sp.solve(sp.Eq(b * S / N - gamma, 0), S))

# dI/dS e integral primera
dIdS = sp.simplify(f[1] / f[0])
print("\n7) dI/dS =", dIdS, " = ", sp.simplify(-1 + gamma * N / (b * S)))
Sv = sp.Symbol("S", positive=True)
H = sp.Function("H")
integral = sp.integrate(-1 + gamma * N / (b * Sv), Sv)
print("   Integral:  I(S) =", integral, "+ C")
# Verificar que I + S - (gamma N / b) ln S es constante a lo largo de las soluciones
Hexpr = I + S - gamma * N / b * sp.log(S)
dH = sp.simplify(sp.diff(Hexpr, S) * f[0] + sp.diff(Hexpr, I) * f[1])
print("   d/dt [ I + S - (gamma N/beta_ef) ln S ] =", dH)

# R0 y umbral u*
R0u = (1 - u) * beta / gamma
print("\n8) R0(u) =", R0u, "| R0(u)=1 ->  u* =", sp.solve(sp.Eq(R0u, 1), u))
vals = {beta: sp.Rational(3, 10), gamma: sp.Rational(1, 10)}
for uu in (0, sp.Rational(1, 2), sp.Rational(4, 5)):
    print(f"   u = {uu}: R0 = {R0u.subs(vals).subs(u, uu)} = {float(R0u.subs(vals).subs(u, uu))}")
print("   u* con parametros base =", sp.solve(sp.Eq(R0u.subs(vals), 1), u))

# Umbral exacto para decrecimiento inicial con S0 = 990, N = 1000 (refinamiento)
S0v = 990
u_exacto = sp.solve(sp.Eq(b.subs(beta, sp.Rational(3, 10)) * S0v / 1000 - sp.Rational(1, 10), 0), u)[0]
print(f"   Umbral exacto de decrecimiento inicial con S(0)=990: u = {u_exacto} = {float(u_exacto):.6f}")

# Pico teórico con parámetros base
Sst = sp.Rational(1, 10) * 1000 / (sp.Rational(3, 10))
Ipico = 10 + 990 - Sst + Sst * sp.log(Sst / 990)
print(f"\n9) Pico teórico (u=0): S* = {Sst} = {float(Sst):.4f}; I_max = {float(Ipico):.6f}")
Sst2 = sp.Rational(1, 10) * 1000 / (sp.Rational(15, 100))
Ipico2 = 10 + 990 - Sst2 + Sst2 * sp.log(Sst2 / 990)
print(f"   Pico teórico (u=0.5): S* = {float(Sst2):.4f}; I_max = {float(Ipico2):.6f}")
print("\nVerificación simbólica terminada.")
