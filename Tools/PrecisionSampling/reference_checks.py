#!/usr/bin/env python3
"""Independent small-system references. SciPy is a qualification dependency only.

No protein chemistry or production trajectories are implied by these fixtures.
"""
import json
import math
from scipy.integrate import quad


def eckart_probability(e: float, v1: float, v2: float, hnu: float) -> float:
    """Direct hyperbolic-cosine form, independent of Swift's log-domain evaluator."""
    if e <= max(0.0, v1-v2):
        return 0.0
    alpha1, alpha2 = 2*math.pi*v1/hnu, 2*math.pi*v2/hnu
    denominator = 1/math.sqrt(alpha1)+1/math.sqrt(alpha2)
    a = 2*math.sqrt(alpha1*e/v1)/denominator
    b = 2*math.sqrt(alpha1*(e-v1+v2)/v1)/denominator
    d2 = alpha1*alpha2-math.pi**2/4
    d = math.cosh(2*math.sqrt(d2)) if d2 >= 0 else math.cos(2*math.sqrt(-d2))
    return (math.cosh(a+b)-math.cosh(a-b))/(math.cosh(a+b)+d)


def main() -> None:
    hnu = 6.62607015e-34*299_792_458*100*1000*6.02214076e23/1000
    rt = 8.31446261815324*300/1000
    # Integrate E directly with adaptive Gauss-Kronrod, not the t=sqrt(E/RT)
    # composite Simpson quadrature in the native implementation.
    integral, error = quad(lambda e: eckart_probability(e,40,60,hnu)*math.exp((40-e)/rt)/rt,
                           0.0, 400.0, epsabs=1e-11, epsrel=1e-11,
                           points=[20,40,60], limit=500)
    print(json.dumps({"temperature_K":300,"forward_barrier_kJ_mol":40,
                      "reverse_barrier_kJ_mol":60,"wavenumber_cm-1":1000,
                      "Eckart_factor":integral,"log_factor":math.log(integral),
                      "quadrature_absolute_error":error,
                      "upper_tail_bound":math.exp((40-400)/rt)}, indent=2))

if __name__ == '__main__':
    main()
