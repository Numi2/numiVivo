# Discrete block error bounds

For a C4 function f on the integer interval [a,b], let H match f and f' at both
endpoints. Hermite interpolation gives

    f(k)-H(k) = f''''(xi_k) (k-a)^2 (k-b)^2 / 24.

Summing the nonnegative factors at integer k, with n=b-a, gives
sum j^2(n-j)^2=(n^5-n)/30. Thus the absolute discrete-sum error is bounded by
M4*(n^5-n)/720. Summing H analytically gives the formula in PROTOCOL.md.
No continuous integral is substituted for the count distribution. A one-point
block is evaluated directly and a two-point block is the exact endpoint sum.

For NB size r and q=mu/(mu+r), the continuous extension of the mass satisfies

    L'(x) = psi(x+r)-psi(x+1)+log(q), L=log P.

From the [polygamma series](https://dlmf.nist.gov/5.15), for m>=1,
|psi^(m)(z)| <= (m-1)!/z^m + m!/z^(m+1). Applying this to the difference over
an interval of length |r-1| yields the derivative bounds in PROTOCOL.md. L' is
monotone because trigamma is decreasing; its endpoint magnitudes bound |L'|.
If both slopes have the same sign, the larger endpoint mass bounds the block.
For r<=1, log mass is convex and the same endpoint bound applies. For a concave
block that crosses the continuous mode, use the discrete modal mass times a
bound for the half-integer distance to a nearest integer. The derivative bound
on that enlarged interval controls the exponential multiplier. Small margins
in code accommodate ordinary rounding; they are not an interval-arithmetic
proof. All blocks begin at counts >=32; the short prefix is summed directly.

If Ln bounds |L^(n)| and p bounds the mass, then derivatives of exp(L) have
bounds p, p L1, p(L1^2+L2), p(L1^3+3 L1 L2+L3), and
p(L1^4+6 L1^2 L2+3 L2^2+4 L1 L3+L4).

D is convex and D' increasing. Its endpoint values bound |D| and |D'|, and,
for n>=2, |D^(n)(x)|=2(n-2)![x^(1-n)-(x+r)^(1-n)] decreases with x.
The product rule therefore bounds the fourth derivatives of P D and P D^2.
These bounds control subdivision, with a shared absolute error budget allocated
in proportion to each disjoint block's number of integer counts.

If approximate retained mass/first/second sums are S0,S1,S2 with absolute
summation errors e0,e1,e2, set m=S1/S0 and t=S2/S0. The mean error is at most
(e1+|m|e0)/(S0-e0), and raw-second-moment error at most
(e2+|t|e0)/(S0-e0), provided S0>e0. Propagate variance subtraction using
Evar <= Esecond+(2|m|+Emean)Emean. Add the geometric omitted-tail bounds to each
raw error before applying the same formula for the combined result. The returned
summation and truncation contributions add to these combined normalized bounds.
All three combined checks must pass before returning a value.

PMF evaluation uses saturated NB log likelihood plus -D/2, expressed with
[Stirling corrections](https://dlmf.nist.gov/5.11), to avoid subtracting large
log-Gamma values. The derivative uses stable power differences via expm1.
Those numerical special-function evaluations remain ordinary FP64; the analytic
interpolation/tail bounds do not certify their floating-point roundoff.
