# Reproduction material for the two-contact filtered-PID paper (IJC 266601489)

Assembled 2026-07-16. This package contains every script that exists in the current
working environment, recovered verbatim from the session transcripts and re-verified
against the paper's anchor numbers before packaging. It does NOT contain the 15 July
symbolic/certification package, which was delivered in earlier sessions and lives on
your machine; those files are named below so you can locate them.

## Coverage of the five requests

| # | Request | Status here | Where the rest is |
|---|---------|-------------|-------------------|
| 1 | Certified full-loop search (Kp,Ki,Kd,N) | NOT in this environment | `order2_certify_v2_20260715_1158.ipynb` + `order2_certified_v2_fine.json` (your machine) |
| 2 | Certified optima behind Table 2 | Values reproduced inline below (from the paper); N per row only in the JSON | `order2_certified_v2_fine.json` |
| 3 | Boundary tracker (Table 1, Fig 2, six-to-seven count) | Original tracker NOT here; an independent tracer that reproduces the curve and locates the corner IS here (verified) | `flat_boundary_tracker_20260715_1409.py` (your machine) |
| 4 | Competitor gains and simulator (MRDP, AMIGO, CHR) | COMPLETE, all scripts here and verified | - |
| 5 | Method-of-steps recursion and C1-C3 / corner solver | Symbolic recursion and C1-C3 solver NOT here; numerical contact-residual solvers that regenerate (8) and (11) ARE here (verified) | `steps_method_of_steps_20260715_1530.py`, `solve_contact_system_20260715_1530.py`, `spectrum/decompose_20260715_1450.py`, `verify_two_contact_rule_20260715_1712.m` (your machine) |

## Files in this package (all re-run 2026-07-16; anchor results noted)

### Item 4: competitors (Tables 2 and 3)
- `huba_qrdp_compare.py` : Huba-Vrancic 2018 QRDP tuning formulas (their eqs 17-18, PDF-garbled
  terms corrected). VERIFIED: a=0 anchor gives IAEs=2.1547, IAEi=4.7626 exactly. Its `simulate()`
  prefilter is the FIRST reading and is wrong (prints IAEs=39.998 self-check); superseded by
  fix_prefilter.py. Use this file for the tuning formulas only.
- `fix_prefilter.py` : corrected 2DOF prefilter cancelling the DOUBLE dominant pole
  so = -(6+A-S)/2: c = 1/(Ti*TD*so^2), b = -2/(Ti*so). Euler delay-buffer + Pade cross-check.
- `final_numbers.py` : produces the Table 3 MRDP row. VERIFIED: Ts/L = 3.747 (T/L=0.7), 3.990 (1.0),
  4.460 (2.0), 4.838 (5.0), OS=0, IAEs 1.68-2.05. NOTE: its "TWO-CONTACT RULE" block uses the
  flawed Euler realization and prints ~4.2 L; it is NOT authoritative for the rule. Use verify_rule.py.
- `rule_iae.py`, `xcheck.py`, `conv.py` : IAEs of the rule (1.638), cross-checks, step-size convergence.
- `verify_rule.py` : rule (10) on the full un-collapsed loop, Pade delay. VERIFIED: Ts/L = 3.0957
  (Pade-14) to 3.0967 (Pade-18) converging to 3.0985, plant-invariant at T/L = 0.5, 1, 5.
- `chr_amigo.py` : CHR 0% series-form PID as reduced loop 0.6(0.5Ls+1)e^{-Ls}/(Ls). VERIFIED Ts=5.00 L, OS=0.
- `amigo_fix.py` : AMIGO PID with published gains, beta = 1 if L>T else 0, gamma = 0, Tf = Td/10;
  plus CHR series vs parallel. VERIFIED: AMIGO Ts/L 3.574 (0.5), 8.360 (1.0), 10.763 (2.0), 14.355 (5.0),
  OS 1.87 / 3.50 / 5.19 / 4.29 %. CHR parallel OS 8.03 / 5.54 / 3.12 / 1.07 %.

### Items 3 and 5 (numerical): two-contact curve, band dependence, corner
- `trace2.py` : partial-fraction realization of the reduced loop (no step-derivative delta error),
  contact residuals as windowed slope minima, fsolve continuation along the C1-C2 curve in N,
  per-band settling. VERIFIED anchor: Ts(2%)=3.0983 (paper 3.0985), contacts 3.7654 / 4.2767
  grazing at +2e-5 / -7e-6, OS=0. Its w2 watch window [2.05,2.99] MISSES the t=3+ layer trough;
  that is why a first pass mislocated the corner at 22.57. Kept for the curve; use corner2.py for the corner.
- `corner2.py` : trough window spanning the interval boundary [2.5,3.16]; three-contact solve.
  VERIFIED: k=0.619234, tau=0.353259, N=21.351; trough at t=3.0043; Ts(2%)=3.1068, Ts(5%)=2.6303.
- `digits.py` : step-size convergence for the corner (N=21.339 at dt=1e-4, 21.351 at 5e-5) and the
  1%-band optimum (N=15.6, k=0.60898, tau=0.3532, Ts1=3.3435 at both steps). Contact times at
  corner: t0=3.0043, t1=3.6961, t2=4.2030.

Dependencies: numpy, scipy, python-control (`pip install control`).

## Item 2: certified optima behind Table 2 (K = L = 1)

Values as printed in the submitted paper. The certified N per row is not in Table 2 and must be
read from `order2_certified_v2_fine.json`; the paper reports the certified filter constant as
loose in [13, 21] on the cancelling branch and "twice as fast" as the rule at T/L = 0.2.

| T/L | Kp cert | Ki cert | Kd cert | Ts cert /L |
|-----|---------|---------|---------|------------|
| 0.2 | 0.305 | 0.600 | 0.027 | 3.100 |
| 0.3 | 0.364 | 0.595 | 0.042 | 3.155 |
| 0.5 | 0.518 | 0.607 | 0.109 | 3.095 |
| 0.7 | 0.633 | 0.600 | 0.150 | 3.120 |
| 1.0 | 0.817 | 0.608 | 0.211 | 3.125 |
| 1.5 | 1.123 | 0.608 | 0.318 | 3.115 |
| 2.0 | 1.431 | 0.609 | 0.423 | 3.125 |
| 3.0 | 2.008 | 0.603 | 0.594 | 3.180 |
| 4.0 | 2.662 | 0.610 | 0.858 | 3.125 |
| 5.0 | 3.233 | 0.601 | 1.065 | 3.135 |

Rule (10) for comparison: Ki = 0.6105, Kp = 0.6105(T + 0.3525), Kd = 0.2152 T, N = 16.35, Ts = 3.0985 L.

## Item 3: Table 1 tracker points (dt-corrected, from the paper)

| N | k analytic | Ts/L analytic | tau analytic | k tracker | Ts/L tracker |
|---|-----------|---------------|--------------|-----------|--------------|
| 12.0 | 0.599576 | 3.11082 | 0.35921 | 0.5998 | 3.1106 |
| 15.0 | 0.607567 | 3.09938 | 0.35363 | 0.6075 | 3.1009 |
| 17.6 | 0.612929 | 3.09910 | 0.35203 | 0.6128 | 3.0993 |
| 20.0 | 0.617095 | 3.10313 | 0.35237 | 0.6166 | 3.1035 |

The six-to-seven contact-count step is a tolerance-based observation of the original tracker
(it registers the t=3+ layer trough as a near-zero minimum before it strictly binds). The strict
corner from corner2.py is N = 21.35; the two are consistent (trough depth at N=21 is 5.9e-3).

## Item 4: competitor gains per T/L (as used)

MRDP filtered PID (Huba-Vrancic 2018 QRDP; ideal PID Kc(1+1/(Ti s)+TD s), prefilter
Fp = (c Ti TD s^2 + b Ti s + 1)/(Ti TD s^2 + Ti s + 1) with b, c cancelling the double pole so;
simulated with N = 200 i.e. the ideal Qn -> 0 limit, which is the fastest case for that design):

| T/L | A=L/T | Kc | Ti | TD | Kp | Ki | Kd | b | c | so |
|-----|-------|----|----|----|----|----|----|---|---|----|
| 0.5 | 2.00 | 0.8120 | 0.7500 | 0.1667 | 0.8120 | 1.0827 | 0.1353 | 1.3333 | 2.0000 | -2.0000 |
| 0.7 | 1.43 | 0.7475 | 0.9512 | 0.1855 | 0.7475 | 0.7858 | 0.1386 | 1.1422 | 1.6729 | -1.8407 |
| 1.0 | 1.00 | 0.7211 | 1.2360 | 0.2040 | 0.7211 | 0.5834 | 0.1471 | 0.9534 | 1.3772 | -1.6972 |
| 2.0 | 0.50 | 0.7252 | 1.9259 | 0.2308 | 0.7252 | 0.3765 | 0.1673 | 0.6923 | 1.0000 | -1.5000 |
| 5.0 | 0.20 | 0.7524 | 2.7745 | 0.2495 | 0.7524 | 0.2712 | 0.1877 | 0.5281 | 0.7754 | -1.3651 |

AMIGO filtered PID (Astrom-Hagglund 2006; Kc = (0.2+0.45T/L)/K, Ti = L(0.4L+0.8T)/(L+0.1T),
Td = 0.5LT/(0.3L+T), beta = 1 if L>T else 0, gamma = 0, derivative filter Tf = Td/10):

| T/L | Kc | Ti | Td | Tf | beta |
|-----|----|----|----|----|------|
| 0.5 | 0.425 | 0.762 | 0.312 | 0.0312 | 1 |
| 1.0 | 0.650 | 1.091 | 0.385 | 0.0385 | 0 |
| 2.0 | 1.100 | 1.667 | 0.435 | 0.0435 | 0 |
| 5.0 | 2.450 | 2.933 | 0.472 | 0.0472 | 0 |

CHR 0% setpoint PID: Kc = 0.6T/(KL), Ti = T, Td = 0.5L; series (interacting) form cancels
exactly and gives Ts = 5.00 L, OS = 0; parallel form does not cancel (OS 1.1-8.0%).

## Item 1: what is known about the certified search (for the appendix)

Recorded in the paper and the 15 July handover memo, and reproducible only from the notebook:
- Objective: minimum 2%-band settling time of the unit setpoint step subject to y'(t) >= 0 for all t,
  over (Kp, Ki, Kd, N), plant K e^{-Ls}/(Ts+1) with K = L = 1, T/L in {0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5}.
- Certificates were rebuilt "rate-consistently" in the v2 certification; the memo records Table 5
  of that stage as dt = 0.025 certificates, and the fine grid as `order2_certified_v2_fine.json`.
- Cancellation detection: the optimum factors as Kd(s + 1/T)(s + 1/tau2) with the slow zero within
  0.5% of 1/T for T/L >= 0.7; the wall Ts ~ 3.13 L is the lag-independent settling across that range.
- Parameter ranges, grid steps, simulator step, and the y' tolerance used in the search are NOT
  recorded in any document available here; read them from `order2_certify_v2_20260715_1158.ipynb`.
  Do not quote values for them from memory.
