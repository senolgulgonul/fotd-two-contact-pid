# fotd-two-contact-pid

Reproduction material for

> S. Gulgonul, "Monotonic, minimum-settling-time filtered-PID tuning for FOTD plants:
> a two-contact characterization" (manuscript, 2026).

The paper characterizes the fastest strictly monotone setpoint step response that a
filtered PID controller can produce for the first-order-plus-dead-time (FOTD) plant. With
the plant pole cancelled, the loop is the delay-normalized object

    L(s) = k (tau s + 1) e^{-s} / ( s (1 + s/N) )

and the optimum is pinned by two simultaneous boundary contacts of the output slope plus
one stationarity condition, giving

    k = 0.6105/L,  tau = 0.3525 L,  N = 16.35/L,  Ts = 3.0985 L (2% band).

Un-collapsed, the tuning rule for K e^{-Ls}/(Ts+1) is

    Ki = 0.6105/(K L),  Kp = 0.6105 (T + 0.3525 L)/(K L),  Kd = 0.2152 T/K,  N = 16.35/L.

## Requirements

    pip install -r requirements.txt

(numpy, scipy, python-control). The notebook runs on CPU in roughly 20 to 40 minutes.
The single MATLAB file needs the Control System Toolbox.

## Files and what they produce

### Certified full-loop search (Appendix A, Table 2)

| file | produces |
|---|---|
| `order2_certify_v2_20260715_1158.ipynb` | Zero-tolerance search over (Kp, Ki, Kd, N) at h = 0.025 L with PI-seeded cross-plant sweeps, fine re-verification at h = 0.005 L, cancellation check (slow zero against 1/T), and the settling wall. |
| `order2_certified_v2_fine.json` | Output of the notebook: certified gains, filter constant, fine-step settling, PI ceiling, and TV(u) for T/L in {0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5}. These are the certified columns of Table 2. |

### Interval calculus and the characterization system (Sections 3 to 5, Eq. 8)

| file | produces |
|---|---|
| `steps_method_of_steps_20260715_1530.py` | Method-of-steps recursion: the slope on each unit interval as a polynomial plus one boundary-layer exponential, symbolic through [4, 5). Eqs. (3) and (4). |
| `solve_contact_system_20260715_1530.py` | Solves the C1 to C3 system (Eqs. 5 to 7) for (k, tau, N) and the contact times t1, t2. Eq. (8). |
| `spectrum_20260715_1450.py`, `decompose_20260715_1450.py` | Closed-loop spectrum of the reduced loop and the modal decomposition of the slope into the smooth part and the delay-echo ladder (Section 4). |
| `verify_two_contact_rule_20260715_1712.m` | MATLAB check of the rule on the full un-collapsed loop with `stepinfo`. |

### Two-contact curve, band dependence, corner (Section 6, Table 1, Figure 2, Section 9)

| file | produces |
|---|---|
| `flat_boundary_tracker_20260715_1409.py` | Original boundary tracker: the dt-corrected tracker points of Table 1 and Figure 2, and the six-to-seven contact-count step near N = 21. |
| `trace2_20260716_1740.py` | Partial-fraction realization of the reduced loop, contact residuals as windowed slope minima, fsolve continuation along the C1 to C2 curve in N, per-band settling. Anchor: Ts(2%) = 3.0983, contacts 3.7654 and 4.2767. |
| `corner2_20260716_1740.py` | Three-contact corner with the layer trough included: k = 0.6192, tau = 0.3533, N = 21.35, Ts(2%) = 3.1068, Ts(5%) = 2.6303. Eq. (11). |
| `digits_20260716_1740.py` | Step-size convergence of the corner and of the 1%-band optimum (N = 15.6, Ts(1%) = 3.3435). |

### Rule verification and comparison designs (Section 7, Tables 2 to 4, Appendix B)

| file | produces |
|---|---|
| `verify_rule_20260716_1740.py` | Rule on the full loop with a Pade delay: Ts/L = 3.0985, zero overshoot, plant-invariant at T/L = 0.5, 1, 5. |
| `rule_iae_20260716_1740.py` | Setpoint IAE of the rule, 1.638 L. |
| `chr_amigo_20260716_1740.py` | CHR 0% setpoint PID in series form as the reduced loop 0.6 (0.5 L s + 1) e^{-Ls}/(Ls): Ts = 5.00 L, no overshoot. |
| `amigo_fix_20260716_1740.py` | AMIGO filtered PID with published gains, setpoint weights and Td/10 filter; CHR series against parallel form. |
| `huba_qrdp_compare_20260716_1740.py` | Huba and Vrancic 2018 quadruple-real-dominant-pole tuning formulas; use for the formulas only, its own `simulate()` prefilter is superseded by `fix_prefilter`. |
| `fix_prefilter_20260716_1740.py` | Corrected 2DOF prefilter cancelling the double dominant pole. |
| `final_numbers_20260716_1740.py` | MRDP row of Table 3: Ts/L = 3.75 (T/L = 0.7), 3.99 (1), 4.46 (2), 4.84 (5). Its "two-contact rule" block uses an Euler realization and is not authoritative for the rule; use `verify_rule`. |
| `xcheck_20260716_1740.py`, `conv_20260716_1740.py` | Cross-checks of the MRDP realization and step-size convergence. |

## Notes

- `corner2` and `digits` import the realization from `trace2` by reading that file; run them from this directory.
- All settling times are to the 2% band unless stated, with K = L = 1.
- The robustness margins of Table 4 are evaluated on the exact delay from the loop transfer
  functions given in the paper and need no script beyond a frequency sweep.
