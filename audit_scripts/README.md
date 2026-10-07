Independent verification scripts used for the LXY-repair audit (Python 3, sympy, numpy).

- e7lie.py        E7 Chevalley basis built from the manuscript's stated sign rule; checks Jacobi, invariant form, grading.
- e7_table.py     EVII trace table T(e), Kazhdan weights, orbit dimensions for all 22 characteristics.
- orbit20.py      Orbit-20 triple, p^f basis, triple centralizer and its center, the component sigma_e, slice invariants.
- cubic.py        E6-invariant cubic norms restricted to the orbit-20 slice; branch stabilizer dimensions.
- orbit15.py      Orbit-15 weight-zero conormal Gram matrix.
- radial.py N     Full radial operator and stabilizer kernel at orbit N=17 or 18, matrices M17/M18 (mod prime 33554393).
- semisimple.py   EVII semisimple-slice table dimensions.
- pleasant_e7.py  Pleasantness of inner involutions of E7 via the Weyl-orbit test (confirms (E7,D6+A1) pleasant).
- cii_sp12.py     The Sp12/(Sp6 x Sp6) p-distinguished, non-swap-stable orbit (Target 2).
- e7a7.py         Appendix A E7/A7 lower-bound table.
- nice_delta.py   Numerical spot checks of Appendix A delta formulas (AIII, CI, BDI).
