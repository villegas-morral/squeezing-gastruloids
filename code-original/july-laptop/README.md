# Your notebook from the laptop (July 2026)

`squeezing-gastruloids.ipynb` was committed on Oct 8 from your laptop, where it had been left uncommitted since
July (`code/code-old-july-anothermachine/squeezing-gastruloids.ipynb` in commit 2772c15). It is a version of
`../squeezing-gastruloids.ipynb` (your main notebook, same functions). Kept as it is.

How it differs from `../squeezing-gastruloids.ipynb`:

- **Three named parameter sets** in the parameters cell, chosen with `params_tuned = params_pol`:
  - `params_pol` (UNCONFINED POLARISATION): `params_do` as it is
  - `params_pat` (CONFINED PATTERNS): `params_do` with `f_att` halved
  - `params_fus` (FUSION): Sotiris' `mu`, `f_att` and `f_rep` (`f_rep` = 3·[3 1 1; 1 2 1; 1 1 2])

  where `params_do` has `mu` = 2 everywhere, `f_rep` = 3 everywhere except CC 5.1, and
  `f_att` = [6 6 2.4; 6 6 1.6; 2.4 1.6 4.8].
- **The initial condition** is seeded with `initialize_diff!(com; g_on = false, cluster = true, center = 0.8,
  a0 = 0.15, b0 = 0.75)`. Inside a core of radius 0.8 × the aggregate radius, the cells are shuffled:
  - 75 % become B
  - 15 % become A
  - the remaining 10 % stay C

  So the core is salt-and-pepper, with C cells mixed in. Growth is then switched back on
  (`g_on = true, cluster = false`).
- `n_cells` 600 and `layers = 1`, confined (instead of 500, unconfined, `layers = 2`).

The seeding with C mixed into the core is what `initialize_diff!(...; c_core)` now reproduces in the B/C model
(see `reports/REPORT.md`, Oct 8).
