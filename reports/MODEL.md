# The model: physics, equations and parameters

*The agent-based model in `code-claude/squeezing-gastruloids-claude.ipynb`, built on the original notebook
(`code-original/`) with [CellBasedModels.jl](https://github.com/dsb-lab/CellBasedModels.jl).
The parameter values actually used are tracked in [REPORT.md](REPORT.md).*

---

## 1. What is simulated

Every cell is a sphere of radius $r = 1$ (the unit of length). There are two cell types:

| type | `cell_state` | biology | colour in plots |
|---|---|---|---|
| **B** | 2 | T/Bra+ | green |
| **C** | 3 | T/Bra− | orange |

Type 1 (A) still exists in the code, so the 3×3 parameter matrices keep their original shape, but no A cell is
ever created.

A run has four phases:

1. **Growth.** One cell grows by division to $N_0$ cells. If the run is confined, it grows between the plates.
2. **Relaxation.** The mechanics switch to the differentiation parameters and the aggregate relaxes.
3. **Seeding.** A core holding 60 % of the cells becomes B; the rest is C.
4. **Differentiation.** Cells keep dividing, B cells turn into C, and B and C sort out mechanically.

---

## 2. Mechanics

### 2.1 Equation of motion (overdamped)

Cells have no inertia: the velocity of cell $i$ follows from the forces on it,

$$
\mathbf v_i \;=\; \frac{\mathbf F_i + \alpha \sum_{j \in \mathcal N_i} \mathbf v_j}{\lambda\, n_i}
\qquad (n_i \ge 2), \qquad\qquad
\mathbf v_i = \frac{\mathbf F_i}{\lambda} \qquad (n_i < 2),
$$

and $\dot{\mathbf x}_i = \mathbf v_i$.

* $\mathbf F_i$ is the sum of the pair forces (§2.2) plus the plate force (§4).
* $\mathcal N_i$ are the neighbours closer than `range` $\times\, 2r$, and $n_i$ is their number.
* $\lambda$ (`lambda`) is a friction **per neighbour**, so a crowded cell is slower.
* $\alpha$ (`alpha_ecm`) couples a cell to its neighbours' motion. With $\alpha = 0.9\,\lambda$, neighbours drag
  each other along strongly. Whole regions of tissue therefore flow together easily, while relative motion (one
  cell sliding past another, which is what sorting needs) stays limited by the friction.

The equations are integrated with Heun's method (`CBMIntegrators.Heun`), with time step `dt`.

### 2.2 Pair force

Two cells $i, j$ at distance $d$ interact if $d < \mu_{ij}\, 2r$. The force on $i$ is along the line joining
them, with magnitude

$$
F(d) \;=\; f\left(\frac{2r}{d}-1\right)\left(\mu_{ij}\,\frac{2r}{d}-1\right),
\qquad
f=\begin{cases} f^{\text{rep}}_{ij} & d < 2r \\ f^{\text{att}}_{ij} & 2r \le d < \mu_{ij}\,2r \end{cases}
$$

* $F > 0$ (**repulsion**) for $d < 2r$, when the cells overlap. $f^{\text{rep}}$ (`f_rep`) sets how hard
  cells are to squash: their stiffness.
* $F < 0$ (**attraction**) for $2r < d < \mu\,2r$. $f^{\text{att}}$ (`f_att`) sets how strongly cells stick.
* $F = 0$ at contact ($d = 2r$) and beyond the range $\mu\,2r$. $\mu$ (`mu`) is therefore the reach of the
  adhesion, in units of a cell diameter.

All three are 3×3 matrices indexed by the two cell types, so every pair of types (BB, CC, BC) can stick and
push differently.

### 2.3 Adhesion energy, tensions, and "glue vs sphere"

Pulling two attached cells apart costs the work done against the attractive branch:

$$
W_{ij} = \int_{2r}^{\mu_{ij}2r} |F|\,\mathrm d d
= f^{\text{att}}_{ij}\; 2r\; A(\mu_{ij}),
\qquad
A(\mu) = (1+\mu)\ln\mu - 2(\mu-1).
$$

$f^{\text{rep}}$ does not appear. For example, $A(2) = 0.0794$ and $A(2.14) = 0.1089$.

Following the Differential Adhesion Hypothesis, the work sets effective tensions, up to a common geometric
factor:

$$
\sigma_i = \tfrac12 W_{ii} \quad\text{(surface tension of a pure } i \text{ tissue)},
\qquad
\gamma_{BC} = \tfrac12\left(W_{BB}+W_{CC}\right) - W_{BC} \quad\text{(B–C interfacial tension)}.
$$

When a B aggregate meets a C aggregate (the fusion assay, §7), these tensions decide the final shape:

| condition | outcome |
|---|---|
| $\gamma_{BC} \ge \sigma_B + \sigma_C$ | they stay apart |
| $\lvert\sigma_B - \sigma_C\rvert \ge \gamma_{BC}$ | one engulfs the other: a single round aggregate |
| otherwise | partial engulfment: **two lobes glued together**, the "glue" we need for BC |

Within the glue regime, a larger $\gamma_{BC}$ relative to $\sigma_B + \sigma_C$ means sharper lobes.

The same tensions drive **sorting** inside an aggregate. B and C unmix because $\gamma_{BC} > 0$. Each tissue
also wants to minimise its own exposed surface, which is what pushes a central B core out to the side
(polarisation).

`tension_table(f_att, mu)` in the notebook prints all of this for a parameter set.

**The numbers for the chosen sets** ($r = 1$):

| | $W_{BB}$ | $W_{CC}$ | $W_{BC}$ | $\sigma_B$ | $\sigma_C$ | $\gamma_{BC}$ | $\gamma_{BC} / (\sigma_B+\sigma_C)$ |
|---|---|---|---|---|---|---|---|
| Oct 7 (`final_oct7`, $\mu_{BB}$ = 2) | 1.02 | 0.65 | 0.25 | 0.51 | 0.33 | 0.58 | 0.70 |
| **Oct 8 (`final`, $\mu_{BB}$ = 1.9)** | **0.79** | 0.65 | 0.25 | **0.39** | 0.33 | **0.47** | 0.65 |

Both sets are in the glue regime:
- $\lvert\sigma_B-\sigma_C\rvert$ = 0.18 / 0.07 is below $\gamma_{BC}$
- $\gamma_{BC}$ is below $\sigma_B+\sigma_C$

**Why the shorter B–B range gives patterns in thick slabs** (stage 7 in REPORT.md):
- **Less cohesive B.** Shortening the B–B attraction range from $2\cdot 2r$ to $1.9\cdot 2r$ cuts the B–B
  adhesion energy by a quarter.
- **Thin slabs.** In a slab of 1–2 layers, the C cells made inside the B core cannot get out, because cells
  cannot pass each other there. They split the core whatever its cohesion.
- **Thick slabs.** From 3 layers on, cells pass over and under each other, so with the Oct 7 cohesion the new
  C cells were squeezed out and the core stayed one blob. With the weaker B cohesion of $\mu_{BB}$ = 1.9, a
  big core hollows out into a ring and then islands.
- **Small cores** still stay compact.
- **Too weak.** At $\mu_{BB}$ = 1.85 and below, $W_{BB}$ falls to about $W_{CC}$ and B dissolves into C:
  salt-and-pepper, even in small aggregates.

### 2.4 What each matrix entry does

Only the **BB and CC entries** are tuned. The **BC entries are frozen** at the FUSION values
($\mu_{BC} = 2$, $f^{\text{att}}_{BC} = 1.6$, $f^{\text{rep}}_{BC} = 3$).

| entry | raising it… |
|---|---|
| $f^{\text{att}}_{BB}$, $f^{\text{att}}_{CC}$ | makes B (C) more cohesive. It raises $\sigma$ and $\gamma_{BC}$, so sorting is faster and the patterns more compact. BC glue gets sharper. |
| $\mu_{BB}$, $\mu_{CC}$ | extends the reach of the adhesion. It raises $W$ through $A(\mu)$ and adds interacting neighbours (smoother, more cohesive tissue). |
| $f^{\text{rep}}_{BB}$, $f^{\text{rep}}_{CC}$ | makes the cells stiffer. It leaves the tensions unchanged but makes rearrangements harder (towards jamming). Sorting slows down, so patterns stay looser. |

---

## 3. Growth: division

Each cell divides when the time passes its `t_div`. It is replaced by two daughters placed $\pm\,$`olap`$\cdot r$
apart along a random direction. Each daughter gets a new `t_div` drawn uniformly in
$[\tau(1-\sigma_{\text{div}}),\ \tau(1+\sigma_{\text{div}})]$, where $\tau$ = `tau_div` depends on the cell
type.

Division only happens while `g_on` is true.

---

## 4. Confinement: the plates

Plates sit at $z = 0$ and $z = L$, with $L = 2.1\,r \times$ `layers`, so a cell centre should stay in
$[r,\ L - r]$. The code has three plate models (`wall_type`), compared in REPORT.md (stage 2):

| `wall_type` | model | force / constraint |
|---|---|---|
| 1 | constant (original) | $F_z = \pm$ `rep_wall` while a cell overlaps a plate |
| 2 | Hookean (July) | $F_z = \pm$ `rep_wall` $\cdot\,(\text{overlap}/r)$ |
| 3 | rigid | a cell touching a plate cannot move further into it; an overshoot of at most one time step is undone in `tau_wall` |

**Optional plate friction** (`lambda_wall`, default 0 = frictionless). A cell within $0.05$ of a plate feels
an extra in-plane drag:

$$
v_{x,y} = \frac{F_{x,y} + \alpha \sum_j v_{j,(x,y)}}{\lambda\, n_i + \lambda_{\text{wall}}}.
$$

Because the neighbour coupling $\alpha = 0.9\lambda$ leaves collective motion with a very small effective drag
($\approx 0.1\,\lambda\, n_i$), even a moderate $\lambda_{\text{wall}}$ slows the *collective* migration of a
confined domain much more than local sorting:

- In a monolayer ($n_i \approx 6$) with $\lambda_{\text{wall}} = 1.5$, sliding the whole B domain becomes
  about $(0.6 + 1.5)/0.6 \approx 3.5$ times slower.
- One cell rearranging relative to its neighbours becomes only $(6 + 1.5)/6 = 1.25$ times slower.
- It acts on every cell at 1–2 layers, half of them at 4 layers, and fewer in thicker slabs. So it holds a
  small squeezed B core in place without stopping a thick slab from polarising.

**Used: $\lambda_{\text{wall}} = 1.5$** (`LAMBDA_WALL`), from the relaxation on; growth is frictionless.
Physically it stands for the cells adhering to, or rubbing on, the confining surfaces. Why this value
(stage 5):
- **Without friction** a small B core in a thin slab often slides to the rim before the readout (2 layers,
  N0 150: 1 of 3 runs).
- **With 1.5** it stays central in 3 of 3 runs.
- **Stronger friction slows sorting** enough to loosen the monolayer patterns towards salt-and-pepper: S at
  1 layer, N0 2400 is 0.49 without friction, 0.44 at 1.5 and 0.39 at 3; at 10 (stage 4) the monolayers were
  nearly salt-and-pepper.

---

## 5. Differentiation

While `d_on` is true, every B cell turns into C at a constant rate $q$ (probability $q\,\mathrm dt$ per step).
The original A → B step (rate $p$) and the feedback from A cells (strength $k$) disappeared together with the A
cells.

**The B fraction over time.** Write $g_B = \ln 2/\tau_B$ and $g_C = \ln 2/\tau_C$ for the growth rates. Then

$$
\frac{\mathrm d f_B}{\mathrm dt} \;\approx\; -f_B\,\Big[\,q - (g_B - g_C)(1-f_B)\Big].
$$

With $q = 0.004667$, $\tau_B = 300$ and $\tau_C = 750$, $f_B$ falls from 0.6 to 0.2 in about 290 time units.

**Time units.** Oriola et al. 2026 set the time unit to $T_0 = R_0/V_0 = 4$ min ($R_0 = 5\ \mu$m, the cell
radius; $V_0 \approx 1.25\ \mu$m/min, the instantaneous cell speed in gastruloids). So **1 h = 15 model time
units**, the same `t_scale = 4/60` as in Sotiris' notebook. With it, $1/p \approx 7.1$ h and
$1/q \approx 14.3$ h ($q = 0.070$ h⁻¹). The 60 % → 20 % drop takes ≈ 19 h in the model (≈ 21 h in the
saved runs), against 24 h (72 → 96 h) in the experiments.

Caveat (November 2026): measured over one 4-min frame, model cells move only about 0.05 µm/min, ~25× less
than $V_0$. Their motion comes from divisions alone (section 7, diffusion).

---

## 6. Initial condition

1. **Growth** from one cell, with `parameters_tfm`:
   - fast division (`tau_div` = 5)
   - heavy daughter overlap (`olap` = 0.75)
   - uniform mechanics
   - between the plates if confined

   All cells are C during this phase.
2. **Relaxation:** switch to `parameters_do` and the chosen mechanics; relax for 100 time units.
3. **Seeding** (`initialize_diff!`):
   - Exactly 60 % of the cells become B: the innermost ones by distance from the aggregate centre.
   - The distance is in-plane when confined, so the core is a cylinder spanning the slab. Unconfined, it is
     the 3D distance, so the core is a sphere.
   - With `width` > 0 the edge is diluted: a cell at normalised distance $\rho$ is B with probability
     $1/(1+e^{(\rho-\rho_0)/w})$, with $\rho_0$ chosen so the total is still 60 %.
   - `f_in` asks for the same dilution by its result: the fraction of the B cells that stay inside the sharp
     core (the 60 % innermost cells). `f_in = 0.75` is the diluted start from the meeting, *"75% of the b cells
     … inside and 25% outside"*. The width that gives it is found by bisection. In practice that is
     74–75 % inside, with the core itself 75 % B and the outer shell about 37 % B, densest near the core.
   - `c_core` (Oct 8) mixes C *into* the core, the other kind of salt-and-pepper. This is how the July runs
     were seeded: 75 % B, 15 % A and 10 % C, shuffled inside a core of radius 0.8 R.
     - The core is enlarged to the round(0.6 N / (1 − c_core)) innermost cells.
     - A random `c_core` fraction of them stays C, so exactly 60 % of all cells are still B.
     - `c_core = 0.1` is the July mix with A counted as B.
   - `fB` sets the B fraction; 0.6 is the default.
4. **A short relaxation** (25 units), then differentiation with growth on.

---

## 7. Readouts

* **Clusters**, segmented the way the experimental images are:
  - A B cell with at least 3 other B cells within 1.3 cell diameters is a *core* cell.
  - Core cells linked within that distance, plus the B cells touching them, form a cluster.
  - Clusters of ≥ 10 cells are **islands**. All other B cells are *isolated* (the white points in the
    experimental segmentations).
* **Looseness** $S$: the mean fraction of B among a B cell's neighbours. Compact islands give $S \approx$
  0.8–0.9; random mixing gives $S \approx f_B$.
* **Offset**: 3D distance of the B centroid from the aggregate centre, in units of the aggregate radius $R$
  (the 98th percentile of the cells' distance from the centre). Polarisation pushes it towards 1.
  `offset_max` is the same for the largest island alone.
* **Spread**: the mean distance of the B cells from the centre, over $R$. A compact central blob gives
  ≈ 0.1–0.25; rings and islands give 0.4–0.6.
* **Classes** (`classify_pattern`), tested in this order:

  | class | criterion |
  |---|---|
  | NONE | no B left, or ≥ 70 % of B isolated |
  | POLARISED | offset ≥ 0.4 |
  | CENTRAL | spread < 0.30 in a slab (< 0.45 unconfined), or a centred cluster that fills the middle (`fill`, below) |
  | RING | the largest island holds ≥ half of B and is centred (offset_max < 0.2) |
  | ISLANDS | at least 2 islands |
  | CENTRAL (again) | at most one island and offset < 0.25: a single loose cluster that has not moved (very small aggregates) |
  | OTHER | anything else, e.g. a cluster drifting off-centre |

  The first CENTRAL test is spread < 0.30 in a slab (< 0.45 unconfined), **or** a largest cluster holding
  ≥ half of B with `fill` ≥ 0.55. `fill` is the B fraction of the cells within 0.7 × the radius that a
  compact central disc (ball, unconfined) of the current B would have.

* **The experiment-like readout** (`exp_clusters`, `exp_polarisation`, `classify_exp`; November 2026, used
  in the figures since then). It is the experimentalists' quantification in model units, with
  $r = R_0 = 5\ \mu$m so that 13 µm = 2.6 r:
  - A B (T+) cell is a *core* cell if at least 40 % of all its neighbours within 2.6 r are B.
  - Core cells within 2.6 r are linked, directly or through other core cells.
  - A non-core B cell joins the cluster of its nearest core cell within 2.6 r.
  - Clusters of at least `nmin` cells are kept. The experiments use 20; the figures use 10, because our
    aggregates are 10–100× smaller.
  - $N_\text{eff} = 1/\sum_i p_i^2$, with $p_i$ the share of the clustered cells in cluster $i$.
  - $P = |\mathbf r - \mathbf r_B| / R$, with $\mathbf r$, $\mathbf r_B$ the centres of mass of all cells
    and of the B cells, and $R$ half the extent of the cells along their first principal axis (PC1).

  The classes, in order:

  | class | criterion |
  |---|---|
  | NONE | no cluster kept |
  | POLARISED | P ≥ 0.3 |
  | ISLANDS | N_eff ≥ 1.5 |
  | CENTRAL | otherwise; an unbroken ring is central, as in the experiments' medium confinement |
* **The phase space** groups RING and ISLANDS as **PATTERN** (with the October readout a ring counts as a
  pattern; it breaks into islands later). The outcome is read at the **readout fB = 0.2**: the first saved
  frame with 20 % B left, ≈ 19–21 h after seeding, the model's 96 h. The class at fB = 0.3 and 0.12 is kept
  too, to show how stable the outcome is.
* **Diffusion** (`msd_2d`, `fit_diffusion`; November 2026). The model has no noise, so cells move only
  because neighbours divide and push and because differentiation changes the adhesion. The in-plane MSD of
  cells is followed by id between divisions, in growing monolayers, with the growth flow removed (an affine
  fit per lag). Fitting MSD = 4Dτ + c at 2–8 h gives D ≈ 0.004 r²/unit ≈ 1.5 µm²/h. The Damköhler number
  is Da = q ℓ²/D: ≈ 4–8 for ℓ = one cell diameter, and 500–3800 for ℓ = the pancake radius at squeezing.
* **Fusion assay**:
  1. Two aggregates of 150 cells each are grown and relaxed separately, each as one pure type.
  2. They are placed touching, with no division and no differentiation, for 24 h.
  3. The metric is COM separation / $(R_{g,1}+R_{g,2})$: about 1.6 at contact, about 0.5 once fully rounded
     into one sphere.

---

## 8. Parameters

| symbol | code | value | meaning | source |
|---|---|---|---|---|
| $r$ | `r` | 1 | cell radius (unit of length) | original |
| $\lambda$ | `lambda` | 1 | friction per neighbour | original |
| $\alpha$ | `alpha_ecm` | $0.9\,\lambda$ | velocity coupling to neighbours | original |
| range | `range` | 1 | neighbour-counting radius, in units of $2r$ | Oriola paper (`parameters_do`) |
| $\tau_B,\ \tau_C$ | `tau_div[2:3]` | 300, 750 | mean division times | Oriola paper |
| $\sigma_{\text{div}}$ | `sigma_div` | 0.5 | spread of division times | original |
| olap | `olap` | 0.3 | daughter separation, in units of $r$ | Oriola paper |
| $q$ | `q` | 0.004667 | B → C rate (fixed) | Oriola paper |
| $\mu_{BC}$, $f^{\text{att}}_{BC}$, $f^{\text{rep}}_{BC}$ | `mu[2,3]`, `f_att[2,3]`, `f_rep[2,3]` | 2, 1.6, 3 | BC interaction (**frozen**) | FUSION set |
| $\mu_{BB}$ | `mu[2,2]` | **1.9** | B–B adhesion reach: sets B cohesion, and so whether a big core breaks up in thick slabs | this work, Oct 8 (FUSION / July / Oct 7: 2) |
| $\mu_{CC}$ | `mu[3,3]` | 2.14 | C–C adhesion reach | FUSION / July (unchanged) |
| $f^{\text{att}}_{BB}$, $f^{\text{att}}_{CC}$ | `f_att[2,2]`, `f_att[3,3]` | 6.4, 3.0 | adhesion strength | FUSION / July (unchanged) |
| $f^{\text{rep}}_{BB}$ | `f_rep[2,2]` | **4.5** | B stiffness: sets how loose the patterns are | this work (July 3, FUSION 6) |
| $f^{\text{rep}}_{CC}$ | `f_rep[3,3]` | **6** | C stiffness: keeps a small squeezed B core in place | FUSION (July 5.1) |
| $\lambda_{\text{wall}}$ | `lambda_wall` | **1.5** (`LAMBDA_WALL`) | friction with the plates | this work |
| $L$ | `L` | $2.1\,r\,\times$ layers | plate gap | original |
| plates | `wall_type` | 3 (rigid), `tau_wall` 0.15 | keep the slab at its nominal thickness | this work (stage 2) |
| `dt` | `dt` | 0.015 | time step | original |

The whole set is `PSETS["final"]` in the notebook (`"july_frB4.5_frC6_muB1.9@lw1.5"` in the stage 5 and 7
scripts). It is the July/FUSION mechanics with only $f^{\text{rep}}_{BB}$, $f^{\text{rep}}_{CC}$ and
$\mu_{BB}$ changed, plus plate friction. The set chosen on Oct 7 (with $\mu_{BB}$ = 2, patterns only in
1–2 layers) is kept as `PSETS["final_oct7"]`. The A row and column are kept for shape and never used.

Growth phase (`parameters_tfm`): `tau_div` = 5, `olap` = 0.75, uniform `mu` = 2, `f_rep` = 2.5, `f_att` = 1,
`dt` = 0.002.
