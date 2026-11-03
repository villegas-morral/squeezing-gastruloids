# Appendix: why the model is set up the way it is

The modelling choices that a reader of the report might question, each with the finding that settled it.
The full investigation is in `REPORT.md` (living report) and in the dated answers in `../claude-answers/`;
only the key points are here.

---

### A1. Rigid plates, not force walls

The original constant-force plates let cells through: at 1 layer, 93 % of the cells sat outside the slab and
the "monolayer" was 1.7–2 cells thick. Springs (the July fix) swell under the in-plane crowding: a 2-layer slab
ends at 2.4–2.6 layers, even with a stiff spring at a third of the time step. Rigid plates (a cell touching a
plate cannot move further into it) keep 1.00 / 2.01 / 4.00 layers at the normal time step, as the constant
thickness requires. When the thickness was right, every plate model gave the same pattern.
Figure: `../claude-answers/img/02-walls-thickness.png`.

### A2. A sharp B core, not a diluted one

B is seeded as the innermost 60 % of the cells. With a diluted edge (25 % of the B cells outside the core), a
2-layer N0 150 aggregate polarised in 20 of 21 runs across 7 parameter sets, where the experiments show a
central cluster. With a sharp core it stays central.
Figure: `../claude-answers/img/03-core-and-friction.png`.

### A3. Plate friction 1.5

Without friction, a small B core in a thin slab slides through C to the rim before the readout (in-plane
polarisation; 2 layers, N0 150). An extra in-plane drag λ_wall = 1.5 on cells touching a plate keeps it central
(3 of 3 runs). It slows collective sliding about 3.5× but single-cell rearrangement only 1.25×. More friction
makes monolayers salt-and-pepper (looseness S: 0.49 without, 0.44 at 1.5, 0.39 at 3).
Figure: `../claude-answers/img/03-friction-monolayers.png`.

### A4. Stiffness sets how loose the islands are, not adhesion

Weakening the attraction loosens the islands, but it also stops unconfined polarisation and makes the BC
doublet round up into a sphere. Raising the repulsion of B–B (f_rep 3 → 4.5) and C–C (5.1 → 6) slows
rearrangements without touching the tensions, so the islands stay sparse (S ≈ 0.4–0.6 instead of the July
0.8–0.9) and polarisation and BC glue survive.
Figure: `../claude-answers/img/03-looseness.png`.

### A5. The B–B range μ_BB = 1.9 (patterns in slabs of 3+ layers)

With μ_BB = 2, every slab of 3 or more layers ended as one compact central cluster, whatever we tried: B 65 %,
B scattered outside the core, C mixed into it (the July seeding), spring plates, friction 5, the July
mechanics, λ × 2–4. What matters is the cohesion of B: the adhesion energy W = f_att·2r·A(μ),
A(μ) = (1+μ)ln μ − 2(μ−1). μ_BB 2 → 1.9 lowers W_BB from 1.02 to 0.79, still above W_CC = 0.65: a large core in a
thick slab then hollows into a ring and breaks into islands, while small cores stay compact. At μ_BB = 1.85
(W_BB ≈ W_CC), B simply dissolves into C. Keeping W_BB = 1.02 with a shorter range keeps the core compact,
so it is the energy and not the range.
Figures: `../claude-answers/img/04-ring-3layers.png`, `../claude-answers/img/04-dissolves.png`.

### A6. BC glue, not a sphere

Two aggregates glue without rounding into one sphere when |σ_B − σ_C| < γ_BC < σ_B + σ_C, with σ_X = W_XX/2 and
γ_BC = (W_BB + W_CC)/2 − W_BC (`MODEL.md`, section 2.3). With the final set, after 24 h, the separation of the
two centres over (Rg1 + Rg2) is 0.54 (BB) and 0.56 (CC): fused into one sphere. For BC it is 0.91: two lobes
glued together. The BC entries were never changed (f_att 1.6, f_rep 3, μ 2).

### A7. When the outcome is read: 20 % B

The outcome is read at the first frame with 20 % B cells, the model's equivalent of 96 h (T+ falls from 60 to
20 % between 72 and 96 h in the experiments). It is also recorded at 30 % and 12 % to show how stable it is:
rings break into islands later, and large aggregates polarise late.

### A8. The time unit: 4 minutes

Oriola et al. 2026 set T0 = R0/V0 = 4 min, with R0 = 5 µm and V0 ≈ 1.25 µm/min the instantaneous cell speed
in gastruloids. So 15 time units = 1 h, and the same 4 min is in Sotiris' notebook (`t_scale = 4/60`).
Then q = 0.004667 per unit = 0.070 per hour (a B cell stays B for 14.3 h on average). With τ_B = 20 h and
τ_C = 50 h, the B fraction falls from 60 to 20 % in 19.4 h, against 24 h in the experiments.

### A9. Three ways to measure polarisation

- **Ours (October):** offset = |r_B − r| / R, with R the 98th-percentile distance of the cells from the
  centre (≈ the radius); polarised if ≥ 0.4.
- **The experiments:** P = |r − r_c| / R, with R half the length of the aggregate along its longest axis
  (PC1); polarised if P ≥ 0.3, a cutoff read from the histogram of control gastruloids.
- **Oriola et al. 2026:** the centre of mass of the largest B cluster (≥ 30 cells, link 1.5 r) against that
  of the aggregate, over the radius of gyration of the largest C cluster.

For round aggregates the first two use nearly the same R, so the difference is mainly the cutoff. Since
November the figures use the experiments' definition (A10).

### A10. The readout since November: as in the experiments

The experiments cluster T+ cells with a DBSCAN-like rule: within 13 µm (= 2.6 r here), a T+ cell with at
least 40 % T+ neighbours is a core cell; core cells within 13 µm are linked; other T+ cells join the cluster of
their nearest core cell; clusters of at least 20 cells are kept; N_eff = 1/Σp_i². We classify a run as
**polarised** if P ≥ 0.3, **pattern** if N_eff ≥ 1.5 (B split into clusters of comparable size), otherwise
**central**. An unbroken ring is therefore central, as in the experiments' medium confinement. The one
change: clusters of 10 cells instead of 20, because our aggregates are 10–100× smaller (with 20, monolayers of
N0 ≤ 300 have no cluster at all).

### A11. Cell diffusion and the Damköhler number

There is no noise in the model: cells move only because neighbours divide and push, and because
differentiation changes the adhesion. In a growing monolayer, with the growth flow removed, cells diffuse
with D ≈ 0.004 r² per time unit ≈ 1.5 µm²/h (B and C alike, N0 600 and 2400). A cell needs about 17 h to move
by one diameter, longer than a B cell lives (14 h). So Da = q ℓ²/D is about 4–8 for ℓ = one cell diameter and
500–3800 for ℓ = the pancake radius. Differentiation is much faster than mixing, so the B pattern is set where
the cells already are.

Over one frame (4 min) a model cell moves about 0.04 r, i.e. about 0.05 µm/min: some 25 times slower than the
V0 ≈ 1.25 µm/min that sets T0 (A8). So the model's cells are much less motile than real ones at this time unit.
Figure: `../claude-answers/img/05-msd.png`.
