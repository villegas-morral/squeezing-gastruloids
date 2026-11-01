# Squeezing gastruloids: living report

*Updated as the work goes. Last update: 2026-10-08 early morning: stage 7 (patterns in thick slabs: new set with `mu` BB 1.9), the four phase-space plots and the presentation PDF.*

## Goal (meeting, October 2026)

One parameter set for a model with two cell types, **B** (T/Bra+) and **C** (T/Bra−). B differentiates into C.

* **Initial condition.** A B-only core holding ≈ 60 % of all cells, with C outside. There are no A cells, so
  the A-based differentiation feedback is dropped. A diluted edge, as in the continuum-model slides, is
  welcome if it does not break the simulations; otherwise we use a sharp core.
* **Squeezed (confined) behaviour.** The phase space, N0 = cells at squeezing vs layers = plate gap, should
  show three regions:

  ```
  layers
    12 | P  P  P  P  P
     8 | P  P  P  P  C
     4 | C  C  C  I  I
     2 | C  C  I  I  I
     1 | C  I  I  I  I
       +---------------- N0
        150 300 600 1200 2400
  P polarised   C central cluster (stays, differentiates away)   I sparse islands
  ```
* **The islands must be sparse.** That means many irregular, loose B domains plus isolated B cells, as in
  `experiments/pattern-examples/` and your first simulations (`code-original/plots/pancake_3layers_n2000_bottom.png`).
  The compact round islands of the July set do not qualify.
* **Same parameters for every check:**
  - unconfined polarisation, which only has to happen, not fast or perfectly
  - BC fusion that "glues" instead of forming a sphere
* **Final plots** of every behaviour are 3D. The data behind each one (positions and cell types) is saved,
  and the way each was made is documented in `FIGURES.md`.

## Constraints

These come from the July chat and the October 6 answers. Quotes are your words.

| constraint | source |
|---|---|
| Tune only `f_att`, `f_rep`, `mu`, and *"mainly do not touch the [2,3] and [3,2] elements of those matrices"*. The BC entries stay at the FUSION values: `mu` = 2, `f_att` = 1.6, `f_rep` = 3. | 2026-07-25; re-confirmed 2026-10-06 |
| *"you cannot change p,q"*. Only `q` = 0.004667 remains, since `p` is gone with A. | 2026-07-26 |
| No random active forces. `lambda`, `alpha_ecm` and `tau_div` stay as in your parameters. | 2026-10-06 |
| The thickness must stay constant: *"If layers=2, then the aggregate should be of thickness 2 the whole time."* | 2026-07-26 |
| BC fusion: *"they dont form a sphere, they sort of like glue"*. There is no numerical target; the Sotiris figure is no longer used. | 2026-07-26; 2026-10-06 |
| Polarisation is only a check: *"It does not matter if polarisation is not perfect or fast."* | 2026-07-26 |
| The continuum (Cahn–Hilliard) model is a different model and is not a guide for this one. | 2026-10-06 |
| The readout time is not fixed. Bra+ going from 60 % to 20 % (72 → 96 h) is a reference, and we look at the results to choose. | 2026-10-06 |

## Decisions log

| date | decision | why |
|---|---|---|
| 10-06 | Repo split into your code and mine. You then renamed the folders: `code-original/` (yours, untouched), `code-claude/` (mine; superseded July code in `code-claude/july-2026/`), `experiments/` (pattern examples and presentations). All reports go in `reports/`: this file, plus the July reports in `reports/july-2026/`. | You asked for it, and the July report had become messy. |
| 10-06 | The new notebook starts from the July notebook, which is your notebook plus the July changes. Each July piece is kept, revised or re-tested; nothing is kept or dropped wholesale. | Your code is the backbone, and the July work improved parts of it. |
| 10-06 | Rigid plates are the default (`wall = :rigid` in `run_confined`, `wall_settings` and the Initialization cell). The other models stay available. | The only model that keeps the nominal thickness (stage 2). Confirmed 10-07 ("whatever works best"). |
| 10-06 | The classifier is based on spread (central vs ring vs islands) instead of cluster connectivity. | Loose rings were being called "central" (stage 3). |
| 10-06 | Weakening `f_att` was dropped; stiffness (`f_rep` BB/CC) is the tuning direction. Plate friction is kept as an option, off. | Stage 4: weaker adhesion kills polarisation and BC glue; friction makes monolayers salt-and-pepper. |
| 10-06 | **Pending your choice:** sparseness level (`f_rep` BB 4.5 to 6; FUSION = 6), rings as patterns, readout time, sharp vs diluted core. | See `claude-answers/02-model-plates-screen.md`. |
| 10-07 | Your answers: (1) *"the best one is loose islands july_frB4.5, but something a little bit more loose would also be ok"*; (2) rings count as patterns; (3) the readout at fB = 0.2 is fine if it catches the pattern; (4) plates: whatever works best; (5) core: whatever works best, and if equal a diluted one, *"75% of the b cells … inside and 25% outside aprox"*. | `claude-answers/02-model-plates-screen.md` |
| 10-07 | Phase-space classes: CENTRAL, PATTERN (= RING or ISLANDS), POLARISED, OTHER, read at fB = 0.2. The class at fB = 0.3 and 0.12 is recorded too, to show how stable it is. | Answers (2) and (3). |
| 10-07 | New option `f_in` in `initialize_diff!`, `run_confined` and `run_unconfined`: the fraction of the B cells inside the sharp core. 0.75 is your diluted start; the logistic width that gives it is found by bisection. | Answer (5). |
| 10-07 | Stage 5a, tuning: the candidates between july_frB4.5 and FUSION (`f_rep` BB 4.5, 5, 5.5 or 6, with CC 5.1 or 6), with both cores and 3 seeds, at the points that separate the three regions. | Answer (1). One seed was not enough to tell the candidates apart at the small-aggregate points. |
| 10-07 | **Sharp core, not diluted.** At 2 layers, N0 150, the diluted core (75 % inside) polarised in 20 of 21 runs across all 7 candidates by fB 0.15. With a sharp core, FUSION stayed central to the readout in 1 run of 2 (the other was undecided). The diluted runs of the other candidates were dropped from the queue; FUSION and july_frB5_frC6 keep one diluted seed for the record. | Your rule: dilution only if it works as well. The B cells scattered outside the core seed a B patch at the rim, and the core moves there. |
| 10-07 | **Classifier: a loose central cluster is CENTRAL.** New statistic `fill` = B fraction of the cells within 0.7 × the radius that a compact central disc of the current B would have. CENTRAL now also covers `f_max` ≥ 0.5 with `fill` ≥ 0.55. All runs are re-scored from their trajectories (`scripts/stage5_common.jl`, cached). | At the readout a perfectly compact central disc of 20 % B has spread (2/3)√0.2 ≈ 0.30, exactly the old threshold, so slightly loose central clusters came out as RING or OTHER. Checked by eye on 12 runs: central clusters have fill 0.65–0.85, rings 0.24–0.41. |
| 10-07 | **Classifier: a single loose cluster that has not moved is CENTRAL** (at most one island, offset < 0.25), checked after RING and ISLANDS. Re-scored again (cache version `2026-10-07-fill-single`). | At 1 layer, N0 75, about 30 B cells are left at the readout. They form one loose cluster near the centre, which is neither a pattern nor polarised, and was coming out as OTHER. |
| 10-07 | **Gentle plate friction re-tested** (wave B, `lambda_wall` 1.5 and 3; it was 10 in stage 4). |
| 10-07 | Friction 1.5–3 keeps the small 2-layer aggregate central when C is stiff (FUSION: B offset at the readout 0.07–0.24 over 6 runs), but not with july_frB4.5's softer C (CC 5.1; offsets 0.29, 0.59, 0.58). **Next candidate: july_frB4.5_frC6** (B as in july_frB4.5, whose looseness you liked; C as in FUSION), at friction 0, 3 and 5 (wave C). | Stiff C resists the sideways slide of the B domain, and friction slows it further. Looseness is set mainly by B (S at 2 layers, N0 2400: BB 4.5 → 0.59–0.63, BB 5 → 0.33–0.49, BB 6 → 0.37). | Without friction, a small B domain in a thin slab slides through C to the rim (in-plane polarisation) before it can break up. That happens for every stiffness tried (2 layers, N0 150: polarised by fB 0.15 in most runs). You said plates are "whatever works best". |
| 10-07 | Wave C, 2 layers N0 150, B offset at the readout: july_frB4.5_frC6 without friction 0.25 / 0.60 (one seed polarised); friction 3: 0.17 / 0.29; friction 5: 0.03 / 0.17 / 0.28. **Provisional choice: july_frB4.5_frC6 with plate friction 3, sharp core.** The phase space (`results/stage5/phase_A`) started with it while the pattern points of wave C were still running. | Friction 3 is enough at that point, and the gentler friction is less likely to make the monolayers salt-and-pepper (stage 4: friction 10 did). |
| 10-07 | **Final choice: friction 1.5, not 3.** The first phase-space runs with friction 3 showed a nearly salt-and-pepper big monolayer (1 layer, N0 4800, S 0.33). A quick test (wave D, `tune_d.csv`) gave: at 1 layer N0 2400, S 0.49 without friction, 0.44 at 1.5, 0.39 at 3; at 1 layer N0 1500, friction 1.5 gave islands with S 0.42–0.46; at 2 layers N0 150, friction 1.5 kept B central in 3 of 3 runs (offsets 0.12, 0.26, 0.11). The phase space was restarted with friction 1.5 (`results/stage5/phase_B`). The friction-3 runs already finished (15) stay in `phase_A`, and the unfinished ones were stopped. | Friction 1.5 is as good as 3 at keeping small cores central but loosens the monolayers much less, which matters for "patterns like july_frB4.5". |

| 10-08 | **Your July laptop commit merged and sorted** (`git merge origin/dev`, commit 2772c15):<br>- the LaTeX report (`report-july-tex/`) went into `reports/july-2026/report/`, next to the same images (on Oct 8 these went into its `figures/` subfolder, with the paths in `main.tex` updated); byte-identical duplicates were dropped, and `main.tex`, `refs.bib`, the PDF and the two new images were moved in<br>- your notebook from the laptop went to `code-original/july-laptop/`<br>- your edited copy of my July notebook went to `code-claude/july-2026/squeezing-gastruloids_claude-july_laptop.ipynb`<br>Each folder has a README saying what the files are. | Your request (task 1). Your code and mine stay apart, as in the rest of the repo. |
| 10-08 | **What the July runs did differently**, found by reading those notebooks (details in stage 7 below). The main one is that July seeded the core **salt-and-pepper, with C mixed in**: inside a cylinder of radius 0.8 R, 75 % B, 15 % A and 10 % C, shuffled. New option `c_core` in `initialize_diff!`, `run_confined` and `run_unconfined`: the fraction of the core left C, mixed in at random, keeping B at 60 %. The batch runner also takes `fB0.65` (B fraction) and `@hookean` (the July spring plates). | Your question (task 2): July had patterns at 3 and 4 layers; we have none above 2. |
| 10-08 | **New final set: `mu` BB 2 → 1.9** (`PSETS["final"]`; the Oct 7 set is kept as `final_oct7`). None of your suggestions changed the thick-slab blob: B 65 %, salt-and-pepper starts (B outside the core or C inside it, as in July), spring plates, or stronger friction. Neither did the July mechanics without A, or slower mechanics (λ). Weaker B–B cohesion did: at `mu` BB 1.9 a big core in a slab of 3+ layers opens into a ring and then islands, while small cores stay central, monolayers pattern, the free aggregate polarises and BC glues. The phase space was rerun (`results/stage5/phase_D`). | Your task 2, and your answer "yes, rerun" if a variant patterns at 3+ layers. Stage 7 below. |
| 10-08 | **The phase space is drawn four ways** (`replot_figures`): x = log₁₀ N0 or N0; y = plate gap (layers) or pancake diameter at squeezing (measured on each run, column `D`). | Your request (task 3). |
| 10-08 | **LaTeX report** `reports/october-2026/main.tex`, written like the July one. Its parameter table has the July report's values next to the current ones, with the changes in red. | Your request. |
| 10-08 | **Presentation figure** `figures/presentation.pdf` (`scripts/stage7_presentation.jl`), laid out like your sketch, in its colours, to finish in Inkscape. | Your request (task 4). |

## The July work, piece by piece

| piece | verdict | notes |
|---|---|---|
| `cellEdge` from the current `mu`; `show_progress` flags | keep | plain fixes |
| `save_frames` / `load_frames` | keep | JLD2 cannot reload a `Community` |
| `fA` / `push_fA!` feedback | drop | there are no A cells; differentiation is B → C at rate `q` |
| centroid-centred cylinder seeding | revise | B-only, 60 % of all cells, optional diluted edge |
| Hookean plates + `wall_settings` | re-test | soft at ≥ 2 layers (≈ 2.6 / 3.5 / 4.5 effective layers) and needs dt = 0.005 at 1 layer; compared with alternatives at checkpoint 2 |
| whole-trajectory island scoring | revise | it called loose patterns "dispersed", which are now the target; add 3 classes, offset, looseness, isolated cells |
| fusion assay | simplify | BC doublet, plus BB and CC for context |
| `run_confined` / `run_unconfined` / `phase_scan` | adapt | new initial condition and metrics |
| export, figure and report scripts | revise | the July report generator is retired; this file is written by hand |
| D/h law, old timesteps, two-category phase plot | drop | refuted or superseded |

## Stages

- [x] **1. Reorganise the repository** (checkpoint 1)
- [x] **2. Model and initial condition:** B → C only; B = 60 %, sharp vs diluted; wall comparison (checkpoint 2)
- [x] **3. Readout and classifier**, calibrated by eye against `experiments/pattern-examples/` (checkpoint 3)
- [x] **4. Parameter screen:** B/C entries of `f_att`, `f_rep`, `mu`; BC frozen (checkpoint 4)
- [x] **5. Phase space:** rough pass, then refined (stage 5a tuning + 5b phase space)
- [x] **6. Final 3D figures with saved data**, plus `FIGURES.md` and `MODEL.md`
- [x] **7. Patterns in slabs of 3+ layers** (Oct 8): new set, phase space rerun, four phase-space plots, presentation PDF

## Results

### Stage 2: model, plates, initial condition (2026-10-06)

**Model changes**, made in `code-claude/squeezing-gastruloids-claude.ipynb` with the library re-exported. Every
function was smoke-tested (`code-claude/scripts/smoke_test.jl`).

* Differentiation is B → C at the constant rate `q`. Removed: `p`, `k`, `fA`, `push_fA!`, `ni_a`, `r_ab`,
  `r_bc`.
* The aggregate grows as C. The original grew type 1 (A); growth mechanics are identical for every type, so
  this only affects the relaxation, which now runs under CC mechanics.
* `initialize_diff!` seeds B so that **exactly 60 %** of all cells are B. The core is the innermost cells: a
  cylinder through the slab when confined, a sphere otherwise. An optional `width` dilutes its edge with a
  smooth profile while keeping exactly 60 %.
* `differentiate!` can stop at a given B fraction (`fB_stop`). By default it runs the whole `tf`.
* The model has a `wall_type` switch: 1 = constant force (original), 2 = Hookean (July), 3 = rigid.

**The plate models compared.** July set, sharp core, seed 2345, about 22 h (`scripts/stage2_walls.jl`):

| wall | 1 layer (N0 600) | 2 layers (N0 1500) | 4 layers (N0 1500) | dt | cost (1 / 2 / 4 layers) |
|---|---|---|---|---|---|
| constant force (original) | **1.66** | **2.59** | 4.01 | 0.015 | 66 / 205 / 188 s |
| Hookean (July) | 1.01 | **2.60** | 4.20 | 0.005 at 1 layer | 180 / 186 / 238 s |
| Hookean, stiff (150) | 1.01 | **2.38** | (running) | 0.005 | 177 / 500 / – s |
| **rigid** | **1.00** | **2.01** | **4.00** | 0.015 | 53 / 198 / 181 s |

*Effective layers at the end of the run, computed as (central-98 % z-spread + 2r) / 2.1r.*

The pattern outcome was the same whenever the thickness was right:
- 1 layer: ISLANDS for the Hookean, stiff Hookean and rigid walls
- 2 and 4 layers: one central cluster for every wall

**Decision: rigid plates** (`wall = :rigid`), reversible with one keyword.

- **It is the only model that keeps the thickness at the nominal value**, which was your requirement. Any
  spring settles where the in-plane crowding balances it, so it swells. Even a stiff one at dt 0.005 leaves a
  2-layer slab at 2.4 layers, at 2.5× the cost.
- **It runs at the normal time step**, which matters for the hundreds of runs the phase space needs.

*Flagged for you: you said rigid walls might not be the way. If there is a physical reason to prefer a soft
wall, the numbers above show what it costs.*

**The readout for the B/C model.**
- **Clusters**, segmented like the experimental images:
  - a *core* B cell has ≥ 3 B neighbours within 1.3 diameters
  - core cells linked within that distance, plus the B cells touching them, form a cluster
  - clusters of ≥ 10 cells are islands; every other B cell is *isolated*
- **Looseness S:** the mean fraction of B among a B cell's neighbours.
- **Offset:** distance of a centroid from the aggregate centre.
- **Classes:** POLARISED / CENTRAL / ISLANDS / NONE.

Reference: the July compact islands score **S ≈ 0.80–0.87, with 0–1 % isolated B cells**.

### Stage 3: the readout, revised after looking at the first runs (2026-10-06)

The first classifier called a loose, ring-shaped B network "CENTRAL" because all its cells are connected. By
eye it is clearly a pattern. The revised `classify_pattern` checks, in order:

| class | criterion |
|---|---|
| POLARISED | the B centroid has shifted ≥ 0.4 R (3D shift, so a thick slab can also polarise towards a plate) |
| CENTRAL | B is gathered in the middle: `spread` < 0.30 in a slab, < 0.45 in 3D |
| RING | B spread out but one connected network centred on the middle |
| ISLANDS | B spread out in several islands |
| NONE / OTHER | no B left / B mostly isolated, or anything else to inspect |

`spread` is the mean distance of the B cells from the centre over R. A compact central blob scores ≈ 0.1;
rings and islands score 0.4–0.6.

The thresholds are provisional and need your eye: see the galleries in `claude-answers/`.

### Stage 4: parameter screen, wave 1 (2026-10-06)

Rigid plates, sharp core, seed 2345 (`scripts/stage4_screen.jl`). Only the BB and CC entries vary; the BC
entries stay frozen. Looseness S and classes are taken at fB = 0.2.

| candidate | BB / CC change | 1 layer, N0 200 | 1 layer, N0 1500 | 2 layers, N0 2500 | polarises? | BC at 24 h |
|---|---|---|---|---|---|---|
| July | — | compact cluster(s) | 9 islands, S 0.68 | **compact central blob** (S 0.92) | yes (0.52 → 0.85) | 0.78 |
| july_frB4.5 | `f_rep` BB 3 → 4.5 | loose central cluster | 17 islands, S 0.52 | breaking ring → islands | yes (0.67) | 0.91 |
| july_frB4.5_frC6 | `f_rep` BB 4.5, CC 6 | loose ring | 17 islands, S 0.50 | ring → ~7 islands | yes (0.54) | 0.86 |
| **FUSION** (Sotiris, unchanged) | `f_rep` BB 6, CC 6 | loose ring | 12 islands, S 0.44 | loose ring → scattered islands | yes (0.61) | **0.94** |
| CONFINED PATTERNS (+ FUSION BC) | weaker `f_att`, `mu` | dispersed | very sparse, 26 % isolated | ring | **no** (0.07) | **0.55 (sphere)** |
| july_fa0.7 | `f_att` BB, CC × 0.7 | dispersed | 17 islands, S 0.45 | central | **no** (0.07) | **0.53 (sphere)** |

What the screen shows:

1. **Weakening B or C adhesion is a dead end.** With the BC entries frozen it makes B and C too miscible. The
   free aggregate stops polarising and the BC doublet rounds into a sphere. `july_fa0.5` and `cp_frB4.5`
   (further in that direction) were dropped without running.
2. **Stiffer B and C cells (`f_rep`) are the sparseness knob.** As `f_rep` goes up, patterns get looser
   (S 0.68 → 0.44), and polarisation and BC glue survive. BC glue even improves, 0.78 → 0.94, since stiffer
   aggregates round up more slowly.
3. **At 2 layers the July set keeps one compact central blob, while stiffer sets form a ring that breaks into
   islands.** This happens around fB ≈ 0.1–0.2, so the readout time matters there.
4. **In the small monolayer (N0 200)**, stiffer sets give a loose ring or central network rather than a
   compact cluster. *Open question for you: does that count as "the cluster stays in the centre"?*

Notably, Sotiris' **FUSION set, unchanged**, now does everything qualitatively. The July softening of
`f_rep` was only needed for the old A+B core.

**Wave 2**, around FUSION:

| candidate | change from FUSION | 1 layer, N0 200 | 1 layer, N0 1500 | 2 layers, N0 2500 | polarises? | BC |
|---|---|---|---|---|---|---|
| fusion@w0.1 | diluted core (width 0.1) | loose central | 15 islands, S 0.40 | ring, S 0.33 | yes (0.61) | (as FUSION) |
| fusion_frC5.1 | `f_rep` CC 6 → 5.1 | 3 islands | 15 islands, S 0.46 | ring, S 0.40 | yes (0.53) | 0.98 |
| fusion_frB9 | `f_rep` BB 6 → 9 | 3 islands | 15 islands, S 0.41 | ring, S 0.31 | **weak** (0.30) | 1.00 |
| fusion_muB2.2 | `mu` BB 2 → 2.2 | ring | 12 islands, S 0.53 | (running) | | |

Stiffness has an upper limit: with BB at 9, polarisation is too slow. The diluted core changes little.

### A problem the coarse phase maps revealed, and a possible fix: friction with the plates

Coarse phase maps (`scripts/stage4_minimap.jl`; layers 1, 2, 4, 8 × N0 150, 600, 2400) show that **small
aggregates at 2–4 layers polarise**: B slides to the rim. You expect them to stay as a central cluster.
In the model the plates are frictionless, so a confined B domain glides in-plane as easily as a free one.

**Option under test: plate friction.**
- Cells touching a plate get an extra in-plane drag `lambda_wall` (new model parameter, default 0, so
  every result above is unchanged).
- Why it targets the right thing: with `alpha_ecm = 0.9 λ`, collective motion has a very small effective
  drag, about $0.1\,\lambda\, n_i$. A plate drag of the same order slows *collective* migration
  (polarisation) roughly tenfold, but local sorting only about twofold.
- It acts on every cell at 1–2 layers, half of them at 4 layers, and a quarter at 8. That would give "few
  layers: the cluster stays; many layers: behaves as unconfined".

*This changes the confinement mechanism, so it is your call. Results go to the answer file.*

**Coarse maps at fB = 0.2** (one seed; P polarised, C central, R ring, I islands; `→` = what it becomes by
fB 0.12):

| layers | FUSION 150 | 600 | 2400 | july_frB4.5 150 | 600 | 2400 | july_frB4.5 + friction 10 150 | 600 | 2400 |
|---|---|---|---|---|---|---|---|---|---|
| 8 | P | C→P | C | P | P | C→P | P | P | C→P |
| 4 | P | C | C | P | →P | C | P | P | C |
| 2 | C | R→I | R | P | C | R→I | C | C | R |
| 1 | I | I | I | R | I | I | I | R | R→I |

* **FUSION, unchanged, reproduces the three regions:**
  - patterns in thin slabs above a threshold N0 that grows with layers: < 150 at 1 layer, 150–600 at 2,
    > 2400 at 4
  - polarisation in slabs thick relative to the aggregate
  - central clusters in between
* **Plate friction (lambda_wall 10)** keeps small, thin aggregates central for longer, but also slows sorting,
  so the monolayer becomes almost salt-and-pepper (S 0.30–0.39). It is **not needed** with FUSION. Kept in the
  model as an option, default 0.
* **Open questions are in `claude-answers/02-model-plates-screen.md`:**
  - which sparseness to pick (`f_rep` BB 4.5 vs 6)
  - whether rings count as patterns
  - the readout time
  - rigid plates without friction
  - the diluted core

### Stage 5a: tuning with 3 seeds (2026-10-07)

**Setup.**
- Scripts: `scripts/stage5_run.jl` (queue workers), `scripts/stage5_pool.sh` (keeps 24 running), and
  `scripts/stage5_tune_report.jl`.
- Job lists in `results/stage5/`:
  - `tune.csv`, `tune_fast.csv`: wave A (stiffness, both cores)
  - `tune_fric.csv`: wave B (friction)
  - `tune_c.csv`: wave C (the final candidate)
- Rigid plates, 3 seeds; every run continues to fB = 0.10.
- Classes come from the re-scored trajectories (current classifier).
- The diluted-core runs of the clearly failing sets were dropped from the queue part-way; see the decisions
  log.

**2 layers, N0 150** (should stay central). The class and B offset at the readout, by seed:

| set | core | friction | outcome | B offsets |
|---|---|---|---|---|
| july_frB4.5 (CC 5.1) | sharp | 0 | 1 undecided, 2 polarised | 0.36, 0.68, 0.54 |
| july_frB5 (CC 5.1) | sharp | 0 | central, undecided, polarised | 0.20, 0.39, 0.66 |
| FUSION (CC 6) | sharp | 0 | central, undecided, central | 0.08, 0.29, 0.19 |
| every set tried | diluted | 0 | polarised (20 of 21 runs by fB 0.15) | — |
| july_frB4.5 | sharp | 3 | undecided, polarised, polarised | 0.29, 0.59, 0.58 |
| FUSION | sharp | 1.5 / 3 | central, 6 of 6 | 0.02–0.24 |
| july_frB5_frC6 | sharp | 1.5 / 3 | central in 6 of 6, 2 drifting | 0.09–0.33 |
| july_frB4.5_frC6 | sharp | 0 | central, polarised, central | 0.25, 0.60, 0.30 |
| **july_frB4.5_frC6** | **sharp** | **3** | **central, 3 of 3** | **0.17, 0.29, 0.01** |
| july_frB4.5_frC6 | sharp | 5 | central, 3 of 3 | 0.03, 0.17, 0.28 |

**The other test points**, for july_frB4.5_frC6:

| point | friction 0 | friction 3 | friction 5 |
|---|---|---|---|
| 1 layer, N0 75 | loose central / ring | central-ish (offset ≤ 0.15) | central-ish |
| 1 layer, N0 1500 | islands, S 0.48–0.50 | ring → islands, S 0.40–0.41 | islands, S 0.38–0.39 |
| 2 layers, N0 2400 | ring → islands, S 0.60–0.62 | ring → islands, S 0.51–0.53 | ring, S 0.54 |
| 4 layers, N0 150 | polarised, 3 of 3 | polarised, 3 of 3 | polarised, 3 of 3 |

BC doublet at 24 h: BB 0.53, CC 0.56, BC 0.86 (glue).

**Looseness at 2 layers, N0 2400** (S at fB = 0.2), set mainly by B stiffness:

| BB stiffness | S |
|---|---|
| 4.5 | 0.59–0.63 |
| 5 | 0.33–0.49 |
| 5.5 | 0.33 |
| 6 (FUSION) | 0.29–0.37 |

Stiff C (CC 6) loosens it a little; friction 3 loosens it to 0.51–0.53.

**Friction level (wave D):**

| point | no friction | friction 1.5 | friction 3 |
|---|---|---|---|
| 1 layer, N0 2400, S | 0.49 | 0.44 | 0.39 |
| 1 layer, N0 1500, S | 0.48–0.50 | 0.42–0.46 | 0.40–0.41 |
| 2 layers, N0 150, central | 2 of 3 | **3 of 3** | 3 of 3 |

**Conclusion: july_frB4.5_frC6 + plate friction 1.5, sharp core** (`PSETS["final"]`, `LAMBDA_WALL = 1.5`).
Figures: `claude-answers/img/03-core-and-friction.png`, `03-looseness.png` (`scripts/stage5_answer_figs.jl`).

### Stage 5b: the phase space (2026-10-07)

**Setup.**
- **Set and runs:** `july_frB4.5_frC6@lw1.5` (the final set), sharp core, rigid plates, 3 seeds per point,
  each run to fB = 0.10.
- **Grid:** layers {1, 2, 3, 4, 6, 8, 12} × N0 {75, 150, 300, 600, 1200, 2400}, plus N0 4800 at 1–4 layers.
- **Refinement:** points between neighbours whose majority differs (`stage5_phase.jl table` →
  `phase_B_refine.csv`).
- **Location:** `results/stage5/phase_B/`.
- **Comparisons:**
  - the same grid without friction, seed 1 (`phase_A0/`)
  - the friction-3 runs that had finished (`phase_A/`)
- **Unconfined reference:** the "free" row, N0 150–2400, 2 seeds (`tune/`).

**Result**, majority over seeds at fB = 0.2. There are 210 runs; every seed is listed in
`results/stage5/phase_B/phase_table_*.md`.

| layers \ N0 | 75 | 106 | 150 | 178 | 212 | 252 | 300 | 424 | 504 | 600 | 714 | 849 | 1200 | 1427 | 1697 | 2400 | 4800 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | | | P | | | | P | | | P | | | C | | | C | |
| 12 | P | | P | | | | P | P | P | o | C | C | C | | | C | |
| 8 | P | | P | | | | P | P | P | o | C | C | C | | | C | |
| 7 | | | | | | | | | | C | | | | | | | |
| 6 | P | | P | | | | P | P | P | C | | | C | | | C | |
| 5 | | | | | | | P | | | | | | | | | | |
| 4 | P | | P | | P | P | C | | | C | | | C | | | C | C |
| 3 | P | | P | P | C | | C | | | C | | | C | | | C | C |
| 2 | P | C | C | | C | | C | | | C | | | C | PAT | PAT | PAT | PAT |
| 1 | PAT | PAT | PAT | | PAT | | PAT | | | PAT | | | PAT | | | PAT | PAT |

*P polarised, C central cluster, PAT pattern (ring or islands), o drifting (between central and
polarised). The "free" row is the unconfined runs, 2 seeds. Seeds: 3 per point, 5 at 1 layer N0 75 and
150 and at 2 layers N0 75.*

**What it shows:**
1. **The three regions are there.**
   - **Patterns:** thin slabs with enough cells.
   - **Polarised:** a gap wide compared with the aggregate.
   - **Central cluster:** in between, and in all thick slabs with large aggregates.
2. **Patterns.**
   - At 1 layer, patterns form at every N0 tried, down to 75 (the smallest points, with 5 seeds, are mixed with
     central).
   - At 2 layers they need N0 ≈ 1300 or more (1200: central, 1427: ring); the ring breaks into islands by
     fB ≈ 0.12.
   - At 3 or more layers there are none up to N0 4800. B stays one compact central blob (spread 0.12–0.15,
     S ≈ 0.9). Newly made C cells can leave the B domain upwards and downwards, so B never perforates.
3. **Polarisation needs a gap wide compared with the aggregate.**
   - The boundary N0 grows with the gap: ≈ 75–106 at 2 layers, ≈ 180–210 at 3, ≈ 250–300 at 4, ≈ 500–700 at
     6–12 layers.
   - Above N0 ≈ 600 nothing polarises within the readout, not even unconfined: at N0 1200 the free
     aggregate is still central at fB 0.2 (offsets 0.10–0.14) and polarises later. So the central region at
     large N0 and thick gaps is "behaving as unconfined". In this model polarisation takes longer in bigger
     aggregates than the B phase lasts.
4. **The readout time matters at the boundaries.** 106 of 210 runs change class between fB 0.3 and 0.12:
   - rings → islands
   - drifting → polarised
   - small central clusters → islands at 1 layer

   `figures/phase_space_late.png` shows the same runs at fB 0.12.
5. **The plates hold the thickness:** in slabs of 1–3 layers, layers_eff / nominal = 1.003 ± 0.001 (93 runs).

**Against your sketch.** The structure is the same. Differences:
- The pattern threshold rises much more steeply with thickness: no patterns at 3–4 layers.
- The smallest monolayers already pattern.
- Large thick aggregates stay central, because they also do unconfined in this model.

**The friction-free comparison** (`phase_A0`: the same grid, seed 1, plus the 3-seed tuning points) gives
nearly the same map: the same regions, with the boundaries within one grid step. Friction only matters at
the edges:
- 2 layers, N0 150: 2 of 3 central without it, 3 of 3 with it.
- 1 layer, N0 75: central in 3 of 3 without it, mixed with it.
- Monolayer islands: S ≈ 0.49 without, ≈ 0.44 with.

So friction is fine tuning, not essential (`claude-answers/03-tuning-and-phase-space.md`).

### Stage 7: patterns in slabs of 3+ layers (2026-10-08)

**The problem.** In the final set every slab of 3 or more layers ends as one compact central blob, at every
N0 up to 4800: the largest B cluster holds 100 % of B, with looseness S ≈ 0.9. This is not the classifier
being strict: there is only one cluster. The July runs had islands at 3 layers (N0 ≥ 2100) and 4 layers
(N0 ≥ 3200).

**What the July runs did differently**, from the July notebooks (your laptop copy included) and the July
report:

| | July | now |
|---|---|---|
| cell types | A → B → C, rates slowed by the A fraction (`p/(1+k fA)`, `q/(1+k fA)`) | B → C at `q` |
| core | a cylinder of radius 0.8 R, **shuffled: 75 % B, 15 % A, 10 % C** (A: stiff, barely sticks to A, sticks to B) | the innermost 60 %, pure B |
| `f_rep` BB / CC | 3 / 5.1 | 4.5 / 6 |
| plates | springs: nominal 3 layers was ≈ 3.5, 4 was ≈ 4.5; no friction | rigid, exact thickness; friction 1.5 |
| pattern call | ≥ 2 islands still present at the end of the run (fAB ≈ 0.04); 2–3 islands at 3–4 layers, largest 40–53 % | class at fB 0.2; one dominant cluster is central |

**Screen 1** (`results/stage5/thick/`, seed 1, 3 layers N0 2400; `scripts/stage5_run.jl` jobs
`thick.csv`, `thick2.csv`). Every change from your list and from the table above:

| change | 3 layers, N0 2400 at fB 0.2 | S |
|---|---|---|
| B 65 % (`fB0.65`) | central | 0.88 |
| B diluted outwards, 90 % / 80 % inside (`fin0.9`, `fin0.8`) | central | 0.88 / 0.86 |
| B 65 % + diluted (`fB0.65+fin0.85`) | central | 0.88 |
| **C mixed into the core, the July seeding**: 10 / 20 / 30 % (`cc0.1`, `cc0.2`, `cc0.3`) | central | 0.89 / 0.91 / 0.89 |
| B 65 % + 15 % C in the core | central (also at 4 layers) | 0.89 |
| spring plates (`@hookean`) | central | 0.91 |
| plate friction 5 | central | 0.89 |
| the July mechanics (BB 3, CC 5.1), sharp / C in the core / + springs | central | 0.90–0.91 |
| `f_att` CC 3 → 3.6 | central | 0.88 |
| `f_att` BB 6.4 → 5.6 | central, but looser | 0.72 |
| **`mu` BB 2 → 1.8** | **ring, from fB 0.3 to 0.12** | **0.30** |
| λ × 2 / × 3 / × 4 (slower mechanics; a test only) | central / central / central (ring by fB 0.15) | 0.84 / 0.66 / 0.52 |

- The mixed core starts salt-and-pepper, but B re-sorts into one blob within a few hours.
- In a slab of 3+ layers the cells can pass over and under each other, so the new C cells are squeezed out of
  the B domain before they can split it.
- In 1–2 layers they cannot, and the domain breaks up.
- **λ is not a way out.** At λ × 2 and × 3 the unconfined N0 600 aggregate no longer polarises by the
  readout (central, offset 0.06–0.09).
- **What does work is B–B cohesion.**
  - B–B adhesion energy (`mu` 2 → 1.8): W_BB 1.02 → 0.59.
    - The core turns porous by fB 0.4 and opens into a loose network of B by fB 0.3–0.2.
    - That network looks like the sparse experimental patterns.
    - It is close to salt-and-pepper by fB 0.12.
  - `f_att` BB 5.6 (W_BB 0.89) only loosens the blob.
**Screen 2: the B–B range** (`thick4.csv`, `thick5.csv`; seed 1 unless noted; class at fB 0.2, S in brackets):

| `mu` BB (W_BB) | 3 layers N0 2400 | 3 layers N0 1200 | 2 layers N0 150 / 300 | 1 layer N0 600 | unconfined N0 600 | BC fusion (BB, CC) |
|---|---|---|---|---|---|---|
| 2 (1.02), Oct 7 | central (0.89) | central | central | islands | polarised | 0.86 (0.53, 0.56) |
| **1.9 (0.79)** | **ring (0.45)**, islands by fB 0.12 | central (0.61) | **central (0.47–0.55, 2 seeds) / central (0.57)** | islands, 6 (0.43) | **polarised** (offset 0.66) | **0.91** (0.54, 0.56) |
| 1.85 (0.68) | ring (0.35) | | ring (0.38) / ring (0.36): B dissolves | islands, 7 | polarised | 0.94 |
| 1.8 (0.59) | ring (0.30) | ring (0.31) | ring (0.36) / ring (0.33): B dissolves | islands, 6 | | 0.91 |
| 1.8, `f_att` BB 11.1 (1.02) | central (0.91) | | central (0.70) / central (0.83) | islands, 9 | polarised (0.45) | 0.82 |

- **At `mu` BB ≤ 1.85, B dissolves.** W_BB is then about W_CC (0.65), so B is no longer more cohesive than C.
  Small aggregates turn salt-and-pepper instead of keeping a central cluster. The classifier calls that
  "ring", because the loose B cells link into one cluster, but it is mixing, not a pattern.
- **At 1.9, B is still more cohesive than C.**
  - Small cores stay central.
  - A big core in a thick slab hollows out into a ring, which breaks into islands by fB 0.12.
- **Every other behaviour still holds:**
  - Monolayers give sparse, irregular islands, which look closer to the experimental examples than before.
  - The unconfined aggregate polarises.
  - BC glues: 0.91, nearer Sotiris' 0.95 than the Oct 7 set's 0.86.
- **It is the energy, not the range.** With `f_att` BB 11.1, the shorter range keeps the Oct 7 adhesion
  energy, and everything is as in the Oct 7 set: the 3-layer slab is one compact blob (S 0.91), and small
  slabs are compact.
- **λ is not needed, and does not help.** Slower mechanics (λ × 2–4) keep the thick blob central at the
  readout and stop the unconfined aggregate from polarising in time.

**Decision (10-08, 03:05): `mu` BB 1.9 becomes the final set** (`PSETS["final"]`;
`"july_frB4.5_frC6_muB1.9@lw1.5"` in the scripts). It changes one more BB entry and nothing else:
- the BC entries, `q`, λ, α_ecm and `tau_div` are untouched
- sharp core, rigid plates, plate friction 1.5

The Oct 7 set is kept as `PSETS["final_oct7"]`, with its figures in `code-claude/figures-oct7/`.

**Phase space of the new set** (`results/stage5/phase_D/`):
- **Grid:** layers {1, 2, 3, 4, 6, 8, 12} × N0 {75, …, 2400}, plus N0 4800 at 1–4 layers.
- **Seeds:** 2 per point, 1 for N0 4800.
- **Free row:** unconfined N0 150–2400, 2 seeds.
- **Reused:** the screen runs of the set, copied from `results/stage5/thick/` with their rows (same seeds,
  so the same simulations).

**Results** (`figures/phase_space*.png`; the class at fB 0.2):

| layers | polarised up to N0 | central | patterns from N0 |
|---|---|---|---|
| 1 | — | 75–150 | ≈ 160 (islands; 178 is 2 of 3) |
| 2 | 75 (one seed "other") | 89–212 | ≈ 250 (rings, islands by fB 0.12; 300 is 2 of 3) |
| 3 | 150–178 | 212–1200 | ≈ 1300 (rings; 1200 central and 1427 ring, 3 of 3 each) |
| 4 | 212 | 300–3394 | between 3394 and 4800 (ring at 4800; 3394 central in 2 of 2) |
| 5 | 300 | 600 | (not run above 600) |
| 6 | 600 | 714–2400 | none up to 2400 (4800 not run) |
| 7 | 849 | | (one point) |
| 8, 12 | 1009 | 1200–2400 | none up to 2400 |
| unconfined | 600 | 849 (mixed: 2 other, 1 central)–4800 | none up to 4800 (N0 4800: central, S 0.81, B offset 0.01) |

Seeds: 2 per point, 3 at the points next to a threshold, 1 at N0 4800; 178 runs in all.

- **The three regions of your sketch, with patterns at 1–4 layers.**
  - The pattern threshold rises steeply with thickness.
  - The smallest monolayers are now central clusters (N0 75–150), as in the sketch. With the Oct 7 set they
    patterned.
- **The pancake diameter at squeezing explains most of the threshold.** The pattern threshold spans
  N0 ≈ 170 to ≈ 4000 over 1–4 layers, a factor of 25. In diameter it spans only 8–17 cell diameters
  (≈ 80–170 µm if r = 5 µm). Last central → first pattern, in cell diameters (`phase_space.csv`, column D):

  | layers | last central | first pattern |
  |---|---|---|
  | 1 | D 12.4 (N0 150) | D 13.4 (N0 178, 2 of 3) |
  | 2 | D 7.8 (N0 212) | D 9.0 (N0 300, 2 of 3) |
  | 3 | D 11.6 | D 12.5 |
  | 4 | D 14.4 | D 16.8 |

  It shows in `phase_space_diameter.png`.
  - **The threshold diameter still rises with thickness** from 2 to 4 layers, by about 3 cell diameters
    per layer. The monolayer needs a wider pancake than the bilayer.
  - **Why the threshold in N0 rises so steeply:** at a fixed diameter, N0 grows in proportion to the number
    of layers.
  - **Why slabs of 6+ layers never pattern in the grid:** their N0 ≤ 2400 pancakes stay at D ≤ 10.
  - **Compared with July:** the July report found no single cut in D/h either.
- **Polarisation** works as before: a gap wide compared with the aggregate, and N0 up to ≈ 600–1000.
- **Thickness:** the plates hold it, with layers_eff / nominal = 1.002 ± 0.001.
- **The late readout** (fB 0.12, `phase_space_late.png`), compared with fB 0.2:
  - Polarisation spreads, to N0 ≈ 1200 at 6–12 layers and 849 unconfined. Big aggregates polarise late.
  - The smallest monolayer patterns (N0 212–300) shrink back into a central cluster as B runs out, so the
    1-layer threshold moves to ≈ 400.
  - The patterns at 2–4 layers stay, and the 2-layer rings turn into islands.

## Caveats

* **Friction is a modelling choice, and a small one.** Plate friction (λ_wall = 1.5) makes small squeezed
  aggregates stay central reliably. Without it they polarise in-plane in about 1 run in 3, and the rest of
  the phase space barely changes (`phase_A0`). It stands for adhesion to the confining surfaces, which we
  have not checked against the experimental setup.
* **Monolayers come out looser than at 2 layers**, because friction slows sorting most where every cell
  touches a plate. At 1 layer, N0 1500–2400, S is ≈ 0.44 with friction 1.5, against 0.49 without.
* **The diluted core is not used.** It makes small squeezed aggregates polarise.
* ~~**No patterns at 3 or more layers**~~ (Oct 7 set). **Solved on Oct 8** by weaker B–B cohesion (`mu` BB 1.9):
  patterns at 3 layers from N0 ≈ 1300 and at 4 layers at N0 4800. Slabs of 6+ layers do not pattern up to
  N0 2400 (their pancakes stay narrow); larger N0 there was not run.
* **`mu` BB 1.9 sits close to the edge.** At 1.85, B dissolves into C (W_BB ≈ W_CC). The margin is about
  0.05 in `mu` BB, or 13 % in B–B adhesion energy (0.79 → 0.68).
* **2–3 seeds per point** in the Oct 8 phase space (3 next to the thresholds, 1 at N0 4800); the Oct 7 one had
  3–5.
* **The central-to-ring transition is gradual, so the thresholds depend on the classifier cut.** As N0 grows
  past the threshold, the central cluster hollows out. Its `fill` falls smoothly:
  - ≈ 0.65 for clearly central
  - 0.50–0.54 for "opening" clusters (2 layers N0 300–357, 3 layers N0 1427), which by eye are in between
  - 0.33–0.48 just above
  - ≈ 0.22 for clear rings

  The cut is CENTRAL at fill ≥ 0.55, set on Oct 7 and not changed. With a stricter 0.45 the thresholds move
  one refinement step up: ≈ 400 at 2 layers, ≈ 1700–2400 at 3 layers. Everything else stays.
* **Big aggregates do not polarise within the B phase**, not even unconfined. At N0 1200 polarisation comes
  only after fB 0.2; at N0 2400 the B offset never exceeds 0.08. The polarised region therefore stops at
  N0 ≈ 600 (Oct 8 set: 600–1000), at any thickness; it grows to ≈ 1200 by fB 0.12.
* **Small monolayers (N0 75–150) are noisy:** only about 30 B cells are left at the readout. Five seeds there
  split between pattern and central cluster.

