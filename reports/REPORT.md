# Squeezing gastruloids: living report

*Updated as the work goes. Last update: 2026-10-06, stage 1 (reorganisation) done.*

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
- [ ] **2. Model and initial condition:** B → C only; B = 60 %, sharp vs diluted; wall comparison (checkpoint 2)
- [ ] **3. Readout and classifier**, calibrated by eye against `experiments/pattern-examples/` (checkpoint 3)
- [ ] **4. Parameter screen:** B/C entries of `f_att`, `f_rep`, `mu`; BC frozen (checkpoint 4)
- [ ] **5. Phase space:** rough pass, then refined (checkpoint 5)
- [ ] **6. Final 3D figures with saved data**, plus `FIGURES.md` and `MODEL.md`

## Results

*(filled in as the stages complete)*

## Caveats

*(filled in as they come up)*
