# Ring 2 — multi-family substrate hardening

Exercising the broadened lens battery (Ring-2 admitted lenses + the model-free tripwire)
on substrate families OUTSIDE the catalog's agent/network corpus, to find where the
lenses break out-of-distribution. Each run = a new generative family + the full battery.

## Run 1 — Gray-Scott reaction-diffusion (continuous FIELD substrate)

- Model: `epc/models/reaction_diffusion.py` (Gray-Scott; (F,k) walks the Pearson atlas:
  solitons / mitosis-chaos / stripes / spots / worms / high-F).
- Scripts: `analysis/ring2/_rd_lens_exercise.py` (battery sweep), `_rd_diag.py`
  (structure-factor + field-stats diagnostic), `_tripwire_null_diag.py` (tripwire internals).
  Run via the venv: `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. /home/matthewhmaxwell/epc-venv/bin/python …`.

### Findings

1. **structure_factor works on RD fields** — patterned regimes sk_peak 262–1112 vs a
   synthetic UNIFORM-null 9.2. Clean separation. (A radial S(k) prototype gives the same
   verdict plus the characteristic wavenumber k0≈8 — a possible future enhancement.)
   - **Dispatch lesson:** `structure_factor` checks `positions` BEFORE `field`. The RD model
     originally emitted both (field + spot centroids), so the lens used the sparse centroids
     and under-read (sk_peak 1–8). A substrate must advertise ONE primary representation —
     RD now emits `field` only. General rule for new substrates: don't double-advertise.

2. **Emergence indicator is correct on RD** — all tested regimes fire em 0.70–0.77, kind
   "field", classified, not tripped. None of the tested (F,k) actually die (v.std > 0
   everywhere; "high-F" F=0.06 has the *highest* std 0.134), so the firing is correct —
   the earlier "death/null" label was wrong, not a false positive.

3. **PH-on-centroids cannot see field loops** — stripes/worms return n/a (connected → too
   few centroids). A labyrinth's holes live in the FIELD, not the point set. Confirms the
   deferred persistent-homology field path (sublevel-set filtration on the field).

4. **★ Tripwire false-positive on a temporally-iid uniform-noise field.** The model-free
   bridge flagged a flat field + iid per-frame noise as complex+unclassified (`tripped=Y`,
   em 0.11 correctly low). Mechanism (`_tripwire_null_diag.py`): the macro-selector picks
   the noise-driven `std` macro → C=0.2245 > C_THR 0.16. Robustness boundary:
   - static dead field (`flat_zero_var`) → C=0 (correct, NOT complex);
   - corpus nulls `null_noise` C=0.07, `null_walk` C=0.13 (correctly below threshold);
   - RD_stripes (real pattern) C=0.33 on the `mean` macro (correctly complex).
   The gap = the macro-MAXIMIZING selector over-reports on noise-driven secondary macros,
   an OOD null (field substrate) the stage-0 calibration never saw.
   - **Principled fix — DONE (this pass).** Surrogate structure-gate on the C-path:
     `structure_score = mean(H_shuffle) - H_obs` (permutation-entropy deficit vs time-shuffled
     surrogates of the mean macro). The shuffle shares the identical finite-size sparsity, so
     the bias cancels — iid noise ~0, genuine temporal structure >0. The C-path now requires
     `C > C_THR AND structure_score > STRUCT_THR` (0.12); the psi-path is untouched (sound on
     noise). Calibrated (`_tripwire_surrogate_calib.py`): null struct ≤ 0.059, emergent-with-
     C>thr struct ≥ 0.184 — a clean gap. **Re-validated: exactly one system flips** — the
     uniform-noise field (is_complex True→False) — while the stage-0 baseline is UNCHANGED
     (0/3 nulls trip, 17/17 classified, every complex corpus system still complex) and
     RD_stripes stays complex. (Follow-up, tested: the same surrogate principle does NOT transfer
     to the deferred `recurrence` lens — FT/shuffle/long-range rescues all empirically overlap
     null_walk with genuine emergents in `_recurrence_surrogate_test.py`; recurrence needs a
     proper delay embedding, not a surrogate on the sorted-centered trajectory.)

### Outcome
The broadened battery behaves correctly on the RD field family (field lenses + emergence
indicator). One bridge robustness gap identified AND FIXED (the tripwire surrogate gate, above).

## Run 2 — Kuramoto 2D oscillator lattice (phase substrate)

- Model: `epc/models/kuramoto_lattice.py` (4-neighbour Kuramoto; regimes incoherent /
  global-sync / spiral / travelling plane-wave). Script: `analysis/ring2/_osc_lens_exercise.py`.
  Chosen to stress `directed_info_flow` (the newest admitted lens, which RD couldn't reach).

### Findings

1. **Emergence indicator + tripwire behave correctly OOD** — sync (em 0.85) and spiral
   (0.57) flagged emergent; incoherent (0.12) and the imposed plane-wave twisted state
   (0.01) not; NOTHING false-trips. Order parameter sanity: sync r=0.46, incoherent/waves r~0.

2. **`directed_info_flow` mean_te detects coupling magnitude OOD** — incoherent 0.116 vs
   coupled sync/spiral/plane 0.46–0.54.

3. **★ directionality is correctly ~0 on the lattice — a real-substrate NEGATIVE control.**
   The hypothesis (waves = directional flow) was WRONG, and instructively so: a phase-locked
   travelling/sync wave is informationally REDUNDANT, not directed — every oscillator is a
   deterministic rotator at a fixed offset, so the source's history adds no predictive info
   about the follower beyond its own past → genuine TE ≈ 0 both ways, even with mean_te high.
   Ordering components along the propagation axis (a row transect) did not change this
   (plane_wave dir_row 0.154 ≈ sync 0.138; incoherent's 0.33 is ratio-noise below the
   coupling floor). So `directed_info_flow` distinguishes DIRECTED TRANSFER (a drive cascade —
   its validated positive) from mere coupling / synchrony / spatial propagation. The lattice
   is now a real-substrate negative for the directionality axis, complementing the synthetic
   symmetric-mesh negative — this STRENGTHENS the lens's admission. Scope recorded in the lens
   docstring: mean_te = coupling magnitude; directionality = drive, not propagation.

### Outcome
Both substrate runs behave correctly. directed_info_flow's directionality axis is now
validated against two independent negatives (synthetic mesh + Kuramoto sync/wave) and one
positive (cascade).

## Run 3 — Lenia (autonomous moving creatures; ASAL's own search space)

- Model: `epc/models/lenia.py` (Bert Chan 2019; canonical Orbium seed + params R=13/T=10/
  mu=0.15/sigma=0.015, single Gaussian-ring kernel). Verified the Orbium GLIDES (mass stable
  73.9→72.9, centroid path 178 over 400 steps). Script: `analysis/ring2/_lenia_lens_exercise.py`.
  Chosen because Lenia is the qualitative regime RD/Kuramoto lacked (self-propelled localized
  creatures) AND it is ASAL's substrate — running our interpretable battery on it is the
  differentiation demo (named/validated axes vs one opaque foundation-model descriptor).

### Findings (regimes: single glider / 5-glider / random soup)

1. **The emergence indicator recognises Lenia creatures** — em 0.57–0.75, kind "field",
   across all three regimes; nothing false-trips. A lone gliding Orbium gets em 0.75 (the
   localized dense creature fires the clustering/autocorr channel), correctly CLASSIFIED.
   The hypothesised "autonomous-locomotion blind spot → novelty lead" did NOT materialise —
   the named field channel already covers it. Honest, reassuring negative.
2. **structure_factor works on Lenia fields** — sk_peak 116–126 (localized creatures have a
   clear characteristic scale).
3. **PH-on-positions stays low** (h1_max 0.046–0.068) — Orbium crescents carry no loops;
   consistent with RD, reconfirms the deferred field-sublevel-set PH path for field topology.
4. directed_info_flow not exercised here — Lenia is a field/grid substrate, so it has the
   same spatial-blindness as the Kuramoto wave (a raveled subsample loses creature motion
   direction); a genuine directed-TE positive needs creature-tracking (centroid identities),
   which is left for later.
5. Sanity: single glider moves (path 258, mass stable); 5-glider collides and grows (mass
   2.9x); soup partially self-organises (final/initial mass 0.52). All clean, no crashes.

### Outcome
Across THREE substrate families (RD field / Kuramoto phase / Lenia creatures), the broadened
battery behaves correctly OOD: it recognises genuine emergence, stays quiet on nulls, and the
one gap found (the tripwire field false-positive) was fixed and re-validated. Lenia adds the
ASAL-substrate demonstration. **Field-sublevel-set PH now CLOSED** (persistent_homology gained a
field path: persistent enclosed-hole `field_loop_area`, gaussian-denoised — RD-worms 0.23 /
stripes 0.38 / ring 0.21 / swiss 0.13 vs blob/spots 0 / noise 0.009, gap +0.12; raw count is
noise-confounded so area is the discriminator, the h1_max lesson again; wired into the descriptor,
PH fire-count 8→15, no regression). Remaining lens follow-up: a genuine directed-TE positive via
creature-interaction tracking (optional — the lens is already admitted). The admitted lenses are
now wired into the live descriptor (the Stage-2 novelty-search on-ramp).
