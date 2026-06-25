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
   - **Principled fix (next unit, NOT yet done):** surrogate/shuffle correction — flag
     complex only if C exceeds a TIME-SHUFFLED surrogate's C by a margin (iid noise: ΔC≈0;
     genuine temporal structure: ΔC large). Same tool the deferred `recurrence` lens needs.
     This changes a validated core (`model_free_complexity`), so it must re-validate against
     the probe corpus (keep 0/3 nulls tripping, 17/17 classified) before landing.

### Outcome
The broadened battery behaves correctly on the RD field family (field lenses + emergence
indicator). One bridge robustness gap identified and specified. Next substrate families to
exercise (each hits different lenses): Lenia (moving creatures → positions → PH +
directed_info_flow), coupled-oscillator lattice (→ directed_info_flow + spectral).
