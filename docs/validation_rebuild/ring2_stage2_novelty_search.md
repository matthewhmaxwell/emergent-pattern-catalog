# Ring 2 — Stage-2 novelty search (the audited frontier)

The genuine-novelty hunt: run the live `ring2_descriptor` over generative-family parameter
sweeps, cluster each observation by its lens fingerprint, and let the model-free tripwire
surface the COMPLEX-and-UNCLASSIFIED quadrant (structure no admitted lens can name). Survivors
are vetted (artifact / mundane check, then the literature gate). Honest prior (Stage-1): genuine
novelty is rare; the deliverable is an AUDITED, DISCLOSED FRONTIER, not a manufactured claim.

Harness: `analysis/ring2/stage2_novelty_search.py <N> <family>`; results JSON per family;
lead vetting `analysis/ring2/_stage2_vet_leads.py`. Run via the project venv.

## Run 1 — particle_life (N=250)
- **0 tripwire leads / 250.** Two fingerprint clusters: `{structure_factor + persistent_homology
  + fractal_dimension | positions}` (158) and `{… | orientation}` (92) — all named-classified.
- Significance: Stage-1's tripwire-ONLY baseline produced 4/300 mundane leads on this same
  family. The broadened battery + hardened tripwire give **0/250** — the upgrade eliminated the
  false leads (higher precision, the entire point of the lens-comprehensiveness thesis).

## Run 2 — Lenia (N=100) — ASAL's own search space
- 7 fingerprint clusters (field-kind, heavy-tail(power-law), temporal(spectral-peak),
  synergy) — Lenia is genuinely more diverse than particle_life, and the battery resolves it
  into interpretable sub-regions.
- **4/100 tripwire leads → ALL vetted MUNDANE.** Each (configs #0/#19/#29/#34) is a DEAD field:
  mass collapses to 0 by mid-run, final field all-zero. The death TRANSIENT is temporally
  complex (C 0.19–0.24, struct 0.16–0.22 pass the hardened gate) but the final state is trivial
  (em 0.06). A genuine false-lead mode, not novelty.

## Instrument hardening surfaced by Stage-2
- **Dead-state gate** added to the tripwire (`model_free_complexity`): a system that COLLAPSES
  to a trivial constant final state (death / dissipation) has a complex transient but no
  SUSTAINED structure → not a lead. Detected via spatial variance across components dropping to
  ~0 in the late window (late < 5% of early). Re-validated: stage-0 UNCHANGED (0/3 nulls, 17/17
  classified); all 4 Lenia death-transient leads drop; Lenia re-run → 0 leads.
- This is the second OOD tripwire hardening from substrate exposure (the first: the iid-noise
  FIELD surrogate gate). Both make the novelty search higher-precision.

## Honest conclusion
**0 genuine-novelty leads this pass** — the prior held. The deliverable is the AUDITED FRONTIER:
the lens battery characterizes the explored generative space (particle_life + Lenia) with
interpretable, validated axes, and NOTHING complex-and-unclassified survived vetting. That is a
legitimate, disclosed result — the differentiation from ASAL is precisely that our frontier is
interpretable and audited, not an opaque descriptor's cluster map.

## Next
Broaden coverage: more generative families (reaction-diffusion (F,k) grid; wider Lenia / 
particle-life parameter ranges; coupled-oscillator regimes) and larger N. Any future
complex+unclassified survivor goes through the literature-novelty gate before any claim.
FM-as-lens reconsidered only if it earns the bar (ring2_strategy.md).
