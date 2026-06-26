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

## Broadened pass — 4 families (run 2, seed_base 3000)
Added two families the search had not seen (reaction_diffusion, kuramoto) + wider
ranges on lenia/particle_life. Harness: `stage2_novelty_search.py <N> <family> <seed_base>`.
- **reaction_diffusion (120): 0 leads.** ~all fractal+PH+SF | field — Turing-class fields, named-classified.
- **lenia (100, wider + soup): 0 leads.** 6 clusters; dead-state gate holding (no death-transient leads).
- **particle_life (120, wider): 0 leads.** positions / orientation.
- **kuramoto (100): 7 leads → vetted = TWISTED / PHASE-WAVE states, NOT novelty.** All 7 have global
  coherence r≈0 but HIGH local order (0.75–0.95), init=plane (6) / spiral (1): neighbours phase-aligned
  but the phase winds across the lattice so global coherence cancels. Sustained (dead-state gate leaves
  them — they don't die). They trip because the emergence indicator's phase channel is GLOBAL-coherence-
  keyed (Kuramoto r), so a locally-ordered/globally-incoherent state reads em≈0.06 (unclassified) while
  its drifting macro is model-free-complex. **This is the P33 active-nematic lesson (local vs global order)
  re-surfacing for oscillators.** Honest: these are known Kuramoto twisted/spiral states, IMPOSED by the
  IC and persisting — not novel-to-science, not self-organized. Vetting `analysis/ring2/_stage2_vet_kuramoto.py`.
  -> Surfaced a genuine CHARACTERIZATION blind spot: the phase channel needs a LOCAL-order component.

## Honest conclusion (across both passes — 690 configs, 4 families)
**0 genuine-novelty leads.** The prior held. Every complex+unclassified lead vetted to a known/
mundane explanation: Lenia death transients (dead-state gate) and Kuramoto twisted/phase-wave
states (imposed, known physics). The deliverable is the AUDITED FRONTIER — the lens battery
characterizes the explored generative space (particle_life / Lenia / RD / Kuramoto) with
interpretable, validated axes, and nothing survived vetting. That auditability IS the
differentiation from ASAL's opaque descriptor: we can say exactly what we searched and why
nothing qualified. The search ALSO did its proper second job — it surfaced two instrument gaps,
both now characterized (one fixed: dead-state gate; one identified: the global-only phase channel).

## Next
1. **Close the phase blind spot** (the P33 lesson for oscillators): add a LOCAL phase-order
   component to the emergence indicator so locally-coherent/globally-incoherent twisted & wave
   states are classified — a net-widening (like the catalog's rounds 2-5), with full T2c +
   blind-spot re-validation. This touches the validated catalog core, so it is its own careful unit.
2. Then broaden further (RD (F,k) grid at finer resolution; more oscillator regimes; larger N).
   Any future complex+unclassified survivor goes through the literature-novelty gate before any
   claim. FM-as-lens reconsidered only if it earns the bar.
