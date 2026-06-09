# EPC Instrument / Out-of-Distribution (OOD) Readiness Roadmap

**Purpose:** Ensure the Emergent Pattern Catalog ships as a **usable measurement instrument**
— a detector battery you can point at *any* system and ask "what emergent behaviors does this
exhibit, with what calibrated confidence?" — not merely a validated taxonomy of 32 known models.

This roadmap exists because the catalog's value engine is the **outward** use: pointing the
battery at multi-agent AI systems (LLM agent swarms, RL collectives), understudied ABMs, and
real time-series. The dim4 panel only proves detectors discriminate *within our native model
family*; using them on alien substrates is out-of-distribution and needs explicit support.

---

## The OOD problem statement

Each detector currently:
- consumes its own model's specific state-history format, and
- has screening/confirmation/definitive thresholds tuned on native positives + the negative panel.

Pointing it at an alien system (different ABM, LLM agents, real data) fails on two axes:
1. **Interface** — the alien system doesn't speak the detector's input format.
2. **Calibration** — raw tier/confidence means nothing off-distribution; we can't trust a score.
3. **Coverage** — the battery must be able to say "strong emergent signal, matches NONE of the 32"
   (essential for both honest measurement and for discovery of new patterns).

---

## Track 1 — Bake into Milestone B (cheap now, expensive to retrofit)

Add to EVERY remaining new-pattern *implement* brief (Waves 2–4) and apply to Wave 1 patterns
in the Milestone C validation pass:

- **T1a — Canonical observation contract.** Each detector must read from a stable, documented
  "observation bundle" schema (e.g. agents × time × {position, velocity, discrete-state}, scalar
  fields, network adjacency, opinion/scalar vectors) via an adapter — NOT directly from the model
  object. Formalize the existing adapter pattern into one public input contract (`docs/observation_schema.md`).
- **T1b — Cross-model generalization test.** Each detector must fire on **≥1 independent
  implementation** of its phenomenon that it was NOT tuned on (e.g. P5 flocking detector also fires
  on a Couzin/Reynolds flock, not just Vicsek). Extend `tests/test_cross_model.py`. This is the
  minimum evidence that a detector recognizes the *phenomenon*, not its *implementation*.

## Track 2 — Milestone C: instrument layer (after the catalog is complete, ~32/32)

- **T2a — Calibration layer.** Characterize each detector's score distribution on native positives
  vs. the negative panel; emit a **calibrated score/probability** (conformal-style: "matches P5 at
  X% relative to the catalog null distribution") instead of a raw tier.
- **T2b — Novelty / "none-of-the-above".** A battery-level abstention signal: generic emergence
  indicators (order parameters rising, structure forming, entropy dropping) trip while NO specific
  detector reaches confirmation → flag as "emergent but unclassified." This is the same primitive
  that powers the **discovery** use (search rule-space for strong-but-unclassified signatures).
- **T2c — OOD validation suite.** A held-out battery run on systems no detector was built on
  (independent ABMs + at least one real dataset), reporting precision/recall of the battery as an
  instrument. This is what lets us *claim* generalization in the paper.
- **T2d — External demonstration (the headline).** Point the full battery at a **multi-agent
  LLM/AI system** (and/or a canonical real dataset) and report its emergent-pattern profile. One
  well-executed external application is the proof that turns "taxonomy" into "instrument."
- **T2e — Public API/CLI.** `epc-detect <observation-bundle> → calibrated emergent-pattern profile`.
  The reusable-tool deliverable; pairs with the Docker/Zenodo/web-demo companion outputs.

---

## Sequencing

1. **Now → Milestone B (Waves 2–4):** enforce T1a + T1b in each implement brief.
2. **At each wave gate:** confirm new detectors satisfy the observation contract + cross-model test.
3. **Milestone C (post-32/32):** T2a → T2b → T2c → T2d → T2e, in that order (each depends on the
   prior). T2d (external demo) is the highest-impact single deliverable and defines the paper's reach.

## Why this ordering protects the project

- The instrument framing is what makes the catalog *used and cited*, not just read.
- Building the input contract + cross-model test per pattern now costs ~minutes/pattern; retrofitting
  32 detectors to a shared contract later costs a full refactor.
- Calibration/novelty/demo genuinely *require* the complete battery, so deferring them to Milestone C
  is correct, not procrastination.
- "None-of-the-above" (T2b) is the shared primitive behind both the instrument and the discovery use —
  build it once, get both.
