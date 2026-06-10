# Validation Rebuild — Phase 1: Honest Measurement

**Branch:** `validation-rebuild` · **Date:** 2026-06-10 · **Scope:** measure, don't fix.

Phase 1 does **not** repair any pattern's validation. It instruments the
negative-discrimination panels so we can state, per pattern and without
self-deception, *how much of the "validation" is real*. Two instruments were
added, both in `epc/phase2a/panel.py`.

---

## Instrument 1 — `rejection_stage`

For every **negative** substrate the panel runs through a detector, we now
record how it was handled:

- **`FIRED`** — the detector fired on a negative → a **false positive** (the
  discriminating metric did *not* reject it).
- **`GUARD:<reason>`** — rejected by a **prerequisite/guard** that runs *before*
  the discriminating metric (e.g. a type-constancy check, a run-length floor, a
  substrate-format check). The headline metric never ran.
- **`METRIC`** — reached the discriminating metric and it did the rejection.
  This is the only stage that represents the detector's *advertised* mechanism
  doing genuine work.

## Instrument 2 — `cohens_d` audit

The panel's effect-size gate is `cohens_d_positive_vs_panel >= 1.0`. It was
computed over `_score = result.confidence` — the **discrete tier-confidence**
(≈0.9 for everything that fires, 0.0 for everything rejected). When positives
all fire and negatives all reject, both arrays are constant → `cohens_d`
returned `float("inf")` → `inf >= 1.0` → **auto-PASS**, with the verdict logic
even rationalising it as a "degenerate-perfect detector". This is a **vacuous
effect size**: there is no within-group spread to standardise against.

---

## Finding 1 — Guard-dependence (catalog-mate negatives)

| pattern | mates | FIRED (fp) | GUARD-gated | METRIC | verdict |
|--|--|--|--|--|--|
| p1  | 7 | 0 | **7** | 0 | GUARD-GATED |
| p2  | 5 | 0 | **5** | 0 | GUARD-GATED |
| p3  | 2 | 0 | 0 | 2 | ok (metric) |
| p4  | 2 | 0 | 0 | 2 | ok (metric) |
| p5  | 5 | 0 | 0 | 5 | ok (metric) |
| p6  | 5 | 0 | 0 | 5 | ok (metric) |
| p7  | 5 | 0 | 0 | 5 | ok (metric) |
| p8  | 3 | 0 | **3** | 0 | GUARD-GATED |
| p9  | 3 | 0 | 0 | 3 | ok (metric) |
| p10 | 3 | 0 | **3** | 0 | GUARD-GATED |
| p11 | 7 | 0 | 0 | 7 | ok (metric) |
| p12 | 7 | 0 | 0 | 7 | ok (metric) |
| p13 | 7 | 0 | 0 | 7 | ok (metric) |
| p14 | 7 | 0 | 0 | 7 | ok (metric) |
| p15 | 7 | 0 | 0 | 7 | ok (metric) |
| p16 | 2 | 0 | 0 | 2 | ok (metric) |
| p17 | 5 | 0 | 0 | 5 | ok (metric) |
| p18 | 3 | 0 | 0 | 3 | ok (metric) |
| p19 | 5 | 0 | 0 | 5 | ok (metric) |
| p20 | 2 | 0 | 0 | 2 | ok (metric) |
| p21 | 3 | 0 | 0 | 3 | ok (metric) |
| p22 | 7 | 0 | 0 | 7 | ok (metric) |
| p23 | 2 | 0 | 0 | 2 | ok (metric) |
| p24 | 2 | 0 | 0 | 2 | ok (metric) |
| p25 | 2 | 0 | 0 | 2 | ok (metric) |
| p26 | 2 | 0 | 0 | 2 | ok (metric) |
| p27 | 7 | 0 | 0 | 7 | ok (metric) |
| p28 | — | — | — | — | **NO PANEL** |
| p29 | 2 | 0 | 0 | 2 | ok (metric) |
| p30 | 2 | 0 | 0 | 2 | ok (metric) |
| p31 | — | — | — | — | **NO PANEL** |
| p32 | 2 | 0 | 0 | 2 | ok (metric) |

**Guard reasons doing the gating:** `type_constancy_failed` ×6,
`run_too_short` ×5, `substrate_mismatch` ×3, `no_coexistence` ×2,
`multi_cluster_failed` ×1, `pos_vel_ac_below_floor` ×1.

**Read-out:**

- **4 patterns are guard-gated** (P1, P2, P8, P10): *100%* of their catalog-mate
  negatives are killed by a prerequisite check before the advertised metric
  runs. For these the headline metric is **decorative**. This confirms the code
  audit's P1 result directly — with the guards bypassed, Moran's I
  false-positives on 3/7 of P1's catalog mates (P11=0.344, P12=0.058, P15=0.269).
- **2 patterns have NO PANEL** (P28, P31).
- **26 patterns reach the metric** and it rejects the mates with zero false
  positives. This is **necessary but not sufficient**: it does *not* verify (a)
  that the catalog mates are real same-substrate lookalikes that genuinely
  stress the metric, or (b) that the metric is the canonical phenomenon metric
  rather than a proxy. Those are Phase 2 / Phase 4.
- Two of the guard reasons — `substrate_mismatch` and `run_too_short` — mean the
  "catalog mate" was excluded on a **technicality** (wrong substrate type / too
  short) rather than discriminated. Those panels need *real* same-substrate
  mates built in Phase 2.

---

## Finding 2 — The effect-size gate cannot fail

Pre-fix scan of all 30 panels with a panel:

- **30/30 = PASS or PASS-with-weakness. 0 FAIL, 0 PARTIAL.**
- **15/30 report `cohens_d = Infinity`** (degenerate-constant scores → vacuous).
- The other 15 report *finite* d (2.6–14.5) — but computed over **discrete
  tier-confidence buckets**, not the canonical continuous metric. "Less
  obviously broken," still not faithful.
- **No panel computes an honest positive-metric-vs-negative-metric separation.**
  The PASS criterion is non-discriminating: it passes everything.

---

## Fix applied in Phase 1 (`epc/phase2a/panel.py`)

This is the *only* code repair Phase 1 makes — it stops the system from
*reporting* a lie, without yet claiming any new validation:

1. **`cohens_d`** — degenerate-constant inputs now return **`NaN`** (was
   `±inf`). NaN is the honest value: a standardized effect size is undefined
   without within-group spread.
2. **`overall_verdict`** — new label **`TNR-PASS-EFFECT-UNDEFINED`**: negatives
   are correctly rejected (TNR passes its gate) but the standardized effect size
   is uncomputable over the current scores. Deliberately distinct from **PASS**
   (the effect-size criterion is *unmet-because-undefined*, pending a
   continuous-metric recompute) and from **FAIL** (the detection separation is
   real).
3. **`summary`** — additive honesty fields: `effect_size_undefined` (bool),
   `effect_size_basis: "tier_confidence"`, and an `effect_size_note` stating the
   effect size is not yet over the canonical metric.

**Verified** on P3 (the "cleanest" pattern): was `cohens_d: Infinity → PASS`,
now `cohens_d: nan`, `effect_size_undefined: true`,
`verdict: TNR-PASS-EFFECT-UNDEFINED`.

**Post-fix verdict distribution (all 32, re-run 2026-06-10 16:06 UTC):**
`d = Infinity` count is now **0/32** (was 15). Verdicts:

- **15 × `TNR-PASS-EFFECT-UNDEFINED`** — the degenerate panels; effect size is
  now honestly undefined (P1, P3, P5, P6, P8, P11, P12, P13, P16, P18, P20, P22,
  P24, P25, P26).
- **12 × `PASS` + 3 × `PASS-with-weakness`** — finite d (2.6–14.5), but still
  computed over discrete tier-confidence (`effect_size_basis: "tier_confidence"`),
  **not** the canonical metric. These escaped the `undefined` flag only because
  some negatives drew *intermediate* confidence, giving non-zero spread — the d
  is defined but not yet faithful. None is a canonical-metric effect size;
  faithfulness is still Phase 2 work.
- **2 × `NO PANEL`** (P28, P31).

---

## Honest state after Phase 1

- The **TNR machinery is sound** — negatives do get rejected (0.95–1.0).
- The **effect-size validation is not faithful for any pattern** — 15 vacuous
  (Infinity→NaN), 15 over discrete confidence tiers, 0 over the canonical metric.
- **4 patterns are guard-gated** (metric decorative): P1, P2, P8, P10.
- **2 patterns have no panel**: P28, P31.
- **No pattern has yet demonstrated faithful, guard-independent, continuous-metric
  discrimination.** That is the work of Phases 2–4.

## Triage into later phases

- **Phase 2 (rebuild panels):** every pattern needs a faithful effect size
  recomputed over the continuous discriminating metric *and* real same-substrate
  catalog mates. Highest priority: the 4 guard-gated patterns and the
  `substrate_mismatch` / `run_too_short` panels.
- **Phase 3 (downgrade/pull):** P28, P31 (no panel) plus the non-reproducing /
  ignore-input patterns from the code audit (P7, P9, P15, P17, P30, P32).
- **Phase 4 (faithfulness):** detectors whose metric is a proxy (P21, P24, P13)
  → switch to the canonical metric or disclose the proxy.
