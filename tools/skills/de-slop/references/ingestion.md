# Ingestion flow (curation mode)

How a piece of slop becomes a memorialized rule. This runs only when something
is flagged for capture. It is off the hot path: drafting and rewriting do not
pay this overhead.

## Trigger

- The user pastes text marked `slop:` (or says "this is slop" + paste).
- I catch my own slop mid-conversation and flag it for capture.
- The user points at a published example and says "ingest this".
- A model self-report during a skill revision (the 2026-09 batch). These file
  as provisional and are promoted or retired on first contact with real work.

Before step 1, run `scripts/slop_scan.py` on the pasted span. It names the
mechanical pattern (dash, contrast shape, tricolon, register hit, label colon)
and gives the anatomy a starting point. It does not name the mechanism; that
is step 2, and it is the part that matters.

## The six steps

For each flagged span, produce:

1. **Anatomy.** Quote the exact span. Name the phrasing and the grammatical move
   (copula avoidance, present-participle tail, tricolon, negative parallelism,
   false agency, etc.).
2. **Mechanism.** Why the model generates it. Tag to a cause:
   - `reward-tuning`: hedging, sycophancy, both-sidesing, over-balancing.
   - `repetition-penalty`: synonym cycling, elegant variation.
   - `instruction-tuning`: list reflex, rule of three, header stacking,
     signposting.
   - `pretraining-register`: formal-corpus words (delve, tapestry, robust).
   - `uncertainty-conditioning`: cutoff disclaimers, over-qualification.
   - `assistant-persona`: "let me know if", "I hope this helps", "great
     question".
   - `displacement`: a new tell created by a crackdown on an old one (see the
     anti-em-dash entry in the living corpus, and the displacement table in
     current-register.md).
   - `register-migration`: a whole vocabulary tier trained out and refilled
     from a different corpus (the post-delve register: load-bearing,
     non-trivial, the delta). Distinct from displacement, which is one rule
     pushing one impulse into one new form.
   - `self-edit`: a tell that appears only when the model edits its own draft
     under these rules (flattened rhythm, register swap, a new closer).
   Naming the mechanism is the differentiator. It is what lets us predict the
   next tell before it is common.
3. **Context.** Where it is actually fine versus where it is a tell. Not
   everything is a universal ban.
4. **Strength.** Assign a tier, never a blunt "never" unless it is a true
   universal (em dash, the fatal pattern):
   - Tier 1: always replace.
   - Tier 2: flag in clusters (2+ in a paragraph).
   - Tier 3: flag by density.
   - context-dependent: fine in some profiles, a tell in others.
5. **Replacement.** The rewrite, plus the general rule behind it so it
   generalizes past this one example.
6. **Dedup + file.** Check the living corpus and patterns.md. If it is already
   covered, merge as a variant. If new, write a fresh entry.

## The protect-list cross-check (mandatory before filing)

Run every proposed entry against `protect-list.md` (canonical: your voice spec).
If the pattern collides with one of your protected signatures (for example, you
pasted something that happens to contain a phrase on your protect list), STOP.
Surface the collision and ask before filing. Do not let a slop entry accidentally
ban your own voice.

## File-step friction (your dial)

- **Deliberate flag** (the user pastes `slop:`): show the proposed entry, file on
  a one-word ok.
- **Self-caught** (I notice my own slop mid-conversation): auto-file and show a
  revertible diff, so it stays frictionless.
- A protect-list collision overrides both: always ask first.

You own this dial. Tighten or loosen it any time.

## Entry template (written into living-corpus.md)

```
### <short-name-of-tell>

- Added: YYYY-MM-DD  |  Tier: 1 | 2 | 3 | context-dependent
- Mechanism: <one of the tags above>
- Context: <where it is fine vs a tell>
- Before: "<the AI version>"
- After: "<the human version>"
- Rule: <the generalized directive>
- Source: <the writer in the wild | my own output | published example URL | model self-report (provisional)>
```
