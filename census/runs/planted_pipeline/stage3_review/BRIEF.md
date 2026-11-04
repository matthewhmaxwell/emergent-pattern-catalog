# Brief for the Stage 3 (b) review agents, planted pipeline test (written 2026-11-04, before any Stage 3 card or reply existed)

PLANTED_PIPELINE_PLAN step 6 (b): "one Claude research agent per card checks the cited sources with the brief of the
census source review (VERIFICATION_PLAN Stage 3 (b), V4), told to count only sources dated before the target paper's
first version, without being told why. Labels DOCUMENTED, CLOSELY RELATED, UNDOCUMENTED CANDIDATE."

Each agent (claude-opus-5-5, Claude Code, web search and page fetching) gets the census common brief
(census/runs/verify/stage3_review/BRIEFS.md, "Common brief (groups A to D)"), unchanged except for the input and output
paths and the group description, plus one added rule:

  DATE RULE. Count only sources whose first public version (first arXiv version, or first publication if there is no
  preprint) is dated before {CUTOFF}. A source dated on or after that day does not count toward any label; list it
  under citation_problems with its date.

{CUTOFF} is the target paper's v1 date from census/controls/planted/N_truth.json. The agent is not told why the rule
exists, what the card's origin is, or anything about planted behaviors. Input: stage3_review/inputs/<Sxx>.md (the
Stage 3 (a) card and reply, card named by a neutral id). Output: stage3_review/outputs/<Sxx>.json with the census
keys. final_labels.json maps class_id -> {"final_label": ...} for the report.
