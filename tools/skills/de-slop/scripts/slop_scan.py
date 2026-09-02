#!/usr/bin/env python3
"""slop_scan: mechanical pre-scan and post-check for the de-slop skill.

Counts the things a model miscounts by eye: dashes, colon rhythm, sentence
length bands, sentence openers, paragraph fragmentation, tricolon candidates,
binary-contrast shapes, tiered vocabulary hits, chatbot and tool artifacts,
markdown on plain-text surfaces. The compare mode checks a rewrite against its
original for fact drift, invented numbers, edit size, register migration, and
the flattened-rhythm failure.

The scan is a counter, not a judge. Every flag is a candidate for the human or
model pass that follows. A tricolon can be a real list of three. A colon can be
the right mark. Protected signature moves on a byline will trip flags at dose;
ignore those per the protect list.

Standard library only. Usage:

    python3 slop_scan.py draft.md
    python3 slop_scan.py draft.md --profile client-deliverable
    python3 slop_scan.py draft.md --surface email --json
    python3 slop_scan.py --compare original.md rewrite.md
    cat draft.txt | python3 slop_scan.py -

Profiles: auto (default), blog, linkedin, technical-blog, investor-email, docs,
casual, client-deliverable, research-paper, email. Surface overrides how
markdown and quotes are judged: markdown (default), plain (email, chat, text),
typeset (Word, PowerPoint: curly quotes are fine).
"""
from __future__ import annotations

import argparse
import difflib
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------
# Vocabulary. Regex fragments, matched case-insensitively on word boundaries.
# Sources: patterns.md tiers (upstream lineage), current-register.md (2026-09
# self-report), Liang et al. 2024 for the scientific-paper words.
# ---------------------------------------------------------------------------

TIER1 = [
    r"delve(?:s|d)?", r"deep[- ]dive", r"dive into", r"unpack(?:s|ed|ing)?",
    r"landscape", r"realm", r"tapestry", r"paradigm", r"embark(?:s|ed|ing)?",
    r"beacon", r"nestled", r"testament to", r"robust(?:ly|ness)?",
    r"comprehensive(?:ly)?", r"cutting[- ]edge", r"leverag(?:e|es|ed|ing)",
    r"utiliz(?:e|es|ed|ing|ation)", r"pivotal", r"underscor(?:e|es|ed|ing)",
    r"meticulous(?:ly)?", r"seamless(?:ly)?", r"game[- ]chang(?:er|ing)",
    r"hits? different(?:ly)?", r"watershed", r"the future looks bright",
    r"only time will tell", r"vibrant", r"thriving", r"bustling",
    r"showcas(?:e|es|ed|ing)", r"intricate", r"intricacies", r"ever[- ]evolving",
    r"daunting", r"holistic(?:ally)?", r"actionable", r"impactful", r"learnings",
    r"thought leader(?:ship)?", r"best practices", r"at its core",
    r"synerg(?:y|ies|istic)", r"interplay", r"in order to",
    r"due to the fact that", r"serves as", r"boasts?", r"commence(?:s|d)?",
    r"ascertain", r"endeavou?r", r"symphony", r"garner(?:s|ed|ing)?", r"adept",
    r"commendable", r"noteworthy", r"plays? a (?:crucial|key|vital|pivotal|central) role",
]

TIER2 = [
    r"harness(?:es|ed|ing)?", r"navigat(?:e|es|ing)", r"foster(?:s|ed|ing)?",
    r"elevat(?:e|es|ed|ing)", r"unleash(?:es|ed|ing)?", r"streamlin(?:e|es|ed|ing)",
    r"empower(?:s|ed|ing|ment)?", r"bolster(?:s|ed|ing)?", r"spearhead(?:s|ed|ing)?",
    r"resonat(?:e|es|ed|ing)", r"revolutioniz(?:e|es|ed|ing)",
    r"facilitat(?:e|es|ed|ing)", r"underpin(?:s|ned|ning|nings)?", r"nuanced?",
    r"crucial(?:ly)?", r"multifaceted", r"ecosystem", r"myriad", r"plethora",
    r"encompass(?:es|ed|ing)?", r"catalyz(?:e|es|ed|ing)", r"reimagin(?:e|es|ed|ing)",
    r"galvaniz(?:e|es|ed|ing)", r"augment(?:s|ed|ing)?", r"cultivat(?:e|es|ed|ing)",
    r"illuminat(?:e|es|ed|ing)", r"elucidat(?:e|es|ed|ing)", r"juxtapos(?:e|es|ed|ing)",
    r"transformative", r"transformation", r"paradigm[- ]shifting", r"cornerstone",
    r"paramount", r"poised to", r"burgeoning", r"nascent", r"quintessential",
    r"overarching",
]

TIER3 = [
    r"significant(?:ly)?", r"innovat(?:ive|ion|ions)", r"effective(?:ly)?",
    r"dynamics?", r"scalab(?:le|ility)", r"compelling", r"unprecedented",
    r"exceptional(?:ly)?", r"remarkabl[ey]", r"sophisticated", r"instrumental",
    r"world[- ]class", r"state[- ]of[- ]the[- ]art", r"best[- ]in[- ]class",
]

# The register newer models moved into once the Tier 1 words were trained out.
# Each is fine alone. Two or more distinct hits in a paragraph is the tell.
REGISTER = [
    r"load[- ]bearing", r"non[- ]?trivial(?:ly)?", r"orthogonal", r"first[- ]order",
    r"second[- ]order", r"at the margin", r"the delta", r"high[- ]signal",
    r"low[- ]signal", r"high[- ]leverage", r"meaningfully", r"materially",
    r"substantively", r"directionally", r"pressure[- ]test(?:s|ed|ing)?",
    r"stress[- ]test(?:s|ed|ing)?", r"sanity[- ]check(?:s|ed|ing)?", r"gut[- ]check",
    r"tease(?:s|d)? apart", r"disentangl(?:e|es|ed|ing)", r"the ask\b",
    r"the unlock\b", r"the lift\b", r"through[- ]line", r"north star",
    r"table stakes", r"long pole", r"forcing function", r"reach(?:es|ed|ing)? for",
    r"lean(?:s|ed|ing)? (?:on|into)", r"the move\b", r"footguns?", r"sharp edges",
    r"happy path", r"paper(?:s|ed|ing)? over", r"the shape of", r"the contours? of",
    r"the texture of", r"zoom(?:ing)? out", r"at a high level", r"in the weeds",
    r"failure modes?", r"degenerate case", r"opinionated", r"the whole (?:game|point|ballgame)",
    r"crisp(?:ly|ness)?", r"compounds?\b", r"compounding", r"flywheel",
    r"maps? (?:onto|cleanly)", r"mental model", r"reframe", r"the lens\b",
    r"at scale", r"bak(?:e|es|ed|ing) in", r"surface area", r"gated on",
    r"legib(?:le|ility)", r"durable", r"hygiene", r"path of least resistance",
    r"primitives?\b", r"scaffold(?:ing)?", r"the plumbing", r"does (?:real|actual|specific|the) work",
    r"earns? (?:its|their) (?:place|keep)", r"pulls? (?:its|their) weight",
    r"carr(?:y|ies) (?:its|the|real) weight", r"sits? in tension", r"there'?s a tension",
    r"holds? the tension", r"double[- ]edged", r"cuts both ways",
    r"there'?s a version of this", r"if anything\b", r"it turns out", r"turns out\b",
    r"there'?s a reason\b", r"it'?s no accident", r"this is where .{1,30} comes in",
    r"this is where it gets", r"some of this is", r"part of the answer",
    r"I'?d (?:argue|push back|gently push back|note|add|flag|offer)", r"fair enough",
    r"that'?s fair", r"quietly\b", r"the (?:honest|short|long|boring|uncomfortable|simple|real|hard) version",
]

# Announcing the move instead of making it. Tier 1: replace or cut.
ANNOUNCERS = [
    r"here'?s the thing", r"the thing is\b", r"let'?s (?:dive|explore|unpack|break|take a look|get into|start|begin)",
    r"without further ado", r"here'?s what you need to know", r"in today'?s (?:world|landscape|fast[- ]paced|digital)",
    r"it is important to note", r"it'?s (?:worth|important) (?:noting|to note|mentioning|flagging|calling out|saying)",
    r"worth (?:noting|flagging|calling out|mentioning)", r"it should be noted", r"needless to say",
    r"in conclusion", r"in summary", r"to sum up", r"all in all", r"at the end of the day",
    r"the bottom line", r"bottom line:", r"the upshot", r"the punchline", r"the (?:key )?takeaway",
    r"tl;?dr", r"which is to say", r"put simply", r"put differently", r"said differently",
    r"in other words", r"that is to say", r"to put it another way", r"simply put",
    r"^(?:translation|read|meaning|the catch|the kicker|the twist|the fix|the result|the problem|the upshot|the short version|concretely|in practice|specifically|net[- ]net):",
    r"think of it (?:as|this way)", r"picture this", r"imagine (?:a|this)",
]

CHATBOT = [
    r"I hope this helps", r"certainly!", r"(?:great|good|excellent|fantastic) question",
    r"feel free to (?:reach out|ask|contact)", r"let me know if (?:you|there|that|this)",
    r"(?:please )?don'?t hesitate", r"I'?d be happy to", r"happy to help",
    r"as an AI", r"as of my last (?:update|training)", r"as of my knowledge cutoff",
    r"while specific details are limited", r"you'?re absolutely right",
    r"that'?s a (?:really |very )?(?:insightful|great|good|excellent) (?:observation|point|question)",
    r"I hope this (?:email|message|note) finds you well", r"I hope you'?re (?:doing )?well",
    r"I hope you had a (?:great|good|wonderful)", r"I wanted to (?:reach out|touch base|follow up|circle back)",
    r"just (?:checking|touching base|circling back)", r"good instinct", r"right call\b",
    r"fair push\b", r"^(?:yes|no|short answer)[.:]", r"happy to (?:draft|expand|dig|adjust|revise)[^.]{0,40}if (?:useful|helpful|that helps)",
    r"if (?:useful|helpful), I can", r"want me to\b",
]

PERFORMED_HONESTY = [
    r"to be (?:honest|frank|candid|clear|fair)", r"honestly", r"the honest (?:answer|truth)",
    r"let me be clear", r"candidly", r"frankly", r"truth be told", r"the real truth",
    r"if I'?m being honest", r"I'?ll be honest", r"(?:in )?full transparency",
    r"for clarity", r"to be precise",
]

PADDING = [
    r"ensur(?:e|es|ed|ing)", r"conversely", r"in essence", r"essentially", r"fundamentally",
    r"moving forward", r"going forward", r"in terms of", r"the reality is",
    r"it goes without saying", r"ultimately",
]

CONFIDENCE_ADVERBS = [
    r"notably", r"interestingly", r"importantly", r"surprisingly", r"undoubtedly",
    r"crucially", r"remarkably", r"arguably", r"admittedly",
]

HEDGE_WORDS = [
    r"typically", r"often", r"generally", r"usually", r"potentially", r"perhaps",
    r"somewhat", r"relatively", r"fairly", r"largely", r"mostly", r"tends? to",
    r"may\b", r"might\b", r"could\b", r"in some cases", r"to some extent", r"arguably",
]

TAIL_QUALIFIERS = [
    r"at least for now", r"for what it'?s worth", r"as far as I can tell",
    r"in my experience", r"for the most part", r"more or less", r"at least in theory",
    r"or so I'?d argue", r"if you squint",
]

SINCERITY = [r"actual(?:ly)?", r"real(?:ly)?", r"genuine(?:ly)?", r"truly", r"literally"]

TOOL_ARTIFACTS = [
    r"oaicite", r"turn\d+(?:search|view|file)\d*", r"contentReference", r"utm_source=",
    r"【", r"】", r"\[citation needed\]", r"\{\{", r"\}\}", r"\bcite(?:turn|_)",
    r"<\|", r"\|>",
]

CONTRAST = [
    (r"\bnot (?:just|only|merely|simply) [^.;:!?]{1,80}?\b(?:but|it'?s|it is|rather)\b", "not just X, but Y"),
    (r"\b(?:isn'?t|is not|aren'?t|are not|wasn'?t|was not) (?:about|just|only|merely|simply) [^.;:!?]{1,60}[.;:!?]\s*(?:it'?s|it is|that'?s|they'?re|it was) (?:about )?", "isn't about X. It's about Y"),
    (r"\b(?:isn'?t|is not|aren'?t|are not|wasn'?t|was not) (?:a|an|the) [\w-]+(?: [\w-]+)? (?:problem|question|issue|matter|story|thing|challenge)\b", "isn't a X problem"),
    (r"\bless about [^.;:!?]{1,60}\bmore about\b", "less about X, more about Y"),
    (r"\bmore than (?:just )?(?:a|an) [\w-]+\b", "more than just a X"),
    (r"\bgoes beyond\b", "goes beyond X"),
    (r"\bstops? being [^.;:!?]{1,60}\bstarts? being\b", "stops being X and starts being Y"),
    (r"\bthe (?:question|answer|problem|issue|point|goal|challenge|risk|real [\w-]+) (?:isn'?t|is not|was never|is no longer)\b", "the question isn't X"),
    (r"\bnot because [^.;:!?]{1,80}\b(?:but )?because\b", "not because X, because Y"),
    (r"\bdoesn'?t mean [^.;:!?]{1,60}\b(?:it means|but)\b", "doesn't mean X, it means Y"),
    (r",\s*not (?:a|an|the|just|because|about)?\s?[\w-]+(?: [\w-]+){0,3}[.!?]", "Y, not X."),
    (r"\bit'?s (?:a|an) [\w-]+ (?:problem|question|issue), not (?:a|an)\b", "it's a X problem, not a Y problem"),
    (r"\bnot [^.;:!?]{1,50}?[,;]\s*(?:but|rather|instead)\b", "not X, but Y"),
]

FALSE_AGENCY = [
    r"\bthe (?:data|numbers|evidence|research|model|framework|process|system|culture|market|decision|strategy|roadmap|plan|document|rulebook) (?:tells? us|shows us|suggests?|ensures?|decides?|emerges?|rewards?|demands?|wants?|knows?|recognizes?|understands?|becomes?)\b",
]

LABEL_COLON = re.compile(
    r"^(?:the (?:catch|problem|result|upshot|fix|kicker|twist|issue|short version|takeaway|answer|reason|difference|good news|bad news)|translation|read|meaning|bottom line|net|concretely|in practice|specifically|put differently|said differently|in other words|the thing)\s*:",
    re.I,
)
WHICH_TAIL = re.compile(r",\s*which is (?:exactly |precisely )?(?:the |the whole )?(?:point|problem|issue|catch|reason|thing)\b", re.I)
RESTATERS = re.compile(r"\b(?:in other words|put simply|put differently|said differently|that is to say|which is to say|to put it another way|simply put|stated differently)\b", re.I)
ENUMERATOR_OPENER = re.compile(r"^(?:two|three|four|a few|a couple of|some|several) (?:things|notes|points|thoughts|observations|caveats|reactions)[.:]", re.I)
ONE_WORD_OPENER = re.compile(r"^(?:yes|no|correct|right|exactly|agreed|short answer)[.:!]", re.I)
IMPERATIVE_OPENER = re.compile(r"^(?:consider|take|look at|notice|think about|picture|imagine)\b", re.I)
CLOSER_PHRASES = re.compile(r"^(?:in conclusion|in summary|to sum up|ultimately|at the end of the day|the bottom line|all in all|in short|in the end)\b", re.I)

FURNITURE = re.compile(
    r"^(?:(?:subject|to|from|cc|bcc|re|fwd?|date)\s*:.*"
    r"|(?:hi|hello|dear|hey|good (?:morning|afternoon|evening))\b[^.!?]{0,40}[,:]?"
    r"|(?:thanks|thank you|many thanks|best|best regards|kind regards|regards|cheers|sincerely|talk soon|warmly|all the best)[,.!]?(?:\s+[A-Z][\w'’.-]*){0,3}"
    r"|(?:[A-Z][\w'’.-]*(?:\s+[A-Z][\w'’.-]*){0,3}),?)$", re.I)
OPENER_WORDS = {"the", "this", "it", "in"}
SUBORDINATORS = re.compile(r"\b(?:because|although|though|which|when|while|since|unless|whereas|until|after|before|if|so that|even if|as long as)\b", re.I)

ABBREVIATIONS = [
    "e.g.", "i.e.", "etc.", "vs.", "Dr.", "Mr.", "Mrs.", "Ms.", "Prof.", "Inc.", "Ltd.",
    "St.", "No.", "Fig.", "Eq.", "cf.", "al.", "U.S.", "U.K.", "a.m.", "p.m.", "approx.",
    "Jan.", "Feb.", "Mar.", "Apr.", "Jun.", "Jul.", "Aug.", "Sep.", "Sept.", "Oct.", "Nov.", "Dec.",
]

PROFILE_EXEMPT = {
    # words with legitimate technical meaning on these surfaces (patterns.md matrix)
    "technical-blog": {"robust", "comprehensive", "seamless", "ecosystem", "leverag", "facilitat", "underpin", "streamlin", "dynamics", "scalab", "primitives", "failure modes", "happy path", "footguns", "orthogonal"},
    "docs": {"robust", "comprehensive", "seamless", "ecosystem", "leverag", "facilitat", "underpin", "streamlin", "dynamics", "scalab", "primitives", "failure modes", "happy path", "footguns", "utiliz", "ensur"},
    "research-paper": {"robust", "comprehensive", "significant", "dynamics", "scalab", "facilitat", "underpin", "nuanced", "paradigm", "novel", "orthogonal", "non-trivial", "ensur", "may", "might", "could"},
}

# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------

BULLET_RE = re.compile(r"^\s*(?:[-*+•▪►]|\d+[.)]|[✅✔☑❌➡️→]|\(?[a-z]\))\s+")
HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")
FENCE_RE = re.compile(r"^\s*(```|~~~)")


def strip_code(lines: list[str]) -> tuple[list[str], int]:
    out, fenced, n = [], False, 0
    for ln in lines:
        if FENCE_RE.match(ln):
            fenced = not fenced
            n += 1
            out.append("")
            continue
        out.append("" if fenced else ln)
    return out, n // 2


def split_sentences(text: str) -> list[str]:
    t = text
    for a in ABBREVIATIONS:
        t = t.replace(a, a.replace(".", "\x00"))
    t = re.sub(r"(?<=\d)\.(?=\d)", "\x00", t)
    t = re.sub(r"\b([A-Z])\.(?=\s?[A-Z]\.)", lambda m: m.group(1) + "\x00", t)
    parts = re.split(r"(?<=[.!?][\"'”’)\]])\s+(?=[\"'“‘(\[]?[A-Z0-9])|(?<=[.!?])\s+(?=[\"'“‘(\[]?[A-Z0-9])", t)
    return [p.replace("\x00", ".").strip() for p in parts if p.strip()]


def word_count(s: str) -> int:
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


def first_word(s: str) -> str:
    m = re.match(r"[\"'“‘(\[*_]*([A-Za-z][A-Za-z'’-]*)", s)
    return m.group(1).lower() if m else ""


def segment(lines: list[str]) -> list[dict]:
    """Split lines into blocks: heading, list, table, prose. Each block carries its
    1-based start line, text, and (for prose) sentences with line numbers."""
    blocks, cur, start = [], [], None
    for i, ln in enumerate(lines, 1):
        if ln.strip() == "":
            if cur:
                blocks.append((start, cur))
                cur = []
            continue
        if not cur:
            start = i
        cur.append(ln)
    if cur:
        blocks.append((start, cur))

    out = []
    for start, blines in blocks:
        first = blines[0]
        if HEADING_RE.match(first) and len(blines) == 1:
            m = HEADING_RE.match(first)
            out.append({"kind": "heading", "line": start, "level": len(m.group(1)), "text": m.group(2).strip()})
            continue
        if all(ln.strip().startswith("|") for ln in blines):
            out.append({"kind": "table", "line": start, "text": "\n".join(blines)})
            continue
        if all(BULLET_RE.match(ln) or ln.startswith("  ") for ln in blines) and BULLET_RE.match(first):
            items = []
            for j, ln in enumerate(blines):
                if BULLET_RE.match(ln):
                    items.append({"line": start + j, "text": BULLET_RE.sub("", ln).strip(), "block": len(out)})
                elif items:
                    items[-1]["text"] += " " + ln.strip()
            out.append({"kind": "list", "line": start, "items": items})
            continue
        # prose (possibly with an inline heading line first)
        text_lines = list(blines)
        if HEADING_RE.match(first):
            m = HEADING_RE.match(first)
            out.append({"kind": "heading", "line": start, "level": len(m.group(1)), "text": m.group(2).strip()})
            text_lines = blines[1:]
            start += 1
            if not text_lines:
                continue
        joined = " ".join(s.strip() for s in text_lines)
        sents, pos = [], 0
        line_offsets = []
        acc = 0
        for ln in text_lines:
            line_offsets.append(acc)
            acc += len(ln.strip()) + 1
        for s in split_sentences(joined):
            idx = joined.find(s, pos)
            if idx < 0:
                idx = pos
            pos = idx + len(s)
            line_no = start + max(k for k, off in enumerate(line_offsets) if off <= idx)
            sents.append({"text": s, "line": line_no, "words": word_count(s), "first": first_word(s), "block": len(out)})
        out.append({"kind": "prose", "line": start, "text": joined, "sentences": sents})
    return out


# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------

def compile_list(frags: list[str]) -> list[tuple[re.Pattern, str]]:
    out = []
    for f in frags:
        if f.startswith("^"):
            pat = re.compile(f, re.I | re.M)
        else:
            pat = re.compile(r"(?<![\w-])(?:" + f + r")(?![\w-])", re.I)
        out.append((pat, f))
    return out


def find_all(patterns, lines: list[str], exempt: set[str] | None = None) -> list[dict]:
    hits = []
    for pat, frag in patterns:
        if exempt and any(frag.startswith(e) or e in frag for e in exempt):
            continue
        for i, ln in enumerate(lines, 1):
            for m in pat.finditer(ln):
                hits.append({"line": i, "match": m.group(0), "pattern": frag, "snippet": snippet(ln, m.start(), m.end())})
    hits.sort(key=lambda h: (h["line"], h["snippet"]))
    return hits


def snippet(line: str, a: int, b: int, width: int = 34) -> str:
    lo, hi = max(0, a - width), min(len(line), b + width)
    s = line[lo:hi].strip()
    return ("…" if lo > 0 else "") + s + ("…" if hi < len(line) else "")


def paragraph_of(blocks, line_no: int) -> int:
    idx = 0
    for k, b in enumerate(blocks):
        if b["line"] <= line_no:
            idx = k
    return idx


def detect_profile(text: str, words: int) -> str:
    low = text.lower()
    if words < 300 and re.search(r"#\w+", text):
        return "linkedin"
    if "```" in text or re.search(r"\b(?:def |import |function\(|SELECT |curl )", text):
        return "technical-blog"
    if re.search(r"^(?:hi|hello|dear|hey)\b", low, re.M) and re.search(r"\b(?:raise|round|investors?|runway|arr|mrr)\b", low):
        return "investor-email"
    if re.search(r"^(?:hi|hello|dear|hey) [A-Z]", text, re.M) or re.search(r"^(?:subject|to|from):", low, re.M):
        return "email"
    if re.search(r"\b(?:we (?:propose|show|present|report|find)|et al\.|abstract|related work|our (?:method|approach|results))\b", low):
        return "research-paper"
    if re.search(r"\b(?:steering committee|governance|deliverable|workstream|stakeholders?|engagement|pre-read|memo)\b", low):
        return "client-deliverable"
    if re.search(r"^\s*(?:step \d|\d+\.\s)", low, re.M) and re.search(r"\b(?:parameter|config|install|run)\b", low):
        return "docs"
    return "blog"


def scan(text: str, profile: str = "auto", surface: str = "markdown") -> dict:
    raw_lines = text.replace("\r\n", "\n").split("\n")
    lines, n_code = strip_code(raw_lines)
    body = "\n".join(lines)
    blocks = segment(lines)
    for b in blocks:
        if b["kind"] == "prose" and len(b["sentences"]) <= 1 and FURNITURE.match(b["text"].strip()) and word_count(b["text"]) <= 6:
            b["kind"] = "furniture"
    prose = [b for b in blocks if b["kind"] == "prose"]
    sents = [s for b in prose for s in b["sentences"]]
    list_sents = []
    for b in blocks:
        if b["kind"] == "list":
            for it in b["items"]:
                for t in split_sentences(re.sub(r"\*\*", "", it["text"])):
                    list_sents.append({"text": t, "line": it["line"], "words": word_count(t), "first": first_word(t), "block": it["block"]})
    all_sents = sents + list_sents
    total_words = sum(word_count(ln) for ln in lines)
    if profile == "auto":
        profile = detect_profile(body, total_words)
    if profile in ("email", "casual") and surface == "markdown":
        surface = "plain"
    exempt = PROFILE_EXEMPT.get(profile, set())

    F: list[dict] = []

    def flag(sev, kind, line, span, note=""):
        F.append({"severity": sev, "kind": kind, "line": line, "span": span, "note": note})

    # --- artifacts and chatbot voice (P0)
    for h in find_all(compile_list(TOOL_ARTIFACTS), lines):
        flag("P0", "tool artifact", h["line"], h["snippet"], "leaked machinery; strip")
    for h in find_all(compile_list(CHATBOT), lines):
        flag("P0", "chatbot artifact", h["line"], h["snippet"], "assistant voice; cut")

    # --- dashes
    for i, ln in enumerate(lines, 1):
        if re.fullmatch(r"\s*\|?[\s:|-]+\|?\s*", ln):
            continue  # markdown table separator row
        ln_nc = re.sub(r"`[^`]*`", lambda m: " " * len(m.group(0)), ln)  # ignore inline code
        for m in re.finditer(r"—|–|(?<=\w)--(?=\w)|\s--\s", ln_nc):
            flag("P1", "dash", i, snippet(ln, m.start(), m.end()), "em/en dash or double hyphen; recast with a period, comma, or parentheses")
        for m in re.finditer(r"(?<=[A-Za-z,;])\s-\s(?=[A-Za-z])", ln_nc):
            flag("P2", "spaced hyphen as dash", i, snippet(ln, m.start(), m.end()), "hyphen doing dash work")

    # --- binary contrast
    contrast_keys: list[tuple[int, str]] = []
    for pat, label in CONTRAST:
        rx = re.compile(pat, re.I)
        for i, ln in enumerate(lines, 1):
            for m in rx.finditer(ln):
                sev = "P0" if profile == "client-deliverable" else "P1"
                flag(sev, "binary contrast", i, snippet(ln, m.start(), m.end()), label)
                contrast_keys.append((i, m.group(0).lower()))
    # cross-sentence "X isn't Y. It's Z." shape, within one paragraph or list
    for k in range(len(all_sents) - 1):
        a, b = all_sents[k], all_sents[k + 1]
        if a.get("block") != b.get("block"):
            continue
        if re.search(r"\b(?:not|never|no longer)\b|n't\b", a["text"], re.I) and a["words"] <= 20 \
                and re.match(r"(?:It'?s|It is|That'?s|They'?re|They are|They have|They do|We do|You do|Instead|Rather|What (?:it|they|we) (?:is|are|have|do))\b", b["text"]) and b["words"] <= 14:
            pair = (a["text"] + " " + b["text"]).lower()
            if any(ln in (a["line"], b["line"]) and key in pair for ln, key in contrast_keys):
                continue
            sev = "P0" if profile == "client-deliverable" else "P1"
            flag(sev, "binary contrast", a["line"], a["text"] + " " + b["text"], "negation then short 'It's Y' reversal")

    # --- vocabulary tiers
    for h in find_all(compile_list(TIER1), lines, exempt):
        flag("P1", "tier-1 word", h["line"], h["match"], f"replace ({h['snippet']})")
    t2 = find_all(compile_list(TIER2), lines, exempt)
    by_para = Counter(paragraph_of(blocks, h["line"]) for h in t2)
    for h in t2:
        if by_para[paragraph_of(blocks, h["line"])] >= 2:
            flag("P1", "tier-2 cluster", h["line"], h["match"], "2+ tier-2 words in one paragraph")
        else:
            flag("info", "tier-2 word", h["line"], h["match"], "fine alone; a tell in clusters")
    t3 = find_all(compile_list(TIER3), lines, exempt)
    if total_words and len(t3) / max(total_words, 1) >= 0.03 and len(t3) >= 3:
        for h in t3:
            flag("P2", "tier-3 density", h["line"], h["match"], "leaning on intensifiers instead of specifics")

    # --- current register (2026-09): flag in clusters, report otherwise
    reg = find_all(compile_list(REGISTER), lines, exempt)
    reg_para = Counter(paragraph_of(blocks, h["line"]) for h in reg)
    for h in reg:
        if reg_para[paragraph_of(blocks, h["line"])] >= 2:
            flag("P1", "register cluster", h["line"], h["match"], "2+ current-model idioms in one paragraph; see current-register.md")
        else:
            flag("info", "register word", h["line"], h["match"], "current-model idiom; fine alone")

    # --- announcers, honesty markers, padding
    for h in find_all(compile_list(ANNOUNCERS), lines):
        flag("P1", "announcer", h["line"], h["snippet"], "announces the move instead of making it")
    for h in find_all(compile_list(PERFORMED_HONESTY), lines):
        flag("P1", "performed honesty", h["line"], h["snippet"], "delete the marker, keep the claim")
    pad = find_all(compile_list(PADDING), lines, exempt)
    ens = [h for h in pad if h["match"].lower().startswith("ensur")]
    for h in pad:
        if h in ens and len(ens) < 2:
            flag("info", "padding", h["line"], h["match"], "one 'ensure' is fine; two is a tell")
        else:
            flag("P2", "padding", h["line"], h["match"], "cut or commit")
    for h in find_all(compile_list(TAIL_QUALIFIERS), lines):
        flag("P2", "tail qualifier", h["line"], h["snippet"], "appended hedge; cut or move the hedge into the claim")
    for h in find_all(compile_list(FALSE_AGENCY), lines):
        flag("P2", "false agency", h["line"], h["snippet"], "name the human")

    # --- densities
    conf = find_all(compile_list(CONFIDENCE_ADVERBS), lines)
    if total_words and (len(conf) >= 3 and len(conf) / total_words * 500 >= 1.5):
        for h in conf:
            flag("P2", "confidence adverb density", h["line"], h["match"], "let the fact carry its weight")
    sinc = find_all(compile_list(SINCERITY), lines)
    sinc_rate = len(sinc) / total_words * 1000 if total_words else 0
    if len(sinc) >= 3 and sinc_rate >= 6:
        for h in sinc:
            flag("P2", "sincerity intensifier", h["line"], h["match"], f"{sinc_rate:.1f} per 1000 words; performed sincerity")
    hedge = find_all(compile_list(HEDGE_WORDS), lines, exempt)
    hedge_para = Counter(paragraph_of(blocks, h["line"]) for h in hedge)
    thresh = 6 if profile == "research-paper" else 4
    for pidx, n in hedge_para.items():
        if n >= thresh:
            flag("P2", "hedge stack", blocks[pidx]["line"], f"{n} hedges in one paragraph", "commit or quantify")

    # --- colons and semicolons
    colon_hits = []
    for i, ln in enumerate(lines, 1):
        if HEADING_RE.match(ln) or BULLET_RE.match(ln) or ln.strip().startswith("|"):
            continue
        stripped = re.sub(r"https?://\S+", "", ln)
        stripped = re.sub(r"\d:\d", "", stripped)
        stripped = re.sub(r"\*\*[^*]{1,40}:\*\*", "", stripped)
        for m in re.finditer(r"[:;]", stripped):
            at_end = m.end() >= len(stripped.rstrip())
            colon_hits.append({"line": i, "mark": m.group(0), "list_intro": at_end, "snippet": snippet(stripped, m.start(), m.end())})
    n_sent = len(sents)
    inline_marks = [c for c in colon_hits if not c["list_intro"]]
    rate = (len(inline_marks) / n_sent * 100) if n_sent else 0
    if n_sent >= 4 and len(inline_marks) >= 3 and rate >= 10:
        for c in inline_marks:
            flag("P1", "colon rhythm", c["line"], c["snippet"], f"{rate:.0f} colons/semicolons per 100 sentences; most become periods")
    elif n_sent >= 4 and len(inline_marks) >= 2 and rate >= 5:
        for c in inline_marks:
            flag("P2", "colon rhythm", c["line"], c["snippet"], f"{rate:.0f} per 100 sentences; check for displaced dash rhythm")
    for i, ln in enumerate(lines, 1):
        if LABEL_COLON.match(ln.strip()):
            flag("P1", "label colon", i, snippet(ln, 0, 30), "label-colon sentence; the dash rhythm moved into a label")
    for s in all_sents:
        if LABEL_COLON.match(s["text"]):
            flag("P1", "label colon", s["line"], s["text"][:60], "label-colon sentence")
        if WHICH_TAIL.search(s["text"]):
            flag("P1", "which-tail verdict", s["line"], s["text"][-70:], "trailing ', which is the point' verdict; cut or make it a sentence with a subject")
    rest = [(i, m.group(0)) for i, ln in enumerate(lines, 1) for m in RESTATERS.finditer(ln)]
    if len(rest) >= 2:
        for i, m in rest:
            flag("P1", "restatement chain", i, m, "saying it again more simply; keep the clearer one")
    elif len(rest) == 1:
        flag("P2", "restater", rest[0][0], rest[0][1], "one restatement; cut unless the first version was needed")

    # --- rhythm
    lengths = [s["words"] for s in sents]
    rhythm = {"sentences": n_sent, "mean": None, "median": None, "stdev": None, "cv": None,
              "band_17_23_longest_run": 0, "seesaw_runs": 0, "short_runs": 0}
    if n_sent >= 2:
        mean = statistics.mean(lengths)
        sd = statistics.pstdev(lengths)
        rhythm.update({"mean": round(mean, 1), "median": statistics.median(lengths), "stdev": round(sd, 1),
                       "cv": round(sd / mean, 2) if mean else None})
        run, best, best_end = 0, 0, 0
        for k, L in enumerate(lengths):
            run = run + 1 if 17 <= L <= 23 else 0
            if run > best:
                best, best_end = run, k
        rhythm["band_17_23_longest_run"] = best
        if best >= 3:
            s0 = sents[best_end - best + 1]
            flag("P1", "17-23 band run", s0["line"], f"{best} consecutive sentences of 17-23 words starting: {s0['text'][:60]}", "cadence uniformity; vary irregularly")
        # seesaw: alternating short (<=8) and long (>=25) for 4+
        alt = 0
        for k in range(1, n_sent):
            a, b = lengths[k - 1], lengths[k]
            if (a <= 8 and b >= 25) or (a >= 25 and b <= 8):
                alt += 1
                if alt >= 3:
                    rhythm["seesaw_runs"] += 1
                    flag("P1", "bimodal seesaw", sents[k]["line"], f"fragment/long alternation around: {sents[k]['text'][:50]}", "mechanical burstiness; two mediums, a fragment, a long one")
                    alt = 0
            else:
                alt = 0
        # runs of very short declaratives
        run = 0
        for k, L in enumerate(lengths):
            run = run + 1 if L <= 5 else 0
            if run == 3:
                rhythm["short_runs"] += 1
                flag("P1", "fragment triad", sents[k]["line"], " ".join(s["text"] for s in sents[k - 2:k + 1]), "three clipped sentences in a row; slogan rhythm")
        if n_sent >= 6 and rhythm["cv"] is not None and rhythm["cv"] < 0.35:
            flag("P2", "uniform cadence", sents[0]["line"], f"cv {rhythm['cv']} over {n_sent} sentences", "metronome; vary length irregularly")
        subs = sum(len(SUBORDINATORS.findall(s["text"])) for s in sents)
        rhythm["subordinators_per_sentence"] = round(subs / n_sent, 2)
        if n_sent >= 6 and mean < 11 and subs / n_sent < 0.35:
            flag("P2" if profile in ("linkedin", "casual") else "P1", "flattened rhythm", sents[0]["line"], f"mean {mean:.1f} words, {subs} subordinators in {n_sent} sentences", "telegraphic; connect related thoughts with because/which/although")

    # --- openers
    openers = Counter(s["first"] for s in sents if s["first"])
    rhythm["openers"] = openers.most_common(5)
    for b in prose:
        ss = b["sentences"]
        if len(ss) >= 3:
            n_open = sum(1 for s in ss if s["first"] in OPENER_WORDS)
            if n_open / len(ss) > 0.5:
                flag("P1", "opener repetition", b["line"], f"{n_open}/{len(ss)} sentences open with The/This/It/In", "recast openers")
        for k in range(2, len(ss)):
            if ss[k]["first"] and ss[k]["first"] == ss[k - 1]["first"] == ss[k - 2]["first"]:
                flag("P2", "same opener x3", ss[k]["line"], ss[k]["first"], "three consecutive sentences open the same way")
    if sents:
        s0 = sents[0]["text"]
        if ONE_WORD_OPENER.match(s0):
            flag("P2", "one-word opener", sents[0]["line"], s0[:40], "'Yes.' / 'Short answer:' directness tic")
        if ENUMERATOR_OPENER.match(s0):
            flag("P2", "enumerator opener", sents[0]["line"], s0[:40], "'Two things.' opener; start with the first thing")
    for b in prose:
        if b["sentences"] and IMPERATIVE_OPENER.match(b["sentences"][0]["text"]) and len(b["sentences"]) > 1:
            flag("P2", "imperative opener", b["line"], b["sentences"][0]["text"][:50], "'Consider X.' staging; state the point")

    # --- paragraphs
    para_sizes = [len(b["sentences"]) for b in prose]
    rhythm["paragraphs"] = len(prose)
    rhythm["sentences_per_paragraph"] = para_sizes
    if len(prose) >= 5:
        med = statistics.median(para_sizes)
        if med <= 2 and max(para_sizes) <= 3:
            flag("P1", "over-fragmented paragraphs", prose[0]["line"], f"{len(prose)} paragraphs, median {med} sentences, none over 3", "let an argument run to 5-8 sentences when it needs it")
        if max(para_sizes) - min(para_sizes) <= 1 and len(prose) >= 4:
            flag("P2", "uniform paragraphs", prose[0]["line"], f"all {len(prose)} paragraphs are {min(para_sizes)}-{max(para_sizes)} sentences", "paragraph-length metronome")
    if prose:
        last = prose[-1]
        ls = last["sentences"]
        if ls and CLOSER_PHRASES.match(ls[0]["text"]):
            flag("P1", "summary closer", last["line"], ls[0]["text"][:60], "cut the wrap-up; end on the last real point")
        hashtags_only = bool(re.fullmatch(r"(?:#\w+\s*)+", last["text"].strip()))
        signoff = bool(re.match(r"^(?:best|thanks|thank you|regards|cheers|sincerely|talk soon|warmly)\b", last["text"].strip(), re.I))
        if len(ls) == 1 and ls[0]["words"] <= 12 and len(prose) >= 2 and not hashtags_only and not signoff and profile not in ("email", "investor-email"):
            flag("P2", "aphoristic closer", last["line"], ls[0]["text"], "one-line verdict paragraph at the end; a common new-model tell")

    # --- tricolon candidates
    W = r"[\w'’-]+(?: [\w'’-]+)?"
    tri = re.compile(r"(?=(?<![\w-])(" + W + r"), (" + W + r"),? (?:and|or) (" + W + r")(?![\w-]))")
    asyn = re.compile(r"(?<![\w-])([\w-]+ [\w-]+), ([\w-]+ [\w-]+), ([\w-]+ [\w-]+)(?=[.!?;]|$)")
    for i, ln in enumerate(lines, 1):
        seen = set()
        for m in tri.finditer(ln):
            key = m.start(2)
            if key in seen:
                continue
            seen.add(key)
            sev = "P1" if profile == "client-deliverable" else "P2"
            flag(sev, "tricolon candidate", i, f"{m.group(1)}, {m.group(2)}, and {m.group(3)}", "drop test: could a member go unnoticed?")
        for m in asyn.finditer(ln):
            sev = "P1" if profile == "client-deliverable" else "P2"
            flag(sev, "tricolon candidate", i, m.group(0), "asyndetic triad; check the real count")

    # --- lists
    for b in blocks:
        if b["kind"] != "list" or len(b["items"]) < 3:
            continue
        items = [it["text"] for it in b["items"]]
        bold_lead = sum(1 for t in items if re.match(r"\*\*[^*]{1,40}:?\*\*:?", t))
        if bold_lead >= 3:
            flag("P1", "inline-header list", b["line"], items[0][:60], "bold label + colon on every bullet; write the point")
        firsts = Counter(first_word(re.sub(r"\*\*[^*]*\*\*:?\s*", "", t)) for t in items)
        common_first, n_first = firsts.most_common(1)[0]
        if common_first and n_first >= 3 and n_first / len(items) >= 0.75:
            flag("P1", "templated list", b["line"], f"{n_first}/{len(items)} items open with '{common_first}'", "same syntactic template; vary or write prose")
        for conn in ("enabling", "ensuring", "allowing", "which", "driving", "delivering", "helping", "so that"):
            if sum(1 for t in items if re.search(r"\b" + conn + r"\b", t, re.I)) >= 3:
                flag("P1", "templated list", b["line"], f"'{conn}' in 3+ items", "repeated connective; the template is doing the writing")
                break
        lens = [word_count(t) for t in items]
        if len(lens) >= 3 and max(lens) - min(lens) <= max(2, int(0.15 * max(lens))):
            flag("P2", "uniform list items", b["line"], f"item lengths {lens}", "every item the same size")
        if len(items) == 3 and all(word_count(t) <= 6 for t in items):
            flag("P2", "three-item slogan list", b["line"], " / ".join(items), "check whether three is the real count")

    # --- headings and formatting
    for b in blocks:
        if b["kind"] != "heading":
            continue
        words = [w for w in re.findall(r"[A-Za-z][\w'-]*", b["text"])]
        small = {"a", "an", "the", "of", "in", "on", "for", "and", "or", "to", "with", "at", "by", "from", "vs", "is"}
        caps = [w for w in words if w[0].isupper() and w.lower() not in small]
        if len(words) >= 4 and len(caps) >= 0.75 * len([w for w in words if w.lower() not in small]) and b["level"] >= 1:
            sev = "P1" if profile == "client-deliverable" else "P2"
            flag(sev, "title-case heading", b["line"], b["text"], "sentence case")
        if re.search(r"—|–|:", b["text"]):
            flag("P2", "punctuated heading", b["line"], b["text"], "dash or colon inside a heading")
        if re.search(r"[\U0001F300-\U0001FAFF\u2600-\u27BF]", b["text"]):
            flag("P1", "emoji heading", b["line"], b["text"], "remove")
    bolds = re.findall(r"\*\*[^*\n]+\*\*", body)
    if len(bolds) > max(1, total_words // 150) and profile not in ("linkedin", "casual"):
        flag("P2", "bold overuse", 0, f"{len(bolds)} bold spans in {total_words} words", "one per major section at most")
    emoji = [(i, m.group(0)) for i, ln in enumerate(lines, 1) for m in re.finditer(r"[\U0001F300-\U0001FAFF\u2600-\u27BF\u2B50\u2705]", ln)]
    if emoji and profile not in ("linkedin", "casual"):
        flag("P2", "emoji", emoji[0][0], f"{len(emoji)} emoji", "remove outside social surfaces")
    elif len(emoji) > 3:
        flag("P2", "emoji", emoji[0][0], f"{len(emoji)} emoji", "LinkedIn tolerates one or two, end of line")
    curly = len(re.findall(r"[“”‘’]", body))
    if curly and surface in ("plain", "markdown"):
        flag("P2", "curly quotes", 0, f"{curly} curly quotes", "straight quotes on plain-text and markdown surfaces; fine in Word or PowerPoint")
    if surface == "plain":
        # Only markdown syntax counts: ASCII headers, bold, backticks, and -/*/1. bullets.
        # Emoji bullets are convention on social surfaces and are judged by the emoji check.
        MD_BULLET = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s+")
        md = [(i, ln[:40]) for i, ln in enumerate(lines, 1)
              if HEADING_RE.match(ln) or "**" in ln or "`" in ln
              or (MD_BULLET.match(ln) and profile != "linkedin")]
        for i, s in md:
            flag("P0", "markdown on plain-text surface", i, s, "headers, bold, bullets in an email or message are an instant tell; prose only")

    # --- summary
    counts = Counter(f["severity"] for f in F)
    kinds = Counter((f["severity"], f["kind"]) for f in F)
    dash_count = sum(1 for f in F if f["kind"] == "dash")
    contrast_count = sum(1 for f in F if f["kind"] == "binary contrast")
    return {
        "file": None, "profile": profile, "surface": surface, "words": total_words,
        "sentences": n_sent, "paragraphs": len(prose), "code_blocks_skipped": n_code,
        "counts": {"P0": counts["P0"], "P1": counts["P1"], "P2": counts["P2"], "info": counts["info"]},
        "dashes": dash_count, "binary_contrasts": contrast_count,
        "colons_semicolons_inline": len(inline_marks), "colon_rate_per_100_sentences": round(rate, 1),
        "rhythm": rhythm,
        "kinds": {f"{sev} {kind}": n for (sev, kind), n in sorted(kinds.items())},
        "findings": sorted(F, key=lambda f: ({"P0": 0, "P1": 1, "P2": 2, "info": 3}[f["severity"]], f["line"])),
    }


# ---------------------------------------------------------------------------
# Compare mode
# ---------------------------------------------------------------------------

NUM_RE = re.compile(r"(?<![\w.])(?:\$|€|£)?\d[\d,]*(?:\.\d+)?%?(?:\s?(?:M|B|K|k|bn|mn|million|billion|thousand|x|days?|weeks?|months?|years?|sites?|hours?|minutes?|percent|pts?|bps))?(?![\w])")
DATE_RE = re.compile(r"\b(?:Q[1-4]\s?(?:FY)?\d{2,4}|FY\s?\d{2,4}|(?:19|20)\d{2}|(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec)[a-z]*\.? \d{1,2}(?:, \d{4})?|\d{1,2}/\d{1,2}/\d{2,4})\b")
QUALIFIERS = re.compile(r"\b(?:preliminary|provisional|draft|roughly|approximately|estimated|excluding|except|only|at least|at most|up to|no more than|may|might|could|likely|unlikely|not yet|pending|partial|subject to|interim|unaudited|self-reported|in principle|unless)\b", re.I)
ENTITY_RE = re.compile(r"\b(?:[A-Z][\w&'’-]+(?:\s+(?:of|for|and|&|de|the)\s+)?(?:\s+[A-Z][\w&'’-]+)+)\b")
WEEKDAYS_MONTHS = {"monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday", "january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december"}
COMMON_STARTS = {"The", "This", "It", "In", "A", "An", "We", "Our", "They", "These", "Those", "That", "There", "And", "But", "So", "If", "When", "For", "To", "As", "At", "By", "On", "Of", "With", "From"}


def numbers_in(text: str) -> Counter:
    out = Counter()
    for m in NUM_RE.finditer(text):
        tok = re.sub(r"\s+", " ", m.group(0).strip())
        out[tok] += 1
    return out


def core_numbers(nums: Counter) -> set[str]:
    return {re.sub(r"[^\d.]", "", n) for n in nums}


def entities_in(text: str) -> set[str]:
    ents = set()
    for para in re.split(r"\n\s*\n", text):
        for sent in split_sentences(" ".join(para.split())):
            for m in ENTITY_RE.finditer(sent):
                parts = m.group(0).strip().split()
                while parts and (parts[0] in COMMON_STARTS or parts[0].lower() in WEEKDAYS_MONTHS):
                    parts = parts[1:]
                if len(parts) >= 2:
                    ents.add(" ".join(parts))
    return ents


def compare(original: str, rewrite: str, profile: str = "auto", surface: str = "markdown") -> dict:
    a, b = scan(original, profile, surface), scan(rewrite, profile if profile != "auto" else a_profile_guard(original, profile), surface)
    o_nums, r_nums = numbers_in(original), numbers_in(rewrite)
    o_core, r_core = core_numbers(o_nums), core_numbers(r_nums)
    missing_numbers = sorted(n for n in o_nums if re.sub(r"[^\d.]", "", n) not in r_core)
    new_numbers = sorted(n for n in r_nums if re.sub(r"[^\d.]", "", n) not in o_core)
    o_dates, r_dates = set(DATE_RE.findall(original)), set(DATE_RE.findall(rewrite))
    o_ents, r_ents = entities_in(original), entities_in(rewrite)
    r_low = rewrite.lower()
    missing_entities = sorted(e for e in o_ents if e.lower() not in r_low and not all(w.lower() in r_low for w in e.split()))
    o_qual = Counter(q.lower() for q in QUALIFIERS.findall(original))
    r_qual = Counter(q.lower() for q in QUALIFIERS.findall(rewrite))
    dropped_qualifiers = sorted(q for q, n in o_qual.items() if r_qual[q] < n)
    ow, rw = original.split(), rewrite.split()
    sm = difflib.SequenceMatcher(a=ow, b=rw, autojunk=False)
    similarity = sm.ratio()
    o_reg = {f["span"].lower() for f in a["findings"] if f["kind"] in ("register word", "register cluster")}
    r_reg = {f["span"].lower() for f in b["findings"] if f["kind"] in ("register word", "register cluster")}
    swapped_in = sorted(r_reg - o_reg)
    o_mean, r_mean = a["rhythm"].get("mean") or 0, b["rhythm"].get("mean") or 0
    o_sub, r_sub = a["rhythm"].get("subordinators_per_sentence") or 0, b["rhythm"].get("subordinators_per_sentence") or 0
    flattened = bool(o_mean and r_mean and r_mean < 0.67 * o_mean) or bool(o_sub and r_sub < 0.5 * o_sub and b["sentences"] >= 6)
    new_closer = False
    r_prose = [ln for ln in rewrite.strip().split("\n\n") if ln.strip()]
    if len(r_prose) >= 2:
        last = r_prose[-1].strip()
        if FURNITURE.match(" ".join(last.split())):
            last = ""
        if last and word_count(last) <= 14 and difflib.SequenceMatcher(a=last.lower(), b=original.lower(), autojunk=False).find_longest_match(0, len(last), 0, len(original)).size < 0.6 * len(last):
            new_closer = True
    verdicts = []
    if missing_numbers:
        verdicts.append(f"P0 fact drift: numbers from the original are missing or changed: {', '.join(missing_numbers)}")
    if new_numbers:
        verdicts.append(f"P0 invented specificity: numbers not in the original: {', '.join(new_numbers)}")
    if o_dates - r_dates:
        verdicts.append(f"P0 dates missing: {', '.join(sorted(o_dates - r_dates))}")
    if missing_entities:
        verdicts.append(f"P1 named entities missing: {', '.join(missing_entities)}")
    if dropped_qualifiers:
        verdicts.append(f"P1 qualifiers dropped or reduced: {', '.join(dropped_qualifiers)} (check each still holds)")
    if similarity < 0.8:
        verdicts.append(f"P2 edit size: {100 - similarity * 100:.0f}% of words changed; over a fifth of a good draft needs a reason")
    if flattened:
        verdicts.append(f"P1 flattened rhythm: mean sentence length {o_mean} to {r_mean}, subordinators/sentence {o_sub} to {r_sub}")
    if swapped_in:
        verdicts.append(f"P1 register migration: rewrite introduced current-model idioms: {', '.join(swapped_in)}")
    if new_closer:
        verdicts.append("P2 new closer: the rewrite ends on a short paragraph the original did not have")
    for key in ("dashes", "binary_contrasts", "colons_semicolons_inline"):
        if b[key] > a[key]:
            verdicts.append(f"P1 rewrite added {key.replace('_', ' ')}: {a[key]} to {b[key]}")
    return {
        "original": {k: a[k] for k in ("words", "sentences", "paragraphs", "dashes", "binary_contrasts", "colons_semicolons_inline", "counts")},
        "rewrite": {k: b[k] for k in ("words", "sentences", "paragraphs", "dashes", "binary_contrasts", "colons_semicolons_inline", "counts")},
        "similarity": round(similarity, 3), "percent_changed": round(100 - similarity * 100, 1),
        "missing_numbers": missing_numbers, "new_numbers": new_numbers,
        "missing_dates": sorted(o_dates - r_dates), "missing_entities": missing_entities,
        "dropped_qualifiers": dropped_qualifiers, "register_swapped_in": swapped_in,
        "flattened_rhythm": flattened, "new_closer": new_closer,
        "mean_sentence_length": {"original": o_mean, "rewrite": r_mean},
        "verdicts": verdicts,
        "rewrite_scan": b,
    }


def a_profile_guard(text: str, profile: str) -> str:
    return detect_profile(text, word_count(text)) if profile == "auto" else profile


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def render(result: dict, show_info: bool = False) -> str:
    r = result
    out = [f"slop_scan  {r.get('file') or '-'}  profile={r['profile']} surface={r['surface']}  "
           f"{r['words']} words, {r['sentences']} sentences, {r['paragraphs']} paragraphs"]
    out.append(f"  P0 {r['counts']['P0']}   P1 {r['counts']['P1']}   P2 {r['counts']['P2']}   "
               f"dashes {r['dashes']}   contrasts {r['binary_contrasts']}   "
               f"colons/semicolons {r['colons_semicolons_inline']} ({r['colon_rate_per_100_sentences']}/100 sentences)")
    rh = r["rhythm"]
    if rh.get("mean") is not None:
        out.append(f"  rhythm: mean {rh['mean']} words, sd {rh['stdev']}, cv {rh['cv']}; longest 17-23 run {rh['band_17_23_longest_run']}; "
                   f"subordinators/sentence {rh.get('subordinators_per_sentence')}; paragraphs {rh.get('sentences_per_paragraph')}")
        out.append(f"  openers: {', '.join(f'{w} x{n}' for w, n in rh.get('openers', []))}")
    for sev in ("P0", "P1", "P2", "info"):
        items = [f for f in r["findings"] if f["severity"] == sev]
        if not items or (sev == "info" and not show_info):
            continue
        label = {"P0": "P0 credibility", "P1": "P1 obvious", "P2": "P2 polish", "info": "info (fine alone)"}[sev]
        out.append(f"\n{label}")
        for f in items:
            loc = f"L{f['line']}" if f["line"] else "doc"
            note = f"  ({f['note']})" if f["note"] else ""
            out.append(f"  - {f['kind']}  {loc}: {f['span']}{note}")
    if not any(f["severity"] in ("P0", "P1") for f in r["findings"]):
        out.append("\nNo P0 or P1 findings. Read it once more for the things a counter cannot see (noun-swap test, facts, voice), then deliver.")
    return "\n".join(out)


def render_compare(c: dict) -> str:
    o, r = c["original"], c["rewrite"]
    out = ["slop_scan compare",
           f"  original: {o['words']} words, {o['sentences']} sentences, {o['paragraphs']} paragraphs, dashes {o['dashes']}, contrasts {o['binary_contrasts']}, colons {o['colons_semicolons_inline']}, P0/P1/P2 {o['counts']['P0']}/{o['counts']['P1']}/{o['counts']['P2']}",
           f"  rewrite:  {r['words']} words, {r['sentences']} sentences, {r['paragraphs']} paragraphs, dashes {r['dashes']}, contrasts {r['binary_contrasts']}, colons {r['colons_semicolons_inline']}, P0/P1/P2 {r['counts']['P0']}/{r['counts']['P1']}/{r['counts']['P2']}",
           f"  words changed: {c['percent_changed']}%   mean sentence length {c['mean_sentence_length']['original']} to {c['mean_sentence_length']['rewrite']}"]
    if c["verdicts"]:
        out.append("\nchecks")
        for v in c["verdicts"]:
            out.append(f"  - {v}")
    else:
        out.append("\nchecks: no fact drift, no invented numbers, no dropped qualifiers, no register migration, rhythm intact.")
    out.append("\nremaining findings in the rewrite")
    rs = c["rewrite_scan"]
    rem = [f for f in rs["findings"] if f["severity"] in ("P0", "P1", "P2")]
    if rem:
        for f in rem:
            out.append(f"  - {f['severity']} {f['kind']}  L{f['line']}: {f['span']}")
    else:
        out.append("  none")
    return "\n".join(out)


def read_input(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    return Path(path).read_text(encoding="utf-8", errors="replace")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="de-slop mechanical scan")
    p.add_argument("file", nargs="?", help="text or markdown file, or - for stdin")
    p.add_argument("--compare", nargs=2, metavar=("ORIGINAL", "REWRITE"), help="check a rewrite against its original")
    p.add_argument("--profile", default="auto", help="auto, blog, linkedin, technical-blog, investor-email, docs, casual, client-deliverable, research-paper, email")
    p.add_argument("--surface", default="markdown", choices=["markdown", "plain", "typeset"], help="how to judge markdown and quotes")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--info", action="store_true", help="also list info-level hits (single tier-2 and register words)")
    args = p.parse_args(argv)

    if args.compare:
        o, r = (read_input(x) for x in args.compare)
        c = compare(o, r, args.profile, args.surface)
        if args.json:
            print(json.dumps(c, indent=2, ensure_ascii=False))
        else:
            print(render_compare(c))
        return 0
    if not args.file:
        p.print_help()
        return 2
    text = read_input(args.file)
    res = scan(text, args.profile, args.surface)
    res["file"] = args.file
    if args.json:
        print(json.dumps(res, indent=2, ensure_ascii=False))
    else:
        print(render(res, args.info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
