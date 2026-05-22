# EPC Sprint Orchestrator

Source of truth for the orchestrator deployed at /opt/epc-orchestrator on
ai-server (root@100.90.46.93). See SETUP.md for installation instructions.

Files in this directory:
- epc_orchestrator.py      — main script (Python, no deps outside stdlib + anthropic SDK)
- config.example.json      — config template (real values filled in at deploy time;
                              live config at /opt/epc-orchestrator/config.json on the VPS
                              is gitignored)
- epc-orchestrator.service — systemd unit
- SETUP.md                 — full installation guide

Sprint 36 deploy patches (vs the v1 bundle):
- Queue-priority brief reading (epc/phase2a/briefs/sprint_NN_brief.md)
- Chain-aware pre-flight pytest skip (EPC_SKIP_PREFLIGHT_PYTEST env)
- Dynamic poll interval (idle=300s, active-chain=60s)
- state.json atomic write + 5-backup rotation
- Pre-chain git tag (pre-chain-YYYYMMDD-HHMMSS) for rollback
- Disk-free gate (<10% blocks new chain link)
- Brief checksum logging
- Queued briefs trump escalate_if_chat_led (they are pre-approved)
