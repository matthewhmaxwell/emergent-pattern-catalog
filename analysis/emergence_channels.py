"""Back-compat shim — the estimators now live in epc.phase2a.info_channels (core
layer, so the emergence indicator can use them without an analysis→core dependency)."""
from epc.phase2a.info_channels import (  # noqa: F401
    _discretize,
    _mi_bits,
    _shannon_bits,
    micro_macro,
    mpr_complexity,
    psi_ce,
    psi_ce_best,
)
