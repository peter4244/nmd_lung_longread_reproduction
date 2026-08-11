"""
Unified significance-asterisk convention for the NMD paper.

Single source of truth — import from here in every panel script and every
legend / README that mentions star thresholds. Mirror in `p_to_stars.R` for
R-side code (e.g., verify_pass2 helpers).

Convention (4 tiers):

    ****  p < 10⁻¹⁰
    ***   p < 10⁻⁴
    **    p < 10⁻³
    *     p < 0.05
    n.s.  otherwise

Rationale: p-values in this paper span 17+ orders of magnitude. The classical
3-tier 0.001 / 0.01 / 0.05 scheme saturates everything at ***. The 4-tier
scheme keeps the moderate tiers distinct while adding **** for the genuinely
extreme p-values (Fig 3 Panel C SE-prevalence ≤ 10⁻⁸¹, Fig 4 Panel C 5'UTR
length ≤ 10⁻¹⁷, Fig 5 Panel G uORF attention ≤ 10⁻²⁵, etc.).

Canonical legend string is exported as STARS_LEGEND so any legend that
includes the threshold key can reference it verbatim — same string everywhere.
"""

from __future__ import annotations

import math


# Unified pasteable legend string. Use this verbatim in every figure legend
# / README / methods text that needs to disclose the asterisk thresholds.
STARS_LEGEND = "Stars: *p < 0.05, **p < 10⁻³, ***p < 10⁻⁴, ****p < 10⁻¹⁰, n.s. otherwise."

# Threshold ladder (descending p). Exposed for code paths that need to test
# tier identity directly.
THRESHOLDS = (
    (1e-10, "****"),
    (1e-4,  "***"),
    (1e-3,  "**"),
    (0.05,  "*"),
)


def p_to_stars(p: float, ns: str = "n.s.") -> str:
    """Return the asterisk symbol for a p-value under the project convention.

    Parameters
    ----------
    p : float
        Two-sided p-value. NaN / None / negative inputs all yield `ns`.
    ns : str
        String to return when p is not significant (default "n.s." — matches
        the manuscript legend convention).
    """
    if p is None or (isinstance(p, float) and math.isnan(p)) or p < 0:
        return ns
    for thresh, mark in THRESHOLDS:
        if p < thresh:
            return mark
    return ns
