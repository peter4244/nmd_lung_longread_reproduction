# Unified significance-asterisk convention for the NMD paper — R side.
# Mirrors `figures/lib/p_to_stars.py`. Single source of truth across the
# project. Source this from any panel R script or verify_pass2 helper that
# needs the p → asterisk conversion.
#
# Convention (4 tiers):
#
#     ****  p < 10⁻¹⁰
#     ***   p < 10⁻⁴
#     **    p < 10⁻³
#     *     p < 0.05
#     n.s.  otherwise
#
# Rationale: see `p_to_stars.py` docstring.

# Canonical legend string. Paste verbatim into any legend that includes the
# asterisk-threshold key.
STARS_LEGEND <- "Stars: ****p < 10⁻¹⁰, ***p < 10⁻⁴, **p < 10⁻³, *p < 0.05, n.s. otherwise."

# Threshold ladder (descending p).
STARS_THRESHOLDS <- list(
  c(1e-10, "****"),
  c(1e-4,  "***"),
  c(1e-3,  "**"),
  c(0.05,  "*")
)

p_to_stars <- function(p, ns = "n.s.") {
  if (is.null(p) || length(p) == 0L) return(ns)
  if (is.na(p) || p < 0) return(ns)
  for (t in STARS_THRESHOLDS) {
    if (p < as.numeric(t[1])) return(t[2])
  }
  ns
}
