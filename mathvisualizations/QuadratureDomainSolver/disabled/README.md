# Disabled features

Files in this directory are **not loaded by `index.html`**. They contain
work-in-progress features that have been removed from the live app pending
further refinement.

## aqd/

Stage-0 through Stage-2 scaffolding for **Algebraic Quadrature Domains**
(weight ρ = |α|² with α = R' for rational R). The bounded non-singular
solver and accompanying UI were implemented but exposed edge cases
(univalence issues, Newton basin sensitivity for certain α-zero
configurations) that need design work before re-shipping.

When restoring the AQD tab:

1. Move `app/disabled/aqd/` back to `app/aqd/`.
2. In `app/index.html`, re-add the AQD `.tab-btn` and `#controls-aqd`
   panel, the three `<script>` tags (`aqd-common.js`, `aqd-bounded.js`,
   `aqd-ui.js`), and add `aqd: document.getElementById('controls-aqd')`
   back to the tab-switching `panels` map.
3. In `app/node-test.js`, re-enable the AQD loader block and the AQD
   test suite (search for `boundedAQD` to locate where they were).
