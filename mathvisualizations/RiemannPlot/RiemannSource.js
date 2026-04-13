// Utility functions shared across Riemann visualization tools.
// The main conformal map logic lives inline in riemann_plot.html.

function addarrays(a, b)      { return a.map((e, i) => e + b[i]); }
function subtractarrays(a, b) { return a.map((e, i) => e - b[i]); }
function arraymult(arr, m)    { return arr.map(e => e * m); }
