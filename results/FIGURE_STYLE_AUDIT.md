# Figure style audit

Audit date: 2026-11-13

## Palette

The manuscript figures use one six-color palette:

- `#66BC98`
- `#AAD09D`
- `#E3EA96`
- `#FDC089`
- `#E26844`
- `#8A233F`

White is used for panel backgrounds and separators. Neutral gray is limited to
axes, zero lines, text, and explicitly uncovered categories.

## Semantic mapping

- Tier 1 and focal evidence: `#8A233F`.
- Tier 2 and active or positive direction: `#E26844`.
- Tier 3, healthy, control, and negative direction: `#66BC98`.
- Intermediate or high-percentile categories: `#FDC089`.
- Low-to-moderate supporting categories: `#AAD09D`.
- Light categorical contrast: `#E3EA96`.

Continuous maps use green for negative values, a near-white center, and
peach-to-orange-to-burgundy for positive values. Coverage maps use a pale
background and green for available evidence. Concordance maps use pale peach to
burgundy.

Color is not the only encoding: line markers, marker shapes, point size,
position, labels, bar orientation, heatmap values, and panel grouping remain
available.

## Scope

The palette is defined centrally in `scripts/figure_palette.py` and used by the
scripts that generate Figure 1 through Figure 12. Figure 12 applies the same
palette to MR, regional coupling, SuSiE-colocalization, QTL coverage,
signature projection and virtual-knockdown panels.
