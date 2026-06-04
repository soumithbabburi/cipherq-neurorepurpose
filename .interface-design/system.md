# Design System — Clinical Cream & Botanical

The design system for **NeuroRepurpose · CipherQ (RepurposeIQ)**, a clinical
drug-repurposing intelligence platform. Built with the `interface-design` skill.

## Direction

- **Personality:** Warm, credible, calm. "Opening a well-typeset monograph in a
  pharmacology library" — not a cold scientific console.
- **Foundation:** Warm cream (light). Surfaces are aged-paper, text is forest ink.
- **Depth:** Subtle soft shadows + warm hairline borders. Shadows do the lifting;
  borders only separate. One strategy, committed.
- **Why this direction:** The product gives an existing molecule a *second life*.
  Its world is the medicinal-chemistry bench — amber reagent glass, kraft lab
  notebooks, botanical natural-product roots, aged pharmacopeia paper. The palette
  comes FROM that world, not applied TO a template.

## Intent

- **Who:** computational biologist / pharma BD analyst, between ChEMBL, the
  literature, and a pitch deck. Deciding whether an existing drug deserves a
  second life against a new disease.
- **Must:** find credible candidates, understand the *why* (evidence chain),
  validate with PK/docking, build a business case. Trust the number.
- **Feel:** quiet confidence. Warm enough to sit with for an hour; precise enough
  to stake a drug program on.

## Tokens

### Surfaces (warm, lightness-stacked)
```
--bg-canvas:  #f6f5f0   soft cream — page AND sidebar (same, hairline divides)
--bg-sunken:  #ecebe3   inset: inputs, table headers, tracks (slightly darker)
--bg-card:    #fffdf9   paper card (barely-warm white, never stark #fff)
--bg-raised:  #ffffff   popovers/dropdowns — one level up, carried by shadow
```

### Borders (warm ink, low alpha — disappear until needed)
```
--border-subtle: rgba(31,42,36,0.06)
--border:        rgba(31,42,36,0.10)
--border-strong: rgba(31,42,36,0.16)
--ring:          rgba(47,111,91,0.28)   focus
```

### Text (forest ink, four honest levels)
```
--text-primary:   #1f2a24
--text-secondary: #45524a
--text-muted:     #6b776e
--text-faint:     #9aa39c
```

### Brand + accent (one accent, used with meaning)
```
--accent:       #2f6f5b   botanical green — actions, active nav, links/info
--accent-deep:  #245546   hover
--accent-soft:  rgba(47,111,91,0.10)
--gold:         #b8860b   warm gold — score seal, sparing emphasis only
--gold-soft:    rgba(184,134,11,0.12)
```

### Semantic — warm earthy clinical scale (phase coding)
```
phase-4 FDA:      #2f6f5b green       (approved = the brand color earns it)
phase-3:          #b5552f terracotta  (advanced, warm rust)
phase-2:          #97751a gold-olive
phase-1:          #8a6a4f warm brown
phase-0 preclin:  #6b776e stone
warning:          #b5552f clay
error:            #b3261e warm brick (never cold rose)
success:          #2f6f5b
```

### Score-seal ring color by band
```
>=0.65 #2f6f5b green · >=0.40 #b8860b gold · >=0.20 #a86b2d amber · <0.20 #9aa39c stone
```

## Typography
```
Display/Headlines: Newsreader (serif) — warm, editorial. weights 400/500/600 + italic
Body / UI:         Inter
Data / mono:       JetBrains Mono (tabular: ChEMBL IDs, SMILES, PK values)
```
- Headlines carry personality (serif). Labels are Inter small-caps,
  `0.08em` tracking. Any measured value or identifier is mono.

## Spacing & shape
```
base 4px · scale 4,8,12,16,24,32,48
radius: sm 8 (inputs/buttons) · md 12 (cards) · lg 16 (focal/modals) · pill 999
shadow-sm: 0 1px 2px rgba(31,42,36,.05), 0 2px 6px -2px rgba(31,42,36,.08)
shadow-md: 0 6px 24px -8px rgba(31,42,36,.14)
```

## Signature patterns

1. **Dossier card** — each compound/candidate is a page in a bound case file:
   left **evidence ribbon** (green vertical bar), serif drug name, mono ChEMBL id,
   a `▸ why` reasoning line, and a **confidence seal** on the right.
2. **Confidence seal** — circular SVG medallion; ring fills proportional to score,
   color by band (green/gold/amber/stone); value in the center set in Newsreader.
3. **Monograph section header** — serif title + thin rule + Inter small-caps kicker.
4. **Notebook-spine sidebar** — same cream as canvas, hairline border; active item
   carries a green left ribbon + soft green tint (not a different background world).

## Decisions

| Decision | Rationale | Date |
|----------|-----------|------|
| Warm cream foundation (light) | Product = second life of a molecule; world is the apothecary bench, not a console | 2026-06-03 |
| Botanical green single accent | Natural-product pharma roots; one accent used with meaning beats five | 2026-06-03 |
| Newsreader serif headlines | "Well-typeset monograph" feel; strongest warm signal, survives the swap test | 2026-06-03 |
| Subtle shadows (not borders-only) | Approachable depth for a tool you sit with; borders only separate | 2026-06-03 |
| Score = confidence seal, not % bar | Domain-specific; a seal reads as a verdict, a bare % reads as a KPI tile | 2026-06-03 |
| Sidebar = canvas color + hairline | Skill layering: don't fragment into sidebar-world vs content-world | 2026-06-03 |
