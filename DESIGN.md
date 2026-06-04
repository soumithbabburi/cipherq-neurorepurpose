---
version: alpha
name: NeuroRepurpose · CipherQ
description: >
  A clinical drug-repurposing intelligence platform surfacing compound
  candidates from ChEMBL data, running PBPK pharmacokinetic simulations,
  and coordinating molecular docking. Designed for computational biologists
  and clinical researchers who need data-dense, scientifically credible tooling
  that feels as precise as the science it supports.

colors:
  # ── Brand ──────────────────────────────────────────────────────────────────
  primary: "#7c3aed"
  secondary: "#4f46e5"
  accent: "#0284c7"
  accent-warm: "#ec4899"

  # ── Background hierarchy ────────────────────────────────────────────────────
  background: "#f8fafc"
  surface: "#f1f5f9"
  card: "#ffffff"
  raised: "#e2e8f0"

  # ── Borders ─────────────────────────────────────────────────────────────────
  border: "#e2e8f0"
  border-strong: "#cbd5e1"

  # ── Text hierarchy ───────────────────────────────────────────────────────────
  on-background: "#0f172a"
  on-surface: "#1e293b"
  muted: "#475569"
  dim: "#94a3b8"

  # ── Semantic / clinical ──────────────────────────────────────────────────────
  success: "#059669"
  warning: "#d97706"
  error: "#e11d48"

typography:
  # ── Display & Hero ───────────────────────────────────────────────────────────
  display:
    fontFamily: Inter
    fontSize: 30px
    fontWeight: "900"
    lineHeight: "1.15"
    letterSpacing: -0.048em

  # ── Headlines ────────────────────────────────────────────────────────────────
  headline-lg:
    fontFamily: Inter
    fontSize: 23px
    fontWeight: "800"
    lineHeight: "1.2"
    letterSpacing: -0.04em

  headline-md:
    fontFamily: Inter
    fontSize: 21px
    fontWeight: "800"
    lineHeight: "1.25"
    letterSpacing: -0.035em

  headline-sm:
    fontFamily: Inter
    fontSize: 15px
    fontWeight: "800"
    lineHeight: "1.1"
    letterSpacing: -0.03em

  # ── Stat values (large data numerals) ────────────────────────────────────────
  stat-value:
    fontFamily: Inter
    fontSize: 24px
    fontWeight: "800"
    lineHeight: "1"
    letterSpacing: -0.04em

  score-value:
    fontFamily: Inter
    fontSize: 34px
    fontWeight: "900"
    lineHeight: "1"
    letterSpacing: -0.05em

  # ── Body ─────────────────────────────────────────────────────────────────────
  body-lg:
    fontFamily: Inter
    fontSize: 13px
    fontWeight: "400"
    lineHeight: "1.75"

  body-md:
    fontFamily: Inter
    fontSize: 12px
    fontWeight: "400"
    lineHeight: "1.6"

  body-sm:
    fontFamily: Inter
    fontSize: 12px
    fontWeight: "500"
    lineHeight: "1.6"

  # ── UI chrome ─────────────────────────────────────────────────────────────────
  ui:
    fontFamily: Inter
    fontSize: 12px
    fontWeight: "500"
    lineHeight: "1.6"

  ui-sm:
    fontFamily: Inter
    fontSize: 11px
    fontWeight: "500"
    lineHeight: "1.6"

  # ── Labels (always uppercase in context) ──────────────────────────────────────
  label:
    fontFamily: Inter
    fontSize: 11px
    fontWeight: "600"
    lineHeight: "1.2"
    letterSpacing: 0.08em

  label-sm:
    fontFamily: Inter
    fontSize: 10px
    fontWeight: "600"
    lineHeight: "1"
    letterSpacing: 0.1em

  caption:
    fontFamily: Inter
    fontSize: 9px
    fontWeight: "500"
    lineHeight: "1"
    letterSpacing: 0.11em

  # ── Monospace (identifiers, data, SMILES, PK values) ─────────────────────────
  mono:
    fontFamily: JetBrains Mono
    fontSize: 11px
    fontWeight: "400"
    lineHeight: "1.5"
    letterSpacing: 0em

  mono-sm:
    fontFamily: JetBrains Mono
    fontSize: 10px
    fontWeight: "400"
    lineHeight: "1.5"
    letterSpacing: 0em

rounded:
  sm: 8px
  md: 12px
  lg: 18px
  full: 9999px

spacing:
  xs: 4px
  sm: 8px
  md: 16px
  lg: 24px
  xl: 32px
  2xl: 40px
  3xl: 64px
  card-padding-v: 20px
  card-padding-h: 24px
  main-padding-v: 28px
  main-padding-h: 35px
  sidebar-width: 250px

components:
  # ── Buttons ───────────────────────────────────────────────────────────────────
  button-primary:
    backgroundColor: "{colors.primary}"
    textColor: "#ffffff"
    rounded: "{rounded.sm}"
    padding: 10px
    typography: "{typography.body-sm}"

  button-primary-hover:
    backgroundColor: "{colors.primary}"
    textColor: "#ffffff"

  button-secondary:
    backgroundColor: "{colors.raised}"
    textColor: "{colors.on-surface}"
    rounded: "{rounded.sm}"
    padding: 10px
    typography: "{typography.body-sm}"

  button-secondary-hover:
    textColor: "{colors.primary}"
    borderColor: "{colors.primary}"

  button-ghost:
    backgroundColor: "transparent"
    textColor: "{colors.muted}"
    rounded: "{rounded.sm}"
    padding: 6px

  button-ghost-hover:
    textColor: "{colors.primary}"
    borderColor: "{colors.primary}"

  # ── Cards ─────────────────────────────────────────────────────────────────────
  card:
    backgroundColor: "{colors.card}"
    rounded: "{rounded.md}"
    padding: 20px

  card-hover:
    borderColor: "{colors.border-strong}"

  stat-card:
    backgroundColor: "{colors.card}"
    rounded: "{rounded.md}"
    padding: 20px

  compound-card:
    backgroundColor: "{colors.card}"
    rounded: "{rounded.md}"
    padding: 16px

  compound-card-hover:
    borderColor: "{colors.primary}"
    backgroundColor: "{colors.background}"

  # ── Form controls ─────────────────────────────────────────────────────────────
  input:
    backgroundColor: "{colors.raised}"
    textColor: "{colors.on-background}"
    rounded: "{rounded.sm}"
    padding: 9px
    typography: "{typography.body-md}"

  input-focus:
    borderColor: "{colors.primary}"

  # ── Navigation ────────────────────────────────────────────────────────────────
  sidebar-nav-item:
    backgroundColor: "transparent"
    textColor: "{colors.muted}"
    rounded: "{rounded.sm}"
    padding: 9px
    typography: "{typography.ui}"

  sidebar-nav-item-active:
    backgroundColor: "#3d1a7812"
    textColor: "{colors.primary}"

  sidebar-nav-item-hover:
    backgroundColor: "{colors.raised}"
    textColor: "{colors.on-surface}"

  tab-btn:
    backgroundColor: "transparent"
    textColor: "{colors.muted}"
    typography: "{typography.ui}"

  tab-btn-active:
    textColor: "{colors.primary}"

  # ── Data & clinical labels ────────────────────────────────────────────────────
  phase-pill-fda:
    backgroundColor: "#f0fdf4"
    textColor: "{colors.success}"
    rounded: "{rounded.full}"

  phase-pill-phase3:
    backgroundColor: "#eff6ff"
    textColor: "{colors.accent}"
    rounded: "{rounded.full}"

  phase-pill-phase2:
    backgroundColor: "#f5f3ff"
    textColor: "{colors.primary}"
    rounded: "{rounded.full}"

  phase-pill-phase1:
    backgroundColor: "#fffbeb"
    textColor: "{colors.warning}"
    rounded: "{rounded.full}"

  phase-pill-preclinical:
    backgroundColor: "#f8fafc"
    textColor: "{colors.muted}"
    rounded: "{rounded.full}"

  chip-purple:
    backgroundColor: "#3d1a7808"
    textColor: "{colors.primary}"
    rounded: "{rounded.full}"

  chip-cyan:
    backgroundColor: "#0284c708"
    textColor: "{colors.accent}"
    rounded: "{rounded.full}"

  chip-success:
    backgroundColor: "#05966908"
    textColor: "{colors.success}"
    rounded: "{rounded.full}"

  # ── Alerts ────────────────────────────────────────────────────────────────────
  alert-success:
    backgroundColor: "#05966908"
    textColor: "#047857"

  alert-warning:
    backgroundColor: "#d9770608"
    textColor: "#b45309"

  alert-error:
    backgroundColor: "#e11d4807"
    textColor: "#be123c"

  alert-info:
    backgroundColor: "#0284c708"
    textColor: "{colors.on-surface}"
---

# NeuroRepurpose · CipherQ Design System

## Overview

NeuroRepurpose is a **Clinical Intelligence Platform** — an application where the design must earn the trust of scientists before it earns their attention. The aesthetic is **Professional Scientific SaaS**: clean, precise, and information-forward, with just enough personality to feel alive rather than sterile.

The emotional target is *quiet confidence*. A researcher opening this tool should feel the same way they feel opening a well-typeset journal article: everything is exactly where it should be, nothing is wasted, and the data is never dressed up beyond what it needs to be.

The brand uses a **violet-to-indigo gradient** as its single expressive element — applied to the logo wordmark, primary call-to-action buttons, and 2px header bars on focal cards. Everything else earns its place through function. The gradient signals intelligence and depth without shouting. It appears sparingly: one gradient per major surface, never tiled or used as a background fill.

The platform targets computational biologists, medicinal chemists, and clinical informaticists. These are users who have a high tolerance for density and a low tolerance for decoration. Data tables, score breakdowns, and PK curves are primary content — not buried features. The UI wraps them in the minimum chrome needed to navigate.

## Colors

The palette is built on a single expressive axis from **Ink Slate** (page background) up through pure white (card surfaces), with **Violescent** as the only saturated hue used for interactive elements.

- **Primary — Violescent (#7c3aed):** The singular action color. Used exclusively for interactive affordances: primary buttons, active navigation state, tab indicators, focus rings, and the gradient terminus. It never appears as a background fill on more than one element per viewport.

- **Secondary — Indigo (#4f46e5):** The gradient pair to Violescent. Used only in tandem with Primary as the cool end of the `135deg` diagonal gradient. Also used in chart color scales for secondary data series.

- **Accent — Science Blue (#0284c7):** The informational blue. Used for Phase 3 clinical labels, hyperlinks, the sidebar logo `NeuroRepurpose` wordmark in the Dash view, and `chip-cyan` background. It reads as "information" rather than "action."

- **Background hierarchy (light surface stack):**
  - `#f8fafc` — Page canvas. Off-white with a hint of slate. Never pure white, which would feel stark.
  - `#f1f5f9` — Sidebar and secondary surfaces. Slightly more saturated to create a visible but gentle separation.
  - `#ffffff` — Card foreground. Pure white surfaces float above the background.
  - `#e2e8f0` — Raised elements: table headers, input backgrounds, progress track fills, separators.

- **Borders:** `#e2e8f0` (default, nearly invisible at rest) and `#cbd5e1` (strong, used on hover and for high-emphasis outlines). The difference between resting and hover border color is the primary depth mechanism.

- **Text hierarchy:** Four steps from `#0f172a` (near-black, primary copy) through `#1e293b` (secondary), `#475569` (muted labels and metadata), to `#94a3b8` (dim, placeholders and decorative text).

- **Semantic clinical coding:** Emerald (`#059669`) for FDA-approved and high-confidence scores. Amber (`#d97706`) for early-phase clinical data and caution states. Rose (`#e11d48`) for errors and terminated trials. These colors follow the same convention as traffic signals and are never used decoratively — only to communicate a clinical status.

## Typography

The system uses two typefaces. No exceptions.

- **Inter** handles all UI chrome, body copy, labels, and data display. It was chosen for its exceptional legibility at small sizes (down to 9px captions), its neutral personality that does not compete with data content, and its full weight range from 300 to 900. The heavy weights (800–900) are used for display numerals and hero text; the workhorse weight is 500 for interactive elements and 400 for body copy.

- **JetBrains Mono** handles all structured data: ChEMBL identifiers (`CHEMBL25`), SMILES strings, pharmacokinetic values, coordinate data, and any reference code. Its proportional character widths prevent identifier strings from expanding or collapsing unpredictably in table cells.

The base font size is **14px** (set on `:root`), making `1rem = 14px` throughout. This tighter base keeps data-dense views readable without requiring excessive vertical space.

**Scale intentions:**
- `display` (30px / 900 weight): Hero section title only. The heaviest typographic moment.
- `headline-lg` (23px / 800): Drug name in the analysis header. One per page.
- `stat-value` (24px / 800): Dashboard stat numerals. Their heavy weight carries the page hierarchy.
- `score-value` (34px / 900): The repurposing score in the compound analysis header. The single most important number on the analysis page.
- `body-md` through `body-lg`: Long-form content, descriptions, PubMed abstracts.
- `label` (11px / 600, uppercase): Section labels, table column headers, card titles. Always rendered `text-transform: uppercase` with `0.08em` letter-spacing.
- `mono`: ChEMBL IDs, SMILES, PK values. If text is a structured identifier or a measured value, it uses `mono`.

## Layout & Spacing

The layout follows a **Fixed Sidebar + Fluid Content** model. The sidebar (`250px`) is sticky and occupies the full viewport height. The main content area takes the remaining width and scrolls independently with `overflow-y: auto`.

The spacing scale is built on a **4px unit**, doubling through `xs → sm → md → lg`:

```
xs:  4px    (micro-gaps: between chips, phase pills, tag clusters)
sm:  8px    (intra-component: padding on nav items, small gaps)
md: 16px    (standard gap: between form fields, between card rows)
lg: 24px    (section rhythm: between cards, grid gutters)
xl: 32px    (major sections)
2xl: 40px   (main content top padding)
```

Content grid uses CSS Grid with `repeat(auto-fill, minmax(300px, 1fr))` for compound cards, and fixed `repeat(4, 1fr)` for the four stat cards on the dashboard. Stat cards collapse to a 2-column grid below 1100px and a single column below 900px.

Main content padding is `28px` top/bottom, `35px` left/right — generous enough to prevent the content from touching the viewport edge while leaving room for full-width tables to breathe.

The sidebar collapses to `display: none` below 640px (mobile), where the main content takes full width with `16px` padding.

## Elevation & Depth

This is a **border-first, shadow-minimal** design. Depth is communicated through three mechanisms, in order of prevalence:

1. **Surface layering:** Four distinct background tones create a perceptible stack. The background (`#f8fafc`) is the base, sidebar (`#f1f5f9`) sits slightly above it, cards (`#ffffff`) sit clearly above both. No explicit shadows are required to read this hierarchy.

2. **Border transitions on interaction:** Cards and inputs are defined by a `1px #e2e8f0` border at rest. On hover, the border steps to `#cbd5e1` — this single-step shift communicates interactivity without any shadow animation. Compound cards additionally reveal a `#7c3aed` border (the action color) on hover, signaling clickability.

3. **Gradient accent bars:** The drug analysis header and stat cards carry a `2px` gradient top-edge (`#7c3aed → #4f46e5`). This is not elevation — it is focus: the gradient bar marks the primary information surface on the page. Use it on exactly one card per focal section.

Shadows are reserved for hover states and elevated modals only: `0 4px 20px rgba(15, 23, 42, 0.10)` for cards on hover, `0 8px 28px rgba(124, 58, 237, 0.10)` for compound cards (tinted with the action color). These shadows are triggered by interaction, not by resting state — they communicate that something is about to be selected, not that something is important.

## Shapes

The corner radius vocabulary has three meaningful steps:

- **8px (`rounded.sm`):** Interactive elements — buttons, inputs, nav items, alert banners. The smallest radius that reads as "soft" rather than "sharp."
- **12px (`rounded.md`):** Content containers — cards, compound cards, cytoscape graph wrapper, score breakdown. The default container shape.
- **18px (`rounded.lg`):** Focal containers only — the hero section, the drug analysis header, and the knowledge graph canvas. The largest radius signals a "zone" rather than a "component."
- **9999px (`rounded.full`):** Pills and chips exclusively — phase labels (`FDA Approved`, `Phase 3`, etc.), semantic tags, and icon+text badge combinations. Full-pill rounding visually separates these small labels from the rectangular card world.

No component uses `0px` corners. The design does not employ sharp edges anywhere — the scientific context already carries enough gravity; no need to reinforce it architecturally.

## Components

### Primary Button
The gradient button is the strongest call-to-action signal in the system. Its background is `linear-gradient(135deg, #7c3aed, #4f46e5)` — the only element in the UI that uses both brand colors simultaneously. On hover, `opacity` drops to `0.88` and the element lifts `1px` on the Y axis (`transform: translateY(-1px)`). No color change on hover — the lift alone communicates affordance. Padding is `10px 22px` with `0.83rem` Inter 600 text.

### Secondary Button
Filled with the `raised` surface (`#e2e8f0`), `1px #cbd5e1` border, `#1e293b` text at rest. On hover, text and border color shift to `primary` violet — the border-and-text color change is the hover signal, not a background change. Use for adjacent actions to a primary button that should have similar visual weight but lower priority.

### Card
All cards share: `#ffffff` background, `1px #e2e8f0` border, `12px` radius, `20px` vertical / `24px` horizontal padding. Hover state shifts border to `#cbd5e1`. The `.card-title` class inside cards is always `0.75rem / 600 / uppercase / #475569` — a muted label that organizes the card content without claiming visual hierarchy.

### Stat Card
A card variant for the four dashboard metrics. Adds a `2px` gradient top bar (`var(--accent-gradient)`) and a faint radial glow at the top-right corner in the accent color at 6% opacity. The large numeral uses `stat-value` typography in the accent color. The label below uses `label-sm` in `#475569`.

### Compound Card
The primary discovery surface. Three zones stacked vertically: header (name + ChEMBL ID on left, score on right), tag row (phase pill + mechanism/target chips), and a `2×2` mini-bar grid showing the four score dimensions. On hover: reveals the gradient top-bar, border shifts to `rgba(124,58,237,0.3)`, card lifts `2px`, a subtle violet-tinted shadow appears. The score is color-coded: emerald (≥0.65), violet (≥0.40), sky-blue (≥0.20), muted (< 0.20).

### Sidebar Navigation
The sidebar carries the wordmark at top (gradient clip-text on Inter 800), a navigation list, and a DB connection indicator at bottom. Nav items: `8px` radius, `12px` horizontal padding, `#475569` muted text at rest. Active state: `rgba(124,58,237,0.08)` background, `#7c3aed` text, `2.5px` violet left-edge bar (`border-radius: 0 2px 2px 0`). The left bar is the primary active indicator — it must always be present on the active item.

### Analysis Tabs
A flush tab bar (`border-bottom: 1px #e2e8f0`) with flat tab buttons. Active tab: `#7c3aed` text, `2px #7c3aed` bottom border, `600` weight. Inactive: `#475569` muted, transparent. No background fill on active state — the underline alone marks selection. The tab bar scrolls horizontally without a visible scrollbar on overflow.

### Phase Pills
Color-coded badges that communicate clinical development stage at a glance. All use `9999px` radius and `0.65rem / 700 / 0.04em tracking` text. Color rules are strict and must not be used for other purposes:
- **FDA Approved:** Emerald background (`rgba(5,150,105,0.10)`) + `#059669` text
- **Phase 3:** Sky-blue background (`rgba(2,132,199,0.10)`) + `#0284c7` text
- **Phase 2:** Violet background (`rgba(124,58,237,0.10)`) + `#7c3aed` text
- **Phase 1:** Amber background (`rgba(217,119,6,0.10)`) + `#d97706` text
- **Preclinical:** Slate background (`rgba(100,116,139,0.10)`) + `#475569` text

### Chips (Semantic Tags)
Small full-radius badges for mechanism-of-action text, target gene names, and semantic metadata. Background is `8%` opacity of the accent color, border is `20%` opacity, text is the full accent color. Variants: `chip-purple` (mechanism), `chip-cyan` (target/link), `chip-success` (positive status), `chip-amber` (caution), `chip-rose` (error/termination).

### Data Tables
Tables communicate scientific provenance. Column headers: `0.66rem / 600 / uppercase / #7c3aed / #e2e8f0` background — using the action color (`primary`) for headers distinguishes them from body text and gives the table a strong spine. Data cells: `0.8rem / #1e293b`, `1px #e2e8f0` row dividers. Row hover: background shifts to `#e2e8f0`. Identifiers (ChEMBL IDs, MeSH IDs, PubMed IDs) always use `mono-sm`.

### Inputs & Form Controls
Background: `#e2e8f0` (the `raised` surface). Border: `1px #e2e8f0` at rest (nearly invisible, relying on background contrast). Focus: border shifts to `#7c3aed`, `0 0 0 3px rgba(124,58,237,0.08)` shadow ring. The 3px ring at 8% opacity is the focus indicator — it must remain visible but not aggressively colored. Placeholder text: `#94a3b8`. Select dropdowns use a custom SVG chevron in `#475569`.

### Alerts
Four semantic variants, each with `7–8%` opacity background tint and `25–30%` opacity border in the semantic color. Text uses a darker shade of the semantic color for contrast. Rounded `8px`. Use only to communicate actual status — not for decorative callout boxes.

### Spinner
`18px × 18px` circle, `2px` border in `#e2e8f0`, `border-top-color: #7c3aed`. 0.7s linear rotation. The only loading state affordance in the system — used inside button labels during async operations and as a standalone indicator in tab content panels during data fetch.

## Do's and Don'ts

**Do:**
- Use the violet-to-indigo gradient on exactly one button or header bar per primary content section. Gradient is the loudest design signal; it earns its place only at decision points.
- Render all structured identifiers (ChEMBL IDs, SMILES, MeSH IDs, PubMed IDs, PK numeric values) in `JetBrains Mono`. This is non-negotiable — it signals data precision and prevents confusing similar-looking characters in identifiers.
- Use the five-step text hierarchy (`on-background → on-surface → muted → dim → placeholder`) to establish reading order without changing font size. A label at `#475569` and body text at `#0f172a` are visually distinct even at the same `12px` size.
- Code clinical phase status exclusively through the defined phase-pill color scheme. Never assign a clinical status color to a non-clinical element.
- Allow charts and tables to breathe. Data visualizations (Plotly charts, compound grids) should have generous padding and never be squeezed. If space is tight, remove chrome before compressing data.
- Apply `text-transform: uppercase` and `0.08em+` letter-spacing to all `label` typography. This is the visual contract for category headers — violation makes the interface read as inconsistent.

**Don't:**
- Use `#ffffff` (pure white) as a page background. The slight warmth of `#f8fafc` prevents the stark clinical feeling of an all-white screen.
- Use the gradient as a background fill for panels, sidebars, or section separators. The gradient belongs only to CTAs and the 2px accent bar on focal cards.
- Introduce a third typeface. The Inter + JetBrains Mono dyad handles every typographic need in this system. Additional faces fragment the clinical authority the typography is meant to convey.
- Apply `box-shadow` to resting (non-hovered) card components. Resting shadows make all cards look equally important, destroying the depth hierarchy. Shadows activate on hover only.
- Invent new semantic colors for one-off situations. The success/warning/error triad is sufficient for all status communication. Use `muted` or `dim` for neutral states.
- Mix `rounded.full` (pill) shapes with `rounded.md` (card) shapes in the same component. Phase pills and chips live inside cards but are visually distinct objects — they should never share the same container corner radius.
- Use the primary violet for body text, decorative borders, or background fills. It must appear only on elements that invite interaction or mark active state.
