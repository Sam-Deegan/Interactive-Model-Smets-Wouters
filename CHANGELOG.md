# Changelog

All notable changes to this app. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for a change to the model or its notation, MINOR for new features
(a stage, a worked example, a figure), PATCH for fixes and wording.
Each release is tagged in git as `vX.Y.Z` and shown in the app footer.

## [1.0.4] - 2026-09-28

### App
- Cards, panels, tiles and buttons are square with no shadow: they organise
  the page rather than decorate it.

## [1.0.3] - 2026-09-28

### App
- The QR code returns to the foot of the sidebar, with the name and site
  address, alongside the small one in the title bar.

## [1.0.2] - 2026-09-28

### App
- No figure carries a title or subtitle inside the image; the card header
  and the caption under it name and explain the figure (CONVENTIONS.md 6).
- Figures are drawn on a white ground, so the image sits flat in its card
  instead of showing as a tinted tile.

## [1.0.1] - 2026-09-28

### App
- The In Words tab lays out its three columns at fixed widths, so an
  equation no longer collapses to one term per line beside its note.
- The preset card no longer doubles the word "Stage" in front of a stage
  name that already carries it.

## [1.0.0] - 2026-09-28

First public release as a standalone repository.

### Model
- Smets and Wouters (2007), eqs (1) to (14) with the flexible-price
  economy and seven shocks, in Sims's canonical form, solved in base R by
  a shifted pencil and a real stable basis (no QZ), at the authors'
  posterior mode.
- Two models from one solver: the lecture's version (no trend growth,
  utility separable in hours, `k_2` as eq (8) prints it) and Smets and
  Wouters' own as their code writes it, checked against Dynare to 1e-10.
- The risk premium in the paper's units, with a positive shock a rise in
  the premium; impulse responses at one posterior standard deviation;
  forecast error variance decompositions weighted by each shock's sigma.
- The frictions of Whelan part 11, slide 10, switched off one at a time
  at the deep parameter, together, and with prices and wages flexible.

### App
- Five stages: the frictions, the model assembled, impulse responses,
  variance decompositions, and ours against Smets and Wouters' own.
- Nine worked examples whose numbers are read off the solver at start-up.
- Seven sliders with typed boxes; the other estimated parameters fixed at
  the mode.
- Equations, Notation and In Words tabs in the module's letters, with
  Smets and Wouters' letters and every equation that differs marked at
  stage 5; a What the Lecture Simplifies tab; a Diagnostics tab.
- A scorecard against six results Whelan states in class, at the estimate
  and at the sliders.
- Save PNG and Save PDF buttons under every figure, at the deck's size.
- A verification script with gensys and Dynare fixtures.
