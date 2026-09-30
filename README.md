# Interactive Model: Smets-Wouters

A Shiny app for teaching the Smets and Wouters (2007) model of the US
business cycle: the estimated medium-scale DSGE model, solved at the
authors' posterior mode, with its impulse responses, its variance
decompositions and its frictions switched off one at a time. Built by
[Sam Deegan](https://sam-deegan.com) for ECON42240 Advanced
Macroeconomics, University College Dublin.

**Try it in the browser (nothing to install):**
https://sam-deegan.com/toy-models/smets-wouters/

Current version: **1.0.8** (see [CHANGELOG.md](CHANGELOG.md)). The version
is shown in the app footer; releases are tagged `vX.Y.Z`.

## What it does

The stage selector builds up the lecture's treatment of the model:

| Stage | What is added |
|---|---|
| 1: The Frictions, Switched Off | Output after a policy tightening with each friction (habit, adjustment costs, utilisation costs, indexation, sticky prices) removed in turn; hours after a technology shock, as estimated, with the real frictions off, and with everything off |
| 2: Our Smets-Wouters Model | The lecture's model assembled whole: every equation, one shock's responses, output's variance shares, and the scorecard |
| 3: Impulse Responses | One shock's four panels (output, hours, inflation, the rate), and the three demand shocks on one set of axes |
| 4: What Explains What | Each shock's share of the forecast error variance of output, inflation, the rate, consumption or investment, at four horizons |
| 5: Our Model Against Smets and Wouters | The lecture's model beside Smets and Wouters' own, at the same deep parameters: responses side by side, shares side by side, and every equation that differs marked |

Stages 1 to 4 run *our model*: Smets and Wouters' structure under the
assumptions the module has made so far (no trend growth, utility separable
in hours) and in its letters. Stage 5 sets it beside the model Smets and
Wouters estimated, which the solver keeps as their code writes it.

Each stage opens on a worked example (The Posterior Mode, Take Away Habit,
Let Investment Move, All Frictions Off, Our Model Assembled, A Policy
Shock, Make Prices Flexible, What Explains GDP, Ours Against Theirs).
Seven sliders move the parameters the lecture argues about; the other
estimated parameters stay at the posterior mode. Every slider has a box
beside it for an exact value. The Equations, Notation and In Words tabs
show the model as it stands at the chosen stage; What the Lecture
Simplifies lists the departures from the paper; Diagnostics counts the
Blanchard-Kahn roots and checks the solution on every solve. A scorecard
under the figures holds the model against six results Whelan states in
class, at the estimate and at the sliders. Every figure has a Save PNG
button.

## Run it locally

1. Install [R](https://cran.r-project.org/) (4.1 or later) and, ideally,
   [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install the three packages once:

   ```r
   install.packages(c("shiny", "bslib", "ggplot2"))
   ```

3. Open `app.R` in RStudio and click **Run App**, or from R in this folder:

   ```r
   shiny::runApp()
   ```

Equations are typeset with MathJax from a CDN, so they need an internet
connection; everything else runs offline. The solver uses base R only (no
QZ decomposition), so the app also runs under shinylive.

## Files

```
app.R          the app: settings and text (section B), figures (D),
               interface (E), server (F)
R/model.R      the model: calibration, canonical system, solver, impulse
               responses, variance decompositions, diagnostics, the
               frictions. Sources on its own, so slides can reuse it.
R/toolkit.R    layout and helpers shared with the other toy-model apps
tests/         the verification script, its fixtures and the Dynare files
www/           logo and QR code
README.md      this file
CHANGELOG.md   version history
CONVENTIONS.md how the figures and worked examples are laid out
LICENSE        CC BY-NC-ND 4.0
```

All text on screen (worked examples, guidance, equations, notation, the
caveat) is in sections `B_02` and `B_04` of `app.R`, so it can be edited
without touching the rest. The numbers the worked examples quote are read
off the solver at start-up, never typed.

## Checks

From the repo root:

```
Rscript tests/verify_against_sw2007.R
```

The script prints a pass or fail line per check and exits with a non-zero
status if any fails. It sources `tests/two_models.R`, which holds the
two-model checks (each side's displayed composites and equations,
evaluated on that side's own solution; the stage-5 pills against the rows
that differ between the two solvers). It needs shiny, bslib, ggplot2,
htmltools and grid.

What it compares against:

- `tests/sw2007_gensys_irf.csv`: impulse responses from an ordered-QZ
  gensys (scipy's `ordqz`) at the posterior means of Table 1A/1B; the
  shifted-pencil solver must match to 1e-8.
- `tests/dynare_mode_irf.csv`, `tests/dynare_mode_cvd.csv`: impulse
  responses and the conditional variance decomposition from Dynare 6.0
  running the authors' own model file (Pfeifer's DSGE_mod copy) at the
  values of their mode file, `usmodel_mode.mat`; the Smets-Wouters side
  must match to 1e-10. `tests/sw2007_mode_hessian.csv` carries those
  full-precision values. `tests/dynare/` holds the Dynare model files.
- The results Whelan states in class, read off the charts in part 11:
  what leads GDP, inflation and the fed funds rate at Q1 and Q100, the
  hump in output after a policy shock, hours falling after a technology
  shock, and what each friction does when switched off.
- The app itself: every displayed equation re-evaluated on the solver's
  path, every slider letter found in the notation, every figure at 3:2 or
  wider with no legend and each shock in its own colour, every Save
  button writing the file it promises, every picker guarded with `req()`,
  and every number a story quotes re-derived from the model.

## The model

Smets and Wouters (2007) is the benchmark estimated New Keynesian model:
the real business cycle core of a household, a firm and capital, with
sticky prices and wages (Calvo, with Kimball aggregators and indexation to
past inflation), habit in consumption, investment adjustment costs,
variable capacity utilisation, a fixed cost in production, a Taylor rule
with smoothing, and seven structural shocks, estimated by Bayesian methods
on seven US series from 1966 to 2004. It is the model of Whelan's *MA
Advanced Macroeconomics*, part 11, which presents the log-linearised
equations in composite form and the estimation tables, and it is the
model on which part 12's financial accelerator is built. Every variable is
a log deviation from steady state, in per cent.

The app writes the equations in the module's letters (`n_t` for hours,
`θ_p` and `θ_w` for the Calvo probabilities, `σ` for the elasticity
itself, `x_t` for the output gap, `φ_π`, `φ_x`, `φ_Δx` and `ρ_i` in the
rule, `v` for an innovation). Habit stays `λ`, `i_t` is investment and
`r_t` the policy rate, as in the paper and on Whelan's slides.

```
Production:   y_t = φ_p (α k^s_t + (1−α) n_t + ε^a_t)
Capital use:  k^s_t = k_{t−1} + z_t,   z_t = z_1 r^k_t
Rental rate:  r^k_t = −(k^s_t − n_t) + w_t
Resources:    y_t = c_y c_t + i_y i_t + z_y z_t + ε^g_t,   i_y = δ k_y
Consumption:  c_t = c_1 c_{t−1} + (1−c_1) E_t c_{t+1} − c_3 (r_t − E_t π_{t+1} + ε^b_t)
Investment:   i_t = i_1 i_{t−1} + (1−i_1) E_t i_{t+1} + i_2 q_t + ε^i_t
Tobin's q:    q_t = q_1 E_t q_{t+1} + (1−q_1) E_t r^k_{t+1} − (r_t − E_t π_{t+1} + ε^b_t)
Capital:      k_t = k_1 k_{t−1} + (1−k_1) i_t + k_2 ε^i_t
Price mark-up: μ^p_t = α (k^s_t − n_t) + ε^a_t − w_t
NKPC:         π_t = π_1 π_{t−1} + π_2 E_t π_{t+1} − π_3 μ^p_t + ε^p_t
Wage mark-up: μ^w_t = w_t − ( σ_l n_t + (c_t − λ c_{t−1}) / (σ (1−λ)) )
Wages:        w_t = w_1 w_{t−1} + (1−w_1)(E_t w_{t+1} + E_t π_{t+1}) − w_2 π_t
                    + w_3 π_{t−1} − w_4 μ^w_t + ε^w_t
Policy rule:  r_t = ρ_i r_{t−1} + (1−ρ_i)(φ_π π_t + φ_x x_t) + φ_Δx (x_t − x_{t−1}) + ε^r_t
Potential:    x_t ≡ y_t − y^p_t, y^p_t the same real economy with μ^p_t = μ^w_t = 0
Shocks:       ε^x_t = ρ_x ε^x_{t−1} + v^x_t (five AR(1); spending also takes ρ_ga v^a_t),
              ε^p_t and ε^w_t ARMA(1,1) with moving-average terms μ_p and μ_w
```

**The supply side** is Cobb-Douglas in capital in use and hours, scaled by
one plus the fixed-cost share `φ_p`. Capital in use is last period's
installed stock worked harder or softer, and running it harder costs
output through `z_1 = (1−ψ)/ψ`. The price mark-up `μ^p` is minus real
marginal cost, so the Phillips curve's `−π_3 μ^p_t` is the module's
real-marginal-cost term; indexation `ι_p` supplies the lagged inflation
term, and the Calvo probability `θ_p` flattens the curve through `π_3`.

**The demand side** has habit `λ` putting last period's consumption on the
right of the Euler equation, which is what makes the output response to a
policy shock hump-shaped; the risk premium `ε^b_t` sits inside the
real-rate bracket, so a rise in it lowers consumption and `q` as a rise in
the real rate would. Adjustment costs `φ` put investment's own lag in, so
new capital cannot come on line at once; the investment shock enters both
the investment equation and capital accumulation, with `k_2 = (1−k_1)/i_2`.
Wages are Calvo with indexation `ι_w`, and the wage mark-up is the gap
between the real wage and the marginal rate of substitution, whose
consumption term carries `1/σ` because utility is separable in hours.

**The policy rule** has inertia `ρ_i`, a response to inflation `φ_π`, and
responses to the output gap in level and in change. The gap is against
potential output, defined as output with prices and wages flexible and no
mark-up shocks, so the model carries a second, flexible-price economy.

**Solution.** The composites (`c_1`, `c_3`, `i_1`, `i_2`, `q_1`, `z_1`,
`k_1`, `k_2`, `π_1` to `π_3`, `w_1` to `w_4`) are functions of the deep
parameters in Table 1A; the slides give the composites only, so that
mapping is taken from the paper and shown under Thresholds and
Simplifications. The system is written in Sims's canonical form
`Γ_0 S_t = Γ_1 S_{t−1} + Ψ v_t + Π ζ_t`, with one variable per expectation
and one forecast error each, and solved for `S_t = T S_{t−1} + R v_t`
without a QZ decomposition: the pencil is shifted so `Γ_1 − s Γ_0` is
invertible, a real basis for the stable subspace is built from its
eigenvectors, and the left null space of `Γ_0 B` pins the forecast errors.
The count of unstable roots against forecast errors is the Blanchard-Kahn
condition, and the Diagnostics tab reports it. Impulse responses are
`T^h R e_x σ_x`, one posterior standard deviation of each shock, and the
variance shares weight each shock's squared responses by its own `σ_x^2`.

**What the five stages show with it**

- *1* The frictions. Habit moves the trough of the output response to a
  policy shock from quarter 2 back to impact; adjustment costs are the
  friction doing the most work, and with them near zero the trough is
  about three times as deep. Neither the real frictions nor flexible
  prices and wages alone flip the sign of hours after a technology shock;
  only both together do, which is as close as the model comes to the RBC.
- *2* The model assembled: every equation flagged new, a policy tightening
  taking output to its trough two quarters on, and the scorecard.
- *3* Impulse responses. After a policy tightening output and hours fall,
  inflation falls a little, and the rate is back to zero within seven
  quarters. The three demand shocks move output, hours, inflation and the
  rate the same way, as on Whelan's chart.
- *4* What explains what. In the quarter, output is a demand story
  (exogenous spending, then the risk premium); from ten quarters on the
  wage mark-up and productivity take over, and productivity reaches about
  a quarter of the long-run variance, far from the whole of it as in the
  RBC. The price mark-up leads inflation at Q1; monetary policy leads the
  fed funds rate at Q1.
- *5* Ours against theirs. With trend growth and hours in the Euler
  equation the output trough after a tightening is a little deeper and a
  quarter later, and the Q1 shares of spending and the risk premium are a
  few points higher; the Equations tab marks each equation that differs
  and shows ours beneath it.

**Where it departs from the textbook.** The lecture's model (stages 1 to
4) is Smets and Wouters under four simplifications, each listed on the
What the Lecture Simplifies tab and each visible at stage 5: no trend
growth (`γ = 1`, so `λ/γ`, `β̄γ` and `(1−δ)/γ` become `λ`, `β` and `1−δ`);
utility separable in hours (the paper's `c_2` term drops out and the
marginal rate of substitution carries `1/σ`); the module's letters; and a
potential output that relaxes the wage margin as well as prices, a wider
gap than Lecture 2.3's. Where the printed paper and the authors' code
differ, our model takes the form the derivation gives (`k_2 = (1−k_1)/i_2`,
the paper's eq (8)) and the Smets-Wouters side keeps the code's, because
the posterior was estimated under it. The risk-premium shock is written as
the paper prints eqs (2) and (4): a positive `ε^b_t` is a *rise* in the
premium, so consumption, `q` and output fall. That is the opposite sign to
the authors' code, and so to the chart Whelan reproduces on slide 18 (the
paper's Figure 2), where a positive shock raises output; the app draws
that chart as a fall in the premium and says so on the line and in the
caption, and the shock's standard deviation is Table 1B's `σ_b` divided by
`c_3`, since the code estimates `c_3 ε^b`. Whelan's slides slip in three
places, which the app notes in one line: installed capital where the
paper has capital in use (slides 3 and 7), no expectation on `r^k_{t+1}`
(slide 6), and the wage block's signs and letters (slide 8). The model
runs at the posterior mode, as the paper's variance decomposition does;
its impulse-response figures are posterior means, which no single
parameter vector reproduces, and the app does not tune the model to match
the published charts. It has no financial block (that is part 12), no open
economy and no zero lower bound.

## References

- Smets, F. and Wouters, R. (2007). Shocks and Frictions in US Business
  Cycles: A Bayesian DSGE Approach. *American Economic Review* 97(3),
  586-606.
- Whelan, K. *MA Advanced Macroeconomics*, part 11: The Smets-Wouters
  Model. https://www.karlwhelan.com/ma-advanced-macroeconomics/
- Galí, J. (1999). Technology, Employment, and the Business Cycle: Do
  Technology Shocks Explain Aggregate Fluctuations? *American Economic
  Review* 89(1).
- Sims, C. (2002). Solving Linear Rational Expectations Models.
  *Computational Economics* 20.

## Licence

© Sam Deegan. Released under
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/):
free to use and share for teaching with attribution; not for commercial use
or redistribution in modified form.
