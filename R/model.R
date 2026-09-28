################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Smets-Wouters Model: Solver, Simulator and Decompositions                  ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Sourced automatically by app.R. Can be sourced alone from a lecture
##   .qmd so slide figures come from the same model:
##     source("R/model.R")
##
## Inputs:
##   None. Every function is a pure function of a parameter list.
##
## Outputs:
##   C_01_* parameters, the variable order, the composite-coefficient map
##   C_02_* the canonical system  G0 y_t = G1 y_(t-1) + Psi e_t + Pi eta_t
##   C_03_* the solver
##   C_04_* impulse responses, variance decompositions, simulation
##   C_05_* diagnostics: Blanchard-Kahn counting and structural residuals
##   C_06_* the frictions-off experiment of Whelan part 11, slide 10
##
## The model (Smets and Wouters 2007, eqs (1)-(14); Whelan part 11, slides
## 3-9 present the log-linearised version and slides 12-13 the posterior):
##   y_t   = phi_p (alpha k^s_t + (1 - alpha) l_t + eps^a_t)           (5)
##   k^s_t = k_(t-1) + z_t,   z_t = z_1 r^k_t                       (6, 7)
##   y_t   = c_y c_t + i_y i_t + z_y z_t + eps^g_t                       (1)
##   c_t   = c_1 c_(t-1) + (1 - c_1) E c_(t+1) + c_2 (l_t - E l_(t+1))
##           - c_3 (r_t - E pi_(t+1) + eps^b_t)                          (2)
##   i_t   = i_1 i_(t-1) + (1 - i_1) E i_(t+1) + i_2 q_t + eps^i_t       (3)
##   q_t   = q_1 E q_(t+1) + (1 - q_1) E r^k_(t+1)
##           - (r_t - E pi_(t+1) + eps^b_t)                              (4)
##   k_t   = k_1 k_(t-1) + (1 - k_1) i_t + k_2 eps^i_t                   (8)
##   mu^p_t = alpha (k^s_t - l_t) + eps^a_t - w_t                        (9)
##   pi_t  = pi_1 pi_(t-1) + pi_2 E pi_(t+1) - pi_3 mu^p_t + eps^p_t    (10)
##   r^k_t = -(k^s_t - l_t) + w_t                                       (11)
##   mu^w_t = w_t - (sigma_l l_t + (c_t - lambda/gamma c_(t-1))
##                                  / (1 - lambda/gamma))               (12)
##   w_t   = w_1 w_(t-1) + (1 - w_1)(E w_(t+1) + E pi_(t+1)) - w_2 pi_t
##           + w_3 pi_(t-1) - w_4 mu^w_t + eps^w_t                      (13)
##   r_t   = rho r_(t-1) + (1 - rho)(r_pi pi_t + r_y (y_t - y^p_t))
##           + r_dy [(y_t - y^p_t) - (y_(t-1) - y^p_(t-1))] + eps^r_t   (14)
##   plus the same real economy with flexible prices and wages for y^p_t,
##   five AR(1) shocks, two ARMA(1,1) mark-up shocks, and rho_ga loading the
##   productivity innovation onto spending. C_01_03 maps the deep
##   parameters of Table 1A/1B to the composites; the slides give the
##   composites only, so that mapping comes from the paper.
##
## Two models, one solver. par_lst$model_chr picks the model:
##   "sw"    Smets and Wouters (2007) as their code solves it, checked
##           against Dynare to 1e-10 (tests/verify_against_sw2007.R). Where
##           paper and code differ (k_2 in eq (8), capital in use in eq
##           (11), the risk premium's units) this side follows the code.
##   "ours"  the lecture's version: the same structure with no trend growth
##           (gamma = 1) and utility separable in hours (c_2 = 0, and the
##           marginal rate of substitution carries sigma_c), and k_2 in the
##           form eq (8) prints. Nothing is re-estimated; both run at the
##           posterior mode of Table 1A/1B.
##   The risk premium is written as the paper prints eqs (2) and (4): a
##   positive eps^b_t is a rise in the premium, so consumption, q and
##   output fall. The authors' code carries b_code = -c_3 eps^b_t, so the
##   Table 1B sigma_b is the s.d. of c_3 eps^b_t and C_04_00 divides by c_3.
##
## The solver. Sims's canonical form, solved without an ordered QZ
##   decomposition because base R has none and the app must run under
##   shinylive. The pencil is shifted so (G1 - s G0) is invertible, a real
##   basis for the stable subspace is built from eigen(), and the left null
##   space of G0 B pins the forecast errors; rank(L Pi) = n - ncol(B) is the
##   Blanchard-Kahn condition. Checked against a genuine ordered-QZ gensys
##   (V_04) and against Dynare (V_04b); C_05 re-checks every solve.
##
## Parameter list (par_lst) elements:
##   Table 1A: csadjcost, csigma, chabb, cprobw, csigl, cprobp, cindw, cindp,
##     czcap, cfc, crpi, crr, cry, crdy, constepinf, constebeta, constelab,
##     ctrend, calfa. Table 1B: crho*, cmap, cmaw, cgy, sd*. Calibrated:
##     delta, gy, lambdaw, curvp, curvw. model_chr.
##
## References:
##   Smets, F. and Wouters, R. (2007). Shocks and Frictions in US Business
##     Cycles: A Bayesian DSGE Approach. American Economic Review 97(3),
##     586-606.
##   Whelan, K. MA Advanced Macroeconomics, part 11 (The Smets-Wouters
##     Model); cited inline as [W11 nn], nn the slide.
##   Sims, C. (2002). Solving Linear Rational Expectations Models.
##     Computational Economics 20, for the canonical form.
##   Gali, J. (1999). American Economic Review 89(1), for hours after a
##     technology shock.

#-------------------------------- Script Begin --------------------------------#

################################################################################
## A: Table of Contents ########################################################
################################################################################
# Note: This file is section C; app.R holds sections B, D, E, F and G.
#
#   C: Model
#     C_01  Parameters and coefficients
#     C_02  The canonical system
#     C_03  The solver
#     C_04  Impulse responses, variance decompositions, simulation
#     C_05  Diagnostics
#     C_06  The frictions

################################################################################
## C: The Model ################################################################
################################################################################
# Note: Two calibration lists, the equation builder, the solver, and what is
#   computed from a solution.

#### C_01: Parameters and Coefficients #########################################
# Note: The posterior mode, the variable order, and the deep-to-composite
#   map.

###### C_01_01: Posterior Mode #################################################
# Note: The posterior mode, because the paper's variance decomposition
#   (Figure 1, [W11 15] to [W11 17]) is computed at the mode (p. 598); its
#   impulse responses (Figures 2, 3 and 6, [W11 18] to [W11 20]) are
#   posterior means, which no single parameter vector reproduces exactly.
#   The values are the authors' own mode file, usmodel_mode.mat, to eight
#   decimals; Table 1A/1B's mode column ([W11 12], [W11 13]) is this list
#   truncated to two decimals, and the truncation is not harmless (at 0.95
#   for rho_a productivity's Q40 share of GDP falls from 30.9% to 27.2%).
#   At these values the solver matches Dynare on the authors' model file to
#   1e-10. The five values below the estimated block are calibrated in the
#   paper, not estimated, and Whelan does not show them.

C_01_01_default_lst <- list(

  # Table 1A, structural (mode; the table's two-decimal value in brackets)
  csadjcost  = 5.48819701, # varphi     capital adjustment cost          [5.48]
  csigma     = 1.3951929, # sigma_c    inverse intertemporal elasticity [1.39]
  chabb      = 0.71240064, # lambda     habit persistence (h in Table 1A) [0.71]
  cprobw     = 0.73754132, # xi_w       Calvo wage stickiness            [0.73]
  csigl      = 1.91988384, # sigma_l    labour supply elasticity         [1.92]
  cprobp     = 0.65626626, # xi_p       Calvo price stickiness           [0.65]
  cindw      = 0.59199831, # iota_w     wage indexation                  [0.59]
  cindp      = 0.22835402, # iota_p     price indexation                 [0.22]
  czcap      = 0.54721313, # psi        capacity utilisation cost        [0.54]
  cfc        = 1.61497959, # Phi        1 + fixed cost share             [1.61]
  crpi       = 2.0294674, # r_pi       policy on inflation              [2.03]
  crr        = 0.81532487, # rho        policy smoothing                 [0.81]
  cry        = 0.08468691, # r_y        policy on the output gap         [0.08]
  crdy       = 0.22292571, # r_dy       policy on the change in the gap  [0.22]
  constepinf = 0.81798222, # pi_bar     steady-state inflation, quarterly %
  constebeta = 0.16065411, # 100(1/beta - 1)
  constelab  = -0.10306517, # l_bar     steady-state hours
  ctrend     = 0.43202637, # gamma_bar  trend growth, quarterly %
  calfa      = 0.19280046, # alpha      capital share                    [0.19]

  # Table 1B, shock processes (mode)
  crhoa    = 0.9587741,  crhob    = 0.18243935, crhog  = 0.97616142,
  crhoqs   = 0.70956932, crhoms   = 0.12713148, crhopinf = 0.90380734,
  crhow    = 0.97185377, cmap     = 0.74487185, cmaw   = 0.88814593,
  cgy      = 0.52612122,
  sda      = 0.45178828, sdb      = 0.2424607,  sdg    = 0.52001032,
  sdqs     = 0.45010691, sdms     = 0.23983933, sdpinf = 0.14112385,
  sdw      = 0.2443916,

  # Calibrated in the paper, not estimated
  delta      = 0.025,  # depreciation
  gy         = 0.18,   # steady-state spending share
  lambdaw    = 1.5,    # steady-state wage mark-up
  curvp      = 10,     # Kimball price aggregator curvature
  curvw      = 10,     # Kimball wage aggregator curvature

  # Which model: "sw" is Smets and Wouters' own (see the header)
  model_chr  = "sw")

###### C_01_01b: Our Model #####################################################
# Note: The same deep parameters, the module's assumptions. Nothing is
#   re-estimated: the posterior mode above is run through our model's
#   structure (header). ctrend stays in the list because the Smets-Wouters
#   side needs it; under "ours" C_01_03 ignores it and sets gamma = 1.

C_01_01b_ours_lst <- utils::modifyList(C_01_01_default_lst,
                                       list(model_chr = "ours"))

###### C_01_01c: Which Model a Parameter List Solves ###########################
# Note: TRUE for our model, FALSE for Smets and Wouters' own.

C_01_01c_ours_fn <- function(par_lst) identical(par_lst$model_chr, "ours")

###### C_01_02: Variable Order #################################################
# Note: The order of every matrix row and column in the file. The block
#   after the flexible-price economy holds the seven shocks and the two
#   moving-average terms the mark-up shocks need; the last block holds one
#   variable per expectation, because the canonical form carries E_t x_(t+1)
#   as a variable of its own. Each expectation is defined by x_t =
#   E_(t-1) x_t + eta_t, and C_01_02_pair_lst writes that pairing once.

C_01_02_names_vec <- c(
  # sticky-price economy
  "mup", "zcap", "rk", "k", "inve", "pk", "c", "y", "lab", "pinf", "w", "r",
  "kp",
  # flexible-price shadow economy
  "zcapf", "rkf", "kf", "invef", "pkf", "cf", "yf", "labf", "wf", "rrf",
  "kpf",
  # shocks
  "a", "b", "g", "qs", "ms", "spinf", "sw", "epinfma", "ewma",
  # expectations
  "Ec", "Einve", "Epk", "Epinf", "Ew", "Erk", "Elab",
  "Ecf", "Einvef", "Epkf", "Erkf", "Elabf")

C_01_02_shock_vec <- c("ea", "eb", "eg", "eqs", "em", "epinf", "ew")

C_01_02_shock_lab_vec <- c(
  ea    = "Productivity",     eb     = "Risk Premium",
  eg    = "Exogenous Spending", eqs  = "Investment",
  em    = "Monetary Policy",  epinf  = "Price Mark-Up",
  ew    = "Wage Mark-Up")

C_01_02_eta_vec <- c(
  "Ec", "Einve", "Epk", "Epinf", "Ew", "Erk", "Elab",
  "Ecf", "Einvef", "Epkf", "Erkf", "Elabf")

C_01_02_pair_lst <- list(
  c("Ec", "c"), c("Einve", "inve"), c("Epk", "pk"), c("Epinf", "pinf"),
  c("Ew", "w"), c("Erk", "rk"), c("Elab", "lab"),
  c("Ecf", "cf"), c("Einvef", "invef"), c("Epkf", "pkf"),
  c("Erkf", "rkf"), c("Elabf", "labf"))

###### C_01_03: Steady State and Composite Coefficients ########################
# Note: The mapping the slides do not carry, from Smets and Wouters (2007).
#   The steady-state ratios come first because every composite below is a
#   function of them, and c_y + i_y + g_y = 1 is the check that they are
#   right (C_05_02). Under "ours" gamma is one, so every gamma drops out and
#   the ratios follow; c_2 = 0 because utility is separable in hours; the
#   marginal rate of substitution's weight is sigma_c/(1 - lambda) rather
#   than 1/(1 - lambda/gamma); and k_2 = (1 - k_1)/i_2, the form eq (8)
#   prints, where the authors' code drops the factor (1 + beta gamma^(1 -
#   sigma_c)). Their side keeps the code's form.


C_01_03_derived_fn <- function(par_lst) {

  d_lst <- par_lst

  ours_lgl   <- C_01_01c_ours_fn(par_lst)
  cpie_num   <- 1 + par_lst$constepinf / 100
  # ours has no trend growth, so gamma is one
  cgamma_num <- if (ours_lgl) 1 else 1 + par_lst$ctrend / 100
  cbeta_num  <- 1 / (1 + par_lst$constebeta / 100)
  calfa_num  <- par_lst$calfa
  ctou_num   <- par_lst$delta
  cfc_num    <- par_lst$cfc

  cbetabar_num <- cbeta_num * cgamma_num^(-par_lst$csigma)
  crk_num      <- cbeta_num^(-1) * cgamma_num^par_lst$csigma - (1 - ctou_num)
  cw_num       <- (calfa_num^calfa_num * (1 - calfa_num)^(1 - calfa_num) /
                     (cfc_num * crk_num^calfa_num))^(1 / (1 - calfa_num))
  cikbar_num   <- 1 - (1 - ctou_num) / cgamma_num
  cik_num      <- cikbar_num * cgamma_num
  clk_num      <- ((1 - calfa_num) / calfa_num) * (crk_num / cw_num)
  cky_num      <- cfc_num * clk_num^(calfa_num - 1)
  ciy_num      <- cik_num * cky_num
  ccy_num      <- 1 - par_lst$gy - ciy_num
  crkky_num    <- crk_num * cky_num
  cwhlc_num    <- (1 / par_lst$lambdaw) * (1 - calfa_num) / calfa_num *
    crk_num * cky_num / ccy_num
  hg_num       <- par_lst$chabb / cgamma_num

  d_lst$cpie   <- cpie_num
  d_lst$cgamma <- cgamma_num
  d_lst$cbeta  <- cbeta_num
  d_lst$cbetabar <- cbetabar_num
  d_lst$crk    <- crk_num
  d_lst$cw     <- cw_num
  d_lst$cikbar <- cikbar_num
  d_lst$cky    <- cky_num
  d_lst$ciy    <- ciy_num
  d_lst$ccy    <- ccy_num
  d_lst$crkky  <- crkky_num
  d_lst$cwhlc  <- cwhlc_num
  d_lst$hg     <- hg_num

  # consumption, [W11 5]; ours has no hours term, c_2 = 0
  d_lst$c1 <- hg_num / (1 + hg_num)
  d_lst$c2 <- if (ours_lgl) 0 else ((par_lst$csigma - 1) * cwhlc_num) /
    (par_lst$csigma * (1 + hg_num))
  d_lst$c3 <- (1 - hg_num) / (par_lst$csigma * (1 + hg_num))

  # the marginal rate of substitution's weight on c - (lambda/gamma) c(-1)
  d_lst$wc <- if (ours_lgl) par_lst$csigma / (1 - hg_num) else
    1 / (1 - hg_num)

  # investment and Tobin's q, [W11 6]
  d_lst$i1  <- 1 / (1 + cbetabar_num * cgamma_num)
  d_lst$i2  <- d_lst$i1 / (cgamma_num^2 * par_lst$csadjcost)
  d_lst$q1  <- (1 - ctou_num) / (crk_num + 1 - ctou_num)
  d_lst$q1r <- crk_num / (crk_num + 1 - ctou_num)

  # capital accumulation, eq (8); k_2 differs between the two models
  d_lst$z1 <- (1 - par_lst$czcap) / par_lst$czcap
  d_lst$k1 <- 1 - cikbar_num
  d_lst$k2 <- if (ours_lgl) cikbar_num / d_lst$i2 else
    cikbar_num * cgamma_num^2 * par_lst$csadjcost

  # prices, [W11 7]
  d_lst$p1 <- par_lst$cindp / (1 + cbetabar_num * cgamma_num * par_lst$cindp)
  d_lst$p2 <- cbetabar_num * cgamma_num /
    (1 + cbetabar_num * cgamma_num * par_lst$cindp)
  d_lst$p3 <- (1 / (1 + cbetabar_num * cgamma_num * par_lst$cindp)) *
    ((1 - par_lst$cprobp) *
       (1 - cbetabar_num * cgamma_num * par_lst$cprobp) / par_lst$cprobp) /
    ((cfc_num - 1) * par_lst$curvp + 1)

  # wages, [W11 8]; w2 multiplies current inflation and w3 lagged, as eq
  # (13) numbers them, not as the circulating Dynare replication does
  d_lst$w1 <- 1 / (1 + cbetabar_num * cgamma_num)
  d_lst$w2 <- (1 + cbetabar_num * cgamma_num * par_lst$cindw) /
    (1 + cbetabar_num * cgamma_num)
  d_lst$w3 <- par_lst$cindw / (1 + cbetabar_num * cgamma_num)
  d_lst$w4 <- (1 - par_lst$cprobw) *
    (1 - cbetabar_num * cgamma_num * par_lst$cprobw) /
    ((1 + cbetabar_num * cgamma_num) * par_lst$cprobw) /
    ((par_lst$lambdaw - 1) * par_lst$curvw + 1)

  d_lst
}

#### C_02: The Canonical System ################################################
# Note: One function writes G0, G1, Psi and Pi row by row.

###### C_02_01: Build the Matrices #############################################
# Note: G0 y_t = G1 y_(t-1) + Psi e_t + Pi eta_t. The three helpers write
#   into the current equation's row, so each equation below reads in the
#   order the slide writes it. A term on the right of Whelan's equation
#   enters G0 negated, because everything time-t is collected on the left.


C_02_01_build_fn <- function(d_lst) {

  n_int  <- length(C_01_02_names_vec)
  ix_fn  <- function(nm_chr) match(nm_chr, C_01_02_names_vec)
  jx_fn  <- function(nm_chr) match(nm_chr, C_01_02_shock_vec)
  kx_fn  <- function(nm_chr) match(nm_chr, C_01_02_eta_vec)

  g_zero_mat <- matrix(0, n_int, n_int)
  g_one_mat  <- matrix(0, n_int, n_int)
  psi_mat    <- matrix(0, n_int, length(C_01_02_shock_vec))
  pi_mat     <- matrix(0, n_int, length(C_01_02_eta_vec))
  row_int    <- 0L

  t_fn <- function(nm_chr, val_num = 1) {
    g_zero_mat[row_int, ix_fn(nm_chr)] <<-
      g_zero_mat[row_int, ix_fn(nm_chr)] + val_num
  }
  l_fn <- function(nm_chr, val_num = 1) {
    g_one_mat[row_int, ix_fn(nm_chr)] <<-
      g_one_mat[row_int, ix_fn(nm_chr)] + val_num
  }
  s_fn <- function(nm_chr, val_num = 1) {
    psi_mat[row_int, jx_fn(nm_chr)] <<-
      psi_mat[row_int, jx_fn(nm_chr)] + val_num
  }
  nx_fn <- function() row_int <<- row_int + 1L

  al_num <- d_lst$calfa

  ## sticky-price economy, [W11 3] to [W11 9] ---------------------------------
  # Note: The rows below are the Equations panel (app.R B_04_01), symbol for
  #   symbol; V_16 in the tests re-evaluates every displayed equation on this
  #   solution. Two choices differ from the circulating Dynare replication.
  #   The price mark-up mu^p is the variable, not marginal cost: mc = -mu^p
  #   exactly, so the solution is unchanged and the row reads as [W11 7]. The
  #   risk premium shock is the paper's, in the paper's units: eqs (2) and
  #   (4) put varepsilon^b inside the real-rate bracket, so a positive shock
  #   is a rise in the premium (p. 589). The authors' code writes b_code =
  #   -c_3 varepsilon^b, so Table 1B's sigma_b is the s.d. of b_code and
  #   C_04_00_sigma_fn supplies sigma_b / c_3; every variance share is
  #   unchanged and every response is Dynare's with its sign flipped.

  nx_fn()   # the price mark-up, [W11 7]
  t_fn("mup"); t_fn("k", -al_num); t_fn("lab", al_num); t_fn("a", -1)
  t_fn("w", 1)

  nx_fn()   # capital utilisation
  t_fn("zcap"); t_fn("rk", -d_lst$z1)

  nx_fn()   # rental rate of capital
  t_fn("rk"); t_fn("w", -1); t_fn("lab", -1); t_fn("k", 1)

  nx_fn()   # capital in use
  t_fn("k"); t_fn("zcap", -1); l_fn("kp", 1)

  nx_fn()   # investment
  t_fn("inve"); l_fn("inve", d_lst$i1); t_fn("Einve", -(1 - d_lst$i1))
  t_fn("pk", -d_lst$i2); t_fn("qs", -1)

  nx_fn()   # value of capital
  t_fn("pk"); t_fn("r", 1); t_fn("Epinf", -1); t_fn("b", 1)
  t_fn("Erk", -d_lst$q1r); t_fn("Epk", -d_lst$q1)

  nx_fn()   # consumption
  t_fn("c"); l_fn("c", d_lst$c1); t_fn("Ec", -(1 - d_lst$c1))
  t_fn("lab", -d_lst$c2); t_fn("Elab", d_lst$c2)
  t_fn("r", d_lst$c3); t_fn("Epinf", -d_lst$c3); t_fn("b", d_lst$c3)

  nx_fn()   # resource constraint
  t_fn("y"); t_fn("c", -d_lst$ccy); t_fn("inve", -d_lst$ciy); t_fn("g", -1)
  t_fn("zcap", -d_lst$crkky)

  nx_fn()   # production function
  t_fn("y"); t_fn("k", -d_lst$cfc * al_num)
  t_fn("lab", -d_lst$cfc * (1 - al_num)); t_fn("a", -d_lst$cfc)

  nx_fn()   # price Phillips curve
  t_fn("pinf"); l_fn("pinf", d_lst$p1); t_fn("Epinf", -d_lst$p2)
  t_fn("mup", d_lst$p3); t_fn("spinf", -1)

  nx_fn()   # wage Phillips curve
  t_fn("w"); l_fn("w", d_lst$w1)
  t_fn("Ew", -(1 - d_lst$w1)); t_fn("Epinf", -(1 - d_lst$w1))
  l_fn("pinf", d_lst$w3); t_fn("pinf", d_lst$w2)
  t_fn("lab", -d_lst$w4 * d_lst$csigl)
  t_fn("c", -d_lst$w4 * d_lst$wc)
  l_fn("c", -d_lst$w4 * d_lst$wc * d_lst$hg)
  t_fn("w", d_lst$w4); t_fn("sw", -1)

  nx_fn()   # monetary policy rule
  t_fn("r"); t_fn("pinf", -d_lst$crpi * (1 - d_lst$crr))
  t_fn("y", -d_lst$cry * (1 - d_lst$crr))
  t_fn("yf", d_lst$cry * (1 - d_lst$crr))
  t_fn("y", -d_lst$crdy); t_fn("yf", d_lst$crdy)
  l_fn("y", -d_lst$crdy); l_fn("yf", d_lst$crdy)
  l_fn("r", d_lst$crr); t_fn("ms", -1)

  nx_fn()   # capital accumulation
  # Note: Eq (8) prints k_2 = (1 - (1 - delta)/gamma)(1 + beta gamma^(1 -
  #   sigma_c)) gamma^2 varphi; every Dynare replication in circulation
  #   drops the middle factor, 1.9967 at the mode, so it is not rounding.
  #   The derivation gives the printed form: the raw shock loads on capital
  #   by (1 - k_1), and rescaling it to enter the investment equation with
  #   coefficient one gives k_2 = (1 - k_1)/i_2. Their side keeps the code's
  #   form, because the posterior was estimated under it and the app's
  #   claims rest on internal consistency; ours takes the derived form, since
  #   it is not their estimated model. C_01_03 computes k_2 for each. At the
  #   mode the difference moves the investment shock's impact on output from
  #   +0.301 to +0.312 and its Q1 share of GDP's variance from 13.6% to 14.5%.
  t_fn("kp"); l_fn("kp", d_lst$k1); t_fn("inve", -(1 - d_lst$k1))
  t_fn("qs", -d_lst$k2)

  ## flexible-price shadow economy --------------------------------------------
  # Note: The same model with no price or wage stickiness and no mark-up
  #   shocks. [W11 9] defines potential output as output under fully flexible
  #   prices and wages, and the policy rule reads the gap to it. Drawn only
  #   as y - yf.

  nx_fn()   # the price mark-up is zero: alpha(k^s - l) + a - w = 0
  t_fn("kf", al_num); t_fn("labf", -al_num); t_fn("a", 1); t_fn("wf", -1)

  nx_fn()
  t_fn("zcapf"); t_fn("rkf", -d_lst$z1)

  nx_fn()
  t_fn("rkf"); t_fn("wf", -1); t_fn("labf", -1); t_fn("kf", 1)

  nx_fn()
  t_fn("kf"); t_fn("zcapf", -1); l_fn("kpf", 1)

  nx_fn()
  t_fn("invef"); l_fn("invef", d_lst$i1); t_fn("Einvef", -(1 - d_lst$i1))
  t_fn("pkf", -d_lst$i2); t_fn("qs", -1)

  nx_fn()
  t_fn("pkf"); t_fn("rrf", 1); t_fn("b", 1)
  t_fn("Erkf", -d_lst$q1r); t_fn("Epkf", -d_lst$q1)

  nx_fn()
  t_fn("cf"); l_fn("cf", d_lst$c1); t_fn("Ecf", -(1 - d_lst$c1))
  t_fn("labf", -d_lst$c2); t_fn("Elabf", d_lst$c2)
  t_fn("rrf", d_lst$c3); t_fn("b", d_lst$c3)

  nx_fn()
  t_fn("yf"); t_fn("cf", -d_lst$ccy); t_fn("invef", -d_lst$ciy)
  t_fn("g", -1); t_fn("zcapf", -d_lst$crkky)

  nx_fn()
  t_fn("yf"); t_fn("kf", -d_lst$cfc * al_num)
  t_fn("labf", -d_lst$cfc * (1 - al_num)); t_fn("a", -d_lst$cfc)

  nx_fn()   # the flexible wage equals the marginal rate of substitution
  t_fn("wf"); t_fn("labf", -d_lst$csigl); t_fn("cf", -d_lst$wc)
  l_fn("cf", -d_lst$wc * d_lst$hg)

  nx_fn()
  t_fn("kpf"); l_fn("kpf", d_lst$k1); t_fn("invef", -(1 - d_lst$k1))
  t_fn("qs", -d_lst$k2)

  ## shock processes, [W11 3] to [W11 9] --------------------------------------
  # Note: The two mark-up shocks are ARMA(1,1), which is what the epinfma and
  #   ewma variables are for. rho_ga loads productivity onto spending, which
  #   is Whelan's "net exports may be affected by domestic productivity".

  nx_fn(); t_fn("a"); l_fn("a", d_lst$crhoa); s_fn("ea")
  nx_fn(); t_fn("b"); l_fn("b", d_lst$crhob); s_fn("eb")
  nx_fn(); t_fn("g"); l_fn("g", d_lst$crhog); s_fn("eg"); s_fn("ea", d_lst$cgy)
  nx_fn(); t_fn("qs"); l_fn("qs", d_lst$crhoqs); s_fn("eqs")
  nx_fn(); t_fn("ms"); l_fn("ms", d_lst$crhoms); s_fn("em")
  nx_fn(); t_fn("spinf"); l_fn("spinf", d_lst$crhopinf); s_fn("epinf")
  l_fn("epinfma", -d_lst$cmap)
  nx_fn(); t_fn("epinfma"); s_fn("epinf")
  nx_fn(); t_fn("sw"); l_fn("sw", d_lst$crhow); s_fn("ew")
  l_fn("ewma", -d_lst$cmaw)
  nx_fn(); t_fn("ewma"); s_fn("ew")

  ## expectation definitions --------------------------------------------------
  # Note: x_t = E_(t-1) x_t + eta_t. One forecast error each, and the count of
  #   them is what the Blanchard-Kahn check in C_05_01 compares against.

  for (pair_chr in C_01_02_pair_lst) {
    nx_fn()
    t_fn(pair_chr[2L]); l_fn(pair_chr[1L], 1)
    pi_mat[row_int, kx_fn(pair_chr[1L])] <- 1
  }

  if (row_int != n_int) {
    stop(sprintf("built %d equations for %d variables", row_int, n_int))
  }

  list(g_zero = g_zero_mat, g_one = g_one_mat,
       psi = psi_mat, pi = pi_mat)
}

#### C_03: The Solver ##########################################################
# Note: A pseudo-inverse and the shifted-pencil solve.

###### C_03_01: Pseudo-Inverse #################################################
# Note: base R has no ginv outside MASS, which a shinylive export may lack.


C_03_01_pinv_fn <- function(m_mat, tol_num = 1e-10) {
  svd_lst  <- svd(m_mat)
  keep_lgl <- svd_lst$d > tol_num * max(svd_lst$d[1L], 1)
  if (!any(keep_lgl)) return(t(m_mat) * 0)
  svd_lst$v[, keep_lgl, drop = FALSE] %*%
    diag(1 / svd_lst$d[keep_lgl], nrow = sum(keep_lgl)) %*%
    t(svd_lst$u[, keep_lgl, drop = FALSE])
}

###### C_03_02: Solve ##########################################################
# Note: The method is set out in the header. Everything returned is real:
#   the stable subspace is spanned by a real basis by construction, so no
#   complex arithmetic escapes this function.

C_03_02_solve_fn <- function(par_lst, shift_num = 1.7182818) {

  d_lst   <- C_01_03_derived_fn(par_lst)
  mat_lst <- C_02_01_build_fn(d_lst)
  n_int   <- length(C_01_02_names_vec)

  fail_lst <- list(ok_lgl = FALSE, d_lst = d_lst, mat_lst = mat_lst,
                   problems_chr = character(0))

  # shift the pencil so the inverse exists
  shifted_mat <- mat_lst$g_one - shift_num * mat_lst$g_zero
  step_mat    <- tryCatch(solve(shifted_mat, mat_lst$g_zero),
                          error = function(e) NULL)
  if (is.null(step_mat)) {
    fail_lst$problems_chr <- "The shifted pencil is singular at these values."
    return(fail_lst)
  }

  # eigenvalues of the pencil
  eig_lst <- eigen(step_mat)
  mu_vec  <- eig_lst$values
  lam_vec <- ifelse(Mod(mu_vec) > 1e-12, shift_num + 1 / mu_vec,
                    complex(real = Inf, imaginary = 0))
  stable_int <- which(Mod(lam_vec) < 1 - 1e-09)

  # a real basis for the stable subspace; a cluster of m equal real roots
  # at mu0 can be defective, so its basis is the null space of
  # (M - mu0 I)^k at the lowest k that gives m dimensions, and distinct
  # roots keep their eigenvectors
  cols_lst <- list()
  used_lgl <- rep(FALSE, length(lam_vec))
  real_int <- stable_int[abs(Im(lam_vec[stable_int])) < 1e-10]
  for (j_int in real_int) {
    if (used_lgl[j_int]) next
    grp_int <- real_int[!used_lgl[real_int] &
                          abs(Re(mu_vec[real_int]) - Re(mu_vec[j_int])) <
                          1e-06 * max(1, abs(Re(mu_vec[j_int])))]
    if (length(grp_int) < 2L) next
    mu0_num  <- mean(Re(mu_vec[grp_int]))
    pow_mat  <- diag(n_int)
    base_mat <- step_mat - mu0_num * diag(n_int)
    # the lowest power whose null space has the cluster's dimension
    for (k_int in seq_along(grp_int)) {
      pow_mat <- pow_mat %*% base_mat
      svd_lst <- svd(pow_mat)
      if (sum(svd_lst$d < 1e-08 * svd_lst$d[1L]) >= length(grp_int)) break
    }
    v_mat <- svd_lst$v
    for (c_int in seq_along(grp_int)) {
      cols_lst[[length(cols_lst) + 1L]] <- v_mat[, n_int - c_int + 1L]
    }
    used_lgl[grp_int] <- TRUE
  }
  for (j_int in stable_int) {
    if (used_lgl[j_int]) next
    used_lgl[j_int] <- TRUE
    if (abs(Im(lam_vec[j_int])) < 1e-10) {
      cols_lst[[length(cols_lst) + 1L]] <- Re(eig_lst$vectors[, j_int])
    } else {
      mate_int <- stable_int[!used_lgl[stable_int] &
                               Mod(lam_vec[stable_int] -
                                     Conj(lam_vec[j_int])) < 1e-09]
      if (length(mate_int) > 0L) used_lgl[mate_int[1L]] <- TRUE
      cols_lst[[length(cols_lst) + 1L]] <- Re(eig_lst$vectors[, j_int])
      cols_lst[[length(cols_lst) + 1L]] <- Im(eig_lst$vectors[, j_int])
    }
  }
  if (length(cols_lst) == 0L) {
    fail_lst$problems_chr <- "No stable roots: the model has no solution here."
    return(fail_lst)
  }
  basis_mat <- do.call(cbind, cols_lst)
  ns_int    <- ncol(basis_mat)
  nu_int    <- n_int - ns_int
  neta_int  <- length(C_01_02_eta_vec)

  # the left null space of G0 B, which kills the unstable directions
  k_mat   <- mat_lst$g_zero %*% basis_mat
  ksvd_lst <- svd(k_mat, nu = n_int)
  if (nu_int <= 0L) {
    eta_mat <- matrix(0, neta_int, length(C_01_02_shock_vec))
    rank_int <- 0L
  } else {
    left_mat <- t(ksvd_lst$u[, (ns_int + 1L):n_int, drop = FALSE])
    lpi_mat  <- left_mat %*% mat_lst$pi
    lps_mat  <- left_mat %*% mat_lst$psi
    rank_int <- sum(svd(lpi_mat)$d > 1e-09 * max(svd(lpi_mat)$d[1L], 1))
    eta_mat  <- -C_03_01_pinv_fn(lpi_mat) %*% lps_mat
  }

  # Blanchard and Kahn, counted
  problems_chr <- character(0)
  if (nu_int > neta_int) {
    problems_chr <- c(problems_chr, sprintf(
      paste("%d unstable roots but only %d forecast errors: no stable",
            "solution exists at these values."), nu_int, neta_int))
  }
  if (nu_int < neta_int) {
    problems_chr <- c(problems_chr, sprintf(
      paste("%d unstable roots against %d forecast errors: the solution is",
            "indeterminate, so sunspots are possible."), nu_int, neta_int))
  }
  if (nu_int == neta_int && rank_int < nu_int) {
    problems_chr <- c(problems_chr,
                      "The forecast errors do not span the unstable block.")
  }

  # the transition and shock loadings, in y space
  kp_mat  <- C_03_01_pinv_fn(k_mat)
  bp_mat  <- C_03_01_pinv_fn(basis_mat)
  r_mat   <- kp_mat %*% (mat_lst$g_one %*% basis_mat)
  t_one_mat <- basis_mat %*% r_mat %*% bp_mat
  t_psi_mat <- basis_mat %*% kp_mat %*%
    (mat_lst$psi + mat_lst$pi %*% eta_mat)

  dimnames(t_one_mat) <- list(C_01_02_names_vec, C_01_02_names_vec)
  dimnames(t_psi_mat) <- list(C_01_02_names_vec, C_01_02_shock_vec)

  list(ok_lgl = length(problems_chr) == 0L,
       t_one_mat = t_one_mat, t_psi_mat = t_psi_mat,
       d_lst = d_lst, par_lst = par_lst, mat_lst = mat_lst,
       lam_vec = sort(Mod(lam_vec)),
       ns_int = ns_int, nu_int = nu_int, neta_int = neta_int,
       rank_int = rank_int,
       cond_num = kappa(basis_mat),
       problems_chr = problems_chr)
}

#### C_04: Simulation, Responses and Decompositions ############################
# Note: Everything computed from a solution.

###### C_04_00: Innovation Standard Deviations, in the Equations' Units ########
# Note: Table 1B's sigmas, with one conversion. Six shocks enter the
#   equations as the authors' code enters them, so their Table 1B s.d.
#   applies as printed. The risk premium sits inside the real-rate bracket
#   of eqs (2) and (4), while the code estimates c_3 varepsilon^b, so the
#   s.d. of varepsilon^b itself is sigma_b / c_3, with c_3 taken from the
#   solve's own composites.

C_04_00_sd_vec <- c(ea = "sda", eb = "sdb", eg = "sdg", eqs = "sdqs",
                    em = "sdms", epinf = "sdpinf", ew = "sdw")

C_04_00_sigma_fn <- function(sol_lst) {
  out_vec <- vapply(C_01_02_shock_vec, function(sh_chr) {
    val_num <- sol_lst$par_lst[[C_04_00_sd_vec[[sh_chr]]]]
    if (is.null(val_num) || !is.finite(val_num) || val_num <= 0) {
      stop("C_04_00_sigma_fn: no positive standard deviation for ", sh_chr)
    }
    val_num
  }, numeric(1))
  out_vec[["eb"]] <- out_vec[["eb"]] / sol_lst$d_lst$c3
  out_vec
}

###### C_04_01: Impulse Responses ##############################################
# Note: The default size is the shock's own posterior standard deviation
#   from [W11 13], in the equations' units (C_04_00), so the panels are
#   comparable across shocks. Quarter 0 is the impact period.

C_04_01_irf_fn <- function(sol_lst, shock_chr, n_horizon_int = 20L,
                           size_num = NULL) {

  scale_num <- if (is.null(size_num)) {
    C_04_00_sigma_fn(sol_lst)[[shock_chr]]
  } else size_num

  e_vec <- rep(0, length(C_01_02_shock_vec))
  e_vec[match(shock_chr, C_01_02_shock_vec)] <- scale_num

  out_mat <- matrix(0, n_horizon_int, length(C_01_02_names_vec),
                    dimnames = list(NULL, C_01_02_names_vec))
  y_vec <- as.vector(sol_lst$t_psi_mat %*% e_vec)
  out_mat[1L, ] <- y_vec
  for (t_int in seq_len(n_horizon_int - 1L)) {
    y_vec <- as.vector(sol_lst$t_one_mat %*% y_vec)
    out_mat[t_int + 1L, ] <- y_vec
  }
  as.data.frame(out_mat)
}

###### C_04_02: Forecast Error Variance Decomposition ##########################
# Note: [W11 15] to [W11 17]. The share of the h-step forecast error
#   variance of each variable attributable to each shock; Whelan's horizons
#   are the default. t_psi_mat loads unit-variance innovations, so shock j's
#   contribution at horizon h is sigma_j^2 sum_s psi_j(s)^2; the sigmas in
#   [W11 13] run from 0.14 to 0.52, so they matter. V_06b in the tests
#   checks the Q1 share against the normalised squared impact response.

C_04_02_fevd_fn <- function(sol_lst,
                            horizon_vec = c(1L, 2L, 4L, 10L, 40L, 100L)) {

  n_int <- length(C_01_02_names_vec)
  k_int <- length(C_01_02_shock_vec)
  acc_mat   <- matrix(0, n_int, k_int)
  power_mat <- diag(n_int)
  out_lst   <- list()

  sigma_vec <- C_04_00_sigma_fn(sol_lst)
  load_mat <- sweep(sol_lst$t_psi_mat, 2L, sigma_vec, `*`)

  for (step_int in seq_len(max(horizon_vec))) {
    contrib_mat <- power_mat %*% load_mat
    acc_mat     <- acc_mat + contrib_mat^2
    power_mat   <- power_mat %*% sol_lst$t_one_mat
    if (step_int %in% horizon_vec) {
      total_vec <- rowSums(acc_mat)
      total_vec[total_vec <= 0] <- NA_real_
      share_mat <- acc_mat / total_vec
      dimnames(share_mat) <- list(C_01_02_names_vec, C_01_02_shock_vec)
      out_lst[[as.character(step_int)]] <- share_mat
    }
  }
  out_lst
}

###### C_04_03: Simulate #######################################################
# Note: For anything that needs a path rather than a response. Seeded, so
#   the figure on the slide and the figure in the app are the same figure.

C_04_03_simulate_fn <- function(sol_lst, n_periods_int = 1000L,
                                seed_int = 42L) {

  set.seed(seed_int)
  k_int <- length(C_01_02_shock_vec)
  sd_vec <- C_04_00_sigma_fn(sol_lst)

  out_mat <- matrix(0, n_periods_int, length(C_01_02_names_vec),
                    dimnames = list(NULL, C_01_02_names_vec))
  y_vec <- rep(0, length(C_01_02_names_vec))
  for (t_int in seq_len(n_periods_int)) {
    e_vec <- stats::rnorm(k_int, 0, sd_vec)
    y_vec <- as.vector(sol_lst$t_one_mat %*% y_vec +
                         sol_lst$t_psi_mat %*% e_vec)
    out_mat[t_int, ] <- y_vec
  }
  as.data.frame(out_mat)
}

#### C_05: Diagnostics #########################################################
# Note: Two checks run on every solve.

###### C_05_01: Structural Residual ############################################
# Note: The solution must satisfy the canonical system on its own path;
#   feeding it an arbitrary lagged state is not a valid test, because the
#   solution kills the unstable directions such a state contains. It also
#   checks rational expectations directly: each E variable must equal the
#   model's own forecast.

C_05_01_residual_fn <- function(sol_lst, n_draw_int = 200L, seed_int = 11L) {

  if (!isTRUE(sol_lst$ok_lgl) && is.null(sol_lst$t_one_mat)) {
    return(list(system_num = NA_real_, expect_num = NA_real_))
  }
  set.seed(seed_int)
  m_lst <- sol_lst$mat_lst
  ix_fn <- function(nm_chr) match(nm_chr, C_01_02_names_vec)

  y_vec <- rep(0, length(C_01_02_names_vec))
  worst_sys_num <- 0
  worst_exp_num <- 0

  for (t_int in seq_len(n_draw_int)) {
    e_vec <- stats::rnorm(length(C_01_02_shock_vec))
    lag_vec <- y_vec
    y_vec <- as.vector(sol_lst$t_one_mat %*% lag_vec +
                         sol_lst$t_psi_mat %*% e_vec)
    eta_vec <- vapply(C_01_02_pair_lst, function(pair_chr)
      y_vec[ix_fn(pair_chr[2L])] - lag_vec[ix_fn(pair_chr[1L])], numeric(1))
    res_vec <- m_lst$g_zero %*% y_vec -
      (m_lst$g_one %*% lag_vec + m_lst$psi %*% e_vec + m_lst$pi %*% eta_vec)
    worst_sys_num <- max(worst_sys_num, max(abs(res_vec)))
    fwd_vec <- as.vector(sol_lst$t_one_mat %*% y_vec)
    worst_exp_num <- max(worst_exp_num, max(vapply(
      C_01_02_pair_lst, function(pair_chr)
        abs(y_vec[ix_fn(pair_chr[1L])] - fwd_vec[ix_fn(pair_chr[2L])]),
      numeric(1))))
  }
  list(system_num = worst_sys_num, expect_num = worst_exp_num)
}

###### C_05_02: Steady-State Check #############################################
# Note: The expenditure shares must sum to one, or the composites are built
#   on a steady state that does not exist.

C_05_02_shares_fn <- function(d_lst) {
  sum_num <- d_lst$ccy + d_lst$ciy + d_lst$gy
  list(ccy_num = d_lst$ccy, ciy_num = d_lst$ciy, gy_num = d_lst$gy,
       sum_num = sum_num, ok_lgl = abs(sum_num - 1) < 1e-10)
}

#### C_06: The Frictions #######################################################
# Note: [W11 10], switched off one at a time and together.

###### C_06_01: What [W11 10] Claims ###########################################
# Note: Whelan lists six additions over the RBC and says they "throw sand in
#   the wheels, making variables more sluggish and giving random shocks a
#   more long-lasting effect". Each entry switches one of them off at the
#   deep parameter, so the steady state and every other coefficient are
#   rebuilt around it rather than patched.

C_06_01_friction_lst <- list(
  habit = list(
    label_chr = "Habit Persistence",
    par_lst   = list(chabb = 1e-06),
    note_chr  = paste(
      "Consumption stops looking back. Without it the response to a policy",
      "shock is at its worst on impact, and the hump is gone.")),
  adjcost = list(
    label_chr = "Investment Adjustment Costs",
    par_lst   = list(csadjcost = 0.1),
    note_chr  = paste(
      "Investment can come on line at once. The trough is far deeper and",
      "comes sooner; this is the friction doing the most work.")),
  utilisation = list(
    label_chr = "Capacity Utilisation Costs",
    par_lst   = list(czcap = 0.99),
    note_chr  = paste(
      "Capital in use moves freely with its rental rate. The effect on the",
      "policy response is small - this friction matters for other things.")),
  pindex = list(
    label_chr = "Price Indexation",
    par_lst   = list(cindp = 1e-06),
    note_chr  = paste(
      "Firms that cannot reset do not follow past inflation. Inflation loses",
      "its backward-looking term, which is what [W11 7] added it for.")),
  windex = list(
    label_chr = "Wage Indexation",
    par_lst   = list(cindw = 1e-06),
    note_chr  = paste(
      "The same for wages, [W11 8].")),
  stickyp = list(
    label_chr = "Sticky Prices",
    par_lst   = list(cprobp = 0.05),
    note_chr  = paste(
      "Almost every firm resets every quarter. The output trough shrinks,",
      "but money is far from neutral while wages stay sticky.")))

###### C_06_02: Run One Friction Off ###########################################
# Note: Returns the same shape as a baseline solve so the drawing code can
#   treat the two identically and the ghost layer does the comparison.

C_06_02_off_fn <- function(par_lst, which_chr) {
  if (!which_chr %in% names(C_06_01_friction_lst)) {
    stop("unknown friction: ", which_chr)
  }
  over_lst <- C_06_01_friction_lst[[which_chr]]$par_lst
  C_03_02_solve_fn(utils::modifyList(par_lst, over_lst))
}

###### C_06_03: The Real Frictions Off #########################################
# Note: Habit, adjustment costs, capacity utilisation and both indexations.
#   Prices and wages are still sticky here: this is the real-friction half
#   of [W11 10], not the whole of it.

C_06_03_all_off_fn <- function(par_lst) {
  over_lst <- list(chabb = 1e-06, csadjcost = 0.1, czcap = 0.99,
                   cindp = 1e-06, cindw = 1e-06)
  C_03_02_solve_fn(utils::modifyList(par_lst, over_lst))
}

###### C_06_04: Everything Off, the Closest Approach to Part 7 #################
# Note: The RBC of part 7 sends hours up after a positive technology shock,
#   the failure Gali (1999) named; this model sends them down. Neither half
#   of part 11's additions alone flips the sign: with the real frictions off
#   and prices sticky hours still fall, with flexible prices and wages and
#   the real frictions on they still fall, and only with both off does the
#   sign flip. This is not the RBC: part 7 has a different capital share,
#   no fixed costs, no wage mark-up and one shock. It is the closest this
#   model comes.

C_06_04_rbc_limit_fn <- function(par_lst) {
  over_lst <- list(chabb = 1e-06, csadjcost = 0.1, czcap = 0.99,
                   cindp = 1e-06, cindw = 1e-06,
                   cprobp = 0.05, cprobw = 0.05)
  C_03_02_solve_fn(utils::modifyList(par_lst, over_lst))
}

#--------------------------------- Script End ---------------------------------#
