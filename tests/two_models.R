################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Smets-Wouters Model: The Two-Model Checks                                  ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Sourced by tests/verify_against_sw2007.R (V_16, V_18, V_19) after
##   R/model.R and app.R's builders are loaded. Defines functions only;
##   runs nothing.
##   W_01  each side's displayed composite formulas, evaluated from the
##         deep parameters exactly as the Equations panel prints them
##   W_02  each side's displayed equations, evaluated term by term on that
##         side's own solution path (residuals must be machine zero)
##   W_03  the stage-2.5e pill set, read off the app, and the independent
##         set of equations that differ between the two solvers
##   W_04  the caveat against the pills
##   Solver names: k = k^s (capital in use), kp = k (installed), zcap = z,
##   inve = i, pk = q, mup = mu^p, lab = hours (n_t ours, l_t theirs),
##   a/b/g/qs/ms/spinf/sw = the seven varepsilon processes, yf = y^p.
##
## Inputs:
##   None. Every function takes a parameter list, a solve or the app's
##   equation list.
##
## Outputs:
##   W_* functions and rule lists, used by the verification script.
##
## Packages:
##   None beyond base R.
##
## References:
##   Smets, F. and Wouters, R. (2007). American Economic Review 97(3).

#-------------------------------- Script Begin --------------------------------#

#### W_01: The Displayed Composite Formulas ####################################
# Note: Written from the panel, not read from C_01_03. Each side's formulas
#   are the ones its "Thresholds and Simplifications" items print. Only the
#   steady-state ratios the panel names without a formula (c_y, k_y, z_y,
#   c_whl) are taken from the solver, and i_y is then checked against the
#   resource constraint's own i_y formula.

W_01_01_coef_fn <- function(par_lst, side_chr, d_lst) {
  with(par_lst, {
    beta <- 1 / (1 + constebeta / 100)
    kappa_p <- (cfc - 1) * curvp + 1
    kappa_w <- (lambdaw - 1) * curvw + 1
    if (identical(side_chr, "ours")) {
      sig <- 1 / csigma                                 # sigma = 1/sigma_c
      lam <- chabb
      i1 <- 1 / (1 + beta)
      out <- list(
        lg = lam, bg = beta, gam = 1,
        rkbar = 1 / beta - (1 - delta),
        c1 = lam / (1 + lam), c2 = 0, c3 = sig * (1 - lam) / (1 + lam),
        mrs = 1 / (sig * (1 - lam)),
        i1 = i1, i2 = i1 / csadjcost,
        k1 = 1 - delta, k2 = (1 - (1 - delta)) / (i1 / csadjcost),
        k2_alt = delta * (1 + beta) * csadjcost,
        z1 = (1 - czcap) / czcap,
        q1 = beta * (1 - delta),
        q1_alt = (1 - delta) / ((1 / beta - (1 - delta)) + 1 - delta),
        pi1 = cindp / (1 + beta * cindp),
        pi2 = beta / (1 + beta * cindp),
        pi3 = 1 / (1 + beta * cindp) * (1 - cprobp) * (1 - beta * cprobp) /
          (cprobp * kappa_p),
        w1 = 1 / (1 + beta), w2 = (1 + beta * cindw) / (1 + beta),
        w3 = cindw / (1 + beta),
        w4 = (1 - cprobw) * (1 - beta * cprobw) /
          ((1 + beta) * cprobw * kappa_w),
        iy_formula = delta * d_lst$cky)
    } else {
      gam <- 1 + ctrend / 100
      betabar <- beta * gam^(-csigma)
      rkbar <- beta^(-1) * gam^csigma - (1 - delta)
      lg <- chabb / gam
      bg <- betabar * gam
      i1 <- 1 / (1 + bg)
      out <- list(
        lg = lg, bg = bg, gam = gam, rkbar = rkbar,
        c1 = lg / (1 + lg),
        c2 = (csigma - 1) * d_lst$cwhlc / (csigma * (1 + lg)),
        c3 = (1 - lg) / (csigma * (1 + lg)),
        mrs = 1 / (1 - lg),
        i1 = i1, i2 = i1 / (gam^2 * csadjcost),
        k1 = (1 - delta) / gam,
        k2 = (1 - (1 - delta) / gam) * gam^2 * csadjcost,
        k2_alt = (1 - (1 - delta) / gam) * gam^2 * csadjcost,
        z1 = (1 - czcap) / czcap,
        q1 = (1 - delta) / (rkbar + 1 - delta),
        q1_alt = (1 - delta) / (rkbar + 1 - delta),
        pi1 = cindp / (1 + bg * cindp),
        pi2 = bg / (1 + bg * cindp),
        pi3 = 1 / (1 + bg * cindp) * (1 - cprobp) * (1 - bg * cprobp) /
          (cprobp * kappa_p),
        w1 = 1 / (1 + bg), w2 = (1 + bg * cindw) / (1 + bg),
        w3 = cindw / (1 + bg),
        w4 = (1 - cprobw) * (1 - bg * cprobw) / ((1 + bg) * cprobw * kappa_w),
        iy_formula = (gam - 1 + delta) * d_lst$cky)
    }
    c(out, list(phip = cfc, alpha = calfa, sigl = csigl,
                cy = d_lst$ccy, iy = d_lst$ciy, zy = d_lst$crkky,
                rho = crr, rpi = crpi, ry = cry, rdy = crdy))
  })
}

###### W_01_02: The Displayed Formulas Against the Solver's Composites #########
# Note: Worst gap between each displayed composite and the solver's d_lst,
#   plus the two identities the panel prints as equalities (ours: k_2 =
#   (1 - k_1)/i_2 = delta (1 + beta) varphi and q_1 = beta (1 - delta)), and
#   the resource constraint's i_y formula.

W_01_02_gap_fn <- function(k_lst, d_lst) {
  map_vec <- c(c1 = "c1", c2 = "c2", c3 = "c3", i1 = "i1", i2 = "i2",
               k1 = "k1", k2 = "k2", z1 = "z1", q1 = "q1", pi1 = "p1",
               pi2 = "p2", pi3 = "p3", w1 = "w1", w2 = "w2", w3 = "w3",
               w4 = "w4", mrs = "wc", rkbar = "crk", lg = "hg")
  gap_vec <- vapply(names(map_vec), function(nm)
    abs(k_lst[[nm]] - d_lst[[map_vec[[nm]]]]), numeric(1))
  gap_vec <- c(gap_vec,
               k2_identity = abs(k_lst$k2 - k_lst$k2_alt),
               q1_identity = abs(k_lst$q1 - k_lst$q1_alt),
               iy_formula  = abs(k_lst$iy_formula - d_lst$ciy))
  gap_vec
}

#### W_02: Every Displayed Equation on Its Own Side's Path #####################
# Note: The panel is what the solver solves, or it is a bug, for each side.
#   The model is simulated from its solution with random innovations,
#   E_t x_(t+1) is the solution's own forecast T x_t, and each displayed
#   equation's residual must be zero in every period. The two sides differ
#   exactly where the panel says: c_2 (l - E l) only in theirs, and the
#   marginal rate of substitution's weight k_lst$mrs.

W_02_01_resid_fn <- function(sol_lst, k_lst, n_int = 120L, seed_int = 7L) {
  p_lst <- sol_lst$par_lst
  sd_vec <- C_04_00_sigma_fn(sol_lst)
  # the risk premium's s.d. from the panel's own c_3, not the solver's
  sd_vec[["eb"]] <- p_lst$sdb / k_lst$c3
  set.seed(seed_int)
  e_mat <- matrix(stats::rnorm(n_int * 7L), ncol = 7L) %*% diag(sd_vec)
  colnames(e_mat) <- C_01_02_shock_vec
  x_mat <- matrix(0, n_int, length(C_01_02_names_vec),
                  dimnames = list(NULL, C_01_02_names_vec))
  prev <- rep(0, length(C_01_02_names_vec))
  for (t_int in seq_len(n_int)) {
    prev <- as.vector(sol_lst$t_one_mat %*% prev +
                        sol_lst$t_psi_mat %*% e_mat[t_int, ])
    x_mat[t_int, ] <- prev
  }
  f_mat <- x_mat %*% t(sol_lst$t_one_mat)
  colnames(f_mat) <- C_01_02_names_vec
  rows <- 2:n_int
  X  <- function(v) x_mat[rows, v]
  L  <- function(v) x_mat[rows - 1L, v]
  E  <- function(v) f_mat[rows, v]
  H  <- function(s) e_mat[rows, s]
  HL <- function(s) e_mat[rows - 1L, s]

  with(k_lst, {
    real_fn <- function(sfx, rr_vec) {
      n <- function(v) paste0(v, sfx)
      list(
        "Aggregate Production Function" =
          X(n("y")) - phip * (alpha * X(n("k")) + (1 - alpha) * X(n("lab")) +
                                X("a")),
        "Capital in Use" =
          c(X(n("k")) - L(n("kp")) - X(n("zcap")),
            X(n("zcap")) - z1 * X(n("rk"))),
        "The Marginal Product of Capital" =
          X(n("rk")) + (X(n("k")) - X(n("lab"))) - X(n("w")),
        "The Resource Constraint" =
          X(n("y")) - cy * X(n("c")) - iy * X(n("inve")) - zy * X(n("zcap")) -
          X("g"),
        "Consumption Euler Equation" =
          X(n("c")) - c1 * L(n("c")) - (1 - c1) * E(n("c")) -
          c2 * (X(n("lab")) - E(n("lab"))) + c3 * (rr_vec + X("b")),
        "Investment Equation" =
          X(n("inve")) - i1 * L(n("inve")) - (1 - i1) * E(n("inve")) -
          i2 * X(n("pk")) - X("qs"),
        "Tobin's q" =
          X(n("pk")) - q1 * E(n("pk")) - (1 - q1) * E(n("rk")) +
          (rr_vec + X("b")),
        "Capital Accumulation" =
          X(n("kp")) - k1 * L(n("kp")) - (1 - k1) * X(n("inve")) - k2 * X("qs"))
    }
    mup_fn <- function(sfx) alpha * (X(paste0("k", sfx)) -
                                       X(paste0("lab", sfx))) +
      X("a") - X(paste0("w", sfx))
    muw_fn <- function(sfx) X(paste0("w", sfx)) -
      (sigl * X(paste0("lab", sfx)) +
         mrs * (X(paste0("c", sfx)) - lg * L(paste0("c", sfx))))
    sticky_lst <- real_fn("", X("r") - E("pinf"))
    flex_lst   <- real_fn("f", X("rrf"))
    names(flex_lst) <- paste("Potential Output:", names(flex_lst))
    c(sticky_lst,
      list(
        "The Price Mark-Up" = X("mup") - mup_fn(""),
        "New Keynesian Phillips Curve (NKPC)" =
          X("pinf") - pi1 * L("pinf") - pi2 * E("pinf") + pi3 * X("mup") -
          X("spinf"),
        "The Wage Equation" =
          X("w") - w1 * L("w") - (1 - w1) * (E("w") + E("pinf")) +
          w2 * X("pinf") - w3 * L("pinf") + w4 * muw_fn("") - X("sw"),
        "Monetary Policy (MP) Rule" =
          X("r") - rho * L("r") - (1 - rho) * (rpi * X("pinf") +
                                                 ry * (X("y") - X("yf"))) -
          rdy * ((X("y") - X("yf")) - (L("y") - L("yf"))) - X("ms"),
        "Five AR(1) Shock Processes" = with(p_lst, c(
          X("a") - crhoa * L("a") - H("ea"),
          X("b") - crhob * L("b") - H("eb"),
          X("qs") - crhoqs * L("qs") - H("eqs"),
          X("ms") - crhoms * L("ms") - H("em"),
          X("g") - crhog * L("g") - H("eg") - cgy * H("ea"))),
        "Two ARMA(1,1) Mark-Up Shocks" = with(p_lst, c(
          X("spinf") - crhopinf * L("spinf") - H("epinf") + cmap * HL("epinf"),
          X("sw") - crhow * L("sw") - H("ew") + cmaw * HL("ew"))),
        "Potential Output: mark-ups at zero" = c(mup_fn("f"), muw_fn("f"))),
      flex_lst)
  })
}

###### W_02_02: What Each Side's Display Must and Must Not Say #################
# Note: The link between W_02's R and the panel's LaTeX. W_02 is written
#   term for term from the panel; these substrings pin the terms that
#   distinguish the sides, so an edit to one without the other fails.
#   Whitespace is stripped before matching. W_02_03 picks the version each
#   side shows (ours "2.5a", theirs "2.5e" where it has one) and W_02_04
#   applies the rules.

W_02_02_tex_rules_lst <- list(
  list(label = "Consumption Euler Equation", side = "ours",
       must = "-c_3\\left(r_t-E_t\\pi_{t+1}+\\varepsilon^b_t\\right)",
       never = "c_2"),
  list(label = "Consumption Euler Equation", side = "sw",
       must = "+c_2\\left(l_t-E_tl_{t+1}\\right)", never = "n_t"),
  list(label = "The Wage Mark-Up", side = "ours",
       must = paste0("\\frac{1}{\\sigma(1-\\lambda)}",
                     "\\left(c_t-\\lambdac_{t-1}\\right)"),
       never = "\\gamma"),
  list(label = "The Wage Mark-Up", side = "sw",
       must = "\\frac{1}{1-\\lambda/\\gamma}", never = "n_t"),
  list(label = "Habit in the Consumption Weights", side = "ours",
       must = "c_3=\\frac{\\sigma(1-\\lambda)}{1+\\lambda}", never = "c_2"),
  list(label = "Habit in the Consumption Weights", side = "sw",
       must = "{\\sigma_c(1+\\lambda/\\gamma)}", never = "\\sigma(1"),
  list(label = "Adjustment Costs in the Capital Weights", side = "ours",
       must = "k_2=\\frac{1-k_1}{i_2}=\\delta(1+\\beta)\\varphi",
       never = "\\gamma"),
  list(label = "Adjustment Costs in the Capital Weights", side = "sw",
       must = "k_2=(1-k_1)\\,\\gamma^2\\varphi", never = "\\delta(1+\\beta)"),
  list(label = "The Resource Constraint", side = "ours",
       must = "i_y=\\deltak_y", never = "\\gamma"),
  list(label = "The Resource Constraint", side = "sw",
       must = "i_y=(\\gamma-1+\\delta)k_y", never = "n_t"),
  list(label = "Monetary Policy (MP) Rule", side = "ours",
       must = "\\phi_\\pi\\pi_t+\\phi_xx_t", never = "r_\\pi"),
  list(label = "Monetary Policy (MP) Rule", side = "sw",
       must = "r_\\pi\\pi_t+r_Y(y_t-y^p_t)", never = "\\phi_\\pi"),
  list(label = "Calvo Prices in the Phillips Curve", side = "ours",
       must = "(1-\\theta_p)(1-\\beta\\theta_p)", never = "\\xi_p"),
  list(label = "Calvo Prices in the Phillips Curve", side = "sw",
       must = "(1-\\xi_p)(1-\\bar{\\beta}\\gamma\\xi_p)", never = "\\theta_p"),
  list(label = "Calvo Wages in the Wage Equation", side = "ours",
       must = "(1-\\theta_w)(1-\\beta\\theta_w)", never = "\\xi_w"),
  list(label = "Calvo Wages in the Wage Equation", side = "sw",
       must = "(1-\\xi_w)(1-\\bar{\\beta}\\gamma\\xi_w)", never = "\\theta_w"),
  list(label = "Utilisation and the q Weight", side = "ours",
       must = "=\\beta(1-\\delta)", never = "\\gamma"),
  list(label = "The Steady State Beneath the Composites", side = "ours",
       must = "\\bar{r}^k=\\beta^{-1}-(1-\\delta)", never = "\\gamma"),
  list(label = "The Steady State Beneath the Composites", side = "sw",
       must = "\\beta^{-1}\\gamma^{\\sigma_c}", never = "\\theta"))

W_02_03_tex_fn <- function(eq_lst, label_chr, side_chr) {
  it_lst <- Filter(function(x) identical(x$label, label_chr), eq_lst)[[1L]]
  key_chr <- if (identical(side_chr, "sw") &&
                 "2.5e" %in% names(it_lst$versions)) "2.5e" else "2.5a"
  it_lst$versions[[key_chr]]
}

W_02_04_tex_bad_fn <- function(eq_lst) {
  bad_chr <- character(0)
  for (r_lst in W_02_02_tex_rules_lst) {
    tex_chr <- gsub("\\s+", "", W_02_03_tex_fn(eq_lst, r_lst$label, r_lst$side))
    if (!grepl(gsub("\\s+", "", r_lst$must), tex_chr, fixed = TRUE) ||
        grepl(r_lst$never, tex_chr, fixed = TRUE)) {
      bad_chr <- c(bad_chr, sprintf("%s (%s)", r_lst$label, r_lst$side))
    }
  }
  bad_chr
}

#### W_03: The Pills Against the Two Solvers ###################################
# Note: What the app marks at stage 2.5e, and what the solvers say differs.

###### W_03_01: The Pills the App Shows at Stage 2.5e ##########################
# Note: The same selection F_06_02 makes, run on B_04_01 at "2.5e".

W_03_01_pills_fn <- function(eq_lst, rank_fn, stage_chr = "2.5e") {
  rank_int <- rank_fn(stage_chr)
  out_chr <- character(0)
  for (it_lst in eq_lst) {
    key_int <- vapply(names(it_lst$versions), rank_fn, integer(1))
    if (min(key_int) > rank_int) next
    now_int <- max(key_int[key_int <= rank_int])
    if (now_int != rank_int) next
    out_chr[[it_lst$label]] <- if (now_int == min(key_int)) "new" else "changed"
  }
  out_chr
}

###### W_03_02: Which Solver Rows Differ, and Why ##############################
# Note: Independent of the panel's text. Each canonical row of C_02_01 is
#   mapped to the Equations-panel item it implements; the two models' rows
#   are compared (G0, G1 and Psi, to 1e-12). A row that differs belongs to
#   an item that must be changed, or whose display is common to both and
#   whose difference is carried entirely by composites defined in changed
#   descriptor items: rebuilding our row with their values of exactly those
#   composites (and of the wage mark-up's, for the wage equation) must
#   reproduce their row. An item whose rows are equal and is changed must
#   be a letters-only change. W_03_02_carried_lst names the composites each
#   shared-display item's row is built from, by d_lst name, and the
#   descriptor item that defines them.

W_03_02_row_map_lst <- list(
  "The Price Mark-Up"               = c(1L, 14L),
  "Capital in Use"                  = c(2L, 4L, 15L, 17L),
  "The Marginal Product of Capital" = c(3L, 16L),
  "Investment Equation"             = c(5L, 18L),
  "Tobin's q"                       = c(6L, 19L),
  "Consumption Euler Equation"      = c(7L, 20L),
  "The Resource Constraint"         = c(8L, 21L),
  "Aggregate Production Function"   = c(9L, 22L),
  "New Keynesian Phillips Curve (NKPC)" = 10L,
  "The Wage Equation"               = 11L,
  "The Wage Mark-Up"                = 23L,
  "Monetary Policy (MP) Rule"       = 12L,
  "Capital Accumulation"            = c(13L, 24L),
  "Five AR(1) Shock Processes"      = 25:29,
  "Two ARMA(1,1) Mark-Up Shocks"    = 30:33)

W_03_02_carried_lst <- list(
  "Investment Equation" = list(
    d = c("i1", "i2"), by = "Adjustment Costs in the Capital Weights"),
  "Tobin's q" = list(
    d = c("q1", "q1r"), by = "Utilisation and the q Weight"),
  "New Keynesian Phillips Curve (NKPC)" = list(
    d = c("p1", "p2", "p3"), by = "Calvo Prices in the Phillips Curve"),
  "The Wage Equation" = list(
    d = c("w1", "w2", "w3", "w4", "wc", "hg"),
    by = c("Calvo Wages in the Wage Equation", "The Wage Mark-Up")),
  "Capital Accumulation" = list(
    d = c("k1", "k2"), by = "Adjustment Costs in the Capital Weights"))

W_03_03_row_fn <- function(mat_lst, r_int) {
  c(mat_lst$g_zero[r_int, ], mat_lst$g_one[r_int, ], mat_lst$psi[r_int, ])
}

W_03_04_diff_fn <- function(par_lst) {
  ours_d <- C_01_03_derived_fn(utils::modifyList(par_lst,
                                                 list(model_chr = "ours")))
  sw_d   <- C_01_03_derived_fn(utils::modifyList(par_lst,
                                                 list(model_chr = "sw")))
  ours_m <- C_02_01_build_fn(ours_d)
  sw_m   <- C_02_01_build_fn(sw_d)
  lapply(stats::setNames(names(W_03_02_row_map_lst),
                         names(W_03_02_row_map_lst)),
         function(lab_chr) {
    rows_int <- W_03_02_row_map_lst[[lab_chr]]
    gap_num <- max(vapply(rows_int, function(r)
      max(abs(W_03_03_row_fn(ours_m, r) - W_03_03_row_fn(sw_m, r))),
      numeric(1)))
    carried_lgl <- NA
    if (!is.null(W_03_02_carried_lst[[lab_chr]])) {
      mix_d <- ours_d
      for (nm in W_03_02_carried_lst[[lab_chr]]$d) mix_d[[nm]] <- sw_d[[nm]]
      mix_m <- C_02_01_build_fn(mix_d)
      carried_lgl <- all(vapply(rows_int, function(r)
        max(abs(W_03_03_row_fn(mix_m, r) - W_03_03_row_fn(sw_m, r))) < 1e-12,
        logical(1)))
    }
    list(differs = gap_num > 1e-12, gap = gap_num, carried = carried_lgl)
  })
}

###### W_03_05: The Verdict on the Pills #######################################
# Note: Returns the problems found, empty if none. `letters_chr` is the
#   caveat's letters and potential-output items' pills, the only items that
#   may be changed with equal rows.

W_03_05_check_fn <- function(pill_chr, diff_lst, eq_lst, par_lst,
                             letters_chr) {
  bad_chr <- character(0)
  changed_chr <- names(pill_chr)[pill_chr == "changed"]
  for (lab_chr in names(diff_lst)) {
    d_lst <- diff_lst[[lab_chr]]
    is_changed <- lab_chr %in% changed_chr
    if (d_lst$differs && !is_changed) {
      c_lst <- W_03_02_carried_lst[[lab_chr]]
      if (is.null(c_lst) || !isTRUE(d_lst$carried) ||
          !all(c_lst$by %in% changed_chr)) {
        bad_chr <- c(bad_chr, paste(lab_chr, "differs but is not marked"))
      }
    }
    if (!d_lst$differs && is_changed && !lab_chr %in% letters_chr) {
      bad_chr <- c(bad_chr, paste(lab_chr,
                                  "is marked but its rows are equal"))
    }
  }
  # descriptor items: marked exactly when a composite they define differs
  # between the two solvers, or their letters do
  ours_d <- C_01_03_derived_fn(utils::modifyList(par_lst,
                                                 list(model_chr = "ours")))
  sw_d   <- C_01_03_derived_fn(utils::modifyList(par_lst,
                                                 list(model_chr = "sw")))
  desc_lst <- list(
    "Habit in the Consumption Weights" = c("c1", "c2", "c3"),
    "Adjustment Costs in the Capital Weights" = c("i1", "i2", "k1", "k2"),
    "Utilisation and the q Weight" = c("z1", "q1"),
    "Calvo Prices in the Phillips Curve" = c("p1", "p2", "p3"),
    "Calvo Wages in the Wage Equation" = c("w1", "w2", "w3", "w4"),
    "The Steady State Beneath the Composites" = c("crk", "cbetabar"),
    "Blanchard and Kahn" = character(0))
  for (lab_chr in names(desc_lst)) {
    num_lgl <- length(desc_lst[[lab_chr]]) > 0L &&
      any(vapply(desc_lst[[lab_chr]], function(nm)
        abs(ours_d[[nm]] - sw_d[[nm]]) > 1e-12, logical(1)))
    if (num_lgl != (lab_chr %in% changed_chr) &&
        !(lab_chr %in% changed_chr && lab_chr %in% letters_chr)) {
      bad_chr <- c(bad_chr, sprintf("%s: composites differ %s, marked %s",
                                    lab_chr, num_lgl, lab_chr %in% changed_chr))
    }
  }
  # new items: theirs has them, ours lacks them; the one today is trend
  # growth, absent from ours (gamma = 1) and present in theirs
  new_chr <- names(pill_chr)[pill_chr == "new"]
  if (!identical(new_chr, "Trend Growth") ||
      ours_d$cgamma != 1 || !(sw_d$cgamma > 1)) {
    bad_chr <- c(bad_chr, paste("NEW set is", paste(new_chr, collapse = ", ")))
  }
  bad_chr
}

#### W_04: The Caveat Against the Pills ########################################
# Note: Every caveat item names at least one pill, every named pill is
#   changed or new at stage 2.5e, and every changed or new pill is named by
#   some item. The four items are B_04_07's four, in its order.

W_04_01_check_fn <- function(caveat_lst, pill_chr) {
  bad_chr <- character(0)
  ids_chr <- vapply(caveat_lst, `[[`, "", "id")
  if (!identical(ids_chr, c("trend", "hours", "letters", "potential"))) {
    bad_chr <- c(bad_chr, paste("items are", paste(ids_chr, collapse = " ")))
  }
  named_chr <- unique(unlist(lapply(caveat_lst, `[[`, "pills")))
  for (it_lst in caveat_lst) {
    if (length(it_lst$pills) == 0L) {
      bad_chr <- c(bad_chr, paste(it_lst$id, "names no pill"))
    }
    off_chr <- setdiff(it_lst$pills, names(pill_chr))
    if (length(off_chr) > 0L) {
      bad_chr <- c(bad_chr, sprintf("%s names unmarked %s", it_lst$id,
                                    paste(off_chr, collapse = ", ")))
    }
  }
  orphan_chr <- setdiff(names(pill_chr), named_chr)
  if (length(orphan_chr) > 0L) {
    bad_chr <- c(bad_chr, paste("no caveat item for",
                                paste(orphan_chr, collapse = ", ")))
  }
  bad_chr
}

#--------------------------------- Script End ---------------------------------#
