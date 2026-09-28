################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Smets-Wouters Model: Verification                                          ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage, from the repo root:
##   Rscript tests/verify_against_sw2007.R
##   Two things are checked. The solver is not gensys (base R has no QZ
##   decomposition and the app must run under shinylive), so V_04 holds it
##   to impulse responses from a genuine ordered-QZ gensys and V_04b to
##   Dynare running the authors' own model file. And the model is not on
##   Whelan's slides in full (the composites are, the deep-parameter
##   mapping is not), so V_06 checks that it reproduces the results Whelan
##   states in class. V_02 to V_08, V_14_01 and V_14_03 check the
##   Smets-Wouters side (C_01_01_default_lst); V_16 evaluates each side's
##   displayed equations on that side's solver; V_18 holds our model's
##   composites to the gamma = 1, c_2 = 0 formulas; V_19 holds stage 2.5e's
##   pills to the equations that differ between the two solvers.
##
## Inputs:
##   R/model.R; app.R (from V_09 on) and R/toolkit.R; tests/two_models.R.
##   Fixtures in tests/: sw2007_gensys_irf.csv (ordered-QZ gensys impulse
##   responses at the posterior means), dynare_mode_irf.csv and
##   dynare_mode_cvd.csv (Dynare 6.0 impulse responses and conditional
##   variance decomposition at the authors' mode), sw2007_mode_hessian.csv
##   (the mode file's full-precision values); tests/dynare/ holds the
##   Dynare model files the fixtures came from.
##
## Outputs:
##   A pass/fail line per check and a non-zero exit status if any fails, so
##   this can gate a commit. Nothing is written to the app folder.
##
## Packages:
##   shiny, bslib, ggplot2, htmltools, grid.
##
## References:
##   Smets, F. and Wouters, R. (2007). Shocks and Frictions in US Business
##     Cycles: A Bayesian DSGE Approach. American Economic Review 97(3).
##   Whelan, K. MA Advanced Macroeconomics, part 11; cited as [W11 nn].

#-------------------------------- Script Begin --------------------------------#

options(scipen = 999, digits = 7)

# Note: A null device, so this script leaves nothing behind. ggplotGrob in
#   V_09 and print() in V_11 both need a device, and with none open R drops
#   an Rplots.pdf into the working directory. Every check that wants a real
#   device opens its own.

grDevices::pdf(NULL)

V_01_01_here_dir <- tryCatch({
  args_vec <- commandArgs(trailingOnly = FALSE)
  file_chr <- sub("^--file=", "", args_vec[grepl("^--file=", args_vec)])
  if (length(file_chr) == 1L && nzchar(file_chr)) {
    dirname(normalizePath(file_chr))
  } else "."
}, error = function(e) ".")

source(file.path(V_01_01_here_dir, "..", "R", "model.R"))

V_01_02_fail_int <- 0L

V_01_03_check_fn <- function(label_chr, ok_lgl, detail_chr = "") {
  message(sprintf("  [%s] %-54s %s",
                  if (isTRUE(ok_lgl)) "PASS" else "FAIL",
                  label_chr, detail_chr))
  if (!isTRUE(ok_lgl)) V_01_02_fail_int <<- V_01_02_fail_int + 1L
}

################################################################################
## V: Verification #############################################################
################################################################################
# Note: One section per property; V_01_03_check_fn prints and counts.

message("\nECON42240 Smets-Wouters app: verification\n")

# Note: The Smets-Wouters side. C_01_01_default_lst carries model_chr =
#   "sw"; our model is C_01_01b_ours_lst, the same deep parameters.
V_02_01_par_lst <- C_01_01_default_lst
V_02_01b_ours_lst <- C_01_01b_ours_lst
V_02_02b_ours_sol_lst <- C_03_02_solve_fn(V_02_01b_ours_lst)
V_02_02_sol_lst <- C_03_02_solve_fn(V_02_01_par_lst)

#### V_02: The Steady State ####################################################
# Note: If the expenditure shares do not sum to one then the composites are
#   built on a steady state that does not exist, and nothing downstream is
#   meaningful. This runs first for that reason.

message("V_02  The steady state (Smets-Wouters side)")

V_02_03_share_lst <- C_05_02_shares_fn(V_02_02_sol_lst$d_lst)
V_01_03_check_fn("expenditure shares sum to one",
                 V_02_03_share_lst$ok_lgl,
                 sprintf("c %.4f + i %.4f + g %.4f = %.10f",
                         V_02_03_share_lst$ccy_num, V_02_03_share_lst$ciy_num,
                         V_02_03_share_lst$gy_num, V_02_03_share_lst$sum_num))
V_01_03_check_fn("the rental rate of capital is positive and small",
                 V_02_02_sol_lst$d_lst$crk > 0 &&
                   V_02_02_sol_lst$d_lst$crk < 0.1,
                 sprintf("r^k = %.6f per quarter",
                         V_02_02_sol_lst$d_lst$crk))

#### V_03: Blanchard and Kahn ##################################################
# Note: The count has to balance exactly: one forecast error per unstable
#   root. With 45 equations this is a demanding test of the whole system, and
#   it is the first thing to look at if a change breaks something.

message("\nV_03  Blanchard and Kahn, counted")

V_01_03_check_fn("a unique stable solution exists",
                 isTRUE(V_02_02_sol_lst$ok_lgl),
                 if (length(V_02_02_sol_lst$problems_chr) == 0L) "" else
                   V_02_02_sol_lst$problems_chr[1L])
V_01_03_check_fn("unstable roots equal forecast errors",
                 V_02_02_sol_lst$nu_int == V_02_02_sol_lst$neta_int,
                 sprintf("%d unstable, %d forecast errors, %d stable",
                         V_02_02_sol_lst$nu_int, V_02_02_sol_lst$neta_int,
                         V_02_02_sol_lst$ns_int))
V_01_03_check_fn("the forecast errors span the unstable block",
                 V_02_02_sol_lst$rank_int == V_02_02_sol_lst$nu_int,
                 sprintf("rank %d of %d", V_02_02_sol_lst$rank_int,
                         V_02_02_sol_lst$nu_int))
V_01_03_check_fn("the stable basis is well conditioned",
                 V_02_02_sol_lst$cond_num < 1e06,
                 sprintf("condition number %.0f", V_02_02_sol_lst$cond_num))

#### V_04: Against a Genuine Ordered-QZ gensys #################################
# Note: The check that justifies the solver: 1120 impulse-response values
#   across all seven shocks and eight variables, computed by scipy's ordqz.
#   If the shifted-pencil method ever degrades (a repeated root, a defective
#   eigenvector) this is where it shows.

message("\nV_04  Against an ordered-QZ gensys")

V_04_01_ref_df <- utils::read.csv(
  file.path(V_01_01_here_dir, "sw2007_gensys_irf.csv"),
  stringsAsFactors = FALSE)

# Note: The fixture was computed in the code's sign and units for the risk
#   premium, where b_code = -c_3 varepsilon^b and a positive innovation
#   raises consumption. The solver writes the paper's eqs (2) and (4), in
#   which a positive varepsilon^b is a rise in the premium, with innovation
#   s.d. sigma_b / c_3 (R/model.R C_02_01, C_04_00). A one-s.d. response is
#   therefore the fixture's mirror image, and the eb rows are compared so.
V_04_01_ref_df$value[V_04_01_ref_df$shock == "eb"] <-
  -V_04_01_ref_df$value[V_04_01_ref_df$shock == "eb"]

# Note: The fixture is at the posterior means, built from the replication's
#   own equations by a separate code path, so it checks this file's
#   equation builder as well as its solver. The app runs at the posterior
#   mode (R/model.R C_01_01); V_05 checks that solution on its own path.
#   The means below are Table 1A/1B's mean column, [W11 12] and [W11 13],
#   and are used nowhere else.
V_04_00_mean_lst <- utils::modifyList(C_01_01_default_lst, list(
  csadjcost = 5.74, csigma = 1.38, chabb = 0.71, cprobw = 0.70,
  csigl = 1.83, cprobp = 0.66, cindw = 0.58, cindp = 0.24, czcap = 0.54,
  cfc = 1.60, crpi = 2.04, crr = 0.81, cry = 0.08, crdy = 0.22,
  constepinf = 0.78, constebeta = 0.16, constelab = 0.53, ctrend = 0.43,
  calfa = 0.19, crhoa = 0.95, crhob = 0.22, crhog = 0.97, crhoqs = 0.71,
  crhoms = 0.15, crhopinf = 0.89, crhow = 0.96, cmap = 0.69, cmaw = 0.84,
  cgy = 0.52, sda = 0.45, sdb = 0.23, sdg = 0.53, sdqs = 0.45,
  sdms = 0.24, sdpinf = 0.14, sdw = 0.24))
V_04_00_mean_sol_lst <- C_03_02_solve_fn(V_04_00_mean_lst)

V_04_02_worst_num <- 0
for (shock_chr in unique(V_04_01_ref_df$shock)) {
  irf_df <- C_04_01_irf_fn(V_04_00_mean_sol_lst, shock_chr, 20L)
  sub_df <- V_04_01_ref_df[V_04_01_ref_df$shock == shock_chr, ]
  mine_vec <- mapply(function(v_chr, q_int) irf_df[[v_chr]][q_int + 1L],
                     sub_df$variable, sub_df$quarter)
  V_04_02_worst_num <- max(V_04_02_worst_num,
                           max(abs(mine_vec - sub_df$value)))
}

V_01_03_check_fn("impulse responses match the QZ solution (at the means)",
                 V_04_02_worst_num < 1e-08,
                 sprintf("max |diff| = %.2e over %d values",
                         V_04_02_worst_num, nrow(V_04_01_ref_df)))

#### V_04b: Against Dynare, Running the Authors' Model at Their Mode ###########
# Note: tests/dynare_mode_irf.csv and dynare_mode_cvd.csv come from Dynare
#   6.0 under Octave, running Pfeifer's DSGE_mod copy of Smets and Wouters'
#   own model file with the parameters of their usmodel_mode.mat, i.e. at
#   this app's calibration: stoch_simul(order = 1, irf = 20,
#   conditional_variance_decomposition = [1 2 3 4 5 6 10 11 40 41 100 101]).
#   Dynare's risk premium b is the replication's sign, so its eb rows are
#   compared as their mirror image. The decomposition check also settles
#   the horizon convention: this file's horizon h is Dynare's horizon h.

message("\nV_04b Against Dynare at the authors' mode (Smets-Wouters side)")

V_04b_01_irf_df <- utils::read.csv(file.path(V_01_01_here_dir,
                                             "dynare_mode_irf.csv"))
V_04b_01_irf_df$value[V_04b_01_irf_df$shock == "eb"] <-
  -V_04b_01_irf_df$value[V_04b_01_irf_df$shock == "eb"]
# Note: Dynare ran at usmodel_mode.mat's full-precision values, which
#   tests/sw2007_mode_hessian.csv carries; the app's default rounds them to
#   eight decimals, so the comparison solves at the file's own values.
V_04b_00_h_df   <- utils::read.csv(file.path(V_01_01_here_dir,
                                             "sw2007_mode_hessian.csv"))
V_04b_00_sol_lst <- C_03_02_solve_fn(utils::modifyList(
  V_02_01_par_lst, as.list(stats::setNames(V_04b_00_h_df$mode,
                                           V_04b_00_h_df$name))))
V_04b_02_worst_num <- 0
for (shock_chr in unique(V_04b_01_irf_df$shock)) {
  irf_df <- C_04_01_irf_fn(V_04b_00_sol_lst, shock_chr, 20L)
  sub_df <- V_04b_01_irf_df[V_04b_01_irf_df$shock == shock_chr, ]
  mine_vec <- mapply(function(v_chr, q_int) irf_df[[v_chr]][q_int + 1L],
                     sub_df$variable, sub_df$quarter)
  V_04b_02_worst_num <- max(V_04b_02_worst_num,
                            max(abs(mine_vec - sub_df$value)))
}
V_01_03_check_fn("impulse responses match Dynare at the mode",
                 V_04b_02_worst_num < 1e-10,
                 sprintf("max |diff| = %.2e over %d values",
                         V_04b_02_worst_num, nrow(V_04b_01_irf_df)))

V_04b_03_cvd_df <- utils::read.csv(file.path(V_01_01_here_dir,
                                             "dynare_mode_cvd.csv"))
V_04b_04_h_vec  <- sort(unique(V_04b_03_cvd_df$horizon))
V_04b_05_f_lst  <- C_04_02_fevd_fn(V_04b_00_sol_lst, V_04b_04_h_vec)
V_04b_06_gap_num <- max(mapply(function(v_chr, h_int, s_chr, x_num)
  abs(V_04b_05_f_lst[[as.character(h_int)]][v_chr, s_chr] - x_num),
  V_04b_03_cvd_df$variable, V_04b_03_cvd_df$horizon, V_04b_03_cvd_df$shock,
  V_04b_03_cvd_df$share))
V_01_03_check_fn("variance shares match Dynare, horizon for horizon",
                 V_04b_06_gap_num < 1e-10,
                 sprintf("max |diff| = %.2e at horizons %s", V_04b_06_gap_num,
                         paste(V_04b_04_h_vec, collapse = " ")))

#### V_05: The Solution Satisfies the Model ####################################
# Note: Independent of the reference above: the solution is fed back into the
#   canonical system it came from, on its own path, and the E variables are
#   checked against the model's own forecasts.

message("\nV_05  The solution satisfies the model")

V_05_01_res_lst <- C_05_01_residual_fn(V_02_02_sol_lst)
V_01_03_check_fn("canonical system residual below 1e-09",
                 V_05_01_res_lst$system_num < 1e-09,
                 sprintf("%.2e", V_05_01_res_lst$system_num))
V_01_03_check_fn("rational expectations hold period by period",
                 V_05_01_res_lst$expect_num < 1e-09,
                 sprintf("%.2e", V_05_01_res_lst$expect_num))

#### V_06: What Whelan States in Class #########################################
# Note: [W11 15] to [W11 20]. These are the results the lecture turns on, and
#   they are the only part of this file that checks ECONOMICS rather than
#   arithmetic.

message("\nV_06  The results Whelan states in class (Smets-Wouters side)")

V_06_01_fevd_lst <- C_04_02_fevd_fn(V_02_02_sol_lst)
V_06_02_q1_mat   <- V_06_01_fevd_lst[["1"]]
V_06_03_q100_mat <- V_06_01_fevd_lst[["100"]]

# [W11 16] the price mark-up shock leads inflation on impact. His chart reads
#   about 72% (wage mark-up about 21%).
V_01_03_check_fn("[W11 16] price mark-up leads inflation at Q1",
                 which.max(V_06_02_q1_mat["pinf", ]) ==
                   match("epinf", C_01_02_shock_vec) &&
                   abs(V_06_02_q1_mat["pinf", "epinf"] - 0.72) < 0.08,
                 sprintf("%.1f%%, against about 72%% on the chart",
                         100 * V_06_02_q1_mat["pinf", "epinf"]))

# [W11 17] monetary policy is the biggest single contributor to the fed funds
#   rate at the shortest horizon
V_01_03_check_fn("[W11 17] monetary policy leads fed funds at Q1",
                 which.max(V_06_02_q1_mat["r", ]) ==
                   match("em", C_01_02_shock_vec),
                 sprintf("%.1f%%, the largest share",
                         100 * V_06_02_q1_mat["r", "em"]))

# [W11 15] demand shocks drive output at short horizons, mark-ups at long ones.
#   His Q1 bar: exogenous spending about 37%, risk premium about 27%.
V_06_02b_rank_int <- order(V_06_02_q1_mat["y", ], decreasing = TRUE)
V_01_03_check_fn("[W11 15] spending leads GDP at Q1, risk premium second",
                 identical(V_06_02b_rank_int[1:2],
                           match(c("eg", "eb"), C_01_02_shock_vec)),
                 sprintf("spending %.1f%%, risk premium %.1f%%",
                         100 * V_06_02_q1_mat["y", "eg"],
                         100 * V_06_02_q1_mat["y", "eb"]))
V_01_03_check_fn("[W11 15] wage mark-up leads GDP at Q100",
                 which.max(V_06_03_q100_mat["y", ]) ==
                   match("ew", C_01_02_shock_vec),
                 sprintf("%.1f%%, against %.1f%% at Q1",
                         100 * V_06_03_q100_mat["y", "ew"],
                         100 * V_06_02_q1_mat["y", "ew"]))

# [W11 19] a contractionary policy shock: output and inflation both fall, and
#   the output response is HUMP-SHAPED rather than worst on impact
V_06_04_mp_df <- C_04_01_irf_fn(V_02_02_sol_lst, "em", 25L)
V_06_05_trough_int <- which.min(V_06_04_mp_df$y) - 1L

V_01_03_check_fn("[W11 19] policy tightening lowers output and inflation",
                 min(V_06_04_mp_df$y) < 0 && min(V_06_04_mp_df$pinf) < 0,
                 sprintf("output %.3f, inflation %.3f at the trough",
                         min(V_06_04_mp_df$y), min(V_06_04_mp_df$pinf)))
V_01_03_check_fn("[W11 19] the output response is hump-shaped",
                 V_06_05_trough_int >= 2L,
                 sprintf("trough at quarter %d, not on impact",
                         V_06_05_trough_int))

# [W11 10] / [W7 38] the model puts hours DOWN after a positive technology
#   shock. The RBC of Part 7 puts them up, which is the failure Gali named.
#   This is the single cleanest link between the two apps.
V_06_06_a_df <- C_04_01_irf_fn(V_02_02_sol_lst, "ea", 25L)
V_01_03_check_fn("[W11 10] hours FALL on a technology shock, unlike the RBC",
                 V_06_06_a_df$lab[1L] < 0,
                 sprintf("%+.3f on impact, against the RBC's +1.70",
                         V_06_06_a_df$lab[1L]))

#### V_06b: The Decomposition Carries the Shock Sigmas #########################
# Note: The Q1 shares must weight each shock by its own sigma. Three checks.
#   (1) An exact identity: the one-step forecast error is the impact
#   response, so the Q1 share must equal C_04_01_irf_fn's squared impact
#   response (which scales by sigma on its own path), normalised. (2)
#   Against Whelan's charts read at 200 dpi: GDP Q1 spending ~37, risk
#   premium ~27, productivity ~16, investment ~13; fed funds Q1 monetary
#   policy ~58. (3) Productivity is not "about a tenth" of GDP at long
#   horizons.

message("\nV_06b The decomposition carries the shock sigmas")

V_06b_01_imp_mat <- vapply(C_01_02_shock_vec, function(sh_chr)
  unlist(C_04_01_irf_fn(V_02_02_sol_lst, sh_chr, 1L)[1L, c("y", "pinf", "r")]),
  numeric(3))
V_06b_02_want_mat <- V_06b_01_imp_mat^2 / rowSums(V_06b_01_imp_mat^2)
V_06b_03_gap_num <- max(abs(V_06b_02_want_mat -
                              V_06_02_q1_mat[c("y", "pinf", "r"),
                                             C_01_02_shock_vec]))
V_01_03_check_fn("Q1 shares equal the sigma-scaled impact responses",
                 V_06b_03_gap_num < 1e-10,
                 sprintf("worst gap %.2e", V_06b_03_gap_num))

V_06b_04_ref_vec <- c(eg = 0.37, eb = 0.27, ea = 0.16, eqs = 0.13)
V_06b_05_err_num <- max(abs(V_06_02_q1_mat["y", names(V_06b_04_ref_vec)] -
                              V_06b_04_ref_vec))
V_01_03_check_fn("[W11 15] GDP Q1 shares within 3 points of the chart",
                 V_06b_05_err_num < 0.03,
                 sprintf("worst |model - chart| = %.1f points",
                         100 * V_06b_05_err_num))
V_01_03_check_fn("[W11 17] monetary policy about 58% of fed funds at Q1",
                 abs(V_06_02_q1_mat["r", "em"] - 0.58) < 0.03,
                 sprintf("%.1f%%", 100 * V_06_02_q1_mat["r", "em"]))
V_01_03_check_fn(
  "[W11 15] productivity is a quarter of GDP by Q40, not a tenth",
                 V_06_01_fevd_lst[["40"]]["y", "ea"] > 0.2,
                 sprintf("%.1f%% at Q40",
                         100 * V_06_01_fevd_lst[["40"]]["y", "ea"]))

#### V_07: The Frictions Do What He Says #######################################
# Note: [W11 10] claims the frictions "throw sand in the wheels, making
#   variables more sluggish and giving random shocks a more long-lasting
#   effect". That is a falsifiable claim about this model, so it is tested
#   rather than repeated.

message("\nV_07  [W11 10]: the frictions, switched off")

V_07_01_base_int <- V_06_05_trough_int
V_07_02_base_num <- min(V_06_04_mp_df$y)

for (name_chr in c("habit", "adjcost")) {
  off_lst <- C_06_02_off_fn(V_02_01_par_lst, name_chr)
  if (!isTRUE(off_lst$ok_lgl)) {
    V_01_03_check_fn(sprintf("without %s the model still solves", name_chr),
                     FALSE, off_lst$problems_chr[1L])
    next
  }
  off_df   <- C_04_01_irf_fn(off_lst, "em", 25L)
  off_int  <- which.min(off_df$y) - 1L
  V_01_03_check_fn(
    sprintf("[W11 10] %s is what delays the peak",
            C_06_01_friction_lst[[name_chr]]$label_chr),
    off_int < V_07_01_base_int,
    sprintf("trough moves from quarter %d to %d",
            V_07_01_base_int, off_int))
}

V_07_03_all_lst <- C_06_03_all_off_fn(V_02_01_par_lst)
V_07_04_all_df  <- C_04_01_irf_fn(V_07_03_all_lst, "em", 25L)
V_01_03_check_fn("[W11 10] with the frictions off the shock hits harder",
                 min(V_07_04_all_df$y) < V_07_02_base_num,
                 sprintf("trough %.3f against %.3f with them",
                         min(V_07_04_all_df$y), V_07_02_base_num))
V_01_03_check_fn("[W11 10] with the frictions off the hump goes",
                 (which.min(V_07_04_all_df$y) - 1L) < V_07_01_base_int,
                 sprintf("trough at quarter %d against %d",
                         which.min(V_07_04_all_df$y) - 1L, V_07_01_base_int))

#### V_08: Neither Half Alone Flips the Hours Response #########################
# Note: The bridge figure's caption says neither the real frictions nor the
#   nominal rigidities flip the sign of the hours response on their own,
#   and only both together do. If a recalibration ever makes one of them
#   sufficient, the caption is wrong and this catches it.

message("\nV_08  What it takes to flip the hours response")

V_08_01_impact_fn <- function(over_lst) {
  sol_lst <- C_03_02_solve_fn(utils::modifyList(V_02_01_par_lst, over_lst))
  if (!isTRUE(sol_lst$ok_lgl)) return(NA_real_)
  C_04_01_irf_fn(sol_lst, "ea", 4L)$lab[1L]
}

V_08_02_base_num <- V_06_06_a_df$lab[1L]
V_08_03_real_num <- V_08_01_impact_fn(
  list(chabb = 1e-06, csadjcost = 0.1, czcap = 0.99,
       cindp = 1e-06, cindw = 1e-06))
V_08_04_nom_num  <- V_08_01_impact_fn(list(cprobp = 0.05, cprobw = 0.05))
V_08_05_both_num <- C_04_01_irf_fn(
  C_06_04_rbc_limit_fn(V_02_01_par_lst), "ea", 4L)$lab[1L]

V_01_03_check_fn("the real frictions alone do NOT flip it",
                 V_08_03_real_num < 0,
                 sprintf("%+.3f, still negative (from %+.3f)",
                         V_08_03_real_num, V_08_02_base_num))
V_01_03_check_fn("flexible prices and wages alone do NOT flip it",
                 V_08_04_nom_num < 0,
                 sprintf("%+.3f, still negative", V_08_04_nom_num))
V_01_03_check_fn("both together DO flip it, towards the RBC's sign",
                 V_08_05_both_num > 0,
                 sprintf("%+.3f, now positive", V_08_05_both_num))

#### V_09: The One Builder That Serves a Deck Figure ###########################
# Note: D_04_01_fevd_fn is the figure on deck 2.5's frame "Variance Shares
#   by Shock and Horizon", the only builder on this app with a deck slot.
#   Everything above checks the model; this section checks that what the
#   deck gets off the model matches the slot's own build spec. The builders
#   live in app.R, so R/ is sourced and app.R evaluated down to the user
#   interface and no further; nothing below starts Shiny. ggplot2 is
#   required here and nowhere else in the file: if it is missing the
#   section fails rather than passing quietly.

message("\nV_09  The deck's variance-shares figure")

V_09_01_ok_lgl <- requireNamespace("ggplot2", quietly = TRUE) &&
  requireNamespace("grid", quietly = TRUE)
V_01_03_check_fn("ggplot2 and grid are available to draw the figure",
                 V_09_01_ok_lgl,
                 if (V_09_01_ok_lgl) "" else "install ggplot2 to run V_09")

if (V_09_01_ok_lgl) {

  suppressMessages(library(ggplot2))

  V_09_02_env <- new.env(parent = globalenv())
  local({
    old_chr <- setwd(file.path(V_01_01_here_dir, ".."))
    on.exit(setwd(old_chr))
    for (file_chr in list.files("R", "\\.R$", full.names = TRUE)) {
      sys.source(file_chr, envir = V_09_02_env)
    }
    src_chr <- readLines("app.R", warn = FALSE)
    cut_int <- grep("^## E: User Interface", src_chr)[1L] - 1L
    eval(parse(text = paste(src_chr[seq_len(cut_int)], collapse = "\n")),
         envir = V_09_02_env)
  })

  # V_09_03: the presentation order is the model's order, rearranged
  # Note: B_03_04 is hand-written and a typo in it would silently paint one
  #   shock's share in another shock's colour under another shock's name.
  V_01_03_check_fn(
    "the deck's shock order is the model's seven, rearranged",
    setequal(V_09_02_env$B_03_04_deck_order_vec, C_01_02_shock_vec) &&
      length(V_09_02_env$B_03_04_deck_order_vec) ==
        length(C_01_02_shock_vec),
    paste(V_09_02_env$B_03_04_deck_order_vec, collapse = " "))

  V_09_04_p_obj  <- V_09_02_env$D_04_01_fevd_fn(
    V_02_02_sol_lst, "y", V_09_02_env$B_03_05_deck_horizon_vec)
  V_09_05_df     <- V_09_04_p_obj$data

  # V_09_06: four bars, at the horizons the slot names
  V_01_03_check_fn(
    "the deck figure has four bars, at Q1, Q4, Q10 and Q40",
    identical(levels(V_09_05_df$horizon_cat), c("Q1", "Q4", "Q10", "Q40")),
    paste(levels(V_09_05_df$horizon_cat), collapse = " "))

  # V_09_07: the stack is in the deck's order
  V_01_03_check_fn(
    "the segments stack in the deck's own shock order",
    identical(levels(V_09_05_df$shock_cat),
              unname(C_01_02_shock_lab_vec[
                V_09_02_env$B_03_04_deck_order_vec])),
    paste(levels(V_09_05_df$shock_cat), collapse = ", "))

  # V_09_08: the tablenote's own claim
  # Note: "Shares Within a Horizon Sum to One" is what the frame's tablenote
  #   says, so it is gated here rather than trusted.
  V_09_09_sum_vec <- vapply(split(V_09_05_df$share_num,
                                  V_09_05_df$horizon_cat), sum, numeric(1))
  V_01_03_check_fn(
    "shares within each horizon sum to one",
    all(abs(V_09_09_sum_vec - 1) < 1e-09),
    sprintf("worst |sum - 1| = %.2e", max(abs(V_09_09_sum_vec - 1))))

  # V_09_10: every shock is named somewhere, which is why no legend
  # Note: The property that makes the legend removable. A segment is
  #   labelled only where it is wide enough, so a recalibration that
  #   squeezed one shock below the threshold at every one of the four
  #   horizons would leave that shock unnamed on the slide.
  V_09_11_named_chr <- unique(V_09_05_df$label_chr[
    nzchar(V_09_05_df$label_chr)])
  V_09_12_missing_chr <- setdiff(
    gsub("\n", " ", unname(C_01_02_shock_lab_vec[
      V_09_02_env$B_03_04_deck_order_vec])),
    gsub("\n", " ", V_09_11_named_chr))
  V_01_03_check_fn(
    "every shock is named on at least one bar",
    length(V_09_12_missing_chr) == 0L,
    if (length(V_09_12_missing_chr) == 0L) "all seven named" else
      paste("unnamed:", paste(V_09_12_missing_chr, collapse = ", ")))

  # V_09_13: no legend, which CONVENTIONS.md 6 forbids
  # Note: The layout cells for a guide box exist whether or not a legend is
  #   drawn, so their presence proves nothing; their size does, since an
  #   empty guide box measures zero. ggplot2 >= 3.5.0 lays out four cells
  #   named guide-box-right/left/bottom/top (zeroGrob when empty); 3.4.x has
  #   one cell, "guide-box", only when a legend is drawn, and a right-hand
  #   legend there sits in the panel's own row, whose height converts to
  #   0 mm. So the grob itself is sized, width times height of whatever
  #   sits in any cell called guide-box*, and V_15 proves the measurement
  #   comes back non-zero for a plot that does carry a legend.
  V_09_13b_guide_fn <- function(p_obj) {
    g_obj <- ggplot2::ggplotGrob(p_obj)
    hit_int <- which(grepl("^guide-box", g_obj$layout$name))
    if (length(hit_int) == 0L) return(0)
    sum(vapply(hit_int, function(i_int) {
      gb_obj <- g_obj$grobs[[i_int]]
      if (!inherits(gb_obj, "gtable")) return(0)
      sum(grid::convertWidth(sum(gb_obj$widths), "mm", valueOnly = TRUE)) *
        sum(grid::convertHeight(sum(gb_obj$heights), "mm", valueOnly = TRUE))
    }, numeric(1)))
  }
  V_09_14_guide_num <- V_09_13b_guide_fn(V_09_04_p_obj)
  V_01_03_check_fn(
    "the figure carries no legend",
    V_09_14_guide_num < 1e-06,
    sprintf("guide box measures %.2f mm^2", V_09_14_guide_num))

  # V_09_15: one shock, one colour, across every figure that draws it
  # Note: B_03_02's own note is the claim being tested: the colours are
  #   fixed "so the stacked bars in one figure and the lines in another
  #   name the same thing the same way". A character shock column would let
  #   the discrete scale sort the levels alphabetically while the colours
  #   arrive in the model's order, so the check is on the built plot, not
  #   the constant. The zero line is a layer too, with a colour and a group
  #   of -1, so the layer wanted is the one whose groups cover the shock
  #   levels. levels(factor(x)) rather than levels(x), because a character
  #   column is precisely the fault this check exists to catch, and
  #   levels() of it returns NULL.
  V_09_16_drawn_fn <- function(p_obj, aes_chr) {
    b_lst   <- ggplot2::ggplot_build(p_obj)
    lvl_chr <- levels(factor(b_lst$plot$data$shock_cat))
    hit_int <- which(vapply(b_lst$data, function(d_df)
      aes_chr %in% names(d_df) && "group" %in% names(d_df) &&
        all(seq_along(lvl_chr) %in% d_df$group), logical(1)))
    if (length(hit_int) == 0L) return(stats::setNames(
      rep(NA_character_, length(lvl_chr)), lvl_chr))
    lay_df <- b_lst$data[[hit_int[1L]]]
    stats::setNames(vapply(seq_along(lvl_chr), function(i_int)
      as.character(lay_df[[aes_chr]][match(i_int, lay_df$group)]),
      character(1)), lvl_chr)
  }

  V_09_17_want_vec <- stats::setNames(
    unname(V_09_02_env$B_03_02_shock_col_vec[C_01_02_shock_vec]),
    unname(C_01_02_shock_lab_vec[C_01_02_shock_vec]))

  V_09_18_bad_chr <- character(0)
  for (fig_lst in list(
    list(name = "D_04_01 stacked shares", aes = "fill",
         plot = V_09_04_p_obj),
    list(name = "D_03_02 demand shocks",  aes = "colour",
         plot = V_09_02_env$D_03_02_demand_fn(V_02_02_sol_lst)))) {
    got_vec <- V_09_16_drawn_fn(fig_lst$plot, fig_lst$aes)
    # Note: An empty reading is a failure, not a pass. A check that finds no
    #   series to compare has not compared anything.
    wrong_chr <- if (length(got_vec) == 0L) "no series read" else
      names(got_vec)[
        is.na(got_vec) |
          toupper(got_vec) != toupper(V_09_17_want_vec[names(got_vec)])]
    if (length(wrong_chr) > 0L) {
      V_09_18_bad_chr <- c(V_09_18_bad_chr,
                           sprintf("%s: %s", fig_lst$name,
                                   paste(wrong_chr, collapse = ", ")))
    }
  }

  V_01_03_check_fn(
    "each shock keeps its B_03_02 colour in every figure",
    length(V_09_18_bad_chr) == 0L,
    if (length(V_09_18_bad_chr) == 0L)
      "D_04_01 and D_03_02 agree" else
      paste(V_09_18_bad_chr, collapse = "; "))
}

#### V_10: The Guard on a Picker's Input #######################################
# Note: On the first paint of a stage a plot output can run before the
#   browser has sent a radio's value, so input$shock or input$target is
#   NULL; without req() the builders throw (C_04_01_irf_fn subsets its sd
#   lookup by a zero-length index; the FEVD row slice for a NULL target is
#   0 x 7). This section proves the guard by driving the app: with the
#   stage set and the input unset the output must raise Shiny's own silent
#   exception, the one req() throws, and with the input set it must return
#   a PNG. One testServer run serves the whole section: the solve lives in
#   F_02_02_sol_r, a reactive over the sliders alone, so moving the stage
#   or the shock does not re-solve. V_10_07 reads app.R and asks the
#   general question, so a new unguarded output fails the suite.

message("\nV_10  The guarded outputs, driven headlessly")

# Note: shiny is required here and not above. As with ggplot2 in V_09, a
#   missing package fails the section rather than skipping it.
V_10_00_ok_lgl <- requireNamespace("shiny", quietly = TRUE)
V_01_03_check_fn("shiny is available to drive the app",
                 V_10_00_ok_lgl,
                 if (V_10_00_ok_lgl) "" else "install shiny to run V_10")
suppressMessages(library(shiny))

V_10_01_verdict_fn <- function(expr_fn) {
  tryCatch({
    val_lst <- expr_fn()
    if (is.list(val_lst) && !is.null(val_lst$src) && nzchar(val_lst$src)) {
      "drew"
    } else "empty"
  }, error = function(e_cnd) {
    if (inherits(e_cnd, "shiny.silent.error")) "guarded" else
      paste0("THREW: ", conditionMessage(e_cnd))
  })
}

V_10_02_seen_lst <- local({
  out_lst <- list()
  app_dir <- file.path(V_01_01_here_dir, "..")
  shiny::testServer(app_dir, {

    # stage 2.5c: the shock radio has not reported yet
    session$setInputs(stage = "2.5c")
    out_lst$irf_unset_chr <<- V_10_01_verdict_fn(function() output$plot_irf)
    out_lst$shock_null_lgl <<- is.null(input$shock)
    out_lst$irf_nom_unset_chr <<- V_10_01_verdict_fn(
      function() output$plot_irf_nom)
    session$setInputs(shock = "em")
    out_lst$irf_set_chr <<- V_10_01_verdict_fn(function() output$plot_irf)
    out_lst$irf_nom_set_chr <<- V_10_01_verdict_fn(
      function() output$plot_irf_nom)
    out_lst$demand_chr <<- V_10_01_verdict_fn(function() output$plot_demand)
    out_lst$demand_nom_chr <<- V_10_01_verdict_fn(
      function() output$plot_demand_nom)

    # stage 2.5d: the same for the target radio
    session$setInputs(stage = "2.5d")
    out_lst$fevd_unset_chr <<- local({
      # input$target has to be cleared, because testServer keeps an input
      # once set and the browser does not send one until the radio is
      # painted; NULL is the first-paint state being reproduced
      session$setInputs(target = NULL)
      V_10_01_verdict_fn(function() output$plot_fevd)
    })
    session$setInputs(target = "y")
    out_lst$fevd_set_chr <<- V_10_01_verdict_fn(function() output$plot_fevd)

    # the other four figures, and the three text panels, at every stage
    session$setInputs(stage = "2.5a")
    out_lst$friction_chr <<- V_10_01_verdict_fn(function() output$plot_friction)
    out_lst$ladder_chr   <<- V_10_01_verdict_fn(function() output$plot_ladder)
    out_lst$bridge_chr   <<- V_10_01_verdict_fn(function() output$plot_bridge)

    # stage 2.5b: our model assembled, its own shock picker
    session$setInputs(stage = "2.5b")
    out_lst$s5_unset_chr <<- V_10_01_verdict_fn(function() output$plot_s2_irf)
    session$setInputs(shock2 = "em")
    out_lst$s5_set_chr <<- c(
      V_10_01_verdict_fn(function() output$plot_s2_irf),
      V_10_01_verdict_fn(function() output$plot_s2_fevd))

    # stage 2.5e: its own two pickers, guarded the same way
    session$setInputs(stage = "2.5e")
    out_lst$vs_unset_chr <<- c(
      V_10_01_verdict_fn(function() output$plot_vs_irf),
      V_10_01_verdict_fn(function() output$plot_vs_fevd))
    session$setInputs(shock5 = "em", target5 = "y")
    out_lst$vs_set_chr <<- vapply(c("plot_vs_irf", "plot_vs_irf_nom",
                                    "plot_vs_fevd", "plot_vs_fevd_long"),
      function(id_chr) V_10_01_verdict_fn(function() output[[id_chr]]), "")

    out_lst$panels_lgl <<- vapply(c("2.5a", "2.5b", "2.5c", "2.5d", "2.5e"),
      function(st_chr) {
        session$setInputs(stage = st_chr)
        all(vapply(list(function() output$eq_model,
                        function() output$eq_notation,
                        function() output$eq_explain,
                        function() output$preset_title,
                        function() output$scenario_story,
                        function() output$prompt,
                        function() output$test_ui),
                   function(f_fn) tryCatch({
                     h_chr <- as.character(f_fn())
                     length(h_chr) > 0L && any(nzchar(h_chr))
                   }, error = function(e) FALSE), logical(1)))
      }, logical(1))
  })
  out_lst
})

V_01_03_check_fn(
  "2.5c with the shock unset is GUARDED, not thrown",
  identical(V_10_02_seen_lst$irf_unset_chr, "guarded") &&
    identical(V_10_02_seen_lst$irf_nom_unset_chr, "guarded") &&
    isTRUE(V_10_02_seen_lst$shock_null_lgl),
  sprintf("input$shock NULL, outputs %s and %s",
          V_10_02_seen_lst$irf_unset_chr, V_10_02_seen_lst$irf_nom_unset_chr))
V_01_03_check_fn(
  "2.5c with the shock set draws both halves of the four panels",
  identical(V_10_02_seen_lst$irf_set_chr, "drew") &&
    identical(V_10_02_seen_lst$irf_nom_set_chr, "drew"),
  sprintf("real %s, nominal %s", V_10_02_seen_lst$irf_set_chr,
          V_10_02_seen_lst$irf_nom_set_chr))
V_01_03_check_fn(
  "2.5d with the target unset is GUARDED, not thrown",
  identical(V_10_02_seen_lst$fevd_unset_chr, "guarded"),
  V_10_02_seen_lst$fevd_unset_chr)
V_01_03_check_fn(
  "2.5d with the target set draws the decomposition",
  identical(V_10_02_seen_lst$fevd_set_chr, "drew"),
  V_10_02_seen_lst$fevd_set_chr)
V_01_03_check_fn(
  "the five unguarded figures still draw",
  all(c(V_10_02_seen_lst$friction_chr, V_10_02_seen_lst$ladder_chr,
        V_10_02_seen_lst$demand_chr, V_10_02_seen_lst$demand_nom_chr,
        V_10_02_seen_lst$bridge_chr) == "drew"),
  sprintf("friction %s, ladder %s, demand %s and %s, bridge %s",
          V_10_02_seen_lst$friction_chr, V_10_02_seen_lst$ladder_chr,
          V_10_02_seen_lst$demand_chr, V_10_02_seen_lst$demand_nom_chr,
          V_10_02_seen_lst$bridge_chr))
V_01_03_check_fn(
  "2.5b with the shock unset is GUARDED; set, both figures draw",
  identical(V_10_02_seen_lst$s5_unset_chr, "guarded") &&
    all(V_10_02_seen_lst$s5_set_chr == "drew"),
  sprintf("unset %s; set %s", V_10_02_seen_lst$s5_unset_chr,
          paste(V_10_02_seen_lst$s5_set_chr, collapse = " ")))
V_01_03_check_fn(
  "2.5e with its pickers unset is GUARDED; set, all four figures draw",
  all(V_10_02_seen_lst$vs_unset_chr == "guarded") &&
    all(V_10_02_seen_lst$vs_set_chr == "drew"),
  sprintf("unset %s; set %s",
          paste(V_10_02_seen_lst$vs_unset_chr, collapse = " "),
          paste(V_10_02_seen_lst$vs_set_chr, collapse = " ")))
V_01_03_check_fn(
  "equations, notation, in-words, presets and tests build at every stage",
  all(V_10_02_seen_lst$panels_lgl),
  sprintf("%d of 5 stages", sum(V_10_02_seen_lst$panels_lgl)))

# V_10_03: the target really does select
# Note: The guard stops the throw; this checks the thing the throw was
#   hiding, that two targets give two different decompositions.

if (V_09_01_ok_lgl) {
  V_10_04_y_df    <- V_09_02_env$D_04_01_fevd_fn(V_02_02_sol_lst, "y")$data
  V_10_05_pinf_df <- V_09_02_env$D_04_01_fevd_fn(V_02_02_sol_lst, "pinf")$data
  V_10_06_lead_fn <- function(df_in) as.character(
    df_in$shock_cat[df_in$horizon_cat == "Q1"][
      which.max(df_in$share_num[df_in$horizon_cat == "Q1"])])
  V_01_03_check_fn(
    "the FEVD target selects: output and inflation differ at Q1",
    !isTRUE(all.equal(V_10_04_y_df$share_num, V_10_05_pinf_df$share_num)) &&
      V_10_06_lead_fn(V_10_04_y_df) != V_10_06_lead_fn(V_10_05_pinf_df),
    sprintf("output leads on %s, inflation on %s",
            V_10_06_lead_fn(V_10_04_y_df), V_10_06_lead_fn(V_10_05_pinf_df)))
}

# V_10_07: the sweep, so a new unguarded output cannot ship
# Note: Reads app.R and asks the general question. Every input id created
#   by a widget inside a renderUI is collected, and every main-window
#   picker (E_02_03_*) with them; every output block is then read for the
#   input ids it uses and the ones it req()s. An output that reads an id
#   from the first set without req()ing it is the fault, wherever it is.
#   Literal ids only: the preset buttons are built with paste0() and are
#   read by observeEvent(ignoreInit = TRUE), which is not exposed to the
#   fault. The count of ids collected is printed.

V_10_08_widget_chr <- c("radioButtons", "selectInput", "selectizeInput",
                        "sliderInput", "checkboxInput", "checkboxGroupInput",
                        "numericInput", "textInput", "actionButton",
                        "dateInput", "varSelectInput")

V_10_09_walk_fn <- function(node, fn) {
  fn(node)
  if (is.call(node) || is.expression(node)) {
    for (i_int in seq_along(node)) {
      if (!is.null(node[[i_int]])) {
        tryCatch(V_10_09_walk_fn(node[[i_int]], fn), error = function(e) NULL)
      }
    }
  }
  invisible(NULL)
}

V_10_10_collect_fn <- function(node, want_chr) {
  hit_chr <- character(0)
  V_10_09_walk_fn(node, function(n_obj) {
    if (is.call(n_obj) && is.name(n_obj[[1L]]) &&
        as.character(n_obj[[1L]]) %in% want_chr && length(n_obj) > 1L &&
        is.character(n_obj[[2L]])) {
      hit_chr <<- c(hit_chr, n_obj[[2L]])
    }
  })
  unique(hit_chr)
}

V_10_11_dollar_fn <- function(node, obj_chr) {
  hit_chr <- character(0)
  V_10_09_walk_fn(node, function(n_obj) {
    if (is.call(n_obj) && identical(n_obj[[1L]], as.name("$")) &&
        length(n_obj) == 3L && identical(n_obj[[2L]], as.name(obj_chr)) &&
        is.name(n_obj[[3L]])) {
      hit_chr <<- c(hit_chr, as.character(n_obj[[3L]]))
    }
  })
  unique(hit_chr)
}

V_10_12_src_expr <- parse(file.path(V_01_01_here_dir, "..", "app.R"))

# every input id a renderUI creates, and every main-window picker
V_10_13_dynamic_chr <- character(0)
V_10_09_walk_fn(V_10_12_src_expr, function(n_obj) {
  if (is.call(n_obj) && is.name(n_obj[[1L]]) &&
      as.character(n_obj[[1L]]) == "renderUI") {
    V_10_13_dynamic_chr <<- c(V_10_13_dynamic_chr,
                              V_10_10_collect_fn(n_obj, V_10_08_widget_chr))
  }
  if (is.call(n_obj) && length(n_obj) == 3L &&
      as.character(n_obj[[1L]]) %in% c("<-", "=") &&
      is.name(n_obj[[2L]]) &&
      grepl("^E_02_03_", as.character(n_obj[[2L]]))) {
    V_10_13_dynamic_chr <<- c(V_10_13_dynamic_chr,
                              V_10_10_collect_fn(n_obj, V_10_08_widget_chr))
  }
})
V_10_13_dynamic_chr <- unique(V_10_13_dynamic_chr)

# every output block, with the inputs it reads and the ones it req()s
V_10_14_unguarded_chr <- character(0)
V_10_09_walk_fn(V_10_12_src_expr, function(n_obj) {
  if (is.call(n_obj) && length(n_obj) == 3L &&
      as.character(n_obj[[1L]]) %in% c("<-", "=") &&
      is.call(n_obj[[2L]]) && identical(n_obj[[2L]][[1L]], as.name("$")) &&
      identical(n_obj[[2L]][[2L]], as.name("output"))) {
    name_chr <- as.character(n_obj[[2L]][[3L]])
    used_chr <- V_10_11_dollar_fn(n_obj[[3L]], "input")
    req_chr  <- character(0)
    V_10_09_walk_fn(n_obj[[3L]], function(m_obj) {
      if (is.call(m_obj) && is.name(m_obj[[1L]]) &&
          as.character(m_obj[[1L]]) %in% c("req", "validate", "need")) {
        req_chr <<- c(req_chr, V_10_11_dollar_fn(m_obj, "input"))
      }
    })
    bad_chr <- setdiff(intersect(used_chr, V_10_13_dynamic_chr), req_chr)
    if (length(bad_chr) > 0L) {
      V_10_14_unguarded_chr <<- c(V_10_14_unguarded_chr,
        sprintf("output$%s reads %s", name_chr,
                paste(bad_chr, collapse = ", ")))
    }
  }
})

V_01_03_check_fn(
  "no output reads a picker's input without req()",
  length(V_10_14_unguarded_chr) == 0L &&
    all(c("shock", "target", "shock2", "shock5", "target5") %in%
          V_10_13_dynamic_chr),
  if (length(V_10_14_unguarded_chr) == 0L)
    sprintf("%d dynamic id(s) swept: %s", length(V_10_13_dynamic_chr),
            paste(V_10_13_dynamic_chr, collapse = " ")) else
    paste(V_10_14_unguarded_chr, collapse = "; "))

#### V_11: No Figure Squarer Than 3:2 ##########################################
# Note: The rule is about the exported figure, so it is read off the PNG
#   and not on the code. A panel squarer than 3:2 is height-limited inside
#   the deck's own box (beamer 16:9 at 10pt gives \linewidth 398.3pt and
#   \textheight 252.1pt) and shrinks rather than filling the column.
#   D_03_01 and D_03_02 facet four panels at nrow = 2, so each panel comes
#   out near 1:1 inside a wide image; the whole PNG is what is asserted.
#   Bytes 17 to 24 of a PNG are the IHDR width and height, big-endian,
#   which needs no package.

message("\nV_11  Every figure at 3:2 or wider")

V_11_01_dim_fn <- function(path_chr) {
  raw_vec <- readBin(path_chr, "raw", 33L)
  be_fn <- function(i_int) sum(as.integer(raw_vec[i_int]) * 256^(3:0))
  c(w = be_fn(17:20), h = be_fn(21:24))
}

# V_11_02: the size constants themselves
V_11_02_bad_chr <- names(Filter(function(s_lst)
  s_lst$width_mm / s_lst$height_mm < 1.5, V_09_02_env$B_03_03_size_lst))
V_01_03_check_fn(
  "every entry in B_03_03_size_lst is 3:2 or wider",
  length(V_11_02_bad_chr) == 0L,
  paste(vapply(names(V_09_02_env$B_03_03_size_lst), function(nm_chr) {
    s_lst <- V_09_02_env$B_03_03_size_lst[[nm_chr]]
    sprintf("%s %.2f:1", nm_chr, s_lst$width_mm / s_lst$height_mm)
  }, character(1)), collapse = ", "))

V_01_03_check_fn(
  "the export spec is 3:2 or wider",
  V_09_02_env$B_03_06_export_lst$width_px /
    V_09_02_env$B_03_06_export_lst$height_px >= 1.5,
  sprintf("%d x %d px, %.2f:1",
          V_09_02_env$B_03_06_export_lst$width_px,
          V_09_02_env$B_03_06_export_lst$height_px,
          V_09_02_env$B_03_06_export_lst$width_px /
            V_09_02_env$B_03_06_export_lst$height_px))

# V_11_03: nothing forces a square panel from inside a builder
# Note: aspect.ratio and coord_fixed both override whatever size a figure
#   is given, so a size constant cannot fix them.
V_11_03_src_chr <- readLines(file.path(V_01_01_here_dir, "..", "app.R"),
                             warn = FALSE)
V_11_04_forced_int <- grep("aspect\\.ratio|coord_fixed", V_11_03_src_chr)
V_01_03_check_fn(
  "no builder forces a panel shape with aspect.ratio or coord_fixed",
  length(V_11_04_forced_int) == 0L,
  if (length(V_11_04_forced_int) == 0L) "none found" else
    paste("app.R line(s)", paste(V_11_04_forced_int, collapse = " ")))

# V_11_05: every builder, rendered at the deck's own sizes
# Note: The solved model is passed in, not re-solved. D_02_02 and D_05_01
#   solve again by design, as counterfactual figures, and that is the only
#   solving this section does.

V_11_06_builder_lst <- list(
  list(name = "D_02_01 friction",   size = "wide",
       fn = function() V_09_02_env$D_02_01_friction_fn(V_02_02_sol_lst)),
  list(name = "D_02_02 ladder",     size = "wide",
       fn = function() V_09_02_env$D_02_02_ladder_fn(V_02_01_par_lst)),
  list(name = "D_03_01 responses",  size = "panels",
       fn = function() V_09_02_env$D_03_01_irf_fn(V_02_02_sol_lst, "em")),
  list(name = "D_03_02 demand",     size = "panels",
       fn = function() V_09_02_env$D_03_02_demand_fn(V_02_02_sol_lst)),
  list(name = "D_04_01 app shares", size = "wide",
       fn = function() V_09_02_env$D_04_01_fevd_fn(V_02_02_sol_lst, "y")),
  list(name = "D_04_01 deck shares", size = "wide",
       fn = function() V_09_02_env$D_04_01_fevd_fn(
         V_02_02_sol_lst, "y", V_09_02_env$B_03_05_deck_horizon_vec)),
  list(name = "D_05_01 bridge",     size = "wide",
       fn = function() V_09_02_env$D_05_01_bridge_fn(
         V_02_02_sol_lst, V_02_01_par_lst)),
  list(name = "D_06_01 ours vs sw",  size = "panels",
       fn = function() V_09_02_env$D_06_01_vs_irf_fn(
         V_02_02b_ours_sol_lst, V_02_02_sol_lst, "em")),
  list(name = "D_06_02 ours vs sw shares", size = "wide",
       fn = function() V_09_02_env$D_06_02_vs_fevd_fn(
         V_02_02b_ours_sol_lst, V_02_02_sol_lst, "y")),
  list(name = "D_09_01 fallback",   size = "wide",
       fn = function() V_09_02_env$D_09_01_blank_fn("No solution here.")))

V_11_07_dir_chr <- file.path(tempdir(), "sw-ratio-check")
dir.create(V_11_07_dir_chr, showWarnings = FALSE, recursive = TRUE)

V_11_08_report_chr <- character(0)
V_11_09_bad_chr    <- character(0)
for (b_lst in V_11_06_builder_lst) {
  size_lst <- V_09_02_env$B_03_03_size_lst[[b_lst$size]]
  path_chr <- file.path(V_11_07_dir_chr,
                        paste0(gsub("[^A-Za-z0-9]+", "-", b_lst$name), ".png"))
  ok_lgl <- tryCatch({
    grDevices::png(path_chr,
                   width  = size_lst$width_mm,  height = size_lst$height_mm,
                   units  = "mm", res = 200)
    print(b_lst$fn())
    grDevices::dev.off()
    TRUE
  }, error = function(e_cnd) {
    try(grDevices::dev.off(), silent = TRUE)
    V_11_09_bad_chr <<- c(V_11_09_bad_chr,
                          sprintf("%s did not render: %s", b_lst$name,
                                  conditionMessage(e_cnd)))
    FALSE
  })
  if (!ok_lgl) next
  dim_int <- V_11_01_dim_fn(path_chr)
  ratio_num <- dim_int[["w"]] / dim_int[["h"]]
  V_11_08_report_chr <- c(V_11_08_report_chr,
                          sprintf("%s %.2f", b_lst$name, ratio_num))
  if (ratio_num < 1.5) {
    V_11_09_bad_chr <- c(V_11_09_bad_chr,
                         sprintf("%s is %.2f:1", b_lst$name, ratio_num))
  }
  if (file.size(path_chr) <= 0) {
    V_11_09_bad_chr <- c(V_11_09_bad_chr,
                         sprintf("%s wrote an empty file", b_lst$name))
  }
}

V_01_03_check_fn(
  "every rendered builder is 3:2 or wider",
  length(V_11_09_bad_chr) == 0L,
  if (length(V_11_09_bad_chr) == 0L)
    paste(V_11_08_report_chr, collapse = ", ") else
    paste(V_11_09_bad_chr, collapse = "; "))

#### V_12: The Save PNG Buttons ################################################
# Note: One button per figure, driven through testServer, and what it
#   writes is opened and checked: a non-empty file, the deck's declared
#   pixel size, and the name it promised. A handler that silently wrote the
#   browser's size rather than the deck's is exactly the failure this
#   catches, and it is invisible from the code. One testServer run again,
#   for the same reason as V_10.

message("\nV_12  The Save PNG buttons")

# Note: The file is copied out inside the testServer block. A download
#   handler writes into a temporary directory that testServer removes when
#   the block ends, so a path read out of it points at nothing by the time
#   it is checked. Keeping the copy also keeps the name the handler chose.

V_12_00_dir_chr <- file.path(tempdir(), "sw-download-check")
unlink(V_12_00_dir_chr, recursive = TRUE)
dir.create(V_12_00_dir_chr, showWarnings = FALSE, recursive = TRUE)

# Note: The handlers are the toolkit's, T_07_07h: <id>__png and <id>__pdf
#   for every id in B_03_07. Each is driven at its own stage, with the
#   picker it reads set, through V_12_01b_drive_fn, which V_15_11 reuses
#   for the PDFs.

V_12_01a_input_lst <- list(
  plot_friction   = list(stage = "2.5a"),
  plot_ladder     = list(stage = "2.5a"),
  plot_irf        = list(stage = "2.5c", shock = "em"),
  plot_irf_nom    = list(stage = "2.5c", shock = "em"),
  plot_demand     = list(stage = "2.5c"),
  plot_demand_nom = list(stage = "2.5c"),
  plot_fevd       = list(stage = "2.5d", target = "y"),
  plot_bridge     = list(stage = "2.5a"),
  plot_s2_irf     = list(stage = "2.5b", shock2 = "em"),
  plot_s2_fevd    = list(stage = "2.5b"),
  plot_vs_irf     = list(stage = "2.5e", shock5 = "em"),
  plot_vs_irf_nom = list(stage = "2.5e", shock5 = "em"),
  plot_vs_fevd    = list(stage = "2.5e", target5 = "y"),
  plot_vs_fevd_long = list(stage = "2.5e", target5 = "y"))

V_12_01b_drive_fn <- function(ext_chr) {
  out_lst <- list()
  app_dir <- file.path(V_01_01_here_dir, "..")
  keep_fn <- function(path_chr) {
    if (is.null(path_chr) || !file.exists(path_chr)) return(NULL)
    out_chr <- file.path(V_12_00_dir_chr, basename(path_chr))
    file.copy(path_chr, out_chr, overwrite = TRUE)
    out_chr
  }
  shiny::testServer(app_dir, {
    for (id_chr in names(V_12_01a_input_lst)) {
      do.call(session$setInputs, V_12_01a_input_lst[[id_chr]])
      out_lst[[id_chr]] <<- tryCatch(
        keep_fn(output[[paste0(id_chr, "__", ext_chr)]]),
        error = function(e) NULL)
    }
  })
  out_lst
}

V_12_01_file_lst <- V_12_01b_drive_fn("png")

V_12_02_want_int <- c(V_09_02_env$B_03_06_export_lst$width_px,
                      V_09_02_env$B_03_06_export_lst$height_px)
V_12_03_bad_chr  <- character(0)
V_12_04_seen_chr <- character(0)

for (id_chr in names(V_09_02_env$B_03_07_figure_lst)) {
  path_chr <- V_12_01_file_lst[[id_chr]]
  want_chr <- V_09_02_env$B_03_07_file_fn(id_chr)
  if (is.null(path_chr) || !file.exists(path_chr)) {
    V_12_03_bad_chr <- c(V_12_03_bad_chr, sprintf("%s wrote nothing", id_chr))
    next
  }
  if (!identical(basename(path_chr), want_chr)) {
    V_12_03_bad_chr <- c(V_12_03_bad_chr,
                         sprintf("%s is named %s, not %s", id_chr,
                                 basename(path_chr), want_chr))
  }
  if (file.size(path_chr) <= 0) {
    V_12_03_bad_chr <- c(V_12_03_bad_chr, sprintf("%s is empty", id_chr))
    next
  }
  dim_int <- V_11_01_dim_fn(path_chr)
  if (!identical(as.integer(dim_int), as.integer(V_12_02_want_int))) {
    V_12_03_bad_chr <- c(V_12_03_bad_chr,
                         sprintf("%s is %d x %d", id_chr, dim_int[["w"]],
                                 dim_int[["h"]]))
  }
  if (dim_int[["w"]] / dim_int[["h"]] < 1.5) {
    V_12_03_bad_chr <- c(V_12_03_bad_chr,
                         sprintf("%s is %.2f:1", id_chr,
                                 dim_int[["w"]] / dim_int[["h"]]))
  }
  V_12_04_seen_chr <- c(V_12_04_seen_chr,
                        sprintf("%s %.0fkB", want_chr,
                                file.size(path_chr) / 1024))
}

V_01_03_check_fn(
  "every Save PNG button writes a non-empty PNG at the deck's size",
  length(V_12_03_bad_chr) == 0L,
  if (length(V_12_03_bad_chr) == 0L)
    sprintf("%d files, all %d x %d px", length(V_12_04_seen_chr),
            V_12_02_want_int[1L], V_12_02_want_int[2L]) else
    paste(V_12_03_bad_chr, collapse = "; "))

V_01_03_check_fn(
  "every figure on screen has a button, and every button a figure",
  setequal(names(V_09_02_env$B_03_07_figure_lst),
           grep("^plot_", V_10_11_dollar_fn(V_10_12_src_expr, "output"),
                value = TRUE)),
  paste(names(V_09_02_env$B_03_07_figure_lst), collapse = " "))

V_01_03_check_fn(
  "every file name is lower case, hyphenated and carries no date",
  all(grepl("^[a-z0-9]+(-[a-z0-9]+)*\\.png$",
            vapply(names(V_09_02_env$B_03_07_figure_lst),
                   V_09_02_env$B_03_07_file_fn, character(1)))) &&
    !any(grepl("[0-9]{4}", vapply(names(V_09_02_env$B_03_07_figure_lst),
                                  V_09_02_env$B_03_07_file_fn, character(1)))),
  paste(vapply(names(V_09_02_env$B_03_07_figure_lst),
               V_09_02_env$B_03_07_file_fn, character(1)), collapse = ", "))

# V_12_05: the six file names that predate the paired cards
# Note: The export file names are kept across layout changes. These are the
#   six from before the nominal halves were split off, with the stage part
#   in the current stage codes and the figure part unchanged.
V_12_05_old_vec <- c(
  plot_friction = "smets-wouters-2-5a-friction.png",
  plot_ladder   = "smets-wouters-2-5a-friction-ladder.png",
  plot_irf      = "smets-wouters-2-5c-impulse-response.png",
  plot_demand   = "smets-wouters-2-5c-demand-shocks.png",
  plot_fevd     = "smets-wouters-2-5d-variance-shares.png",
  plot_bridge   = "smets-wouters-2-5a-rbc-bridge.png")
V_01_03_check_fn(
  "the six file names from before the restyle are unchanged",
  all(vapply(names(V_12_05_old_vec), function(id_chr)
    identical(V_09_02_env$B_03_07_file_fn(id_chr), V_12_05_old_vec[[id_chr]]),
    logical(1))),
  "")

#### V_12b: The Page Frame #####################################################
# Note: Checked on the page itself. The whole app.R is evaluated (not cut at
#   the user interface, as V_09 does) so the UI object exists; it is then
#   rendered to HTML, which is what the browser gets, and the checks read
#   that HTML: every plotOutput sits in the toolkit's card, T_07_07f, inside
#   its 2:1 holder with a Save PNG and a Save PDF button of its own; the
#   figures the card register holds are exactly B_03_07's; every figure but
#   the two loners is in a pair; no slider label, card header or tile label
#   spells a Greek letter out.

message("\nV_12b The page frame")

V_12b_01_env <- new.env(parent = globalenv())
V_12b_02_html_chr <- tryCatch(local({
  old_chr <- setwd(file.path(V_01_01_here_dir, ".."))
  on.exit(setwd(old_chr))
  for (file_chr in list.files("R", "\\.R$", full.names = TRUE)) {
    sys.source(file_chr, envir = V_12b_01_env)
  }
  src_expr <- parse("app.R")
  # everything but the final shinyApp() call
  for (e_obj in src_expr[-length(src_expr)]) eval(e_obj, envir = V_12b_01_env)
  htmltools::renderTags(V_12b_01_env$E_02_06_page_ui)$html
}), error = function(e_cnd) paste("ERROR:", conditionMessage(e_cnd)))

V_01_03_check_fn(
  "the UI object builds and renders to HTML",
  !startsWith(V_12b_02_html_chr, "ERROR:") &&
    nchar(V_12b_02_html_chr) > 1000L,
  if (startsWith(V_12b_02_html_chr, "ERROR:")) V_12b_02_html_chr else
    sprintf("%d characters", nchar(V_12b_02_html_chr)))

# V_12b_03: every plotOutput is in a T_07_07f card with its buttons
V_12b_03_plot_chr <- regmatches(
  V_12b_02_html_chr,
  gregexpr(paste0('<div class="shiny-plot-output[^"]*" id="[^"]+"|',
                  '<div id="[^"]+" class="shiny-plot-output'),
           V_12b_02_html_chr))[[1L]]
V_12b_03_plot_chr <- sub('.*id="([^"]+)".*', "\\1", V_12b_03_plot_chr)
V_12b_04_bad_chr <- character(0)
for (id_chr in V_12b_03_plot_chr) {
  at_int <- regexpr(sprintf('id="%s"', id_chr), V_12b_02_html_chr,
                    fixed = TRUE)
  before_chr <- substr(V_12b_02_html_chr, max(1L, at_int - 600L), at_int)
  after_chr  <- substr(V_12b_02_html_chr, at_int, at_int + 1500L)
  if (!grepl("fig-card", before_chr, fixed = TRUE)) {
    V_12b_04_bad_chr <- c(V_12b_04_bad_chr,
                          paste(id_chr, "not in a fig-card"))
  }
  if (!grepl('<div class="fig-r21">\\s*<div[^>]*$', before_chr)) {
    V_12b_04_bad_chr <- c(V_12b_04_bad_chr, paste(id_chr, "not held at 2:1"))
  }
  for (ext_chr in c("png", "pdf")) {
    if (!grepl(sprintf('id="%s__%s"', id_chr, ext_chr), after_chr,
               fixed = TRUE)) {
      V_12b_04_bad_chr <- c(V_12b_04_bad_chr,
                            sprintf("%s has no Save %s under it", id_chr,
                                    toupper(ext_chr)))
    }
  }
}
V_01_03_check_fn(
  "every plotOutput sits in a T_07_07f card with PNG and PDF",
  length(V_12b_03_plot_chr) == length(V_09_02_env$B_03_07_figure_lst) &&
    length(V_12b_04_bad_chr) == 0L,
  if (length(V_12b_04_bad_chr) == 0L)
    sprintf("%d plots: %s", length(V_12b_03_plot_chr),
            paste(V_12b_03_plot_chr, collapse = " ")) else
    paste(V_12b_04_bad_chr, collapse = "; "))

V_01_03_check_fn(
  "the card register holds exactly B_03_07's figures",
  setequal(V_12b_01_env$T_07_07e_ids_fn(),
           names(V_09_02_env$B_03_07_figure_lst)) &&
    setequal(V_12b_03_plot_chr, names(V_09_02_env$B_03_07_figure_lst)),
  paste(V_12b_01_env$T_07_07e_ids_fn(), collapse = " "))

# V_12b_05: the pairs
# Note: T_07_07g_pair_fn is a bslib layout_columns. Each figure's card is
#   found inside one, and the pairs are the ones B_03_07's notes promise.
V_12b_05_src_chr <- paste(readLines(file.path(V_01_01_here_dir, "..",
                                              "app.R"), warn = FALSE),
                          collapse = "\n")
V_12b_06_pair_lst <- list(c("plot_friction", "plot_ladder"),
                          c("plot_irf", "plot_irf_nom"),
                          c("plot_demand", "plot_demand_nom"),
                          c("plot_s2_irf", "plot_s2_fevd"),
                          c("plot_vs_irf", "plot_vs_irf_nom"),
                          c("plot_vs_fevd", "plot_vs_fevd_long"))
V_12b_07_ok_lgl <- vapply(V_12b_06_pair_lst, function(p_chr)
  grepl(sprintf(paste0('T_07_07g_pair_fn\\(E_02_01_card_fn\\("%s"\\),\\s*',
                       'E_02_01_card_fn\\("%s"\\)\\)'), p_chr[1L], p_chr[2L]),
        V_12b_05_src_chr), logical(1))
V_01_03_check_fn(
  "related figures sit two to a row",
  all(V_12b_07_ok_lgl),
  paste(vapply(V_12b_06_pair_lst, paste, character(1), collapse = " + "),
        collapse = "; "))

# V_12b_08: no spelled-out Greek where a student reads a symbol
# Note: Entities (&lambda;) and MathJax (\lambda) are what is wanted, so
#   both are stripped first; what is left must not contain a Greek letter's
#   name as a word, nor a spelling such as "xi_p" or "r_pi".
V_12b_08_greek_chr <- paste0(
  "\\b(alpha|beta|gamma|delta|epsilon|varepsilon|zeta|eta|theta|iota|",
  "kappa|lambda|mu|nu|xi|pi|rho|sigma|tau|phi|varphi|chi|psi|omega)\\b|",
  "[a-z]_(p|w|pi|t)\\b")
V_12b_09_strip_fn <- function(x_chr) {
  x_chr <- gsub("&[A-Za-z]+;|&#[0-9]+;", " ", x_chr)
  x_chr <- gsub("\\\\\\(.*?\\\\\\)", " ", x_chr, perl = TRUE)
  gsub("<[^>]+>", " ", x_chr)
}
V_12b_10_text_chr <- c(
  slider = vapply(V_12b_01_env$B_01_01_slider_lst, `[[`, "", "label"),
  help   = unlist(V_12b_01_env$B_01_01b_help_lst),
  card   = vapply(V_12b_01_env$B_03_07_figure_lst, `[[`, "", "title_chr"),
  header = regmatches(V_12b_02_html_chr, gregexpr(
    '<div class="card-header[^"]*">[^<]*', V_12b_02_html_chr))[[1L]],
  label  = regmatches(V_12b_02_html_chr, gregexpr(
    '<div class="ctl-label">.*?</div>', V_12b_02_html_chr, perl = TRUE))[[1L]],
  test   = vapply(V_12b_01_env$B_05_01_test_lst, function(t_lst)
    paste(t_lst$name, t_lst$data, t_lst$src), ""))
V_12b_11_bad_chr <- V_12b_10_text_chr[grepl(
  V_12b_08_greek_chr, V_12b_09_strip_fn(V_12b_10_text_chr),
  ignore.case = TRUE, perl = TRUE)]
V_01_03_check_fn(
  "no slider label, card header or test row spells out a Greek name",
  length(V_12b_11_bad_chr) == 0L &&
    sum(grepl("^slider", names(V_12b_10_text_chr))) == 7L &&
    sum(grepl("^label", names(V_12b_10_text_chr))) == 7L,
  if (length(V_12b_11_bad_chr) == 0L)
    sprintf("%d strings read", length(V_12b_10_text_chr)) else
    paste(V_12b_11_bad_chr, collapse = " | "))

#### V_13: The Equations, Notation and In-Words Structures #####################
# Note: The is-mp-pc contract, checked rather than trusted. An equation
#   with no note, or a version with no note, renders as a blank cell in the
#   In Words panel. A notation symbol that appears in no equation is a
#   symbol the student is asked to learn and never meets, and a symbol in
#   an equation with no notation row is the reverse. The reverse check is
#   over LaTeX command tokens: every Greek letter and named operator is a
#   command, so it covers exactly the class of symbol a reader cannot guess
#   at, while bare Latin letters are indices as often as variables (u and v
#   in the variance share, j and x for a variable and a shock).

message("\nV_13  Equations, notation and in-words")

V_13_01_eq_lst   <- V_09_02_env$B_04_01_equation_lst
V_13_02_note_lst <- V_09_02_env$B_04_03_notation_lst
V_13_03_stage_chr <- unname(V_09_02_env$B_02_01_stage_vec)

V_13_04_flat_fn <- function(x_chr) gsub("[[:space:]]", "", x_chr)

# V_13_05: every item has a label, a version and a note per version
V_13_06_bad_chr <- character(0)
for (it_lst in V_13_01_eq_lst) {
  tag_chr <- if (is.null(it_lst$label)) "<unlabelled>" else it_lst$label
  if (is.null(it_lst$label) || !nzchar(it_lst$label)) {
    V_13_06_bad_chr <- c(V_13_06_bad_chr, "an item has no label")
  }
  if (length(it_lst$versions) < 1L) {
    V_13_06_bad_chr <- c(V_13_06_bad_chr, paste(tag_chr, "has no version"))
  }
  miss_chr <- setdiff(names(it_lst$versions), names(it_lst$notes))
  if (length(miss_chr) > 0L) {
    V_13_06_bad_chr <- c(V_13_06_bad_chr,
                         sprintf("%s has no note at %s", tag_chr,
                                 paste(miss_chr, collapse = " ")))
  }
  if (!it_lst$group %in% names(V_09_02_env$B_04_02_group_vec)) {
    V_13_06_bad_chr <- c(V_13_06_bad_chr,
                         sprintf("%s is in group %s", tag_chr, it_lst$group))
  }
  off_chr <- setdiff(names(it_lst$versions), V_13_03_stage_chr)
  if (length(off_chr) > 0L) {
    V_13_06_bad_chr <- c(V_13_06_bad_chr,
                         sprintf("%s is keyed to %s, which is not a stage",
                                 tag_chr, paste(off_chr, collapse = " ")))
  }
}
V_01_03_check_fn(
  "every equation has a label, a version and a note for each version",
  length(V_13_06_bad_chr) == 0L,
  if (length(V_13_06_bad_chr) == 0L)
    sprintf("%d equations in %d groups", length(V_13_01_eq_lst),
            length(unique(vapply(V_13_01_eq_lst, function(x) x$group,
                                 character(1))))) else
    paste(unique(V_13_06_bad_chr), collapse = "; "))

V_01_03_check_fn(
  "every group in B_04_02 is used, and no label is written twice",
  setequal(unique(vapply(V_13_01_eq_lst, function(x) x$group, character(1))),
           names(V_09_02_env$B_04_02_group_vec)) &&
    !anyDuplicated(vapply(V_13_01_eq_lst, function(x) x$label, character(1))),
  paste(names(V_09_02_env$B_04_02_group_vec), collapse = ", "))

# V_13_07: notation is well formed
V_13_08_sym_chr <- vapply(V_13_02_note_lst, function(x) x$sym, character(1))
V_01_03_check_fn(
  "every notation entry is well formed and no symbol is repeated",
  all(vapply(V_13_02_note_lst, function(x_lst)
    x_lst$grp %in% c("var", "par", "tgt", "shk") &&
      nzchar(x_lst$sym) && nzchar(x_lst$txt) &&
      !grepl("\\.$", x_lst$txt) &&
      identical(substr(x_lst$txt, 1L, 1L),
                tolower(substr(x_lst$txt, 1L, 1L))) &&
      x_lst$from %in% V_13_03_stage_chr, logical(1))) &&
    !anyDuplicated(V_13_08_sym_chr),
  sprintf("%d symbols: %d var, %d par, %d tgt, %d shk",
          length(V_13_02_note_lst),
          sum(vapply(V_13_02_note_lst, function(x) x$grp == "var", logical(1))),
          sum(vapply(V_13_02_note_lst, function(x) x$grp == "par", logical(1))),
          sum(vapply(V_13_02_note_lst, function(x) x$grp == "tgt", logical(1))),
          sum(vapply(V_13_02_note_lst, function(x) x$grp == "shk",
                     logical(1)))))

# V_13_09: every symbol is used, and every command is declared
V_13_09_raw_chr <- unlist(lapply(V_13_01_eq_lst, function(x) x$versions))
V_13_10_tex_chr <- V_13_04_flat_fn(V_13_09_raw_chr)
# Note: The token scan runs on the raw LaTeX, not on the flattened copy.
#   The flattened copy exists so that "y_t - y^p_t" can be found inside
#   "r_y (y_t - y^p_t)"; run the command regex over it and "\alpha k^s_t"
#   becomes "\alphak", which is not a command.
V_13_11_all_chr <- paste(V_13_09_raw_chr, collapse = " ")

V_13_12_unused_chr <- V_13_08_sym_chr[
  !vapply(V_13_08_sym_chr, function(s_chr)
    any(grepl(V_13_04_flat_fn(s_chr), V_13_10_tex_chr, fixed = TRUE)),
    logical(1))]
V_01_03_check_fn(
  "every symbol in B_04_03 appears in at least one equation",
  length(V_13_12_unused_chr) == 0L,
  if (length(V_13_12_unused_chr) == 0L) "all declared symbols are drawn" else
    paste("never drawn:", paste(V_13_12_unused_chr, collapse = ", ")))

V_13_13_skip_chr <- c(
  "left", "right", "frac", "dfrac", "qquad", "quad", "text", "mathrm",
  "sim", "equiv", "to", "cdot", "sum", "dim", "bar", "hat", "tilde",
  "big", "Big", "bigl", "bigr", "operatorname", "times", "ldots", "in",
  "begin", "end", "aligned", "max", "min", "log", "exp")

V_13_14_token_chr <- unique(unlist(
  regmatches(V_13_11_all_chr,
             gregexpr("\\\\[A-Za-z]+", V_13_11_all_chr))))
V_13_14_token_chr <- setdiff(V_13_14_token_chr,
                             paste0("\\", V_13_13_skip_chr))

V_13_15_undeclared_chr <- V_13_14_token_chr[
  !vapply(V_13_14_token_chr, function(t_chr)
    any(grepl(t_chr, V_13_08_sym_chr, fixed = TRUE)), logical(1))]
V_01_03_check_fn(
  "every LaTeX symbol in an equation has a notation entry",
  length(V_13_15_undeclared_chr) == 0L,
  if (length(V_13_15_undeclared_chr) == 0L)
    sprintf("%d command tokens, all declared", length(V_13_14_token_chr)) else
    paste("undeclared:", paste(V_13_15_undeclared_chr, collapse = " ")))

# V_13_16: the LaTeX is well formed, so MathJax cannot fail silently
# Note: An unbalanced brace or a \left with no \right does not raise anything
#   in R. It reaches the browser, MathJax gives up on that one formula, and
#   the panel shows raw backslashes where an equation should be - which is
#   invisible from here and obvious to a student. Cheap to check, so it is.

V_13_16_all_chr <- c(V_13_09_raw_chr,
                     vapply(V_13_02_note_lst, function(x) x$sym, character(1)))
V_13_17_count_fn <- function(s_chr, pat_chr) {
  hit_lst <- gregexpr(pat_chr, s_chr)
  vapply(hit_lst, function(h_int)
    if (h_int[1L] == -1L) 0L else length(h_int), integer(1))
}
# Note: Every brace is counted, escaped or not. An escaped pair, \{ and \},
#   balances on its own, so there is no need for a lookbehind.
V_13_18_bad_tex_chr <- V_13_16_all_chr[
  V_13_17_count_fn(V_13_16_all_chr, "\\{") !=
    V_13_17_count_fn(V_13_16_all_chr, "\\}") |
    V_13_17_count_fn(V_13_16_all_chr, "\\\\left") !=
      V_13_17_count_fn(V_13_16_all_chr, "\\\\right") |
    grepl("\\$", V_13_16_all_chr)]
V_01_03_check_fn(
  "every equation and symbol is balanced LaTeX",
  length(V_13_18_bad_tex_chr) == 0L,
  if (length(V_13_18_bad_tex_chr) == 0L)
    sprintf("%d strings checked", length(V_13_16_all_chr)) else
    substr(V_13_18_bad_tex_chr[1L], 1L, 60L))

# V_13_19: the panels agree with the sliders' letters
# Note: B_01_01's slider labels carry a letter in brackets as HTML and
#   B_04_03 carries the same letters as LaTeX; if one is edited without the
#   other, a student reading a slider cannot find it in an equation. The
#   bracket is matched against the HTML the notation symbol stands for;
#   &phi; is the curly phi MathJax draws for \varphi. The module's letters:
#   theta_p, theta_w, phi_pi and rho_i, with lambda, varphi and iota_p.
V_13_20_want_vec <- c(chabb = "\\lambda", csadjcost = "\\varphi",
                      cprobp = "\\theta_p", cprobw = "\\theta_w",
                      cindp = "\\iota_p", crpi = "\\phi_\\pi",
                      crr = "\\rho_i")
V_13_20b_html_vec <- c(chabb = "&lambda;", csadjcost = "&phi;",
                       cprobp = "&theta;<sub>p</sub>",
                       cprobw = "&theta;<sub>w</sub>",
                       cindp = "&iota;<sub>p</sub>",
                       crpi = "&phi;<sub>&pi;</sub>",
                       crr = "&rho;<sub>i</sub>")
V_13_21_bad_chr <- character(0)
for (nm_chr in names(V_13_20_want_vec)) {
  lab_chr <- V_09_02_env$B_01_01_slider_lst[[nm_chr]]$label
  bare_chr <- V_13_20b_html_vec[[nm_chr]]
  if (!grepl(paste0("(", bare_chr, ")"), lab_chr, fixed = TRUE)) {
    V_13_21_bad_chr <- c(V_13_21_bad_chr,
                         sprintf("%s says %s", nm_chr, lab_chr))
  }
  if (!V_13_20_want_vec[[nm_chr]] %in% V_13_08_sym_chr) {
    V_13_21_bad_chr <- c(V_13_21_bad_chr,
                         sprintf("%s has no notation row",
                                 V_13_20_want_vec[[nm_chr]]))
  }
}
V_01_03_check_fn(
  "every slider's letter is the letter the notation panel uses",
  length(V_13_21_bad_chr) == 0L,
  if (length(V_13_21_bad_chr) == 0L)
    "lambda, varphi, theta_p, theta_w, iota_p, phi_pi, rho_i" else
    paste(V_13_21_bad_chr, collapse = "; "))

#### V_14: What the App Says Agrees With the Model and With Whelan #############
# Note: Each claim in the app's prose, held to the model's own output.

message("\nV_14  What the app says")

# V_14_01: the in-app scorecard passes at the posterior mode
V_14_01_fevd_lst <- C_04_02_fevd_fn(V_02_02_sol_lst)
V_14_01_irf_lst  <- stats::setNames(
  lapply(C_01_02_shock_vec, function(sh_chr)
    C_04_01_irf_fn(V_02_02_sol_lst, sh_chr, 25L)), C_01_02_shock_vec)
V_14_02_fail_chr <- vapply(V_09_02_env$B_05_01_test_lst, function(t_lst)
  if (isTRUE(t_lst$pass(V_02_02_sol_lst, V_14_01_fevd_lst, V_14_01_irf_lst)))
    "" else t_lst$name, character(1))
V_01_03_check_fn("every row of the Tests tab matches at the mode, theirs",
                 !any(nzchar(V_14_02_fail_chr)),
                 if (!any(nzchar(V_14_02_fail_chr)))
                   sprintf("%d rows", length(V_14_02_fail_chr)) else
                   paste(V_14_02_fail_chr[nzchar(V_14_02_fail_chr)],
                         collapse = "; "))
# V_14_02b: the table scores the estimated model and the sliders apart
# Note: After "Let Investment Move" or "Make Prices Flexible" some rows
#   read Does not match at the sliders, which is not the model failing
#   Whelan: the estimated model matches every row. So the table carries two
#   columns. Checked here: the app's own scorer passes every row at the
#   posterior mode; the table renders both columns under headers that say
#   which is which; and, driven through testServer with prices made
#   flexible, the Estimated column still matches every row while At Your
#   Settings does not, and nothing is set in capitals. Our model is
#   reported, not required to pass: the check is that the app scores our
#   model, every row, and prints the count.
V_14_02c_post_lgl <- V_09_02_env$B_05_02_score_fn(
  C_03_02_solve_fn(V_09_02_env$B_01_02_default_lst))
V_01_03_check_fn("the app scores OUR model at the mode, every row",
                 identical(V_09_02_env$B_01_02_default_lst$model_chr,
                           "ours") &&
                   length(V_14_02c_post_lgl) ==
                   length(V_09_02_env$B_05_01_test_lst) &&
                   !anyNA(V_14_02c_post_lgl),
                 sprintf("ours matches %d of %d",
                         sum(V_14_02c_post_lgl %in% TRUE),
                         length(V_14_02c_post_lgl)))

V_14_02d_seen_lst <- local({
  out_lst <- list()
  shiny::testServer(file.path(V_01_01_here_dir, ".."), {
    session$setInputs(stage = "2.5c")
    out_lst$post_chr <<- as.character(output$test_ui$html)
    session$setInputs(preset_flexprice = 1)
    # testServer does not echo an update* back into input, so the slider
    # and its box are set by hand to what the preset sets
    session$setInputs(cprobp = 0.05, cprobp_box = 0.05)
    out_lst$flex_chr <<- as.character(output$test_ui$html)
  })
  out_lst
})

V_14_02e_cells_fn <- function(html_chr) {
  row_chr <- regmatches(html_chr, gregexpr("(?s)<tr>.*?</tr>", html_chr,
                                           perl = TRUE))[[1L]]
  row_chr <- row_chr[!grepl("<th>", row_chr, fixed = TRUE)]
  lapply(row_chr, function(r_chr) {
    cell_chr <- regmatches(r_chr, gregexpr("(?s)<td[^>]*>.*?</td>", r_chr,
                                           perl = TRUE))[[1L]]
    trimws(gsub("<[^>]+>", "", cell_chr))
  })
}
V_14_02f_head_ok_fn <- function(html_chr) {
  all(vapply(c("<th>Estimated</th>", "<th>At Your Settings</th>",
               "<strong>Estimated</strong>",
               "<strong>At Your Settings</strong>"),
             grepl, logical(1), x = html_chr, fixed = TRUE))
}
V_14_02g_flex_lst <- V_14_02e_cells_fn(V_14_02d_seen_lst$flex_chr)
V_14_02h_post_col <- vapply(V_14_02g_flex_lst, `[`, "", 3L)
V_14_02i_now_col  <- vapply(V_14_02g_flex_lst, `[`, "", 4L)

V_01_03_check_fn(
  "the Tests table renders both columns, headed and explained",
  V_14_02f_head_ok_fn(V_14_02d_seen_lst$post_chr) &&
    V_14_02f_head_ok_fn(V_14_02d_seen_lst$flex_chr) &&
    length(V_14_02g_flex_lst) == length(V_09_02_env$B_05_01_test_lst) &&
    all(lengths(V_14_02g_flex_lst) == 5L),
  sprintf("%d rows x %d cells", length(V_14_02g_flex_lst),
          if (length(V_14_02g_flex_lst)) length(V_14_02g_flex_lst[[1L]])
          else 0L))
V_01_03_check_fn(
  "flexible prices: Estimated is unmoved, At Your Settings moves",
  identical(V_14_02h_post_col == "Matches",
            as.vector(V_14_02c_post_lgl %in% TRUE)) &&
    any(V_14_02i_now_col == "Does not match"),
  sprintf("estimated %d of %d match; at your settings %d of %d",
          sum(V_14_02h_post_col == "Matches"), length(V_14_02h_post_col),
          sum(V_14_02i_now_col == "Matches"), length(V_14_02i_now_col)))
V_01_03_check_fn(
  "the verdicts are in sentence case, not capitals",
  length(V_14_02h_post_col) > 0L &&
    !grepl("DOES NOT MATCH|MATCHES", V_14_02d_seen_lst$flex_chr) &&
    all(c(V_14_02h_post_col, V_14_02i_now_col) %in%
          c("Matches", "Does not match", "No solution")),
  paste(unique(c(V_14_02h_post_col, V_14_02i_now_col)), collapse = ", "))

V_01_03_check_fn("the Q1 output row names spending, not the risk premium",
                 any(vapply(V_09_02_env$B_05_01_test_lst, function(t_lst)
                   grepl("^exogenous spending", t_lst$data), logical(1))),
                 "")

# V_14_03: a positive risk premium shock is a rise in the premium
# Note: The solver writes the paper's eqs (2) and (4), so a positive
#   innovation lowers consumption (Whelan's [W11 5] form), and the
#   per-shock figure's title and caption say "Rise". [W11 18] (the paper's
#   Figure 2) plots the replication's sign, output up, i.e. a fall in the
#   premium, so D_03_02, the figure that reproduces it, draws a one-s.d.
#   fall and names the line "Risk Premium Falls". Read off Whelan's chart:
#   output +0.43 on impact.
V_14_03_eb_df <- C_04_01_irf_fn(V_02_02_sol_lst, "eb", 2L)
V_14_04_title_chr <-
  V_09_02_env$D_03_01_irf_fn(V_02_02_sol_lst, "eb")$labels$title
V_01_03_check_fn("a positive risk premium shock lowers c, q and output",
                 V_14_03_eb_df$y[1L] < 0 && V_14_03_eb_df$c[1L] < 0 &&
                   V_14_03_eb_df$pk[1L] < 0 &&
                   grepl("Rise in the Risk Premium", V_14_04_title_chr),
                 sprintf("output %+.3f on impact; title '%s'",
                         V_14_03_eb_df$y[1L], V_14_04_title_chr))
V_14_03b_cap_chr <-
  V_09_02_env$D_03_01_irf_fn(V_02_02_sol_lst, "eb")$labels$caption
V_01_03_check_fn("the per-shock figure says it draws a RISE in the premium",
                 grepl("RISE in the risk premium", V_14_03b_cap_chr,
                       fixed = TRUE),
                 "")
V_14_03c_dem_obj <- V_09_02_env$D_03_02_demand_fn(V_02_02_sol_lst)
V_14_03d_dem_df  <- V_14_03c_dem_obj$data
V_14_03e_rp_num  <- V_14_03d_dem_df$value_num[
  V_14_03d_dem_df$shock_cat == C_01_02_shock_lab_vec[["eb"]] &
    V_14_03d_dem_df$panel_cat == "Output" & V_14_03d_dem_df$quarter_int == 0L]
V_01_03_check_fn("[W11 18] the demand figure draws a FALL in the premium",
                 length(V_14_03e_rp_num) == 1L &&
                   abs(V_14_03e_rp_num - 0.43) < 0.05 &&
                   abs(V_14_03e_rp_num + V_14_03_eb_df$y[1L]) < 1e-12 &&
                   grepl("FALLS", V_14_03c_dem_obj$labels$caption,
                         fixed = TRUE),
                 sprintf("output %+.3f on impact, against +0.43 on his chart",
                         V_14_03e_rp_num))

# V_14_06: the displayed variance share carries the sigmas
V_14_06_tex_chr <- Filter(function(x) identical(x$label,
  "Forecast Error Variance Share"),
  V_09_02_env$B_04_01_equation_lst)[[1L]]$versions[["2.5d"]]
V_01_03_check_fn("the variance-share formula on screen weights by sigma^2",
                 grepl("\\sigma_x^2", V_14_06_tex_chr, fixed = TRUE) &&
                   grepl("\\sigma_v^2", V_14_06_tex_chr, fixed = TRUE), "")

# V_14_07: no story repeats a sigma-less number
V_14_07_story_chr <- paste(vapply(V_09_02_env$B_02_02_example_lst,
                                  function(e) e$story, character(1)),
                           collapse = " ")
V_01_03_check_fn("no worked example repeats the sigma-less shares",
                 !grepl("about half of it|about a tenth|close to neutral",
                        V_14_07_story_chr), "")

#### V_15: The Figures Follow the Toolkit and the Decks ########################
# Note: The theme is the toolkit's, the series colours are the toolkit's in
#   deck order with no two alike, no line figure carries a legend, and the
#   Save PDF buttons write PDFs.

message("\nV_15  Figures on the toolkit")

V_01_03_check_fn("app.R defines no local ggplot theme",
                 !any(grepl("^[^#]*(theme_bw|theme_minimal|theme_classic)\\(",
                            V_11_03_src_chr)), "")

V_15_01_plots_lst <- list(
  friction = V_09_02_env$D_02_01_friction_fn(V_02_02_sol_lst),
  ladder   = V_09_02_env$D_02_02_ladder_fn(V_02_01_par_lst),
  irf      = V_09_02_env$D_03_01_irf_fn(V_02_02_sol_lst, "em"),
  demand   = V_09_02_env$D_03_02_demand_fn(V_02_02_sol_lst),
  fevd     = V_09_04_p_obj,
  bridge   = V_09_02_env$D_05_01_bridge_fn(V_02_02_sol_lst, V_02_01_par_lst))
V_15_02_wash_chr <- V_09_02_env$T_01_01_palette_vec[["wash"]]
V_15_03_bad_chr <- names(Filter(function(p_obj) !identical(
  p_obj$theme$panel.background$fill, V_15_02_wash_chr), V_15_01_plots_lst))
V_01_03_check_fn("every builder draws on the toolkit theme",
                 length(V_15_03_bad_chr) == 0L,
                 if (length(V_15_03_bad_chr) == 0L) "six builders" else
                   paste(V_15_03_bad_chr, collapse = ", "))

V_15_04_line_col_fn <- function(p_obj) {
  b_lst <- ggplot2::ggplot_build(p_obj)
  hit_int <- which(vapply(p_obj$layers, function(l)
    inherits(l$geom, "GeomLine"), logical(1)))
  unique(toupper(b_lst$data[[hit_int[1L]]]$colour))
}
V_15_05_series_vec <- toupper(unname(
  V_09_02_env$T_01_02_series_vec[c("main", "second", "third")]))
V_15_06_bridge_chr <- V_15_04_line_col_fn(V_15_01_plots_lst$bridge)
V_15_07_demand_chr <- V_15_04_line_col_fn(V_15_01_plots_lst$demand)
V_01_03_check_fn("the three-line figures use blue, green, navy, all distinct",
                 setequal(V_15_06_bridge_chr, V_15_05_series_vec) &&
                   setequal(V_15_07_demand_chr, V_15_05_series_vec),
                 sprintf("bridge %s; demand %s",
                         paste(V_15_06_bridge_chr, collapse = " "),
                         paste(V_15_07_demand_chr, collapse = " ")))
V_01_03_check_fn("the one-series figures draw in the main series colour",
                 all(vapply(V_15_01_plots_lst[c("friction", "irf")],
                            function(p_obj)
                              identical(V_15_04_line_col_fn(p_obj),
                                        V_15_05_series_vec[1L]),
                            logical(1))), V_15_05_series_vec[1L])
V_01_03_check_fn("the seven shock colours are all different",
                 !anyDuplicated(toupper(V_09_02_env$B_03_02_shock_col_vec)),
                 paste(V_09_02_env$B_03_02_shock_col_vec, collapse = " "))

# V_15_07b: stage 5 draws ours in the main colour, theirs light blue
V_15_07c_vs_obj <- V_09_02_env$D_06_01_vs_irf_fn(V_02_02b_ours_sol_lst,
                                                 V_02_02_sol_lst, "em")
V_15_07d_vs_df  <- ggplot2::ggplot_build(V_15_07c_vs_obj)$data
V_15_07e_col_chr <- unique(toupper(unlist(lapply(V_15_07d_vs_df, function(d)
  if ("colour" %in% names(d) && nrow(d) > 20L) d$colour))))
V_01_03_check_fn(
  "stage 5: ours in the main colour, theirs in the comparator",
  setequal(V_15_07e_col_chr,
           toupper(unname(V_09_02_env$T_01_02_series_vec[c("main",
                                                           "compare")]))) &&
    identical(unname(V_09_02_env$D_06_00_col_vec),
              unname(V_09_02_env$T_01_02_series_vec[c("main", "compare")])),
  paste(V_15_07e_col_chr, collapse = " "))
V_15_01_plots_lst$vs_irf  <- V_15_07c_vs_obj
V_15_01_plots_lst$vs_fevd <- V_09_02_env$D_06_02_vs_fevd_fn(
  V_02_02b_ours_sol_lst, V_02_02_sol_lst, "y")

# V_15_08: the legend measurement can see a legend
# Note: V_09_13 passing proves nothing unless it could fail. The same
#   reading is taken on a plot that does carry a legend and must come back
#   non-zero, on ggplot2 3.4.x (one "guide-box" cell) and 3.5+ (four
#   "guide-box-<side>" cells) alike.
V_15_08_guide_fn <- V_09_13b_guide_fn
V_15_09_with_num <- V_15_08_guide_fn(
  ggplot2::ggplot(data.frame(x = 1:2, y = 1:2, g = c("a", "b")),
                  ggplot2::aes(x, y, colour = g)) + ggplot2::geom_point())
V_15_10_none_chr <- names(Filter(
  function(p_obj) V_15_08_guide_fn(p_obj) > 1e-06, V_15_01_plots_lst))
V_01_03_check_fn("the legend check sees a real legend, and none here",
                 V_15_09_with_num > 1 && length(V_15_10_none_chr) == 0L,
                 sprintf("ggplot2 %s: test legend %.0f mm^2; %s",
                         as.character(utils::packageVersion("ggplot2")),
                         V_15_09_with_num,
                         if (length(V_15_10_none_chr) == 0L)
                           "no builder has one" else
                           paste("legend on", paste(V_15_10_none_chr,
                                                    collapse = ", "))))

# V_15_11: the Save PDF buttons
V_15_11_pdf_lst <- V_12_01b_drive_fn("pdf")
V_15_12_bad_chr <- names(V_09_02_env$B_03_07_figure_lst)[!vapply(
  names(V_09_02_env$B_03_07_figure_lst), function(id_chr) {
    f_chr <- V_15_11_pdf_lst[[id_chr]]
    !is.null(f_chr) && file.exists(f_chr) &&
      identical(readChar(f_chr, 4L, useBytes = TRUE), "%PDF") &&
      identical(basename(f_chr),
                sub("\\.png$", ".pdf", V_09_02_env$B_03_07_file_fn(id_chr)))
  }, logical(1))]
V_01_03_check_fn("every Save PDF button writes a PDF under the PNG's name",
                 length(V_15_12_bad_chr) == 0L,
                 if (length(V_15_12_bad_chr) == 0L)
                   sprintf("%d files", length(V_15_11_pdf_lst)) else
                   paste(V_15_12_bad_chr, collapse = ", "))

#### V_16: Every Displayed Equation Holds on the Solver's Own Path #############
# Note: The Equations panel is what R/model.R solves, for each side. Our
#   model's display is the "2.5a" version of every item; Smets and Wouters'
#   is "2.5e" where it has one and "2.5a" where the two share an equation.
#   tests/two_models.R writes each side's equations out in R, term for term
#   as the panel prints them, with every composite computed from that
#   side's own descriptor formulas (W_01), not read out of C_01_03; W_02_02
#   pins the LaTeX terms that tell the sides apart, so the R and the panel
#   cannot drift. Each side is simulated from its own solution and every
#   residual must be machine zero.

message("\nV_16  Every displayed equation, evaluated on each side's solution")

sys.source(file.path(V_01_01_here_dir, "two_models.R"), envir = V_09_02_env)

V_16_01_p <- V_02_01_par_lst
for (V_16_00_side_chr in c("ours", "sw")) {
  V_16_02_sol <- if (V_16_00_side_chr == "ours") V_02_02b_ours_sol_lst else
    V_02_02_sol_lst
  V_16_03_k <- V_09_02_env$W_01_01_coef_fn(V_16_02_sol$par_lst,
                                           V_16_00_side_chr, V_16_02_sol$d_lst)
  V_16_03b_gap_vec <- V_09_02_env$W_01_02_gap_fn(V_16_03_k, V_16_02_sol$d_lst)
  V_01_03_check_fn(
    sprintf("%s: every displayed composite formula is the solver's",
            V_16_00_side_chr),
    max(V_16_03b_gap_vec) < 1e-12,
    sprintf("%d composites and identities, worst %.1e",
            length(V_16_03b_gap_vec), max(V_16_03b_gap_vec)))
  V_16_12_res_lst <- V_09_02_env$W_02_01_resid_fn(V_16_02_sol, V_16_03_k)
  V_16_13_worst_vec <- vapply(V_16_12_res_lst, function(r) max(abs(r)),
                              numeric(1))
  V_16_14_bad_chr <- names(V_16_13_worst_vec)[V_16_13_worst_vec > 1e-9]
  V_01_03_check_fn(
    sprintf("%s: every displayed equation holds on its solver's path",
            V_16_00_side_chr),
    length(V_16_14_bad_chr) == 0L,
    if (length(V_16_14_bad_chr) == 0L)
      sprintf("%d equations, worst residual %.1e", length(V_16_13_worst_vec),
              max(V_16_13_worst_vec)) else
      paste(sprintf("%s %.2e", V_16_14_bad_chr,
                    V_16_13_worst_vec[V_16_14_bad_chr]), collapse = "; "))
}
V_16_04_tex_bad_chr <- V_09_02_env$W_02_04_tex_bad_fn(
  V_09_02_env$B_04_01_equation_lst)
V_01_03_check_fn(
  "the panel's LaTeX carries each side's distinguishing terms",
  length(V_16_04_tex_bad_chr) == 0L,
  if (length(V_16_04_tex_bad_chr) == 0L)
    sprintf("%d rules", length(V_09_02_env$W_02_02_tex_rules_lst)) else
    paste(V_16_04_tex_bad_chr, collapse = "; "))

# V_16_15: every model-group row in the panel is covered above
# Note: So a new equation added to the panel without a line in W_02 fails.
V_16_15_model_chr <- vapply(Filter(function(x) x$group %in%
  c("model", "assumption") && !identical(x$label, "The Innovations") &&
    !identical(x$label, "Potential Output") &&
    !identical(x$label, "The Wage Mark-Up") &&
    !identical(x$label, "Trend Growth"),
  V_09_02_env$B_04_01_equation_lst), function(x) x$label, character(1))
V_16_16_miss_chr <- setdiff(V_16_15_model_chr, names(V_16_12_res_lst))
V_01_03_check_fn("every model equation on the panel is evaluated here",
                 length(V_16_16_miss_chr) == 0L,
                 if (length(V_16_16_miss_chr) == 0L)
                   sprintf("%d panel rows", length(V_16_15_model_chr)) else
                   paste("not evaluated:", paste(V_16_16_miss_chr,
                                                 collapse = ", ")))

# V_16_17: the panel's frictions-off values are the solver's
V_16_17_want_lst <- list(habit = list(chabb = 0),
                         adjcost = list(csadjcost = 0.1),
                         utilisation = list(czcap = 0.99),
                         pindex = list(cindp = 0), windex = list(cindw = 0),
                         stickyp = list(cprobp = 0.05))
V_16_18_ok_lgl <- all(vapply(names(V_16_17_want_lst), function(nm)
  isTRUE(all.equal(unlist(C_06_01_friction_lst[[nm]]$par_lst),
                   unlist(V_16_17_want_lst[[nm]]), tolerance = 1e-5)),
  logical(1))) &&
  isTRUE(all.equal(C_06_04_rbc_limit_fn(V_16_01_p)$par_lst[
    c("chabb", "csadjcost", "czcap", "cindp", "cindw", "cprobp", "cprobw")],
    list(chabb = 0, csadjcost = 0.1, czcap = 0.99, cindp = 0, cindw = 0,
         cprobp = 0.05, cprobw = 0.05), tolerance = 1e-5))
V_01_03_check_fn("the Frictions, Switched Off row uses the solver's values",
                 V_16_18_ok_lgl,
                 "lambda = iota = 0, varphi 0.1, psi 0.99, xi 0.05")

#### V_17: The Risk Premium in the Paper's Units, and the Caveat Tab ###########
# Note: The displayed and solved risk premium equations are the paper's
#   eqs (2) and (4) exactly; the innovation s.d. is Table 1B's sigma_b /
#   c_3, because the authors' code estimates c_3 varepsilon^b; and the
#   Equations card carries the caveat tab. The shares are unchanged by
#   construction and the responses are Dynare's with the sign flipped
#   (V_04b); these check the rest.

message("\nV_17  The risk premium in the paper's units, and the caveat tab")

V_17_01_eq_lst <- V_09_02_env$B_04_01_equation_lst
V_17_02_tex_fn <- function(lab_chr) gsub("\\s+", "", Filter(function(x)
  identical(x$label, lab_chr), V_17_01_eq_lst)[[1L]]$versions[["2.5a"]])
V_01_03_check_fn(
  "consumption shows eq (2): varepsilon^b inside c_3's bracket",
  grepl("-c_3\\left(r_t-E_t\\pi_{t+1}+\\varepsilon^b_t\\right)",
        V_17_02_tex_fn("Consumption Euler Equation"), fixed = TRUE), "")
V_01_03_check_fn(
  "Tobin's q shows eq (4): varepsilon^b inside the bracket, no 1/c_3",
  grepl("-\\left(r_t-E_t\\pi_{t+1}+\\varepsilon^b_t\\right)",
        V_17_02_tex_fn("Tobin's q"), fixed = TRUE) &&
    !grepl("frac{1}{c_3}", V_17_02_tex_fn("Tobin's q"), fixed = TRUE), "")

V_17_03_d   <- V_02_02_sol_lst$d_lst
V_17_04_row <- function(lab_chr) {
  # the G0 row whose own variable is lab_chr, with the shock b in it
  m <- V_02_02_sol_lst$mat_lst$g_zero
  r <- which(m[, match(lab_chr, C_01_02_names_vec)] == 1 &
               m[, match("b", C_01_02_names_vec)] != 0)
  m[r[1L], match(c("r", "Epinf", "b"), C_01_02_names_vec)]
}
V_17_05_c_vec <- V_17_04_row("c")
V_17_06_q_vec <- V_17_04_row("pk")
V_01_03_check_fn(
  "the solver's c and q rows load b exactly as the real rate",
  isTRUE(all.equal(unname(V_17_05_c_vec), c(1, -1, 1) * V_17_03_d$c3)) &&
    isTRUE(all.equal(unname(V_17_06_q_vec), c(1, -1, 1))),
  sprintf("c row (r, E pi, b) = %s; q row = %s",
          paste(sprintf("%.4f", V_17_05_c_vec), collapse = ", "),
          paste(sprintf("%.4f", V_17_06_q_vec), collapse = ", ")))

V_17_07_sig_vec <- C_04_00_sigma_fn(V_02_02_sol_lst)
V_01_03_check_fn(
  "sigma(varepsilon^b) = sigma_b / c_3, the rest as Table 1B",
  abs(V_17_07_sig_vec[["eb"]] - V_02_01_par_lst$sdb / V_17_03_d$c3) <
    1e-12 &&
    abs(V_17_07_sig_vec[["em"]] - V_02_01_par_lst$sdms) < 1e-12,
  sprintf("c_3 = %.4f, sigma_b = %.4f, sigma(eps^b) = %.4f",
          V_17_03_d$c3, V_02_01_par_lst$sdb, V_17_07_sig_vec[["eb"]]))

# Note: The note and the notation quote our model's c_3, since stages 2.5a
#   to 2.5d run it; the conversion is the same.
V_17_07b_ours_d <- V_02_02b_ours_sol_lst$d_lst
V_17_07c_sig_vec <- C_04_00_sigma_fn(V_02_02b_ours_sol_lst)
V_01_03_check_fn(
  "ours: sigma(varepsilon^b) = sigma_b / c_3, with our c_3",
  abs(V_17_07c_sig_vec[["eb"]] - V_02_01b_ours_lst$sdb / V_17_07b_ours_d$c3) <
    1e-12,
  sprintf("c_3 = %.4f, sigma(eps^b) = %.4f", V_17_07b_ours_d$c3,
          V_17_07c_sig_vec[["eb"]]))
V_17_08_note_chr <- V_09_02_env$B_04_04_equation_note_chr
V_17_09_nota_chr <- Filter(function(x) identical(x$sym, "\\varepsilon^b_t"),
                           V_09_02_env$B_04_03_notation_lst)[[1L]]$txt
V_17_10_want_chr <- sprintf("= %.2f",
                            V_02_01b_ours_lst$sdb / V_17_07b_ours_d$c3)
V_01_03_check_fn(
  paste("the note and the notation say sigma_b is the s.d. of c_3 eps^b,",
        "and give 0.24/c_3"),
  all(vapply(list(V_17_08_note_chr, V_17_09_nota_chr), function(t_chr)
    grepl("0.24/c<sub>3</sub>", t_chr, fixed = TRUE) &&
      grepl("c<sub>3</sub>&epsilon;<sup>b</sup>", t_chr, fixed = TRUE) &&
      grepl(V_17_10_want_chr, t_chr, fixed = TRUE), logical(1))),
  V_17_10_want_chr)

# V_17_11: the caveat tab
# Note: A short caveat, "What the Lecture Simplifies", sits on the
#   Equations card with its four items and Whelan's slips in one line; no
#   "Differences From the Paper" tab.
V_17_13_src_chr <- readLines(file.path(V_01_01_here_dir, "..", "app.R"),
                             warn = FALSE)
V_17_14_cav_lst <- V_09_02_env$B_04_07_caveat_lst
V_01_03_check_fn(
  "the caveat tab replaces the differences tab on the Equations card",
  any(grepl('nav_panel("What the Lecture Simplifies"', V_17_13_src_chr,
            fixed = TRUE)) &&
    !any(grepl('nav_panel("Differences From the Paper"', V_17_13_src_chr,
               fixed = TRUE)) &&
    identical(V_09_02_env$B_04_07_caveat_title_chr,
              "What the Lecture Simplifies") &&
    grepl("[W11 3, 7]", V_09_02_env$B_04_07_caveat_whelan_chr,
          fixed = TRUE) &&
    grepl("[W11 6]", V_09_02_env$B_04_07_caveat_whelan_chr, fixed = TRUE) &&
    grepl("[W11 8]", V_09_02_env$B_04_07_caveat_whelan_chr, fixed = TRUE),
  paste(vapply(V_17_14_cav_lst, `[[`, "", "title"), collapse = "; "))

#### V_18: Our Model Is Smets and Wouters at gamma = 1 and c_2 = 0 #############
# Note: Our model's composites, read off the solver, against the closed
#   forms with gamma = 1 and utility separable in hours; and the
#   Smets-Wouters side unchanged. Then every number the app's prose quotes
#   is re-derived from the solver.

message("\nV_18  Our model: gamma = 1, c_2 = 0, and the numbers the app quotes")

V_18_01_d <- V_02_02b_ours_sol_lst$d_lst
V_18_02_p <- V_02_01b_ours_lst
V_18_03_beta <- 1 / (1 + V_18_02_p$constebeta / 100)
V_18_04_want_vec <- with(V_18_02_p, c(
  cgamma = 1, cbetabar = V_18_03_beta, c2 = 0,
  c1 = chabb / (1 + chabb),
  c3 = (1 - chabb) / (csigma * (1 + chabb)),
  i1 = 1 / (1 + V_18_03_beta),
  i2 = 1 / ((1 + V_18_03_beta) * csadjcost),
  k1 = 1 - delta,
  k2 = delta * (1 + V_18_03_beta) * csadjcost,
  crk = 1 / V_18_03_beta - (1 - delta),
  q1 = V_18_03_beta * (1 - delta),
  p1 = cindp / (1 + V_18_03_beta * cindp),
  w1 = 1 / (1 + V_18_03_beta),
  hg = chabb, wc = csigma / (1 - chabb)))
V_18_05_gap_vec <- vapply(names(V_18_04_want_vec), function(nm)
  abs(V_18_01_d[[nm]] - V_18_04_want_vec[[nm]]), numeric(1))
V_01_03_check_fn(
  "ours: every composite is the gamma = 1, c_2 = 0 formula",
  max(V_18_05_gap_vec) < 1e-14 &&
    abs(V_18_01_d$ciy - V_18_02_p$delta * V_18_01_d$cky) < 1e-14,
  sprintf(paste("%d composites, worst %.1e;",
                "c_3 = sigma(1-lambda)/(1+lambda) = %.4f"),
          length(V_18_05_gap_vec), max(V_18_05_gap_vec), V_18_01_d$c3))
V_01_03_check_fn(
  "ours: sigma = 1/sigma_c is the same number, and it multiplies",
  abs((1 / V_18_02_p$csigma) * (1 - V_18_02_p$chabb) /
        (1 + V_18_02_p$chabb) - V_18_01_d$c3) < 1e-14,
  sprintf("sigma = %.4f", 1 / V_18_02_p$csigma))
V_01_03_check_fn(
  "ours solves: a unique stable solution, residuals at machine zero",
  isTRUE(V_02_02b_ours_sol_lst$ok_lgl) &&
    C_05_01_residual_fn(V_02_02b_ours_sol_lst)$system_num < 1e-9 &&
    C_05_02_shares_fn(V_18_01_d)$ok_lgl,
  sprintf("%d unstable roots, %d forecast errors",
          V_02_02b_ours_sol_lst$nu_int, V_02_02b_ours_sol_lst$neta_int))
V_01_03_check_fn(
  "theirs is untouched: the default list is model 'sw'",
  identical(C_01_01_default_lst$model_chr, "sw") &&
    identical(V_02_02_sol_lst$d_lst$cgamma,
              1 + C_01_01_default_lst$ctrend / 100),
  "")

# V_18_06: every number in the app's prose comes from the solver
# Note: B_02_01b computes them; this re-derives each with the model's own
#   functions and checks the prose carries the formatted value.
V_18_06_n <- V_09_02_env$B_02_01b_num_lst
V_18_07_em <- C_04_01_irf_fn(V_02_02b_ours_sol_lst, "em", 25L)
V_18_08_f  <- C_04_02_fevd_fn(V_02_02b_ours_sol_lst)
V_18_09_ex <- V_09_02_env$B_02_02_example_lst
V_18_10_pc <- function(x) sprintf("%.0f per cent", 100 * x)
V_18_11_ok_lgl <- c(
  trough = V_18_06_n$trough_num == min(V_18_07_em$y) &&
    grepl(paste("hardest",
                V_09_02_env$B_02_01b_word_fn(which.min(V_18_07_em$y) - 1L),
                "quarters later"), V_18_09_ex$posterior$story, fixed = TRUE),
  adjcost = grepl(sprintf("%.1f times as deep", min(C_04_01_irf_fn(
    C_06_02_off_fn(V_02_01b_ours_lst, "adjcost"), "em", 25L)$y) /
      min(V_18_07_em$y)), V_18_09_ex$fastinv$story, fixed = TRUE),
  q1 = grepl(V_18_10_pc(V_18_08_f[["1"]]["y", "eg"] +
                          V_18_08_f[["1"]]["y", "eb"]),
             V_18_09_ex$explain$story, fixed = TRUE),
  q100 = grepl(V_18_10_pc(V_18_08_f[["100"]]["y", "ea"]),
               V_18_09_ex$explain$story, fixed = TRUE),
  hours = grepl(sprintf("%.2f", V_18_06_n$real_off_num),
                V_18_09_ex$rbc$story, fixed = TRUE) &&
    V_18_06_n$real_off_num < 0 && V_18_06_n$nom_off_num < 0 &&
    V_18_06_n$both_off_num > 0,
  theirs = grepl(sprintf("%.2f in theirs", min(C_04_01_irf_fn(
    V_02_02_sol_lst, "em", 25L)$y)), V_18_09_ex$theirs$story, fixed = TRUE))
V_01_03_check_fn(
  "every number the stories quote is the solver's",
  all(V_18_11_ok_lgl),
  paste(names(V_18_11_ok_lgl)[!V_18_11_ok_lgl], collapse = ", "))

#### V_19: Stage 2.5e's Pills Are the Equations That Differ ####################
# Note: tests/two_models.R W_03 and W_04. The pills are read off the app
#   the way F_06_02 reads them; the differences are read off the two
#   solvers, row by row, independently of the panel's text.

message(paste("\nV_19  Stage 5: the pills against the two solvers,",
              "the caveat against the pills"))

V_19_01_pill_chr <- V_09_02_env$W_03_01_pills_fn(
  V_09_02_env$B_04_01_equation_lst, V_09_02_env$B_02_03_rank_fn)
V_19_02_diff_lst <- V_09_02_env$W_03_04_diff_fn(C_01_01_default_lst)
V_19_03_let_chr <- unique(unlist(lapply(Filter(function(x)
  x$id %in% c("letters", "potential"), V_17_14_cav_lst), `[[`, "pills")))
V_19_04_bad_chr <- V_09_02_env$W_03_05_check_fn(
  V_19_01_pill_chr, V_19_02_diff_lst, V_09_02_env$B_04_01_equation_lst,
  C_01_01_default_lst, V_19_03_let_chr)
V_01_03_check_fn(
  "CHANGED and NEW are exactly the equations that differ",
  length(V_19_04_bad_chr) == 0L,
  if (length(V_19_04_bad_chr) == 0L)
    sprintf("%d changed, %d new; %d of %d solver blocks differ",
            sum(V_19_01_pill_chr == "changed"), sum(V_19_01_pill_chr == "new"),
            sum(vapply(V_19_02_diff_lst, `[[`, TRUE, "differs")),
            length(V_19_02_diff_lst)) else
    paste(V_19_04_bad_chr, collapse = "; "))
V_19_05_bad_chr <- V_09_02_env$W_04_01_check_fn(V_17_14_cav_lst,
                                                V_19_01_pill_chr)
V_01_03_check_fn(
  "every caveat item has its pills, and every pill its caveat item",
  length(V_19_05_bad_chr) == 0L,
  if (length(V_19_05_bad_chr) == 0L)
    paste(vapply(V_17_14_cav_lst, function(x)
      sprintf("%s %d", x$id, length(x$pills)), ""), collapse = ", ") else
    paste(V_19_05_bad_chr, collapse = "; "))
# V_19_05b: stage 2 renders our full system, every item NEW
V_19_05c_s5 <- local({
  out <- list()
  shiny::testServer(file.path(V_01_01_here_dir, ".."), {
    session$setInputs(stage = "2.5b")
    out$html <<- as.character(output$eq_model$html)
  })
  out
})
V_19_05d_n_int <- sum(vapply(V_09_02_env$B_04_01_equation_lst, function(x)
  min(vapply(names(x$versions), V_09_02_env$B_02_03_rank_fn, integer(1))) <=
    V_09_02_env$B_02_03_rank_fn("2.5b"), logical(1)))
V_19_05e_new_int <- lengths(regmatches(V_19_05c_s5$html, gregexpr(
  'class="eq-flag eq-new"', V_19_05c_s5$html))) - 1L
V_01_03_check_fn(
  "stage 2 shows our full system, every item NEW, none from theirs",
  V_19_05e_new_int == V_19_05d_n_int &&
    lengths(regmatches(V_19_05c_s5$html, gregexpr('class="eq-flag eq-changed"',
                                                  V_19_05c_s5$html))) == 1L &&
    !grepl("Trend Growth", V_19_05c_s5$html, fixed = TRUE) &&
    grepl("n_t", V_19_05c_s5$html, fixed = TRUE) &&
    !grepl("l_t", V_19_05c_s5$html, fixed = TRUE),
  sprintf("%d items, %d NEW", V_19_05d_n_int, V_19_05e_new_int))

V_19_06_seen_chr <- local({
  out_chr <- ""
  shiny::testServer(file.path(V_01_01_here_dir, ".."), {
    session$setInputs(stage = "2.5e")
    out_chr <<- as.character(output$eq_model$html)
  })
  out_chr
})
V_01_03_check_fn(
  "stage 5's Equations tab shows each pill, and ours under each changed one",
  lengths(regmatches(V_19_06_seen_chr, gregexpr('class="eq-flag eq-changed"',
                                                V_19_06_seen_chr))) ==
    sum(V_19_01_pill_chr == "changed") + 1L &&
    lengths(regmatches(V_19_06_seen_chr, gregexpr("was, in our model",
                                                  V_19_06_seen_chr))) ==
    sum(V_19_01_pill_chr == "changed"),
  sprintf("%d changed pills, %d was-lines",
          lengths(regmatches(V_19_06_seen_chr, gregexpr(
            'class="eq-flag eq-changed"', V_19_06_seen_chr))) - 1L,
          lengths(regmatches(V_19_06_seen_chr, gregexpr(
            "was, in our model", V_19_06_seen_chr)))))

################################################################################
## V_99: Verdict ###############################################################
################################################################################
# Note: The exit status is what gates a commit.

message("")
if (V_01_02_fail_int == 0L) {
  message("ALL CHECKS PASSED.")
} else {
  message(V_01_02_fail_int, " CHECK(S) FAILED.")
  quit(status = 1L)
}

#--------------------------------- Script End ---------------------------------#
