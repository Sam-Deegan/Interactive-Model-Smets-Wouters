################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Smets-Wouters Model: Interactive Shiny App                                 ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Open app.R in RStudio and click Run App, or from this folder:
##     shiny::runApp()
##   Needs R 4.1 or later with shiny, bslib and ggplot2 installed. A hosted
##   copy runs in the browser at
##   https://sam-deegan.com/toy-models/smets-wouters/
##   The stage selector builds up the lecture's treatment of the model:
##     2.5a  the frictions of Whelan part 11, switched off one at a time
##     2.5b  our Smets-Wouters model, assembled
##     2.5c  impulse responses
##     2.5d  what explains what: variance decompositions
##     2.5e  our model against Smets and Wouters' own
##   Stages 2.5a to 2.5d run OUR MODEL: Smets and Wouters' structure under
##   the module's assumptions (no trend growth, utility separable in hours)
##   and in its letters, at their posterior mode. Stage 2.5e sets it beside
##   Smets and Wouters' own model, which R/model.R keeps as estimated. All
##   text (stories, prompts, equations, notation) lives in B_02 and B_04.
##
## Inputs:
##   R/model.R (the solver) and R/toolkit.R (shared layout and helpers),
##   both sourced automatically by Shiny. www/ holds the two images.
##
## Version:
##   B_03_09_version_chr; history in CHANGELOG.md; git tag vX.Y.Z.
##
## Outputs:
##   None. Save PNG buttons under each figure write through
##   the toolkit's T_02_03c_export_fn.
##
## Packages:
##   shiny, bslib, ggplot2.
##
## References:
##   Smets, F. and Wouters, R. (2007). Shocks and Frictions in US Business
##     Cycles: A Bayesian DSGE Approach. American Economic Review 97(3).
##   Whelan, K. MA Advanced Macroeconomics, part 11 (The Smets-Wouters
##     Model); cited inline as [W11 nn], nn the slide.
##   Gali, J. (1999). Technology, Employment, and the Business Cycle.
##     American Economic Review 89(1), for hours after a technology shock.

#-------------------------------- Script Begin --------------------------------#

################################################################################
## A: Table of Contents ########################################################
################################################################################
# Note: Section C (the model) is R/model.R so the slides can source it alone.
#
#   B: Setup
#     B_00  Packages
#     B_01  Parameters: sliders, help text, the posterior mode
#     B_02  Stages and worked examples
#     B_03  Presentation: palette, sizes, exports, figure names, version
#     B_04  Text: equations, notation, notes, the caveat, guidance, scope
#     B_05  The scorecard
#   C: Model (R/model.R)
#   D: Drawing
#     D_01  Shared furniture
#     D_02  Stage 2.5a, the frictions
#     D_03  Stage 2.5c, impulse responses
#     D_04  Stage 2.5d, variance decompositions
#     D_05  Stage 2.5a, the bridge to part 7
#     D_06  Stage 2.5e, our model against Smets and Wouters
#     D_09  Fallback
#   E: User Interface
#   F: Server
#   G: Run

################################################################################
## B: Setup ####################################################################
################################################################################
# Note: Packages and every soft-coded value.

#### B_00: Packages ############################################################
# Note: Shiny, bslib for the layout, ggplot2 for the figures.

###### B_00_01: Load Packages ##################################################
# Note: All three run under shinylive.

library(shiny)
library(bslib)
library(ggplot2)

#### B_01: Parameters ##########################################################
# Note: Which parameters get a slider, the note under each, and the
#   posterior mode every stage starts from.

###### B_01_01: Which Parameters Get a Slider ##################################
# Note: Seven of the nineteen estimated parameters, the ones the lecture
#   argues about; the rest stay at the posterior mode. Every range brackets
#   the posterior 5-95 interval in [W11 12]. The letter in brackets is the
#   module's, as the Notation panel writes it: theta_p and theta_w for the
#   Calvo probabilities (Lecture 2.2), phi_pi and rho_i for the rule
#   (Lecture 2.3); habit stays lambda. Smets and Wouters' own letters (xi_p,
#   xi_w, r_pi, rho) appear only in stage 2.5e. The labels are HTML with the
#   letters as entities, never literal Greek, because shinylive sources R/
#   under a C locale; &phi; is the curly phi MathJax draws for \varphi.

B_01_01_slider_lst <- list(
  chabb     = list(label = "Habit Persistence (&lambda;)",
                   min = 0, max = 0.95, step = 0.01),
  csadjcost = list(label = "Investment Adjustment Cost (&phi;)",
                   min = 0.1, max = 15, step = 0.1),
  cprobp    = list(label = "Calvo Prices (&theta;<sub>p</sub>)",
                   min = 0.05, max = 0.95, step = 0.01),
  cprobw    = list(label = "Calvo Wages (&theta;<sub>w</sub>)",
                   min = 0.05, max = 0.95, step = 0.01),
  cindp     = list(label = "Price Indexation (&iota;<sub>p</sub>)",
                   min = 0, max = 0.99, step = 0.01),
  crpi      = list(label = paste0("Policy Response to Inflation ",
                                  "(&phi;<sub>&pi;</sub>)"),
                   min = 1.01, max = 3, step = 0.01),
  crr       = list(label = "Policy Smoothing (&rho;<sub>i</sub>)",
                   min = 0, max = 0.98, step = 0.01))

###### B_01_01b: The Grey Note Under Each Slider ###############################
# Note: What the slider does, in a sentence, and the posterior mode as Table
#   1A prints it on [W11 12]. The table truncates the authors' mode file to
#   two decimals, so the box beside the slider (R/model.R C_01_01) can read
#   up to 0.01 above the number quoted here.

B_01_01b_help_lst <- list(
  chabb     = paste(
    "How much last quarter's consumption sets this quarter's. Zero is no",
    "habit. Posterior mode 0.71."),
  csadjcost = paste(
    "How costly it is to change investment quickly. Near 0.1, new capital",
    "comes on line at once. Posterior mode 5.48."),
  cprobp    = paste(
    "The chance a firm cannot reset its price this quarter. Near 0.05,",
    "prices are almost flexible. Posterior mode 0.65."),
  cprobw    = paste(
    "The same for wages. Posterior mode 0.73."),
  cindp     = paste(
    "How far a price not reset follows last quarter's inflation. Zero",
    "leaves the Phillips curve forward-looking. Posterior mode 0.22."),
  crpi      = paste(
    "The long-run response of the policy rate to inflation. The slider",
    "stops just above one, the Taylor principle. Posterior mode 2.03."),
  crr       = paste(
    "How much of last quarter's rate carries over. Posterior mode 0.81."))

###### B_01_02: The Posterior Mode #############################################
# Note: Pulled from the model, so there is one place a calibration can be
#   wrong. Stages 2.5a to 2.5d solve our model at Smets and Wouters'
#   posterior mode (R/model.R C_01_01b), the mode because their variance
#   decomposition (Figure 1, [W11 15] to [W11 17]) is computed at it. The
#   seven sliders' values are read off it, so a preset restates only the
#   parameter it moves.

B_01_02_default_lst <- C_01_01b_ours_lst

B_01_02_slider_mode_lst <- B_01_02_default_lst[names(B_01_01_slider_lst)]

###### B_01_02b: The Smets-Wouters Side of the Bridge ##########################
# Note: The same parameter list, solved as Smets and Wouters' own model. At
#   the default this is C_01_01_default_lst exactly, the calibration the
#   Dynare checks hold to 1e-10.

B_01_02b_sw_fn <- function(par_lst) {
  utils::modifyList(par_lst, list(model_chr = "sw"))
}

#### B_02: Stages and Examples #################################################
# Note: The stage list, the numbers the stories quote, and the worked
#   examples.

###### B_02_01: Stages #########################################################
# Note: Whelan's own order in part 11: the frictions ([W11 10]), then the
#   responses and what the estimated model attributes to what. Stage 2.5b
#   assembles our model whole, with every equation flagged new; stage 2.5e
#   is the bridge to the model Smets and Wouters estimated, side by side
#   with no path between the two. "Smets-Wouters" takes a hyphen: an en dash
#   cannot be translated under the C locale shinylive runs in.

B_02_01_stage_vec <- c(
  "Stage 1: The Frictions, Switched Off"          = "2.5a",
  "Stage 2: Our Smets-Wouters Model"              = "2.5b",
  "Stage 3: Impulse Responses"                    = "2.5c",
  "Stage 4: What Explains What"                   = "2.5d",
  "Stage 5: Our Model Against Smets and Wouters"  = "2.5e")

###### B_02_01b: The Numbers the Stories Quote #################################
# Note: Computed, never typed. Every number a worked example or a note
#   quotes is read off the solver here, at our model's posterior mode and,
#   for stage 2.5e, at Smets and Wouters' own. tests/verify_against_sw2007.R
#   V_18 re-derives each one.

B_02_01b_word_fn <- function(n_int) {
  words_vec <- c("zero", "one", "two", "three", "four", "five", "six",
                 "seven", "eight", "nine", "ten", "eleven", "twelve")
  if (n_int >= 0L && n_int <= 12L) words_vec[[n_int + 1L]] else
    as.character(n_int)
}

B_02_01b_num_lst <- local({
  ours_lst <- C_03_02_solve_fn(B_01_02_default_lst)
  sw_lst   <- C_03_02_solve_fn(B_01_02b_sw_fn(B_01_02_default_lst))
  em_df    <- C_04_01_irf_fn(ours_lst, "em", 25L)
  sw_em_df <- C_04_01_irf_fn(sw_lst, "em", 25L)
  nh_df    <- C_04_01_irf_fn(C_06_02_off_fn(B_01_02_default_lst, "habit"),
                             "em", 25L)
  na_df    <- C_04_01_irf_fn(C_06_02_off_fn(B_01_02_default_lst, "adjcost"),
                             "em", 25L)
  fp_df    <- C_04_01_irf_fn(C_06_02_off_fn(B_01_02_default_lst, "stickyp"),
                             "em", 25L)
  f_lst    <- C_04_02_fevd_fn(ours_lst)
  sw_f_lst <- C_04_02_fevd_fn(sw_lst)
  hours_fn <- function(over_lst) C_04_01_irf_fn(C_03_02_solve_fn(
    utils::modifyList(B_01_02_default_lst, over_lst)), "ea", 2L)$lab[1L]
  list(
    trough_num    = min(em_df$y),
    trough_int    = which.min(em_df$y) - 1L,
    sw_trough_num = min(sw_em_df$y),
    sw_trough_int = which.min(sw_em_df$y) - 1L,
    nohabit_int   = which.min(nh_df$y) - 1L,
    adjcost_ratio = min(na_df$y) / min(em_df$y),
    flex_cut_num  = 1 - min(fp_df$y) / min(em_df$y),
    flex_int      = which.min(fp_df$y) - 1L,
    rate_zero_int = which(em_df$r <= 0)[1L] - 1L,
    hours_ea_num  = C_04_01_irf_fn(ours_lst, "ea", 2L)$lab[1L],
    rbc_hours_num = hours_fn(list(chabb = 0, csadjcost = 0.1, cindp = 0,
                                  cprobp = 0.05, cprobw = 0.05)),
    q40_ew_num    = f_lst[["40"]]["y", "ew"],
    q1_eg_num     = f_lst[["1"]]["y", "eg"],
    q1_eb_num     = f_lst[["1"]]["y", "eb"],
    q100_ea_num   = f_lst[["100"]]["y", "ea"],
    sw_q1_eg_num  = sw_f_lst[["1"]]["y", "eg"],
    sw_q1_eb_num  = sw_f_lst[["1"]]["y", "eb"],
    real_off_num  = hours_fn(list(chabb = 1e-06, csadjcost = 0.1,
                                  czcap = 0.99, cindp = 1e-06,
                                  cindw = 1e-06)),
    nom_off_num   = hours_fn(list(cprobp = 0.05, cprobw = 0.05)),
    both_off_num  = C_04_01_irf_fn(C_06_04_rbc_limit_fn(B_01_02_default_lst),
                                   "ea", 2L)$lab[1L])
})

###### B_02_02: Worked Examples ################################################
# Note: One card per stage in the main window. A preset sets the sliders; it
#   cannot move the stage. Each stage opens on its own first example. Story
#   order and wording follow CONVENTIONS.md 3 to 5.

B_02_02_example_lst <- local({
  n_lst <- B_02_01b_num_lst
  pc_fn <- function(x_num) sprintf("%.0f per cent", 100 * x_num)
  list(

  posterior = list(
    label = "The Posterior Mode",
    stage = "2.5a",
    values = B_01_02_slider_mode_lst,
    story = paste0(
      "Our model at Smets and Wouters' posterior mode from [W11 12]. Start ",
      "here and read the output ",
      "response to a policy tightening: it does not bite hardest on impact, ",
      "it bites hardest ", B_02_01b_word_fn(n_lst$trough_int),
      " quarters later. Nothing in the RBC of Part 7 ",
      "does that. Now drag habit persistence to zero and watch the trough ",
      "walk back to quarter 0. The faded line is where you started.")),

  nohabit = list(
    label = "Take Away Habit",
    stage = "2.5a",
    values = utils::modifyList(B_01_02_slider_mode_lst, list(chabb = 0)),
    story = paste0(
      "Consumption stops looking backwards. The trough moves to quarter ",
      n_lst$nohabit_int, " and ",
      "gets deeper. Whelan's phrase on [W11 10] is that the frictions ",
      "'throw sand in the wheels'; this is the sand, removed. Note what it ",
      "does <em>not</em> change - the long-run response is much the same. ",
      "Habit changes the timing, not the destination.")),

  fastinv = list(
    label = "Let Investment Move",
    stage = "2.5a",
    values = utils::modifyList(B_01_02_slider_mode_lst,
                               list(csadjcost = 0.1)),
    story = paste0(
      "Adjustment costs to nearly nothing, so new capital can come on line ",
      "at once. This is the friction doing the most work in the model: the ",
      "trough is ", sprintf("%.1f", n_lst$adjcost_ratio),
      " times as deep. Whelan's [W11 6] motivates it as 'an ",
      "adjustment cost function that limits the amount of new investment ",
      "that can come on line immediately' - here is what it is worth.")),

  money = list(
    label = "A Policy Shock",
    stage = "2.5c",
    values = B_01_02_slider_mode_lst,
    story = paste0(
      "[W11 19]. Output and hours fall, inflation falls a little, and the ",
      "interest rate is back to zero within ",
      B_02_01b_word_fn(n_lst$rate_zero_int), " quarters. The shock is one ",
      "posterior standard deviation, so the panels are comparable across ",
      "shocks without a note about units.")),

  flexprice = list(
    label = "Make Prices Flexible",
    stage = "2.5c",
    values = utils::modifyList(B_01_02_slider_mode_lst,
                               list(cprobp = 0.05)),
    story = paste0(
      "Almost every firm resets its price every quarter. The output trough ",
      "shrinks by ", pc_fn(n_lst$flex_cut_num),
      if (n_lst$flex_int == n_lst$trough_int) ", in the same quarter," else
        sprintf(" and comes at quarter %d,", n_lst$flex_int),
      " but money is ",
      "far from neutral, because wages are still sticky - which is the ",
      "point of having both.")),

  explain = list(
    label = "What Explains GDP",
    stage = "2.5d",
    values = B_01_02_slider_mode_lst,
    story = paste0(
      "[W11 15]'s question, asked of our model. Read the Q1 bar and the Q100 ",
      "bar as two different ",
      "questions. In the quarter, output (y) is a demand story: exogenous ",
      "spending leads and the risk premium comes second, together ",
      pc_fn(n_lst$q1_eg_num + n_lst$q1_eb_num),
      ". From ten quarters on the wage mark-up and productivity take ",
      "over. Productivity reaches ", pc_fn(n_lst$q100_ea_num),
      ", a large share but far from ",
      "the whole of it, as it is in the RBC.")),

  rbc = list(
    label = "All Frictions Off: Close to the RBC",
    stage = "2.5a",
    values = utils::modifyList(B_01_02_slider_mode_lst, list(
      chabb = 0, csadjcost = 0.1, cindp = 0, cprobp = 0.05, cprobw = 0.05)),
    story = paste0(
      "Every friction a slider reaches, switched off: no habit, investment ",
      "free to move, no price indexation, and prices and wages reset almost ",
      "every quarter. This is <em>close to</em> Lecture 2.1's RBC, not the ",
      "RBC itself: seven shocks, the Taylor rule, capacity utilisation, ",
      "wage indexation and the fixed cost all remain. A positive technology ",
      "shock now moves hours by ", sprintf("%+.2f", n_lst$rbc_hours_num),
      " on impact, ", if (n_lst$rbc_hours_num > 0) "up, the RBC's sign"
      else "still down", ", against ", sprintf("%.2f", n_lst$hours_ea_num),
      " as estimated. The hours figure shows it takes both halves: the real ",
      "frictions off alone give ", sprintf("%.2f", n_lst$real_off_num),
      ", flexible prices and wages alone ", sprintf("%.2f", n_lst$nom_off_num),
      ". Read the scorecard's At Your Settings column: this is what the ",
      "frictions were buying.")),

  assembled = list(
    label = "Our Model, Assembled",
    stage = "2.5b",
    values = B_01_02_slider_mode_lst,
    story = paste0(
      "Every block of the lecture at once, in the module's letters, at ",
      "Smets and Wouters' posterior mode. A policy tightening takes output ",
      "to ", sprintf("%.2f", n_lst$trough_num), " at quarter ",
      n_lst$trough_int, "; a technology shock moves hours by ",
      sprintf("%.2f", n_lst$hours_ea_num), " on impact. In the quarter ",
      "spending and the risk premium explain ",
      pc_fn(n_lst$q1_eg_num + n_lst$q1_eb_num), " of output; at forty ",
      "quarters the wage mark-up explains ", pc_fn(n_lst$q40_ew_num),
      ". The table below holds it against Whelan's six results.")),

  theirs = list(
    label = "Ours Against Theirs",
    stage = "2.5e",
    values = B_01_02_slider_mode_lst,
    story = paste0(
      "The same deep parameters, two structures. Blue is the model the ",
      "lecture teaches; light blue is Smets and Wouters' own, with trend ",
      "growth and hours in the Euler equation. After a policy tightening ",
      "output bottoms at ", sprintf("%.2f", n_lst$trough_num),
      " in ours and ", sprintf("%.2f", n_lst$sw_trough_num),
      " in theirs. In the quarter, spending and the risk premium explain ",
      pc_fn(n_lst$q1_eg_num + n_lst$q1_eb_num), " of output in ours and ",
      pc_fn(n_lst$sw_q1_eg_num + n_lst$sw_q1_eb_num),
      " in theirs, as on Whelan's chart. The Equations tab ",
      "marks every equation that differs.")))
})

###### B_02_03: Stage Order ####################################################
# Note: The equations and notation panels show what the student has reached,
#   which needs the stages ordered, and "2.5a" is not a number. Read off
#   B_02_01 rather than restated.

B_02_03_rank_fn <- function(stage_chr) {
  out_int <- match(stage_chr, unname(B_02_01_stage_vec))
  if (length(out_int) != 1L || is.na(out_int)) 1L else out_int
}

###### B_02_03b: The Assembly Stage ############################################
# Note: Stage 2.5b presents our completed model as a whole, so every item in
#   force there is flagged new, including those first shown at 2.5a. The
#   toolkit's T_06_03 rule cannot produce that (an item is new only at its
#   first key), so F_06_02 applies this one exception. Stage 2.5e's changed
#   pills and "was" lines are unaffected.

B_02_03b_assembly_chr <- "2.5b"

#### B_03: Presentation ########################################################
# Note: Palette, colours, sizes, exports, figure names, panels, version.

###### B_03_01: Palette ########################################################
# Note: The toolkit's, not a local one. Series run blue, green, navy
#   (T_01_02_series_vec main / second / third); a comparator is light blue
#   (compare); annotation is grey (annot). Navy is furniture: text, the
#   trough marker's label, axis titles. Kept as a short alias.

B_03_01_palette_vec <- c(
  navy   = T_01_01_palette_vec[["navy"]],
  muted  = T_01_01_palette_vec[["muted"]],
  zero   = T_01_01_zero_chr,
  ground = T_01_01_palette_vec[["ground"]])

###### B_03_02: Shock Colours ##################################################
# Note: Seven shocks, one colour each, keyed by shock name so the stacked
#   bars in one figure and the lines in another name the same thing the same
#   way. The three demand shocks take the three series colours in [W11 18]'s
#   own order (risk premium, spending, investment), because D_03_02 draws
#   exactly those three. The other four appear only as segments of the
#   stacked shares and take the remaining tokens of the deck's palette. All
#   seven differ, and no two neighbours in B_03_04's stack order are alike.
#   See CONVENTIONS.md 7.

B_03_02_shock_col_vec <- c(
  eb    = T_01_02_series_vec[["main"]],     # blue
  eg    = T_01_02_series_vec[["second"]],   # green
  eqs   = T_01_02_series_vec[["third"]],    # navy
  ea    = T_01_02_series_vec[["compare"]],  # light blue
  epinf = "#A8D6B6",                        # dublin_diverging light green
  ew    = T_01_02_series_vec[["annot"]],    # grey
  em    = "#7FAED4")                        # dublin_sequential mid blue

###### B_03_02b: The Name's Colour Inside Each Segment #########################
# Note: White on the four dark fills, navy on the three light ones, so every
#   name inside a stacked segment can be read.

B_03_02b_shock_text_vec <- c(
  eb = "#FFFFFF", eg = "#FFFFFF", eqs = "#FFFFFF", ew = "#FFFFFF",
  ea = T_01_01_palette_vec[["navy"]], epinf = T_01_01_palette_vec[["navy"]],
  em = T_01_01_palette_vec[["navy"]])

###### B_03_03: Figure Sizes ###################################################
# Note: Millimetre geometry for a figure written for the deck. No figure may
#   be squarer than 3:2, because a squarer panel is height-limited inside a
#   beamer 16:9 box and shrinks rather than filling the column; the rule is
#   about the exported PNG, which is what V_11 measures, not the facets
#   inside it. `panels` is the multi-panel size (D_03_01, D_03_02).

B_03_03_size_lst <- list(
  wide   = list(width_mm = 160, height_mm = 90),
  panels = list(width_mm = 160, height_mm = 100))

###### B_03_04: The Deck's Own Shock Order #####################################
# Note: A presentation order, not the model's. C_01_02_shock_vec is the
#   column order of Psi and of every decomposition matrix and cannot move;
#   this is the order Lecture 2.5 names the shocks in, so a stacked bar
#   reads bottom to top in the order the room heard them. B_03_02 is keyed
#   by shock name, so reordering the stack moves no colour off its shock.

B_03_04_deck_order_vec <- c("ea", "eg", "eb", "eqs", "epinf", "ew", "em")

###### B_03_05: The Deck's Own Horizons ########################################
# Note: Four bars, the horizons the deck's slot asks for. The builder's own
#   default is Whelan's six ([W11 15]: Q1, Q2, Q4, Q10, Q40, Q100), so the
#   horizons are an argument and the deck export passes this vector.

B_03_05_deck_horizon_vec <- c(1L, 4L, 10L, 40L)

###### B_03_05b: The App's Own Horizons ########################################
# Note: Four bars on screen too, and not the deck's four. At half width six
#   bars are too narrow to carry the shocks' names inside them, and the
#   names are what replace a legend. Q1, Q10 and Q100 stay, so stage 2.5d's
#   worked example, which reads the Q1 bar against the Q100 bar, is still
#   true of what is drawn.

B_03_05b_app_horizon_vec <- c(1L, 4L, 10L, 100L)

###### B_03_06: What a Save PNG Button Writes ##################################
# Note: The deck's size, not the browser's. A full-width figure in the deck
#   is 2:1, rendered at 1600 x 800 px at 200 dpi, and every figure here is
#   full width (pair = FALSE). The numbers are what T_02_03c_export_fn
#   writes, restated so the on-screen card can be drawn at the same 2:1 and
#   the button can print the size beside itself; V_12 checks the two agree.
#   The handler calls the same builder the screen calls, at the export size.

B_03_06_export_lst <- list(width_px = 1600L, height_px = 800L)

###### B_03_07: What Each Figure Is Called #####################################
# Note: {app}-{stage}-{figure}.png, lower case, hyphenated, no spaces and no
#   date; the stage's dots become hyphens because a dot in a file name reads
#   as an extension. The card header lives here too, so the id, the file and
#   the title on screen are set in one row. The `_nom` figures are the
#   nominal half (inflation and the rate) of a four-panel response, drawn as
#   a second card beside the real half (output and hours).

B_03_07_figure_lst <- list(
  plot_friction   = list(stage_chr = "2.5a", slug_chr = "friction",
                         title_chr = "Output After a Policy Tightening"),
  plot_ladder     = list(stage_chr = "2.5a", slug_chr = "friction-ladder",
                         title_chr = "Each Friction Switched Off in Turn"),
  plot_irf        = list(stage_chr = "2.5c", slug_chr = "impulse-response",
                         title_chr = "The Response: Output and Hours"),
  plot_irf_nom    = list(stage_chr = "2.5c",
                         slug_chr = "impulse-response-nominal",
                         title_chr = "The Response: Inflation and the Rate"),
  plot_demand     = list(stage_chr = "2.5c", slug_chr = "demand-shocks",
                         title_chr = "Three Demand Shocks: Output and Hours"),
  plot_demand_nom = list(stage_chr = "2.5c", slug_chr = "demand-shocks-nominal",
                         title_chr = paste("Three Demand Shocks: Inflation",
                                           "and the Rate")),
  plot_fevd       = list(stage_chr = "2.5d", slug_chr = "variance-shares",
                         title_chr = "Shares of the Forecast Error Variance"),
  plot_bridge     = list(stage_chr = "2.5a", slug_chr = "rbc-bridge",
                         title_chr = "Hours After a Technology Shock"),
  plot_s2_irf     = list(stage_chr = "2.5b", slug_chr = "our-model-response",
                         title_chr = "Our Model: The Response to One Shock"),
  plot_s2_fevd    = list(stage_chr = "2.5b", slug_chr = "our-model-shares",
                         title_chr = "Our Model: What Explains Output"),
  plot_vs_irf     = list(stage_chr = "2.5e", slug_chr = "ours-vs-sw-response",
                         title_chr = "Ours Against Theirs: Output and Hours"),
  plot_vs_irf_nom = list(stage_chr = "2.5e",
                         slug_chr = "ours-vs-sw-response-nominal",
                         title_chr = paste("Ours Against Theirs: Inflation",
                                           "and the Rate")),
  plot_vs_fevd    = list(stage_chr = "2.5e", slug_chr = "ours-vs-sw-shares",
                         title_chr = "Ours Against Theirs: Shares, Q1 and Q10"),
  plot_vs_fevd_long = list(stage_chr = "2.5e",
                           slug_chr = "ours-vs-sw-shares-long",
                           title_chr = paste("Ours Against Theirs: Shares,",
                                             "Q40 and Q100")))

B_03_07_file_fn <- function(id_chr) {
  spec_lst <- B_03_07_figure_lst[[id_chr]]
  sprintf("smets-wouters-%s-%s.png",
          gsub(".", "-", spec_lst$stage_chr, fixed = TRUE), spec_lst$slug_chr)
}

###### B_03_08: Which Panels Each Half Draws ###################################
# Note: The real side and the nominal side of Whelan's four panels in [W11 19]
#   and [W11 20], two to a figure, so a pair of cards shows all four.

B_03_08_real_vec    <- c("y", "lab")
B_03_08_nominal_vec <- c("pinf", "r")

###### B_03_09: Version ########################################################
# Note: Shown in the footer; CHANGELOG.md has the history.

B_03_09_version_chr <- "1.0.6"

###### B_03_10: Source Repository ##############################################
# Note: The GitHub repo, linked from the footer.

B_03_10_repo_chr <- paste0("https://github.com/Sam-Deegan/",
                        "Interactive-Model-Smets-Wouters")

#### B_04: Text ################################################################
# Note: Everything the student reads: equations, notation, notes, the
#   caveat, guidance and the scope line.

###### B_04_00: c_3 and the Risk Premium's Sigma at the Mode ###################
# Note: For the notation panel and the note under the equations. Computed,
#   so the text cannot drift from the solver (C_04_00_sigma_fn). Our
#   model's, since stages 2.5a to 2.5d run it.

B_04_00_c3_num   <- C_01_03_derived_fn(B_01_02_default_lst)$c3
B_04_00_sigb_num <- B_01_02_default_lst$sdb / B_04_00_c3_num

###### B_04_01: Equations ######################################################
# Note: One list item per equation, not per stage. `versions` is keyed by
#   the stage a version first applies; the panel shows the latest version at
#   or before the stage on screen and marks what is new or changed. `notes`
#   is the in-words, keyed the same way. The "2.5a" version of every model
#   equation is our model in the module's letters (Lectures 2.1 to 2.3: n_t
#   for hours, theta_p and theta_w, sigma the elasticity itself, x_t the
#   output gap, phi_pi, phi_x, phi_Delta-x and rho_i, v an innovation). The
#   "2.5e" version is Smets and Wouters' own, in their letters, and exists
#   only where theirs differs, so stage 2.5e's changed pills are exactly the
#   equations that differ; V_16 and V_19 in the tests hold that to the two
#   solvers. Habit stays lambda, i_t investment and r_t the policy rate.

B_04_01_trend_chr <- sprintf("%.2f", C_01_01_default_lst$ctrend)

B_04_01_equation_lst <- list(

  ## --- The model's equations ------------------------------------------------
  list(
    group = "model", label = "Aggregate Production Function",
    versions = list(
      "2.5a" = paste0("y_t = \\phi_p\\left(\\alpha k^s_t + (1-\\alpha) n_t",
                      " + \\varepsilon^a_t\\right)"),
      "2.5e" = paste0("y_t = \\phi_p\\left(\\alpha k^s_t + (1-\\alpha) l_t",
                      " + \\varepsilon^a_t\\right)")),
    notes = list(
      "2.5a" = paste("Output is Cobb-Douglas in capital in use and hours,",
                     "scaled by one plus the fixed-cost share. Fixed costs",
                     "are the reason the scale factor is there at all."),
      "2.5e" = paste("The same equation in the paper's letter for hours,",
                     "l<sub>t</sub>. Nothing else differs."))),
  list(
    group = "model", label = "Capital in Use",
    versions = list(
      "2.5a" = "k^s_t = k_{t-1} + z_t, \\qquad z_t = z_1 r^k_t"),
    notes = list(
      "2.5a" = paste("The capital a firm actually runs is last period's",
                     "installed stock worked harder or softer. Utilisation",
                     "is chosen, and running capital harder costs output.",
                     "The same in both models."))),
  list(
    group = "model", label = "The Marginal Product of Capital",
    versions = list(
      "2.5a" = "r^k_t = -(k^s_t - n_t) + w_t",
      "2.5e" = "r^k_t = -(k^s_t - l_t) + w_t"),
    notes = list(
      "2.5a" = paste("The rental rate rises with the wage and falls as",
                     "capital in use is added per worker. Capital in use",
                     "k<sup>s</sup><sub>t</sub>, as in the paper's eq (9) and",
                     "the authors' code."),
      "2.5e" = paste("The paper's letter for hours. Its eq (11) prints",
                     "installed capital k<sub>t</sub> here; its eq (9) and",
                     "the authors' code use capital in use, and so do both",
                     "sides of this bridge."))),
  list(
    group = "model", label = "The Resource Constraint",
    versions = list(
      "2.5a" = paste0("y_t = c_y c_t + i_y i_t + z_y z_t + \\varepsilon^g_t,",
                      " \\qquad i_y = \\delta k_y"),
      "2.5e" = paste0("y_t = c_y c_t + i_y i_t + z_y z_t + \\varepsilon^g_t,",
                      " \\qquad i_y = (\\gamma - 1 + \\delta) k_y")),
    notes = list(
      "2.5a" = paste("Output is spent on consumption, investment, the",
                     "resource cost of running capital harder, and exogenous",
                     "spending. With no trend growth, investment only",
                     "replaces depreciation in the steady state. The shares",
                     "sum to one, which the Diagnostics panel asserts on",
                     "every solve."),
      "2.5e" = paste("With trend growth the steady-state capital stock must",
                     "also grow, so investment's share carries",
                     "&gamma; &minus; 1 on top of depreciation."))),
  list(
    group = "model", label = "Consumption Euler Equation",
    versions = list(
      "2.5a" = paste0("c_t = c_1 c_{t-1} + (1-c_1) E_t c_{t+1}",
                      " - c_3\\left(r_t - E_t \\pi_{t+1}",
                      " + \\varepsilon^b_t\\right)"),
      "2.5e" = paste0("c_t = c_1 c_{t-1} + (1-c_1) E_t c_{t+1}",
                      " + c_2\\left(l_t - E_t l_{t+1}\\right)",
                      " - c_3\\left(r_t - E_t \\pi_{t+1}",
                      " + \\varepsilon^b_t\\right)")),
    notes = list(
      "2.5a" = paste("The paper's eq (2) with utility separable in hours, so",
                     "hours drop out. Habit puts last",
                     "period's consumption on the right, which is what makes",
                     "the output response to a policy shock hump-shaped",
                     "rather than worst on impact. The risk premium",
                     "&epsilon;<sup>b</sup><sub>t</sub> sits inside the",
                     "real-rate bracket, so a rise in it lowers consumption",
                     "as a rise in the real rate would."),
      "2.5e" = paste("The paper's eq (2) as printed. Its utility is not",
                     "separable in hours, so expected hours growth enters",
                     "with weight c<sub>2</sub>, which vanishes when",
                     "&sigma;<sub>c</sub> = 1."))),
  list(
    group = "model", label = "Investment Equation",
    versions = list(
      "2.5a" = paste0("i_t = i_1 i_{t-1} + (1-i_1) E_t i_{t+1} + i_2 q_t",
                      " + \\varepsilon^i_t")),
    notes = list(
      "2.5a" = paste("Adjustment costs put investment's own lag in, so new",
                     "capital cannot come on line at once. The trough ladder",
                     "on stage 1 finds this the friction doing the most",
                     "work. The same form in both models; the weights",
                     "differ through their definitions below."))),
  list(
    group = "model", label = "Tobin's q",
    versions = list(
      "2.5a" = paste0("q_t = q_1 E_t q_{t+1} + (1-q_1) E_t r^k_{t+1}",
                      " - \\left(r_t - E_t \\pi_{t+1}",
                      " + \\varepsilon^b_t\\right)")),
    notes = list(
      "2.5a" = paste("The paper's eq (4), as solved. The value of installed",
                     "capital is the discounted stream of expected rentals.",
                     "The same rise in the risk premium that lowers",
                     "consumption lowers q, one for one with the real",
                     "rate."))),
  list(
    group = "model", label = "Capital Accumulation",
    versions = list(
      "2.5a" = "k_t = k_1 k_{t-1} + (1-k_1) i_t + k_2 \\varepsilon^i_t"),
    notes = list(
      "2.5a" = paste("Depreciation sets how much of the stock",
                     "survives. The investment shock enters here as well as in",
                     "the Euler equation, which is what makes it a shock to",
                     "the productivity of new capital rather than to its",
                     "quantity."))),
  list(
    group = "model", label = "The Price Mark-Up",
    versions = list(
      "2.5a" = "\\mu^p_t = \\alpha(k^s_t - n_t) + \\varepsilon^a_t - w_t",
      "2.5e" = "\\mu^p_t = \\alpha(k^s_t - l_t) + \\varepsilon^a_t - w_t"),
    notes = list(
      "2.5a" = paste("The wedge between the marginal product of labour and",
                     "the real wage. It is minus real marginal cost,",
                     "&mu;<sup>p</sup><sub>t</sub> = &minus;(mc<sub>t</sub>",
                     "&minus; p<sub>t</sub>), so the Phillips curve's",
                     "&minus;&pi;<sub>3</sub>&mu;<sup>p</sup><sub>t</sub> is",
                     "Lecture 2.2's real-marginal-cost term."),
      "2.5e" = "The paper's letter for hours; nothing else differs.")),
  list(
    group = "model", label = "New Keynesian Phillips Curve (NKPC)",
    versions = list(
      "2.5a" = paste0("\\pi_t = \\pi_1 \\pi_{t-1} + \\pi_2 E_t \\pi_{t+1}",
                      " - \\pi_3 \\mu^p_t + \\varepsilon^p_t")),
    notes = list(
      "2.5a" = paste("Indexation supplies the lagged inflation term the",
                     "purely forward-looking curve of Lecture 2.3 has no room",
                     "for. That term is what gives inflation its own",
                     "persistence."))),
  list(
    group = "model", label = "The Wage Mark-Up",
    versions = list(
      "2.5a" = paste0("\\mu^w_t = w_t - \\left(\\sigma_l n_t",
                      " + \\frac{1}{\\sigma(1 - \\lambda)}",
                      "\\left(c_t - \\lambda c_{t-1}\\right)",
                      "\\right)"),
      "2.5e" = paste0("\\mu^w_t = w_t - \\left(\\sigma_l l_t",
                      " + \\frac{1}{1 - \\lambda/\\gamma}",
                      "\\left(c_t - \\frac{\\lambda}{\\gamma} c_{t-1}\\right)",
                      "\\right)")),
    notes = list(
      "2.5a" = paste("The gap between the real wage and the marginal rate of",
                     "substitution. With utility separable in hours the",
                     "marginal utility of consumption carries its curvature,",
                     "1/&sigma;, into the marginal rate of substitution."),
      "2.5e" = paste("The paper's eq (12). Its non-separable utility makes",
                     "the consumption term's weight 1/(1 &minus;",
                     "&lambda;/&gamma;), with no &sigma;<sub>c</sub>; the two",
                     "agree only when &sigma;<sub>c</sub> = 1."))),
  list(
    group = "model", label = "The Wage Equation",
    versions = list(
      "2.5a" = paste0("w_t = w_1 w_{t-1} + (1-w_1)\\left(E_t w_{t+1}",
                      " + E_t \\pi_{t+1}\\right) - w_2 \\pi_t + w_3 \\pi_{t-1}",
                      " - w_4 \\mu^w_t + \\varepsilon^w_t")),
    notes = list(
      "2.5a" = paste("Calvo on wages, with indexation to past inflation",
                     "between resets. w<sub>2</sub> multiplies",
                     "<em>current</em> inflation",
                     "and w<sub>3</sub> <em>lagged</em>, which is the",
                     "paper's own numbering; the",
                     "circulating Dynare replication swaps the two labels."))),
  list(
    group = "model", label = "Monetary Policy (MP) Rule",
    versions = list(
      "2.5a" = paste0("r_t = \\rho_i r_{t-1} + (1-\\rho_i)",
                      "\\left(\\phi_\\pi \\pi_t",
                      " + \\phi_x x_t\\right)",
                      " + \\phi_{\\Delta x}\\left(x_t - x_{t-1}\\right)",
                      " + \\varepsilon^r_t"),
      "2.5e" = paste0("r_t = \\rho r_{t-1} + (1-\\rho)\\left(r_\\pi \\pi_t",
                      " + r_Y (y_t - y^p_t)\\right)",
                      " + r_{\\Delta y}\\left[(y_t - y^p_t)",
                      " - (y_{t-1} - y^p_{t-1})\\right]",
                      " + \\varepsilon^r_t")),
    notes = list(
      "2.5a" = paste("Inertia, a response to inflation, and a response to the",
                     "output gap in level and in change: Lecture 2.3's rule",
                     "with smoothing. The gap is measured",
                     "against potential, which is why the model has to carry",
                     "a second economy."),
      "2.5e" = paste("The paper's eq (14) in its own letters: &rho; for",
                     "&rho;<sub>i</sub>, r<sub>&pi;</sub>, r<sub>Y</sub> and",
                     "r<sub>&Delta;y</sub> for the three responses, and the",
                     "gap written out."))),

  ## --- Assumptions ----------------------------------------------------------
  list(
    group = "assumption", label = "Potential Output",
    versions = list(
      "2.5a" = paste0("x_t \\equiv y_t - y^p_t, \\qquad y^p_t:\\;",
                      " \\text{the model less its policy rule,",
                      " with }\\; \\mu^p_t = \\mu^w_t = 0",
                      " \\;\\text{ and }\\; r^f_t \\;\\text{ for }\\;",
                      " r_t - E_t \\pi_{t+1}"),
      "2.5e" = paste0("y^p_t:\\; \\text{the model less its policy rule,",
                      " with }\\; \\mu^p_t = \\mu^w_t = 0",
                      " \\;\\text{ and }\\; r^f_t \\;\\text{ for }\\;",
                      " r_t - E_t \\pi_{t+1}")),
    notes = list(
      "2.5a" = paste("Potential output is the same real economy solved a",
                     "second time with both mark-ups held at their steady",
                     "state, which is what fully flexible prices and wages",
                     "deliver, so the mark-up shocks drop out. Its own real",
                     "rate replaces the policy-set one. Lecture 2.3's",
                     "y<sup>n</sup><sub>t</sub> relaxed price stickiness",
                     "alone; this one relaxes the wage margin too, so",
                     "x<sub>t</sub> here is a wider gap than 2.3's."),
      "2.5e" = paste("The same object; the paper writes the gap out as",
                     "y<sub>t</sub> &minus; y<sup>p</sup><sub>t</sub> and",
                     "gives it no letter."))),
  list(
    group = "assumption", label = "Trend Growth",
    versions = list(
      "2.5e" = paste0("\\gamma = 1 + \\bar{\\gamma}/100, \\qquad",
                      " \\bar{\\gamma} = ", B_04_01_trend_chr,
                      " \\;\\text{ per cent a quarter}")),
    notes = list(
      "2.5e" = paste("Smets and Wouters' economy grows along a deterministic",
                     "trend, and every equation is written in deviations",
                     "from it; our model has no trend at all, which is",
                     "Lectures 2.1 to 2.3's constant steady state. This is",
                     "where every &gamma; in their weights comes from."))),
  list(
    group = "assumption", label = "Five AR(1) Shock Processes",
    versions = list(
      "2.5a" = paste0("\\varepsilon^x_t = \\rho_x \\varepsilon^x_{t-1}",
                      " + v^x_t, \\qquad",
                      "\\varepsilon^g_t = \\rho_g \\varepsilon^g_{t-1}",
                      " + v^g_t + \\rho_{ga} v^a_t"),
      "2.5e" = paste0("\\varepsilon^x_t = \\rho_x \\varepsilon^x_{t-1}",
                      " + \\eta^x_t, \\qquad",
                      "\\varepsilon^g_t = \\rho_g \\varepsilon^g_{t-1}",
                      " + \\eta^g_t + \\rho_{ga}\\eta^a_t")),
    notes = list(
      "2.5a" = paste("Productivity, the risk premium, investment and policy",
                     "are plain AR(1). Exogenous spending also takes the",
                     "productivity innovation directly, which is the paper's",
                     "allowance for net exports. v is the innovation, as in",
                     "Lecture 2.3."),
      "2.5e" = "The paper writes the innovation &eta;.")),
  list(
    group = "assumption", label = "Two ARMA(1,1) Mark-Up Shocks",
    versions = list(
      "2.5a" = paste0("\\varepsilon^p_t = \\rho_p \\varepsilon^p_{t-1}",
                      " + v^p_t - \\mu_p v^p_{t-1}, \\qquad",
                      "\\varepsilon^w_t = \\rho_w \\varepsilon^w_{t-1}",
                      " + v^w_t - \\mu_w v^w_{t-1}"),
      "2.5e" = paste0("\\varepsilon^p_t = \\rho_p \\varepsilon^p_{t-1}",
                      " + \\eta^p_t - \\mu_p \\eta^p_{t-1}, \\qquad",
                      "\\varepsilon^w_t = \\rho_w \\varepsilon^w_{t-1}",
                      " + \\eta^w_t - \\mu_w \\eta^w_{t-1}")),
    notes = list(
      "2.5a" = paste("Both mark-up shocks subtract their own lagged",
                     "innovation. Whelan's reason is that this is an attempt",
                     "to get at temporary price-level shocks. Note the source",
                     "collision: &mu;<sub>p</sub> with a subscript is this",
                     "moving-average coefficient, &mu;<sup>p</sup><sub>t</sub>",
                     "with a superscript is the price",
                     "mark-up, and [W11 7] writes both."),
      "2.5e" = "The paper writes the innovation &eta;.")),
  list(
    group = "assumption", label = "The Innovations",
    versions = list(
      "2.5a" = "v^x_t \\sim N(0, \\sigma_x^2), \\;\\text{ i.i.d.}",
      "2.5e" = "\\eta^x_t \\sim N(0, \\sigma_x^2), \\;\\text{ i.i.d.}"),
    notes = list(
      "2.5a" = paste("Seven innovations, one per shock process. Every impulse",
                     "response on this app is scaled to one posterior",
                     "standard deviation, so the panels are comparable across",
                     "shocks without a note about units. The risk premium's",
                     "&sigma; is in the units of eqs (2) and (4): Table 1B's",
                     "&sigma;<sub>b</sub> divided by c<sub>3</sub>."),
      "2.5e" = paste("The same seven, with the same posterior standard",
                     "deviations, in the paper's letter."))),

  ## --- Solved forms ---------------------------------------------------------
  list(
    group = "solved", label = "The Canonical Form",
    versions = list(
      "2.5a" = paste0("\\Gamma_0 S_t = \\Gamma_1 S_{t-1} + \\Psi v_t",
                      " + \\Pi \\zeta_t")),
    notes = list(
      "2.5a" = paste("Sims's form, which is what the solver is handed. The",
                     "state stacks the sticky-price economy, the shadow",
                     "economy, the seven shocks and one variable per",
                     "expectation, and the forecast errors are what those",
                     "expectations generate. Both models are solved this",
                     "way."))),
  list(
    group = "solved", label = "The Solution",
    versions = list("2.5a" = "S_t = T S_{t-1} + R v_t"),
    notes = list(
      "2.5a" = paste("What the solver returns. Every figure on this app is a",
                     "function of these two matrices and of nothing else."))),
  list(
    group = "solved", label = "Impulse Response at Horizon h",
    versions = list(
      "2.5c" = "\\mathrm{IRF}^x_h = T^{h} R e_x \\sigma_x"),
    notes = list(
      "2.5c" = paste("The path of every variable after one shock, with the",
                     "rest held at zero. Quarter 0 is the impact period, so",
                     "the exponent starts at zero."))),
  list(
    group = "solved", label = "Forecast Error Variance Share",
    versions = list(
      "2.5d" = paste0("s^j_{x,h} = \\frac{\\sigma_x^2\\sum_{u=0}^{h-1}",
                      "\\left[(T^{u}R)_{jx}\\right]^2}",
                      "{\\sum_{v}\\sigma_v^2\\sum_{u=0}^{h-1}",
                      "\\left[(T^{u}R)_{jv}\\right]^2}")),
    notes = list(
      "2.5d" = paste("The share of one variable's h-step forecast error",
                     "variance that one shock accounts for. R loads a",
                     "<em>unit</em>",
                     "innovation, so each shock's squared responses are",
                     "weighted by its own variance; drop the sigmas and a",
                     "shock of 0.14 counts the same as one of 0.52. The",
                     "shares within a horizon sum to one."))),

  ## --- Thresholds and simplifications ---------------------------------------
  list(
    group = "descriptor", label = "Habit in the Consumption Weights",
    versions = list(
      "2.5a" = paste0("c_1 = \\frac{\\lambda}{1 + \\lambda},",
                      " \\quad c_3 = \\frac{\\sigma(1 - \\lambda)}",
                      "{1 + \\lambda}"),
      "2.5e" = paste0("c_1 = \\frac{\\lambda/\\gamma}{1 + \\lambda/\\gamma},",
                      " \\quad c_2 = \\frac{(\\sigma_c - 1)\\,c_{whl}}",
                      "{\\sigma_c (1 + \\lambda/\\gamma)},",
                      " \\quad c_3 = \\frac{1 - \\lambda/\\gamma}",
                      "{\\sigma_c (1 + \\lambda/\\gamma)}")),
    notes = list(
      "2.5a" = paste("Where the habit slider acts. At &lambda; = 0 the weight",
                     "on lagged consumption goes to zero and the hump in the",
                     "output response goes with it. &sigma; is the",
                     "intertemporal elasticity itself, as in Lecture 2.3,",
                     "so it multiplies."),
      "2.5e" = paste("The paper writes &sigma;<sub>c</sub>, the",
                     "<em>inverse</em> elasticity, so it divides:",
                     "&sigma; = 1/&sigma;<sub>c</sub>, the same number. Trend",
                     "growth puts &lambda;/&gamma; where ours has &lambda;,",
                     "and c<sub>2</sub> is the hours weight our model",
                     "lacks."))),
  list(
    group = "descriptor", label = "Adjustment Costs in the Capital Weights",
    versions = list(
      "2.5a" = paste0("i_1 = \\frac{1}{1 + \\beta},",
                      " \\quad i_2 = \\frac{i_1}{\\varphi},",
                      " \\quad k_1 = 1-\\delta,",
                      " \\quad k_2 = \\frac{1 - k_1}{i_2}",
                      " = \\delta(1 + \\beta)\\varphi"),
      "2.5e" = paste0("i_1 = \\frac{1}{1 + \\bar{\\beta}\\gamma},",
                      " \\quad i_2 = \\frac{i_1}{\\gamma^2 \\varphi},",
                      " \\quad k_1 = \\frac{1-\\delta}{\\gamma},",
                      " \\quad k_2 = (1 - k_1)\\,\\gamma^2 \\varphi")),
    notes = list(
      "2.5a" = paste("Where the adjustment-cost slider acts. k<sub>2</sub>",
                     "is the form the derivation gives, as the paper's eq",
                     "(8) prints it: the investment shock is scaled to enter",
                     "the investment equation one for one, so it loads on",
                     "capital by (1 &minus; k<sub>1</sub>)/i<sub>2</sub>."),
      "2.5e" = paste("The authors' code, which produced the paper's",
                     "results, drops eq (8)'s factor (1 +",
                     "&beta;&gamma;<sup>1&minus;&sigma;<sub>c</sub></sup>)",
                     "from k<sub>2</sub>; their side keeps the code's",
                     "form. &beta;&#772;&gamma; is",
                     "&beta;&gamma;<sup>1&minus;&sigma;<sub>c</sub></sup>,",
                     "which is our &beta; when &gamma; = 1."))),
  list(
    group = "descriptor", label = "Utilisation and the q Weight",
    versions = list(
      "2.5a" = paste0("z_1 = \\frac{1-\\psi}{\\psi}, \\qquad",
                      " q_1 = \\frac{1-\\delta}{\\bar{r}^k + 1 - \\delta}",
                      " = \\beta(1-\\delta)"),
      "2.5e" = paste0("z_1 = \\frac{1-\\psi}{\\psi}, \\qquad",
                      " q_1 = \\frac{1-\\delta}{\\bar{r}^k + 1 - \\delta}")),
    notes = list(
      "2.5a" = paste("Neither is on a slider in this app: both are set by",
                     "depreciation, the utilisation cost and the steady-state",
                     "rental rate. With no trend growth",
                     "r&#772;<sup>k</sup> + 1 &minus; &delta; = 1/&beta;, so",
                     "q<sub>1</sub> is simply &beta;(1 &minus; &delta;)."),
      "2.5e" = paste("With trend growth the steady-state rental rate carries",
                     "&gamma;<sup>&sigma;<sub>c</sub></sup>, so q<sub>1</sub>",
                     "does not simplify."))),
  list(
    group = "descriptor", label = "Calvo Prices in the Phillips Curve",
    versions = list(
      "2.5a" = paste0("\\pi_1 = \\frac{\\iota_p}",
                      "{1 + \\beta\\iota_p}, \\quad",
                      "\\pi_2 = \\frac{\\beta}",
                      "{1 + \\beta\\iota_p}, \\quad",
                      "\\pi_3 = \\frac{1}",
                      "{1 + \\beta\\iota_p}\\cdot",
                      "\\frac{(1-\\theta_p)(1 - \\beta\\theta_p)}",
                      "{\\theta_p\\left[(\\phi_p - 1)\\epsilon_p + 1\\right]}"),
      "2.5e" = paste0("\\pi_1 = \\frac{\\iota_p}",
                      "{1 + \\bar{\\beta}\\gamma\\iota_p}, \\quad",
                      "\\pi_2 = \\frac{\\bar{\\beta}\\gamma}",
                      "{1 + \\bar{\\beta}\\gamma\\iota_p}, \\quad",
                      "\\pi_3 = \\frac{1}",
                      "{1 + \\bar{\\beta}\\gamma\\iota_p}\\cdot",
                      "\\frac{(1-\\xi_p)(1 - \\bar{\\beta}\\gamma\\xi_p)}",
                      "{\\xi_p\\left[(\\phi_p - 1)\\epsilon_p + 1\\right]}")),
    notes = list(
      "2.5a" = paste("Where the two price sliders act. A higher Calvo",
                     "probability &theta;<sub>p</sub>, Lecture 2.2's letter,",
                     "flattens the curve through &pi;<sub>3</sub>; higher",
                     "indexation moves weight off expected inflation and onto",
                     "lagged inflation."),
      "2.5e" = paste("The paper's &xi;<sub>p</sub> is our",
                     "&theta;<sub>p</sub>, and trend growth puts",
                     "&beta;&#772;&gamma; where ours has &beta;."))),
  list(
    group = "descriptor", label = "Calvo Wages in the Wage Equation",
    versions = list(
      "2.5a" = paste0("w_1 = \\frac{1}{1+\\beta}, \\;",
                      "w_2 = \\frac{1 + \\beta\\iota_w}",
                      "{1+\\beta}, \\;",
                      "w_3 = \\frac{\\iota_w}{1+\\beta}, \\;",
                      "w_4 = \\frac{(1-\\theta_w)",
                      "(1-\\beta\\theta_w)}",
                      "{(1+\\beta)\\,\\theta_w",
                      "\\left[(\\bar{\\mu}^w - 1)\\epsilon_w + 1\\right]}"),
      "2.5e" = paste0("w_1 = \\frac{1}{1+\\bar{\\beta}\\gamma}, \\;",
                      "w_2 = \\frac{1 + \\bar{\\beta}\\gamma\\iota_w}",
                      "{1+\\bar{\\beta}\\gamma}, \\;",
                      "w_3 = \\frac{\\iota_w}{1+\\bar{\\beta}\\gamma}, \\;",
                      "w_4 = \\frac{(1-\\xi_w)",
                      "(1-\\bar{\\beta}\\gamma\\xi_w)}",
                      "{(1+\\bar{\\beta}\\gamma)\\,\\xi_w",
                      "\\left[(\\bar{\\mu}^w - 1)\\epsilon_w + 1\\right]}")),
    notes = list(
      "2.5a" = paste("The same construction on the wage side. The",
                     "steady-state wage mark-up and the Kimball curvature",
                     "both enter w<sub>4</sub>, so the wage curve is flatter",
                     "than the",
                     "Calvo probability alone would make it."),
      "2.5e" = paste("&xi;<sub>w</sub> for &theta;<sub>w</sub>, and",
                     "&beta;&#772;&gamma; for &beta;."))),
  list(
    group = "descriptor", label = "The Steady State Beneath the Composites",
    versions = list(
      "2.5a" = paste0("\\bar{r}^k = \\beta^{-1} - (1-\\delta),",
                      " \\qquad c_y + i_y + g_y = 1"),
      "2.5e" = paste0("\\bar{\\beta} = \\beta\\gamma^{-\\sigma_c}, \\qquad",
                      " \\bar{r}^k = \\beta^{-1}\\gamma^{\\sigma_c}",
                      " - (1-\\delta), \\qquad c_y + i_y + g_y = 1")),
    notes = list(
      "2.5a" = paste("Every composite above is a function of these. If the",
                     "three expenditure shares stop summing to one the steady",
                     "state does not exist and nothing downstream means",
                     "anything, which is why it is asserted on every solve."),
      "2.5e" = paste("Trend growth raises the steady-state rental rate and",
                     "discounts the future by &beta;&#772; rather than",
                     "&beta;."))),
  list(
    group = "descriptor", label = "Blanchard and Kahn",
    versions = list(
      "2.5a" = paste0("\\text{number of unstable roots of }",
                      "(\\Gamma_0, \\Gamma_1) = \\dim \\zeta_t")),
    notes = list(
      "2.5a" = paste("One forecast error per unstable root, and the errors",
                     "have to span the unstable block. Too few roots and the",
                     "solution is indeterminate; too many and none exists.",
                     "The Diagnostics panel counts both."))),
  list(
    group = "descriptor", label = "The Frictions, Switched Off",
    versions = list(
      "2.5a" = paste0("\\lambda = 0 \\;\\text{ or }\\;",
                      "\\varphi = 0.1 \\;\\text{ or }\\;",
                      "\\psi = 0.99 \\;\\text{ or }\\;",
                      "\\iota_p = 0 \\;\\text{ or }\\;",
                      "\\iota_w = 0 \\;\\text{ or }\\;",
                      "\\theta_p = 0.05 \\;\\text{; all at once: }\\;",
                      "\\lambda = \\iota_p = \\iota_w = 0,\\;",
                      " \\varphi = 0.1,\\; \\psi = 0.99,\\;",
                      "\\theta_p = \\theta_w = 0.05")),
    notes = list(
      "2.5a" = paste0("One deep parameter at a time, with the steady state ",
                     "and ",
                     "every composite rebuilt around it rather than patched; ",
                     "that is what the trough ladder measures. Then all at ",
                     "once, which is as close as this model ",
                      "comes to the RBC of Part 7. Neither half alone flips ",
                      "the sign of the hours response to a technology shock: ",
                      "the real frictions off give ",
                      sprintf("%.2f", B_02_01b_num_lst$real_off_num),
                      ", flexible prices and wages alone give ",
                      sprintf("%.2f", B_02_01b_num_lst$nom_off_num),
                      ", and only both together give ",
                      sprintf("%+.2f", B_02_01b_num_lst$both_off_num), "."))))

###### B_04_02: Equation Group Titles ##########################################
# Note: The four groups of the is-mp-pc app, all four used. "Thresholds and
#   Simplifications" is where the deep-parameter mapping lives, the one
#   block this model has that is not on any slide.

B_04_02_group_vec <- c(
  model      = "Model Equations",
  assumption = "Assumptions",
  solved     = "Solved Forms",
  descriptor = "Thresholds and Simplifications")

###### B_04_03: Notation #######################################################
# Note: One entry per symbol, flat, with `grp` in var / par / tgt / shk and
#   `from` the stage the symbol first appears. The module's letters from
#   stage 2.5a and Smets and Wouters' own from 2.5e, flagged new there.
#   Habit is lambda in both (the paper's equations and Whelan's slides write
#   lambda; only Table 1A writes h); the steady-state wage mark-up is
#   bar-mu^w, since lambda is taken; the state vector is S_t, since x_t is
#   the output gap. B_01_01's slider labels carry the same letters.

B_04_03_notation_lst <- list(

  ## variables
  list(grp = "var", sym = "y_t", txt = "output (log deviation, per cent)",
       from = "2.5a"),
  list(grp = "var", sym = "y^p_t",
       txt = "potential output, prices and wages flexible", from = "2.5a"),
  list(grp = "var", sym = "c_t", txt = "consumption", from = "2.5a"),
  list(grp = "var", sym = "i_t",
       txt = "investment (the paper's letter, not the policy rate)",
       from = "2.5a"),
  list(grp = "var", sym = "k_t", txt = "installed capital", from = "2.5a"),
  list(grp = "var", sym = "k^s_t", txt = "capital in use", from = "2.5a"),
  list(grp = "var", sym = "z_t", txt = "the rate of capacity utilisation",
       from = "2.5a"),
  list(grp = "var", sym = "n_t", txt = "hours, as in Lecture 2.1",
       from = "2.5a"),
  list(grp = "var", sym = "w_t", txt = "the real wage", from = "2.5a"),
  list(grp = "var", sym = "r^k_t",
       txt = "the marginal productivity of capital", from = "2.5a"),
  list(grp = "var", sym = "q_t", txt = "the value of installed capital",
       from = "2.5a"),
  list(grp = "var", sym = "\\pi_t", txt = "inflation", from = "2.5a"),
  list(grp = "var", sym = "r_t",
       txt = "the nominal policy rate (the paper's letter, not i<sub>t</sub>)",
       from = "2.5a"),
  list(grp = "var", sym = "\\mu^p_t",
       txt = paste("the price mark-up over marginal cost, minus real marginal",
                   "cost: &minus;(mc<sub>t</sub> &minus; p<sub>t</sub>)"),
       from = "2.5a"),
  list(grp = "var", sym = "\\mu^w_t", txt = "the wage mark-up", from = "2.5a"),
  list(grp = "var", sym = "r^f_t",
       txt = "the real interest rate with prices and wages flexible",
       from = "2.5a"),
  list(grp = "var", sym = "S_t", txt = "the model's full state vector",
       from = "2.5a"),
  list(grp = "var", sym = "\\zeta_t",
       txt = "the forecast errors, one per expectation", from = "2.5a"),
  list(grp = "var", sym = "l_t",
       txt = "hours, the paper's letter for n<sub>t</sub>", from = "2.5e"),

  ## parameters
  list(grp = "par", sym = "\\alpha", txt = "the capital share", from = "2.5a"),
  list(grp = "par", sym = "\\phi_p", txt = "one plus the fixed-cost share",
       from = "2.5a"),
  list(grp = "par", sym = "\\delta", txt = "the depreciation rate",
       from = "2.5a"),
  list(grp = "par", sym = "\\beta", txt = "the discount factor", from = "2.5a"),
  list(grp = "par", sym = "\\sigma",
       txt = "the intertemporal elasticity of substitution, as in Lecture 2.3",
       from = "2.5a"),
  list(grp = "par", sym = "\\sigma_l",
       txt = "the inverse Frisch elasticity of labour supply", from = "2.5a"),
  list(grp = "par", sym = "\\lambda",
       txt = "the habit parameter (Lecture 2.3's loss weight was also lambda)",
       from = "2.5a"),
  list(grp = "par", sym = "\\varphi", txt = "the investment adjustment cost",
       from = "2.5a"),
  list(grp = "par", sym = "\\psi", txt = "the capacity utilisation cost",
       from = "2.5a"),
  list(grp = "par", sym = "\\theta_p",
       txt = "the probability a price is not reset, as in Lecture 2.2",
       from = "2.5a"),
  list(grp = "par", sym = "\\theta_w",
       txt = "the probability a wage is not reset", from = "2.5a"),
  list(grp = "par", sym = "\\iota_p",
       txt = "price indexation to past inflation", from = "2.5a"),
  list(grp = "par", sym = "\\iota_w", txt = "wage indexation to past inflation",
       from = "2.5a"),
  list(grp = "par", sym = "\\epsilon_p",
       txt = "the Kimball curvature of the price aggregator", from = "2.5a"),
  list(grp = "par", sym = "\\epsilon_w",
       txt = "the Kimball curvature of the wage aggregator", from = "2.5a"),
  list(grp = "par", sym = "\\bar{\\mu}^w",
       txt = "the steady-state wage mark-up", from = "2.5a"),
  list(grp = "par", sym = "\\bar{r}^k",
       txt = "the steady-state rental rate (per quarter)", from = "2.5a"),
  list(grp = "par", sym = "\\rho_i", txt = "the policy rate's own inertia",
       from = "2.5a"),
  list(grp = "par", sym = "\\phi_\\pi",
       txt = "the policy response to inflation", from = "2.5a"),
  list(grp = "par", sym = "\\phi_x",
       txt = "the policy response to the output gap", from = "2.5a"),
  list(grp = "par", sym = "\\phi_{\\Delta x}",
       txt = "the policy response to the change in the gap", from = "2.5a"),
  list(grp = "par", sym = "c_y", txt = "the steady-state consumption share",
       from = "2.5a"),
  list(grp = "par", sym = "i_y", txt = "the steady-state investment share",
       from = "2.5a"),
  list(grp = "par", sym = "g_y", txt = "the steady-state spending share",
       from = "2.5a"),
  list(grp = "par", sym = "z_y",
       txt = "the steady-state utilisation-cost share", from = "2.5a"),
  list(grp = "par", sym = "k_y",
       txt = "the steady-state capital to output ratio", from = "2.5a"),
  list(grp = "par", sym = "c_1",
       txt = "the composite weight on lagged consumption", from = "2.5a"),
  list(grp = "par", sym = "c_3",
       txt = "the composite weight on the real rate", from = "2.5a"),
  list(grp = "par", sym = "i_1",
       txt = "the composite weight on lagged investment", from = "2.5a"),
  list(grp = "par", sym = "i_2", txt = "the composite weight on q",
       from = "2.5a"),
  list(grp = "par", sym = "q_1",
       txt = "the composite weight on expected q", from = "2.5a"),
  list(grp = "par", sym = "z_1",
       txt = "the composite weight on the rental rate in utilisation",
       from = "2.5a"),
  list(grp = "par", sym = "k_1",
       txt = "the composite weight on lagged capital", from = "2.5a"),
  list(grp = "par", sym = "k_2",
       txt = "the composite loading of the investment shock on capital",
       from = "2.5a"),
  list(grp = "par", sym = "\\pi_1",
       txt = "the composite weight on lagged inflation", from = "2.5a"),
  list(grp = "par", sym = "\\pi_2",
       txt = "the composite weight on expected inflation", from = "2.5a"),
  list(grp = "par", sym = "\\pi_3",
       txt = "the composite weight on the price mark-up", from = "2.5a"),
  list(grp = "par", sym = "w_1",
       txt = "the composite weight on the lagged wage", from = "2.5a"),
  list(grp = "par", sym = "w_2",
       txt = "the composite weight on current inflation", from = "2.5a"),
  list(grp = "par", sym = "w_3",
       txt = "the composite weight on lagged inflation in the wage equation",
       from = "2.5a"),
  list(grp = "par", sym = "w_4",
       txt = "the composite weight on the wage mark-up", from = "2.5a"),
  list(grp = "par", sym = "\\Gamma_0",
       txt = "the canonical system's matrix on time t", from = "2.5a"),
  list(grp = "par", sym = "\\Gamma_1",
       txt = "the canonical system's matrix on time t-1", from = "2.5a"),
  list(grp = "par", sym = "\\Psi", txt = "the canonical shock loadings",
       from = "2.5a"),
  list(grp = "par", sym = "\\Pi",
       txt = "the canonical forecast-error loadings", from = "2.5a"),
  list(grp = "par", sym = "T", txt = "the solved transition matrix",
       from = "2.5a"),
  list(grp = "par", sym = "R", txt = "the solved shock loadings",
       from = "2.5a"),
  list(grp = "par", sym = "e_x", txt = "the selection vector for shock x",
       from = "2.5c"),
  list(grp = "par", sym = "\\gamma", txt = "the gross trend growth factor",
       from = "2.5e"),
  list(grp = "par", sym = "\\bar{\\gamma}",
       txt = "trend growth, per cent a quarter", from = "2.5e"),
  list(grp = "par", sym = "\\bar{\\beta}",
       txt = "the growth-adjusted discount factor", from = "2.5e"),
  list(grp = "par", sym = "\\sigma_c",
       txt = "the paper's inverse elasticity, 1/&sigma;", from = "2.5e"),
  list(grp = "par", sym = "\\xi_p",
       txt = "the paper's letter for &theta;<sub>p</sub>", from = "2.5e"),
  list(grp = "par", sym = "\\xi_w",
       txt = "the paper's letter for &theta;<sub>w</sub>", from = "2.5e"),
  list(grp = "par", sym = "\\rho",
       txt = "the paper's letter for &rho;<sub>i</sub>", from = "2.5e"),
  list(grp = "par", sym = "r_\\pi",
       txt = "the paper's letter for &phi;<sub>&pi;</sub>", from = "2.5e"),
  list(grp = "par", sym = "r_Y",
       txt = paste("the paper's letter for &phi;<sub>x</sub> (r<sub>y</sub>",
                   "in Table 1A)"),
       from = "2.5e"),
  list(grp = "par", sym = "r_{\\Delta y}",
       txt = "the paper's letter for &phi;<sub>&Delta;x</sub>", from = "2.5e"),
  list(grp = "par", sym = "c_2",
       txt = "the composite weight on expected hours growth", from = "2.5e"),
  list(grp = "par", sym = "c_{whl}",
       txt = "the steady-state labour-income to consumption ratio",
       from = "2.5e"),

  ## targets, thresholds and shocks
  list(grp = "tgt", sym = "x_t",
       txt = paste("the output gap, y<sub>t</sub> &minus;",
                   "y<sup>p</sup><sub>t</sub>"),
       from = "2.5a"),
  list(grp = "tgt", sym = "h", txt = "the forecast horizon (quarters)",
       from = "2.5c"),
  list(grp = "tgt", sym = "\\mathrm{IRF}^x_h",
       txt = "the response to shock x, h quarters on", from = "2.5c"),
  list(grp = "tgt", sym = "s^j_{x,h}",
       txt = "shock x's share of variable j's h-step forecast error variance",
       from = "2.5d"),
  list(grp = "tgt", sym = "y_t - y^p_t",
       txt = "the output gap as the paper writes it, our x<sub>t</sub>",
       from = "2.5e"),
  list(grp = "shk", sym = "\\varepsilon^a_t", txt = "the productivity shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^b_t",
       txt = paste0("the risk premium shock, inside the real-rate bracket as ",
                    "in the paper's eqs (2) and (4); positive is a rise. ",
                    "Smets and Wouters estimate it in their code's units: ",
                    "Table 1B's &sigma;<sub>b</sub> = 0.24 is the s.d. of ",
                    "c<sub>3</sub>&epsilon;<sup>b</sup>, so here ",
                    "&sigma;(&epsilon;<sup>b</sup>) = 0.24/c<sub>3</sub> = ",
                    sprintf("%.2f", B_04_00_sigb_num)),
       from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^g_t",
       txt = "the exogenous spending shock", from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^i_t", txt = "the investment shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^p_t", txt = "the price mark-up shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^w_t", txt = "the wage mark-up shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^r_t",
       txt = "the monetary policy shock", from = "2.5a"),
  list(grp = "shk", sym = "\\varepsilon^x_t",
       txt = "a shock process, x standing for any of the seven",
       from = "2.5a"),
  list(grp = "shk", sym = "v^x_t",
       txt = "the innovation to shock process x, as in Lecture 2.3",
       from = "2.5a"),
  list(grp = "shk", sym = "\\rho_x",
       txt = "the persistence of shock process x", from = "2.5a"),
  list(grp = "shk", sym = "\\rho_g",
       txt = "the persistence of exogenous spending", from = "2.5a"),
  list(grp = "shk", sym = "\\rho_p",
       txt = "the persistence of the price mark-up shock", from = "2.5a"),
  list(grp = "shk", sym = "\\rho_w",
       txt = "the persistence of the wage mark-up shock", from = "2.5a"),
  list(grp = "shk", sym = "\\rho_{ga}",
       txt = "the productivity innovation's loading on spending",
       from = "2.5a"),
  list(grp = "shk", sym = "\\mu_p",
       txt = "the moving-average term in the price mark-up shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\mu_w",
       txt = "the moving-average term in the wage mark-up shock",
       from = "2.5a"),
  list(grp = "shk", sym = "\\sigma_x",
       txt = paste("the posterior standard deviation of innovation x;",
                   "for the risk premium, &sigma;<sub>b</sub>/c<sub>3</sub>"),
       from = "2.5a"),
  list(grp = "shk", sym = "\\eta^x_t",
       txt = "the paper's letter for the innovation v<sup>x</sup><sub>t</sub>",
       from = "2.5e"))

###### B_04_04: The Note Under the Equations ###################################
# Note: The one caveat a student reading this panel against the slides
#   needs. Two versions: our model's, for stages 2.5a to 2.5d, and the
#   bridge's, for 2.5e.

B_04_04_equation_note_chr <- paste0(
  "The model equations are Smets and Wouters' eqs (1) to (14), in the ",
  "composite form Whelan presents on [W11 3] to [W11 9], under the ",
  "module's assumptions and in its letters: see What the Lecture ",
  "Simplifies. The composites ",
  "(c<sub>1</sub>, c<sub>3</sub>, i<sub>1</sub>, i<sub>2</sub>, ",
  "q<sub>1</sub>, z<sub>1</sub>, k<sub>1</sub>, k<sub>2</sub>, ",
  "&pi;<sub>1</sub> to &pi;<sub>3</sub> and w<sub>1</sub> to w<sub>4</sub>) ",
  "are functions of the deep parameters in [W11 12], and that mapping is on ",
  "no slide. It is under Thresholds and Simplifications above, and in ",
  "R/model.R at C_01_03. The solver computes exactly what is shown. The ",
  "risk premium is the ",
  "paper's, inside the bracket of eqs (2) and (4); Smets and Wouters ",
  "estimate the shock in code units, so Table 1B's &sigma;<sub>b</sub> = ",
  "0.24 is the s.d. of c<sub>3</sub>&epsilon;<sup>b</sup>, and here ",
  "&sigma;(&epsilon;<sup>b</sup>) = 0.24/c<sub>3</sub> = ",
  sprintf("%.2f", B_04_00_sigb_num), " (c<sub>3</sub> = ",
  sprintf("%.3f", B_04_00_c3_num), " at the mode).")

B_04_04b_bridge_note_chr <- paste0(
  "Stage 5 shows Smets and Wouters' own equations, in their letters. ",
  "<strong>Changed</strong> marks each one that differs from ours, and ",
  "the In Words tab shows what it was; <strong>new</strong> marks what ",
  "theirs has and ours lacks. An equation with the same symbols in both, ",
  "such as investment or the Phillips curve, is not marked: its weights ",
  "differ only through their definitions, which are marked under ",
  "Thresholds and Simplifications. Their side keeps their code's ",
  "k<sub>2</sub>, which produced the published results, and is checked ",
  "against Dynare running their model to 1e-12.")

###### B_04_05: Stage Guidance #################################################
# Note: One line under the worked example, telling the reader what to do
#   with the sliders rather than what the figure shows.

B_04_05_guidance_lst <- list(
  "2.5a" = paste(
    "Move one slider at a time and watch where the trough sits. The faded",
    "curve is the loaded example, so the gap between the two is the",
    "friction you just changed. The hours figure below is the walk back",
    "towards Lecture 2.1's RBC; load All Frictions Off to take it."),
  "2.5c" = paste(
    "Choose a shock above the figures. Its two cards are Whelan's four",
    "panels in [W11 19] and [W11 20]: output and hours, then inflation and",
    "the interest rate."),
  "2.5d" = paste(
    "Choose a variable. Each bar is one forecast horizon, and the segments",
    "are the seven shocks. Compare the first bar with the last."),
  "2.5b" = paste(
    "The whole model at once. Choose a shock for the responses; the shares",
    "are output's. Every equation in the Equations tab is flagged new,",
    "because this is where the model is assembled."),
  "2.5e" = paste(
    "Blue is our model, light blue Smets and Wouters' own, at the same",
    "sliders. Choose a shock for the responses and a variable for the",
    "shares; the Equations tab marks every equation that differs."))

###### B_04_06: Scope ##########################################################
# Note: The companion to the RBC app's own scope line. Each app says where
#   it sits, so neither is mistaken for the whole arc.

B_04_06_scope_chr <- paste0(
  "Whelan's Part 11, the full DSGE: fourteen equations plus a flexible-price ",
  "shadow economy, seven shocks, estimated on US data. Stages 1 to 4 run ",
  "the lecture's version, with no trend growth and the module's letters; ",
  "stage 2 assembles it whole and stage 5 sets it beside Smets and ",
  "Wouters' own. Part 7's RBC is a ",
  "separate app; Part 12 builds a financial accelerator on a model of this ",
  "kind and is not built yet. The coefficient mapping comes from the paper, ",
  "not from the slides - see the Equations tab.")

###### B_04_07: What the Lecture Simplifies ####################################
# Note: The caveat: only the departures the lecture makes from Smets and
#   Wouters; the app's comparison is stage 2.5e. Each item names the
#   Equations-tab items (`pills`, by label) where stage 2.5e shows its
#   effect, and V_19 in the tests holds the two to each other. k_2 is not a
#   teaching departure, so it is not an item: where the printed paper and
#   the code differ, our model takes the derived form and their side keeps
#   the code's, and the pill's note says so. Whelan's slide slips are one
#   line at the end.

B_04_07_caveat_title_chr <- "What the Lecture Simplifies"

B_04_07_caveat_lead_chr <- paste(
  "The lecture teaches Smets and Wouters' model under the assumptions and",
  "in the letters this module has used so far. Four departures, and stage",
  "5 shows what each is worth:")

B_04_07_caveat_lst <- list(
  list(id = "trend", title = "No Trend Growth",
       txt = paste(
         "Every steady state is constant, as in Lectures 2.1 to 2.3. The",
         "paper's economy grows at &gamma;, which puts &lambda;/&gamma;,",
         "&beta;&#772;&gamma; and (1 &minus; &delta;)/&gamma; where ours has",
         "&lambda;, &beta; and 1 &minus; &delta;."),
       pills = c("Trend Growth", "The Resource Constraint",
                 "Habit in the Consumption Weights",
                 "Adjustment Costs in the Capital Weights",
                 "Utilisation and the q Weight",
                 "Calvo Prices in the Phillips Curve",
                 "Calvo Wages in the Wage Equation",
                 "The Steady State Beneath the Composites",
                 "The Wage Mark-Up")),
  list(id = "hours", title = "Utility Separable in Hours",
       txt = paste(
         "Hours drop out of the consumption equation (the paper's",
         "c<sub>2</sub> term), and the marginal rate of substitution",
         "carries 1/&sigma;."),
       pills = c("Consumption Euler Equation", "The Wage Mark-Up",
                 "Habit in the Consumption Weights")),
  list(id = "letters", title = "The Module's Letters",
       txt = paste(
         "n<sub>t</sub> for hours, &theta;<sub>p</sub> and",
         "&theta;<sub>w</sub> for Calvo, &sigma; for the elasticity itself",
         "(the paper's 1/&sigma;<sub>c</sub>), &phi;<sub>&pi;</sub>,",
         "&phi;<sub>x</sub>, &phi;<sub>&Delta;x</sub> and",
         "&rho;<sub>i</sub> in the rule, v for an innovation. As in the",
         "paper, i<sub>t</sub> is investment and r<sub>t</sub> the policy",
         "rate, and habit is &lambda;, which was the loss weight in",
         "Lecture 2.3."),
       pills = c("Aggregate Production Function",
                 "The Marginal Product of Capital", "The Price Mark-Up",
                 "Consumption Euler Equation", "The Wage Mark-Up",
                 "Monetary Policy (MP) Rule", "Five AR(1) Shock Processes",
                 "Two ARMA(1,1) Mark-Up Shocks", "The Innovations",
                 "Habit in the Consumption Weights",
                 "Calvo Prices in the Phillips Curve",
                 "Calvo Wages in the Wage Equation")),
  list(id = "potential", title = "Potential Output",
       txt = paste(
         "x<sub>t</sub> is output less output with prices <em>and</em>",
         "wages flexible and no mark-up shocks, a wider gap than Lecture",
         "2.3's y<sub>t</sub> &minus; y<sup>n</sup><sub>t</sub>."),
       pills = c("Potential Output", "Monetary Policy (MP) Rule")))


B_04_07_caveat_whelan_chr <- paste(
  "Whelan's slides slip in three places: installed capital k<sub>t</sub>",
  "where the paper has capital in use [W11 3, 7], no E<sub>t</sub> on",
  "r<sup>k</sup><sub>t+1</sub> [W11 6], and the wage block's signs and",
  "letters [W11 8].")

#### B_05: The Scorecard #######################################################
# Note: The results Whelan states in class, scored on the live model.

###### B_05_01: What the Model Is Held Against #################################
# Note: The same scorecard the test script gates on, shown in the app so a
#   failure is visible in class rather than silent. Two of these are results
#   the model gets right that the RBC gets wrong.

B_05_01_test_lst <- list(
  # [W11 16] reads about 72%: the largest single share and a majority
  list(name = "Price mark-up's share of inflation at Q1",
       data = "the largest single share, about 72%",
       src  = "[W11 16]",
       pass = function(s_lst, f_lst, irf_lst)
         which.max(f_lst[["1"]]["pinf", ]) ==
           match("epinf", C_01_02_shock_vec) &&
           f_lst[["1"]]["pinf", "epinf"] > 0.5),
  list(name = "Monetary policy's share of the fed funds rate at Q1",
       data = "the largest single share",
       src  = "[W11 17]",
       pass = function(s_lst, f_lst, irf_lst)
         which.max(f_lst[["1"]]["r", ]) == match("em", C_01_02_shock_vec)),
  # [W11 15]: spending about 37% at Q1, the risk premium second at about 26%
  list(name = "What drives output in the quarter",
       data = "exogenous spending, then the risk premium",
       src  = "[W11 15]",
       pass = function(s_lst, f_lst, irf_lst)
         identical(order(f_lst[["1"]]["y", ], decreasing = TRUE)[1:2],
                   match(c("eg", "eb"), C_01_02_shock_vec))),
  list(name = "What drives output over a century",
       data = "the wage mark-up",
       src  = "[W11 15]",
       pass = function(s_lst, f_lst, irf_lst)
         which.max(f_lst[["100"]]["y", ]) == match("ew", C_01_02_shock_vec)),
  list(name = "Output after a policy tightening",
       data = "hump-shaped, not worst on impact",
       src  = "[W11 19]",
       pass = function(s_lst, f_lst, irf_lst)
         (which.min(irf_lst$em$y) - 1L) >= 2L),
  list(name = "Hours after a positive technology shock",
       data = "fall",
       src  = "Gal&iacute; (1999); the RBC gets this wrong",
       pass = function(s_lst, f_lst, irf_lst) irf_lst$ea$lab[1L] < 0))

###### B_05_02: Score One Calibration ##########################################
# Note: One logical per row of B_05_01: TRUE matches, FALSE does not, NA when
#   the calibration has no stable solution and nothing can be tested. A row
#   whose test throws counts as not matching.

B_05_02_score_fn <- function(sol_lst) {
  if (!isTRUE(sol_lst$ok_lgl)) {
    return(rep(NA, length(B_05_01_test_lst)))
  }
  fevd_lst <- C_04_02_fevd_fn(sol_lst)
  irf_lst  <- stats::setNames(
    lapply(C_01_02_shock_vec, function(sh_chr)
      C_04_01_irf_fn(sol_lst, sh_chr, 25L)), C_01_02_shock_vec)
  vapply(B_05_01_test_lst, function(t_lst)
    tryCatch(isTRUE(t_lst$pass(sol_lst, fevd_lst, irf_lst)),
             error = function(e) FALSE), logical(1))
}

###### B_05_03: The Table, Two Columns of Verdicts #############################
# Note: "Estimated" is always scored at the posterior mode; "At Your
#   Settings" at the sliders, so a row that stops matching there names the
#   result the friction just changed is responsible for. Stage 2.5e scores
#   Smets and Wouters' model the same two ways beside our Estimated column.
#   Nothing is adjusted to pass.

B_05_03_verdict_fn <- function(pass_lgl) {
  if (is.na(pass_lgl)) {
    return(tags$td(class = "verdict verdict-na", "No solution"))
  }
  tags$td(class = if (pass_lgl) "verdict verdict-yes" else "verdict verdict-no",
          if (pass_lgl) "Matches" else "Does not match")
}

B_05_03_table_fn <- function(col_lst, caption_chr) {
  rows_lst <- lapply(seq_along(B_05_01_test_lst), function(i_int) {
    t_lst <- B_05_01_test_lst[[i_int]]
    do.call(tags$tr, c(
      list(tags$td(t_lst$name), tags$td(t_lst$data)),
      unname(lapply(col_lst, function(v_lgl)
        B_05_03_verdict_fn(v_lgl[[i_int]]))),
      list(tags$td(class = "verdict-src", HTML(t_lst$src)))))
  })
  tagList(
    tags$p(class = "stat-caption", HTML(caption_chr)),
    tags$table(
      class = "table table-sm verdict-table",
      tags$thead(do.call(tags$tr, c(
        list(tags$th("Result"), tags$th("Whelan Says")),
        lapply(names(col_lst), tags$th),
        list(tags$th("Source"))))),
      tags$tbody(rows_lst)))
}

###### B_05_03b: The Caption ###################################################
# Note: Rows our model misses at the estimate while Smets and Wouters'
#   matches are named; none at the mode today, and the sentence says so.

B_05_03b_caption_fn <- function(ours_post_lgl, sw_post_lgl, stage_chr) {
  gap_int <- which(!(ours_post_lgl %in% TRUE) & (sw_post_lgl %in% TRUE))
  gap_chr <- if (length(gap_int) == 0L) {
    paste("At the estimate our model matches every row that Smets and",
          "Wouters' own model matches.")
  } else {
    paste0("Our model misses ", length(gap_int), " result",
           if (length(gap_int) > 1L) "s" else "",
           " that Smets and Wouters' own model matches: ",
           paste(vapply(B_05_01_test_lst[gap_int], `[[`, "", "name"),
                 collapse = "; "), ".")
  }
  if (identical(stage_chr, "2.5e")) {
    return(paste(
      "<strong>Estimated</strong> scores each model at Smets and Wouters'",
      "posterior mode; <strong>At Your Settings</strong> scores their model",
      "at the sliders as they are now.", gap_chr))
  }
  paste(
    "<strong>Estimated</strong> scores our model at Smets and Wouters'",
    "posterior mode. <strong>At Your Settings</strong> scores the same",
    "result at the sliders as they are now, so a row that stops matching",
    "there is the result the parameter you moved is responsible for.",
    gap_chr)
}

################################################################################
## D: Drawing ##################################################################
################################################################################
# Note: Every builder takes a solved model and returns a ggplot. None of them
#   touches a Shiny input, so a deck script can call them directly.

#### D_01: Shared Furniture ####################################################
# Note: The theme alias, the long form of a response, and label helpers.

###### D_01_01: Theme ##########################################################
# Note: The toolkit's theme, T_02_01_theme_fn, kept as a one-line alias so
#   the grid argument reads the same in every builder below.

D_01_01_theme_fn <- function(grid_chr = "h") {
  T_02_01_theme_fn(grid = grid_chr)
}

###### D_01_02: Long Form for a Four-Panel Response ############################
# Note: The four panels Whelan draws in [W11 19] and [W11 20], in his order,
#   set in one place. vars_vec picks which of the four, so the app can draw
#   the real pair and the nominal pair as two figures.

D_01_02_panel_vec <- c(y = "Output", lab = "Hours",
                       pinf = "Inflation", r = "Interest Rate")

D_01_02_long_fn <- function(irf_df, which_chr = "live",
                            vars_vec = names(D_01_02_panel_vec)) {
  keep_vec <- names(D_01_02_panel_vec)[names(D_01_02_panel_vec) %in% vars_vec]
  out_lst <- lapply(keep_vec, function(v_chr) {
    data.frame(quarter_int = seq_len(nrow(irf_df)) - 1L,
               value_num   = irf_df[[v_chr]],
               panel_cat   = factor(D_01_02_panel_vec[[v_chr]],
                                    levels = unname(
                                      D_01_02_panel_vec[keep_vec])),
               which_cat   = which_chr,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out_lst)
}

###### D_01_03: Keep Line Names Apart ##########################################
# Note: Each line is named horizontally at one end. Where two names would
#   touch they are pushed apart by at least gap_num, keeping their order,
#   with the push split so neither drifts far from its own line.

D_01_03_dodge_fn <- function(y_vec, gap_num) {
  ord_int <- order(y_vec)
  z_vec   <- y_vec[ord_int]
  for (pass_int in seq_len(50L)) {
    moved_lgl <- FALSE
    for (i_int in seq_along(z_vec)[-1L]) {
      short_num <- gap_num - (z_vec[i_int] - z_vec[i_int - 1L])
      if (short_num > 1e-12) {
        z_vec[i_int - 1L] <- z_vec[i_int - 1L] - short_num / 2
        z_vec[i_int]      <- z_vec[i_int] + short_num / 2
        moved_lgl <- TRUE
      }
    }
    if (!moved_lgl) break
  }
  out_vec <- y_vec
  out_vec[ord_int] <- z_vec
  out_vec
}

###### D_01_04: Shock Names in Title Case ######################################
# Note: For plot titles only. The model's own labels, C_01_02_shock_lab_vec,
#   are what the segments and lines carry.

D_01_04_title_vec <- c(
  ea = "a Productivity Shock", eb = "a Rise in the Risk Premium",
  eg = "an Exogenous Spending Shock", eqs = "an Investment Shock",
  em = "a Monetary Policy Shock", epinf = "a Price Mark-Up Shock",
  ew = "a Wage Mark-Up Shock")

#### D_02: Stage 2.5a, the Frictions ###########################################
# Note: One friction against the loaded example, and every friction off in
#   turn.

###### D_02_01: One Friction, On and Off #######################################
# Note: [W11 10]. The live line is the sliders; the faded line is the loaded
#   example, the same series at other parameters, so both take the main
#   colour and neither is named (CONVENTIONS.md 8).

D_02_01_friction_fn <- function(sol_lst, ghost_lst = NULL,
                                n_horizon_int = 25L) {

  live_df <- C_04_01_irf_fn(sol_lst, "em", n_horizon_int)
  plot_df <- data.frame(quarter_int = seq_len(n_horizon_int) - 1L,
                        value_num = live_df$y, which_cat = "live")

  if (!is.null(ghost_lst) && isTRUE(ghost_lst$ok_lgl)) {
    ghost_df <- C_04_01_irf_fn(ghost_lst, "em", n_horizon_int)
    plot_df <- rbind(plot_df,
                     data.frame(quarter_int = seq_len(n_horizon_int) - 1L,
                                value_num = ghost_df$y, which_cat = "ghost"))
  }

  trough_int <- which.min(live_df$y) - 1L
  trough_num <- min(live_df$y)

  ggplot(plot_df, aes(quarter_int, value_num, group = which_cat)) +
    T_02_02_zero_fn(v = FALSE) +
    geom_line(aes(alpha = which_cat), colour = T_01_02_series_vec[["main"]],
              linewidth = 0.9) +
    scale_alpha_manual(values = c(live = 1, ghost = T_01_04_ghost_alpha_num),
                       guide = "none") +
    T_02_03_point_fn(trough_int, trough_num, size = 2.6) +
    # label above the trough and to its right, clear of the recovering curve
    annotate("text", x = trough_int + 1.5, y = trough_num,
             label = sprintf("Trough %.2f at quarter %d",
                             trough_num, trough_int),
             hjust = 0, vjust = -0.9, size = 3.6,
             colour = T_01_01_palette_vec[["navy"]]) +
    labs(
         caption = paste("A one-standard-deviation policy shock. Where the",
                         "trough sits is the frictions' doing. Faded: the",
                         "loaded example."),
         x = "Quarters After the Shock",
         y = "Output (%)") +
    D_01_01_theme_fn("h")
}

###### D_02_02: Every Friction, Off in Turn ####################################
# Note: One point per friction: how deep the trough goes and how long it
#   takes to get there. The baseline takes the main colour; each
#   friction-off run is a counterfactual and takes the compare colour. Rows
#   are named on the y axis, so no legend. Short row names leave the panel
#   room at the export's 2:1.

D_02_02_short_vec <- c(habit = "Habit", adjcost = "Adjustment Costs",
                       utilisation = "Utilisation Costs",
                       pindex = "Price Indexation", windex = "Wage Indexation",
                       stickyp = "Sticky Prices")

D_02_02_ladder_fn <- function(par_lst, n_horizon_int = 25L) {

  base_lst <- C_03_02_solve_fn(par_lst)
  rows_lst <- list()

  if (isTRUE(base_lst$ok_lgl)) {
    base_df <- C_04_01_irf_fn(base_lst, "em", n_horizon_int)
    rows_lst[[1L]] <- data.frame(
      label_chr  = "As Estimated",
      trough_num = min(base_df$y),
      peak_int   = which.min(base_df$y) - 1L,
      base_lgl   = TRUE, stringsAsFactors = FALSE)
  }

  for (name_chr in names(C_06_01_friction_lst)) {
    off_lst <- C_06_02_off_fn(par_lst, name_chr)
    if (!isTRUE(off_lst$ok_lgl)) next
    off_df <- C_04_01_irf_fn(off_lst, "em", n_horizon_int)
    rows_lst[[length(rows_lst) + 1L]] <- data.frame(
      label_chr  = paste("No", D_02_02_short_vec[[name_chr]]),
      trough_num = min(off_df$y),
      peak_int   = which.min(off_df$y) - 1L,
      base_lgl   = FALSE, stringsAsFactors = FALSE)
  }

  plot_df <- do.call(rbind, rows_lst)
  if (is.null(plot_df) || nrow(plot_df) == 0L) {
    return(D_09_01_blank_fn("No calibration here solves."))
  }
  plot_df$label_chr <- factor(plot_df$label_chr,
                              levels = plot_df$label_chr[
                                order(plot_df$trough_num)])

  ggplot(plot_df, aes(trough_num, label_chr)) +
    geom_segment(aes(x = 0, xend = trough_num, yend = label_chr),
                 colour = T_01_01_zero_chr, linewidth = 0.5) +
    geom_point(aes(colour = base_lgl), size = 3.2) +
    # quarter label to the right of the point, towards zero
    geom_text(aes(label = sprintf("q%d", peak_int)),
              hjust = -0.6, vjust = -0.9, size = 3.2,
              colour = T_01_02_series_vec[["annot"]]) +
    scale_colour_manual(values = c(`TRUE`  = T_01_02_series_vec[["main"]],
                                   `FALSE` = T_01_02_series_vec[["compare"]]),
                        guide = "none") +
    scale_x_continuous(n.breaks = 4,
                       expand = expansion(mult = c(0.10, 0.04))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
    labs(
         caption = paste("The trough of the output response to a policy",
                         "tightening, with each friction switched off in",
                         "turn; the label is the quarter it happens in."),
         x = "Output at the Trough (%)", y = NULL) +
    D_01_01_theme_fn("v") +
    theme(axis.text.y = element_text(colour = T_01_01_palette_vec[["navy"]]))
}

#### D_03: Stage 2.5c, Impulse Responses #######################################
# Note: One shock's four panels, and the three demand shocks together.

###### D_03_01: Four Panels, One Shock #########################################
# Note: [W11 19] and [W11 20] are exactly this figure for two of the shocks.
#   The shock is a positive innovation of one posterior standard deviation
#   from [W11 13]; live and ghost are the same series, so neither is named.
#   For the risk premium a positive innovation is a rise in the premium, as
#   in the paper's eqs (2) and (4) (R/model.R C_02_01), so output falls and
#   the caption says so. [W11 18] plots the opposite direction, a fall in
#   the premium with output rising; D_03_02 draws that.

D_03_01_cap_fn <- function(shock_chr) {
  if (identical(shock_chr, "eb")) {
    return(paste("A one-standard-deviation RISE in the risk premium, at the",
                 "posterior mode: output falls. The demand-shock figure",
                 "draws a fall, as [W11 18] does. Faded: the loaded",
                 "example."))
  }
  paste("A positive innovation of one posterior standard deviation, at the",
        "posterior mode. Faded: the loaded example.")
}

D_03_01_irf_fn <- function(sol_lst, shock_chr, ghost_lst = NULL,
                           n_horizon_int = 21L,
                           vars_vec = names(D_01_02_panel_vec)) {

  plot_df <- D_01_02_long_fn(
    C_04_01_irf_fn(sol_lst, shock_chr, n_horizon_int), "live", vars_vec)
  if (!is.null(ghost_lst) && isTRUE(ghost_lst$ok_lgl)) {
    plot_df <- rbind(plot_df, D_01_02_long_fn(
      C_04_01_irf_fn(ghost_lst, shock_chr, n_horizon_int), "ghost", vars_vec))
  }
  # Note: Two panels sit side by side, four sit two by two.
  n_row_int <- if (length(unique(plot_df$panel_cat)) <= 2L) 1L else 2L

  ggplot(plot_df, aes(quarter_int, value_num, group = which_cat)) +
    T_02_02_zero_fn(v = FALSE) +
    geom_line(aes(alpha = which_cat), colour = T_01_02_series_vec[["main"]],
              linewidth = 0.8) +
    scale_alpha_manual(values = c(live = 1, ghost = T_01_04_ghost_alpha_num),
                       guide = "none") +
    scale_y_continuous(n.breaks = 4) +
    facet_wrap(~ panel_cat, scales = "free_y", nrow = n_row_int) +
    labs(
         caption = D_03_01_cap_fn(shock_chr),
         x = "Quarters After the Shock", y = "Deviation (%)") +
    D_01_01_theme_fn("h")
}

###### D_03_02: The Three Demand Shocks Together ###############################
# Note: [W11 18] draws risk premium, exogenous spending and investment on
#   one set of axes, the only place in the lecture where shocks are compared
#   rather than shown. The risk premium is drawn falling, as [W11 18] (the
#   paper's Figure 2) draws it: an innovation of -sigma_b, so output, hours,
#   inflation and the rate all rise, as on Whelan's chart; the line's name
#   and the caption both say "falls". Colours are a named vector keyed by
#   shock, never positional (V_09_15 guards it). Names sit at the right-hand
#   end of the Hours panel, where the three paths are still apart at
#   quarter 20, or on the last panel drawn; the x axis runs past quarter 20
#   to make the room.

D_03_02_demand_fn <- function(sol_lst, n_horizon_int = 21L,
                              vars_vec = names(D_01_02_panel_vec)) {

  keep_vec <- c("eb", "eg", "eqs")
  name_vec <- c(eb = "Risk Premium Falls", eg = "Spending", eqs = "Investment")
  sign_vec <- c(eb = -1, eg = 1, eqs = 1)
  var_vec  <- names(D_01_02_panel_vec)[names(D_01_02_panel_vec) %in% vars_vec]
  rows_lst <- lapply(keep_vec, function(sh_chr) {
    sd_num <- C_04_00_sigma_fn(sol_lst)[[sh_chr]]
    irf_df <- C_04_01_irf_fn(sol_lst, sh_chr, n_horizon_int,
                             size_num = sign_vec[[sh_chr]] * sd_num)
    do.call(rbind, lapply(var_vec, function(v_chr)
      data.frame(quarter_int = seq_len(n_horizon_int) - 1L,
                 value_num   = irf_df[[v_chr]],
                 panel_cat   = factor(D_01_02_panel_vec[[v_chr]],
                                      levels = unname(
                                        D_01_02_panel_vec[var_vec])),
                 shock_cat   = C_01_02_shock_lab_vec[[sh_chr]],
                 stringsAsFactors = FALSE)))
  })
  plot_df <- do.call(rbind, rows_lst)
  plot_df$shock_cat <- factor(plot_df$shock_cat,
                              levels = unname(C_01_02_shock_lab_vec[keep_vec]))

  col_vec <- stats::setNames(unname(B_03_02_shock_col_vec[keep_vec]),
                             unname(C_01_02_shock_lab_vec[keep_vec]))

  lab_lvl <- if ("lab" %in% var_vec) D_01_02_panel_vec[["lab"]] else
    utils::tail(levels(plot_df$panel_cat), 1L)
  n_row_int <- if (length(var_vec) <= 2L) 1L else 2L
  end_df  <- plot_df[plot_df$quarter_int == n_horizon_int - 1L &
                       plot_df$panel_cat == lab_lvl, ]
  span_num <- diff(range(plot_df$value_num[plot_df$panel_cat == lab_lvl], 0))
  end_df$label_y   <- D_01_03_dodge_fn(end_df$value_num, 0.24 * span_num)
  end_df$label_chr <- unname(name_vec[keep_vec[
    match(as.character(end_df$shock_cat), C_01_02_shock_lab_vec[keep_vec])]])

  ggplot(plot_df, aes(quarter_int, value_num, colour = shock_cat)) +
    T_02_02_zero_fn(v = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_text(data = end_df, aes(y = label_y, label = label_chr),
              hjust = 0, nudge_x = 0.6, size = 3.2, fontface = "bold") +
    scale_colour_manual(values = col_vec, guide = "none") +
    scale_x_continuous(breaks = seq(0, n_horizon_int - 1L, by = 10L),
                       limits = c(0, n_horizon_int + 17)) +
    scale_y_continuous(n.breaks = 4) +
    coord_cartesian(clip = "off") +
    facet_wrap(~ panel_cat, scales = "free_y", nrow = n_row_int) +
    labs(
         caption = paste("Whelan's [W11 18]: the risk premium FALLS by one",
                         "posterior standard deviation, as he and the paper",
                         "draw it; spending and investment rise by one.",
                         "At the posterior mode."),
         x = "Quarters After the Shock", y = "Deviation (%)") +
    D_01_01_theme_fn("h")
}

#### D_04: Stage 2.5d, Variance Decompositions #################################
# Note: The stacked shares by horizon.

###### D_04_01: Stacked Shares by Horizon ######################################
# Note: [W11 15] to [W11 17]. Horizons are an argument, defaulting to
#   Whelan's six. The stack runs in the deck's shock order (B_03_04),
#   colours are keyed by shock name (B_03_02), and there is no legend: each
#   segment carries its shock's name inside itself, in B_03_02b's colour,
#   and a segment too small for the words carries nothing. label_min_num is
#   the drawing threshold per line of text; the wrap width follows the bar
#   count. Every shock clears the threshold at some horizon, which is what
#   makes the legend unnecessary (V_09). text_size_num is the names' size:
#   2.5 for the deck, larger on a half-width card.

D_04_01_fevd_fn <- function(sol_lst, target_chr = "y",
                            horizon_vec = c(1L, 2L, 4L, 10L, 40L, 100L),
                            label_min_num = 0.055, text_size_num = 2.5) {

  fevd_lst <- C_04_02_fevd_fn(sol_lst, horizon_vec)
  order_vec <- B_03_04_deck_order_vec
  lab_vec   <- unname(C_01_02_shock_lab_vec[order_vec])

  rows_lst <- lapply(names(fevd_lst), function(h_chr) {
    share_vec <- fevd_lst[[h_chr]][target_chr, order_vec]
    data.frame(horizon_cat = factor(paste0("Q", h_chr),
                                    levels = paste0("Q", names(fevd_lst))),
               shock_cat   = factor(lab_vec, levels = lab_vec),
               share_num   = as.numeric(share_vec),
               stringsAsFactors = FALSE)
  })
  plot_df <- do.call(rbind, rows_lst)

  # the name sits at the midpoint of its own segment; geom_col stacks the
  # last level at the bottom, so the cumulative sum runs down the levels
  plot_df$mid_num <- NA_real_
  for (h_chr in levels(plot_df$horizon_cat)) {
    row_int <- which(plot_df$horizon_cat == h_chr)
    up_int  <- row_int[order(plot_df$shock_cat[row_int], decreasing = TRUE)]
    share_vec <- plot_df$share_num[up_int]
    plot_df$mid_num[up_int] <- cumsum(share_vec) - share_vec / 2
  }
  # the name, wrapped, and only where the segment can hold it
  wrap_int <- if (length(horizon_vec) <= 4L) 16L else 12L
  wrap_vec <- vapply(lab_vec, function(nm_chr)
    paste(strwrap(nm_chr, width = wrap_int), collapse = "\n"), character(1))
  lines_vec <- vapply(wrap_vec, function(nm_chr)
    length(strsplit(nm_chr, "\n", fixed = TRUE)[[1L]]), integer(1))

  fits_lgl <- plot_df$share_num >=
    label_min_num * lines_vec[as.character(plot_df$shock_cat)]
  plot_df$label_chr <- ifelse(fits_lgl,
                              wrap_vec[as.character(plot_df$shock_cat)], "")
  text_vec <- stats::setNames(unname(B_03_02b_shock_text_vec[order_vec]),
                              lab_vec)
  plot_df$text_col <- unname(text_vec[as.character(plot_df$shock_cat)])

  title_chr <- c(y = "What Explains Output", pinf = "What Explains Inflation",
                 r = "What Explains the Interest Rate",
                 c = "What Explains Consumption",
                 inve = "What Explains Investment")[[target_chr]]

  ggplot(plot_df, aes(horizon_cat, share_num, fill = shock_cat)) +
    geom_col(width = 0.84) +
    geom_text(aes(y = mid_num, label = label_chr, colour = I(text_col)),
              size = text_size_num, fontface = "bold", lineheight = 0.85) +
    scale_fill_manual(values = unname(B_03_02_shock_col_vec[order_vec]),
                      guide = "none") +
    scale_y_continuous(labels = function(x_num) paste0(100 * x_num, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
         caption = paste("Each bar is one forecast horizon and each segment",
                         "one shock's share of the forecast error variance."),
         x = "Forecast Horizon (Quarters)", y = "Share of Variance (%)") +
    D_01_01_theme_fn("h")
}

#### D_05: Stage 2.5a, the Bridge to Part 7 ####################################
# Note: The figure that joins this app to the RBC app.

###### D_05_01: Hours After a Technology Shock #################################
# Note: Three lines, because two would tell a tidier story than the model
#   supports: at the posterior mode hours fall on impact, agreeing with Gali
#   (1999); with the real frictions off they still fall; only with flexible
#   prices and wages as well does the sign flip towards the RBC's. Three
#   series in the deck's order, all solid, each named at its impact end
#   where the three are apart; no legend. "All Frictions Off" is the real
#   frictions and sticky prices and wages (C_06_04).

D_05_01_bridge_fn <- function(sol_lst, par_lst, n_horizon_int = 21L) {

  level_chr <- c("As Estimated", "Real Frictions Off", "All Frictions Off")

  live_df <- C_04_01_irf_fn(sol_lst, "ea", n_horizon_int)
  plot_df <- data.frame(quarter_int = seq_len(n_horizon_int) - 1L,
                        value_num = live_df$lab,
                        which_cat = level_chr[1L],
                        stringsAsFactors = FALSE)

  add_fn <- function(df_in, sol_in, label_chr) {
    if (!isTRUE(sol_in$ok_lgl)) return(df_in)
    irf_df <- C_04_01_irf_fn(sol_in, "ea", n_horizon_int)
    rbind(df_in, data.frame(quarter_int = seq_len(n_horizon_int) - 1L,
                            value_num = irf_df$lab,
                            which_cat = label_chr,
                            stringsAsFactors = FALSE))
  }
  plot_df <- add_fn(plot_df, C_06_03_all_off_fn(par_lst), level_chr[2L])
  plot_df <- add_fn(plot_df, C_06_04_rbc_limit_fn(par_lst), level_chr[3L])
  plot_df$which_cat <- factor(plot_df$which_cat, levels = level_chr)

  end_df <- plot_df[plot_df$quarter_int == 0L, ]
  span_num <- diff(range(plot_df$value_num, 0))
  end_df$label_y <- D_01_03_dodge_fn(end_df$value_num, 0.09 * span_num)

  ggplot(plot_df, aes(quarter_int, value_num, colour = which_cat)) +
    T_02_02_zero_fn(v = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_text(data = end_df,
              aes(y = label_y, label = as.character(which_cat)),
              hjust = 1, nudge_x = -0.4, size = 3.6, fontface = "bold") +
    scale_colour_manual(values = stats::setNames(
      unname(T_01_02_series_vec[c("main", "second", "third")]), level_chr),
      guide = "none") +
    scale_x_continuous(limits = c(-9.5, n_horizon_int - 1L),
                       breaks = seq(0, n_horizon_int - 1L, by = 4L)) +
    labs(
         caption = paste("A positive technology shock of one posterior",
                         "standard deviation. Gali says hours fall; the RBC",
                         "of Part 7 says they rise. Neither half alone flips",
                         "the sign."),
         x = "Quarters After the Shock", y = "Hours (%)") +
    D_01_01_theme_fn("h")
}

#### D_06: Stage 2.5e, Our Model Against Smets and Wouters #####################
# Note: Side by side only; no path is drawn between the two.

###### D_06_00: The Two Sides' Names and Colours ###############################
# Note: Ours is the main series, theirs the comparator (toolkit T_01_02),
#   and each is named on the figure, never in a legend.

D_06_00_side_vec <- c(ours = "Our Model", sw = "Smets and Wouters")
D_06_00_col_vec  <- stats::setNames(
  unname(T_01_02_series_vec[c("main", "compare")]), unname(D_06_00_side_vec))

###### D_06_01: Responses, Ours and Theirs #####################################
# Note: One shock, the panels vars_vec picks, two lines a panel. The names
#   sit at the right-hand end of the first panel, pushed apart where they
#   would touch; the x axis runs past quarter 20 to make the room.

D_06_01_vs_irf_fn <- function(ours_lst, sw_lst, shock_chr,
                              n_horizon_int = 21L,
                              vars_vec = names(D_01_02_panel_vec)) {

  plot_df <- rbind(
    D_01_02_long_fn(C_04_01_irf_fn(ours_lst, shock_chr, n_horizon_int),
                    D_06_00_side_vec[["ours"]], vars_vec),
    D_01_02_long_fn(C_04_01_irf_fn(sw_lst, shock_chr, n_horizon_int),
                    D_06_00_side_vec[["sw"]], vars_vec))
  plot_df$which_cat <- factor(plot_df$which_cat,
                              levels = unname(D_06_00_side_vec))

  first_lvl <- levels(plot_df$panel_cat)[1L]
  end_df <- plot_df[plot_df$quarter_int == n_horizon_int - 1L &
                      plot_df$panel_cat == first_lvl, ]
  span_num <- diff(range(plot_df$value_num[plot_df$panel_cat == first_lvl],
                         0))
  end_df$label_y <- D_01_03_dodge_fn(end_df$value_num, 0.16 * span_num)
  n_row_int <- if (length(unique(plot_df$panel_cat)) <= 2L) 1L else 2L

  ggplot(plot_df, aes(quarter_int, value_num, colour = which_cat)) +
    T_02_02_zero_fn(v = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_text(data = end_df,
              aes(y = label_y, label = as.character(which_cat)),
              hjust = 0, nudge_x = 0.6, size = 3.2, fontface = "bold") +
    scale_colour_manual(values = D_06_00_col_vec, guide = "none") +
    scale_x_continuous(breaks = seq(0, n_horizon_int - 1L, by = 10L),
                       limits = c(0, n_horizon_int + 17)) +
    scale_y_continuous(n.breaks = 4) +
    coord_cartesian(clip = "off") +
    facet_wrap(~ panel_cat, scales = "free_y", nrow = n_row_int) +
    labs(
         caption = paste("One posterior standard deviation, the same deep",
                         "parameters on both sides. Blue: our model; light",
                         "blue: Smets and Wouters' own."),
         x = "Quarters After the Shock", y = "Deviation (%)") +
    D_01_01_theme_fn("h")
}

###### D_06_02: Variance Shares, Ours and Theirs ###############################
# Note: One variable, two horizons, a pair of bars a shock: ours blue,
#   theirs light blue. Shocks in the deck's order (B_03_04), top to bottom.
#   The two sides are named inside their own bars on the row with the
#   largest share in the first panel, so there is no legend.

D_06_02_vs_fevd_fn <- function(ours_lst, sw_lst, target_chr = "y",
                               horizon_vec = c(1L, 10L)) {

  lab_vec <- unname(C_01_02_shock_lab_vec[B_03_04_deck_order_vec])
  side_lst <- list(ours = ours_lst, sw = sw_lst)
  plot_df <- do.call(rbind, lapply(names(side_lst), function(sd_chr) {
    f_lst <- C_04_02_fevd_fn(side_lst[[sd_chr]], horizon_vec)
    do.call(rbind, lapply(names(f_lst), function(h_chr)
      data.frame(
        horizon_cat = factor(paste0("Q", h_chr),
                             levels = paste0("Q", horizon_vec)),
        shock_cat   = factor(lab_vec, levels = rev(lab_vec)),
        side_cat    = factor(D_06_00_side_vec[[sd_chr]],
                             levels = rev(unname(D_06_00_side_vec))),
        share_num   = as.numeric(f_lst[[h_chr]][target_chr,
                                                B_03_04_deck_order_vec]),
        stringsAsFactors = FALSE)))
  }))

  first_df <- plot_df[plot_df$horizon_cat == levels(plot_df$horizon_cat)[1L], ]
  top_chr  <- as.character(first_df$shock_cat[which.max(first_df$share_num)])
  name_df  <- first_df[as.character(first_df$shock_cat) == top_chr, ]
  name_df$text_col <- ifelse(name_df$side_cat == D_06_00_side_vec[["ours"]],
                             "#FFFFFF", T_01_01_palette_vec[["navy"]])

  title_chr <- c(y = "What Explains Output", pinf = "What Explains Inflation",
                 r = "What Explains the Interest Rate",
                 c = "What Explains Consumption",
                 inve = "What Explains Investment")[[target_chr]]

  ggplot(plot_df, aes(share_num, shock_cat, fill = side_cat)) +
    geom_col(position = position_dodge(width = 0.84), width = 0.8) +
    geom_text(data = name_df,
              aes(x = share_num, label = as.character(side_cat),
                  colour = I(text_col)),
              position = position_dodge(width = 0.84), hjust = 1.05,
              size = 2.6, fontface = "bold") +
    scale_fill_manual(values = D_06_00_col_vec, guide = "none") +
    scale_x_continuous(labels = function(x_num) paste0(100 * x_num, "%"),
                       n.breaks = 4, expand = expansion(mult = c(0, 0.04))) +
    facet_wrap(~ horizon_cat, nrow = 1) +
    labs(
         caption = paste("Each shock's share of the forecast error variance",
                         "at two horizons. Blue: our model; light blue:",
                         "Smets and Wouters' own."),
         x = "Share of Variance (%)", y = NULL) +
    D_01_01_theme_fn("v") +
    theme(axis.text.y = element_text(colour = T_01_01_palette_vec[["navy"]]))
}

#### D_09: Fallback ############################################################
# Note: What a card shows when the model does not solve.

###### D_09_01: A Figure That Says Why There Is No Figure ######################
# Note: A blank panel with an explanation, from the toolkit.

D_09_01_blank_fn <- function(msg_chr) {
  T_02_02_placeholder_fn(msg_chr)
}

################################################################################
## E: User Interface ###########################################################
################################################################################
# Note: The toolkit's page frame: the site's nav bar, a bslib page_sidebar
#   on the toolkit's theme, the stage and the seven sliders in the sidebar,
#   and in the main window the equations card, the worked examples, the
#   stage's note, the readouts, the figures two to a row, the scorecard and
#   the footer.

#### E_01: Sidebar #############################################################
# Note: Control shorthand, the QR code and the sidebar itself.

###### E_01_01: Control Shorthand ##############################################
# Note: A slider, the box beside it for an exact value, the label and the
#   grey note under it, from the toolkit.

E_01_01_ctl_fn <- function(id_chr) {
  T_03_01_control_fn(id_chr, B_01_01_slider_lst, B_01_01b_help_lst,
                     B_01_02_default_lst)
}

###### E_01_02: The QR Code ####################################################
# Note: The toolkit's data URI, so a shinylive export needs no www/ file.

E_01_02_qr_src_chr <- T_07_04_qr_fn()

###### E_01_03: Sidebar ########################################################
# Note: The stage, then the sliders in three groups: the real frictions, the
#   nominal rigidities, and the policy rule. All three groups are open, since
#   every stage uses every slider.

E_01_03_sidebar_lst <- sidebar(
  width = 380,
  radioButtons("stage", "Stage of the Model",
               choices = B_02_01_stage_vec, selected = "2.5a"),
  T_03_05_note_fn(B_04_06_scope_chr),
  accordion(
    open = c("Real Frictions", "Nominal Rigidities", "Monetary Policy"),
    accordion_panel(
      "Real Frictions",
      E_01_01_ctl_fn("chabb"),
      E_01_01_ctl_fn("csadjcost")),
    accordion_panel(
      "Nominal Rigidities",
      E_01_01_ctl_fn("cprobp"),
      E_01_01_ctl_fn("cprobw"),
      E_01_01_ctl_fn("cindp")),
    accordion_panel(
      "Monetary Policy",
      E_01_01_ctl_fn("crpi"),
      E_01_01_ctl_fn("crr"))),
  T_03_05_note_fn(paste(
    "The other estimated parameters stay at their posterior mode.",
    "The Equations tab shows where each slider acts.")),
  actionButton("reset", "Reset Everything",
               class = "btn-outline-secondary btn-sm w-100"),
  T_07_10b_sidebarqr_fn(E_01_02_qr_src_chr))

#### E_02: Main Panel ##########################################################
# Note: Cards, presets, pickers, the caveat and the page.

###### E_02_01: A Figure Card ##################################################
# Note: The toolkit's card, T_07_07f, with its title read off B_03_07 so the
#   id, the file name and the header are set in one place.

E_02_01_card_fn <- function(id_chr) {
  T_07_07f_figcard_fn(id_chr, B_03_07_figure_lst[[id_chr]]$title_chr)
}

###### E_02_01b: A Figure With No Partner ######################################
# Note: Half a row, the same card size as a pair's, so a lone figure is not
#   drawn twice as large as the rest.

E_02_01b_lone_fn <- function(card_obj) {
  layout_columns(col_widths = breakpoints(sm = 12, lg = 6), card_obj)
}

###### E_02_02: Worked-Example Presets #########################################
# Note: The toolkit's preset card: the buttons for the stage on screen and no
#   other, the loaded one marked, and its story underneath.

E_02_02_presets_lst <- T_05_04_presets_fn(
  B_02_02_example_lst, B_02_01_stage_vec, stage_word = "")

###### E_02_03: The Two Pickers ################################################
# Note: The shock and the variable are chosen in the main window, above the
#   figures they drive, because they choose what to draw rather than the
#   model (CONVENTIONS.md 1). Built once with the page; every block that
#   reads one still req()s it, because a browser can run an output before
#   it has sent the radio's first value. Stages 2.5b and 2.5e have their
#   own copies, since an input id may appear once on a page.

E_02_03_shock_ui <- radioButtons(
  "shock", "Choose a Shock", inline = TRUE,
  choices = stats::setNames(C_01_02_shock_vec,
                            C_01_02_shock_lab_vec[C_01_02_shock_vec]),
  selected = "em")

E_02_03_target_ui <- radioButtons(
  "target", "Choose a Variable", inline = TRUE,
  choices = c("Output" = "y", "Inflation" = "pinf",
              "Interest Rate" = "r", "Consumption" = "c",
              "Investment" = "inve"),
  selected = "y")

E_02_03_shock2_ui <- radioButtons(
  "shock2", "Choose a Shock", inline = TRUE,
  choices = stats::setNames(C_01_02_shock_vec,
                            C_01_02_shock_lab_vec[C_01_02_shock_vec]),
  selected = "em")

E_02_03_shock5_ui <- radioButtons(
  "shock5", "Choose a Shock", inline = TRUE,
  choices = stats::setNames(C_01_02_shock_vec,
                            C_01_02_shock_lab_vec[C_01_02_shock_vec]),
  selected = "em")

E_02_03_target5_ui <- radioButtons(
  "target5", "Choose a Variable", inline = TRUE,
  choices = c("Output" = "y", "Inflation" = "pinf",
              "Interest Rate" = "r", "Consumption" = "c",
              "Investment" = "inve"),
  selected = "y")

###### E_02_04: The One Style the Toolkit Lacks ################################
# Note: A table of verdicts. A result that does not match is set in navy
#   and bold; one that cannot be scored is muted.

E_02_04_css_chr <- paste(
  ".verdict-table td, .verdict-table th { font-size: 0.95rem; }",
  ".verdict-table .verdict-no { color: #04204C; font-weight: 700; }",
  ".verdict-table .verdict-na, .verdict-table .verdict-src {",
  "  color: #6C757D; }",
  sep = "\n")

###### E_02_05: What the Lecture Simplifies ####################################
# Note: Static, so it is built once with the page. B_04_07's four items in
#   a short list, each followed by where stage 2.5e shows it, then Whelan's
#   slips in one line.

E_02_05_caveat_ui <- tagList(
  tags$h5(class = "eq-group-title", B_04_07_caveat_title_chr),
  tags$p(HTML(B_04_07_caveat_lead_chr)),
  tags$ul(lapply(B_04_07_caveat_lst, function(it_lst) {
    tags$li(tags$strong(paste0(it_lst$title, ".")), " ", HTML(it_lst$txt))
  })),
  tags$p(class = "stat-caption", HTML(paste(
    "Stage 5, Our Model Against Smets and Wouters, sets the two side by",
    "side and marks every equation that differs."))),
  tags$p(class = "stat-caption", HTML(B_04_07_caveat_whelan_chr)))

###### E_02_06: The Page #######################################################
# Note: The nav bar, then the toolkit's page_sidebar. Two figures to a row;
#   the variance shares and the hours bridge have no partner and take half
#   a row each. Stage 2.5a carries the hours figure whatever preset is
#   loaded; the "All Frictions Off" preset is the one that reads it.

E_02_06_page_ui <- tagList(
  T_07_08b_nav_fn(),
  page_sidebar(
    title        = T_07_09_title_fn(
      HTML("The Smets&ndash;Wouters Model, Solved and Decomposed"),
      E_01_02_qr_src_chr),
    window_title = paste("The Smets-Wouters Model \u00b7",
                         T_07_01_author_chr),
    fillable     = FALSE,
    theme        = T_07_05_theme_fn(),
    sidebar      = E_01_03_sidebar_lst,
    T_07_08_head_fn(),
    tags$head(
      tags$style(HTML(T_05_07_preset_css_chr)),
      tags$style(HTML(E_02_04_css_chr)),
      tags$script(HTML(T_05_05_preset_js_chr))),
    T_07_07j_eqtabs_fn(
      title = textOutput("eq_title", inline = TRUE),
      nav_panel("Equations", uiOutput("eq_model")),
      nav_panel("Notation", uiOutput("eq_notation")),
      nav_panel("In Words", uiOutput("eq_explain")),
      nav_panel("What the Lecture Simplifies", value = "simplifies",
                E_02_05_caveat_ui),
      nav_panel("Diagnostics", uiOutput("diagnostic_ui"))),
    E_02_02_presets_lst,
    uiOutput("prompt"),
    uiOutput("problems"),
    uiOutput("tiles"),
    conditionalPanel(
      "input.stage == '2.5a'",
      T_07_07g_pair_fn(E_02_01_card_fn("plot_friction"),
                       E_02_01_card_fn("plot_ladder")),
      E_02_01b_lone_fn(E_02_01_card_fn("plot_bridge"))),
    conditionalPanel(
      "input.stage == '2.5c'",
      E_02_03_shock_ui,
      T_07_07g_pair_fn(E_02_01_card_fn("plot_irf"),
                       E_02_01_card_fn("plot_irf_nom")),
      T_07_07g_pair_fn(E_02_01_card_fn("plot_demand"),
                       E_02_01_card_fn("plot_demand_nom"))),
    conditionalPanel(
      "input.stage == '2.5d'",
      E_02_03_target_ui,
      E_02_01b_lone_fn(E_02_01_card_fn("plot_fevd"))),
    conditionalPanel(
      "input.stage == '2.5b'",
      E_02_03_shock2_ui,
      T_07_07g_pair_fn(E_02_01_card_fn("plot_s2_irf"),
                       E_02_01_card_fn("plot_s2_fevd"))),
    conditionalPanel(
      "input.stage == '2.5e'",
      E_02_03_shock5_ui,
      T_07_07g_pair_fn(E_02_01_card_fn("plot_vs_irf"),
                       E_02_01_card_fn("plot_vs_irf_nom")),
      E_02_03_target5_ui,
      T_07_07g_pair_fn(E_02_01_card_fn("plot_vs_fevd"),
                       E_02_01_card_fn("plot_vs_fevd_long"))),
    card(
      card_header("What the Model Is Held Against"),
      uiOutput("test_ui")),
    T_07_11_footer_fn(paste0(
      "The model is Smets and Wouters (2007) under the module's ",
      "assumptions, solved at their posterior mode, in the module's ",
      "letters; stage 5 sets it beside their own model. Version ",
      B_03_09_version_chr, "."), repo = B_03_10_repo_chr)))

################################################################################
## F: Server ###################################################################
################################################################################
# Note: Two helpers, then the server function.

#### F_00: Helpers #############################################################
# Note: Objects the server uses that do not need its scope.

###### F_00_01: A Small Helper #################################################
# Note: %||% is in rlang, which is not guaranteed after a shinylive export.

`%||%` <- function(a_val, b_val) if (is.null(a_val)) b_val else a_val

###### F_00_02: Stage 2.5e's Equations Tab, With What Each Item Was ############
# Note: T_06_04's markup, one column, with T_06_06's "was" line (class
#   chg-was) under the equation of every changed item. Only the stage-2.5e
#   panel uses it; the other stages use T_06_04 as it comes.

F_00_02_bridge_model_fn <- function(items_lst, groups_vec) {
  group_fn <- function(grp_chr) {
    rows_lst <- lapply(Filter(function(x) x$group == grp_chr, items_lst),
                       function(x) {
      tags$tr(
        tags$td(class = "eq-label", HTML(x$label), T_06_02_flag_fn(x$status)),
        tags$td(class = "eq-math",
                tags$div(T_06_01_mj_fn(x$tex)),
                if (!is.null(x$was)) {
                  tags$div(class = "chg-was", "was, in our model ",
                           T_06_01_mj_fn(x$was))
                }))
    })
    tags$div(class = "eq-group",
             tags$div(class = "eq-group-title", groups_vec[[grp_chr]]),
             tags$table(class = "eq-table", do.call(tagList, rows_lst)))
  }
  withMathJax(tagList(
    do.call(layout_columns, c(list(col_widths = 12),
                              lapply(names(groups_vec), group_fn))),
    tags$div(
      class = "eq-legend",
      tags$span(class = "eq-flag eq-changed", "changed"),
      " marks an equation of Smets and Wouters' that differs from ours, with",
      " ours underneath; ", tags$span(class = "eq-flag eq-new", "new"),
      " marks one theirs has and ours lacks.")))
}

#### F_01: Server Function #####################################################
# Note: Reactives for the two models, the figures, the exports and the
#   panels.

###### F_01_00: Server #########################################################
# Note: Every output on the page.

F_01_00_server_fn <- function(input, output, session) {

  # --- Figure captions --------------------------------------------------------
  T_07_07d_cap_fn(output)

  # --- Controls ---------------------------------------------------------------
  # Note: The typed box holds the exact value; T_03_04_val_fn reads it before
  #   the slider.
  T_03_02_sync_fn(input, session, B_01_01_slider_lst)
  val_fn <- function(id_chr) T_03_04_val_fn(input, id_chr)

  ###### F_01_01: Assemble the Parameter List ##################################
  # Note: Called once for the live run and once for the ghost, so the two
  #   cannot drift apart. A slider that has not reported yet keeps its mode.

  F_01_01_assemble_fn <- function(values_lst) {
    out_lst <- B_01_02_default_lst
    for (nm_chr in names(B_01_01_slider_lst)) {
      v_num <- values_lst[[nm_chr]]
      if (!is.null(v_num) && !is.na(v_num)) out_lst[[nm_chr]] <- v_num
    }
    out_lst
  }

  ###### F_01_02: Which Example Is Loaded ######################################
  # Note: The toolkit's pattern. Every stage opens on its own first example,
  #   a preset button loads its example, Reset returns the posterior mode and
  #   unloads, and the loaded button is marked through dgPreset.

  F_01_02_loaded_r <- reactiveVal("posterior")

  F_01_02_scn_r <- reactive({
    key_chr <- F_01_02_loaded_r()
    if (is.null(key_chr) || !key_chr %in% names(B_02_02_example_lst)) NULL
    else B_02_02_example_lst[[key_chr]]
  })

  F_01_02_set_loaded_fn <- function(key_chr) {
    F_01_02_loaded_r(if (is.null(key_chr)) "custom" else key_chr)
    session$sendCustomMessage("dgPreset", if (is.null(key_chr)) "" else key_chr)
    invisible(NULL)
  }

  F_01_02_apply_fn <- function(values_lst) {
    for (nm_chr in names(B_01_01_slider_lst)) {
      T_03_03_set_fn(session, B_01_01_slider_lst, nm_chr,
                     values_lst[[nm_chr]])
    }
    invisible(NULL)
  }

  F_01_02_load_fn <- function(key_chr) {
    F_01_02_set_loaded_fn(key_chr)
    F_01_02_apply_fn(utils::modifyList(
      B_01_02_default_lst, B_02_02_example_lst[[key_chr]]$values))
  }

  lapply(names(B_02_02_example_lst), function(key_chr) {
    observeEvent(input[[paste0("preset_", key_chr)]],
                 F_01_02_load_fn(key_chr), ignoreInit = TRUE)
  })

  observeEvent(input$stage, {
    hit_chr <- names(B_02_02_example_lst)[vapply(
      B_02_02_example_lst, function(e_lst) identical(e_lst$stage, input$stage),
      logical(1))]
    if (length(hit_chr) > 0L) F_01_02_load_fn(hit_chr[1L])
    else F_01_02_set_loaded_fn(NULL)
  })

  observeEvent(input$reset, {
    F_01_02_apply_fn(B_01_02_default_lst)
    F_01_02_set_loaded_fn(NULL)
  })

  ###### F_02_01: The Live Model ###############################################
  # Note: A reactiveVal, not a reactive. Moving a slider invalidates the
  #   parameter list twice (the slider, then its box), and a reactiveVal
  #   only fires when its value changes, so the model is solved once per
  #   move. F_02_04 scores B_05_01's results at the sliders.

  F_02_01_par_rv <- reactiveVal(B_01_02_default_lst)
  observe({
    F_02_01_par_rv(F_01_01_assemble_fn(
      stats::setNames(lapply(names(B_01_01_slider_lst), val_fn),
                      names(B_01_01_slider_lst))))
  })

  F_02_01_par_r <- reactive(F_02_01_par_rv())
  F_02_02_sol_r <- reactive(C_03_02_solve_fn(F_02_01_par_r()))
  F_02_04_score_r <- reactive(B_05_02_score_fn(F_02_02_sol_r()))

  ###### F_02_03: The Ghost ####################################################
  # Note: The loaded example's own parameters, solved separately. NULL when
  #   the sliders are already sitting on it, so the figure does not draw the
  #   same line twice.

  F_02_03_ghost_r <- reactive({
    scn_lst <- F_01_02_scn_r()
    if (is.null(scn_lst)) return(NULL)
    ref_lst <- F_01_01_assemble_fn(scn_lst$values)
    if (T_02_03b_ghost_off_fn(F_02_01_par_r(), ref_lst)) return(NULL)
    C_03_02_solve_fn(ref_lst)
  })

  ###### F_02_05: The Smets-Wouters Side #######################################
  # Note: Their model at the same sliders, solved only when stage 2.5e asks
  #   (a reactive is lazy), and scored for its table.

  F_02_05_sw_r <- reactive(C_03_02_solve_fn(B_01_02b_sw_fn(F_02_01_par_r())))
  F_02_06_sw_score_r <- reactive(B_05_02_score_fn(F_02_05_sw_r()))

  ###### F_03_01: Presets, Prompt and Problems #################################
  # Note: The toolkit's preset card, the stage guidance and any solver
  #   problems.

  output$preset_title <- renderUI({
    T_05_06_preset_title_fn(F_01_02_scn_r(), input$stage, B_02_01_stage_vec,
                            stage_word = "Stage")
  })

  output$scenario_story <- renderUI({
    T_05_02_story_fn(F_01_02_scn_r(), B_01_01_slider_lst, B_01_01b_help_lst)
  })

  output$prompt <- renderUI({
    T_07_12_prompt_fn(NULL, input$stage, B_04_05_guidance_lst)
  })

  output$problems <- renderUI({
    T_07_13_problems_fn(F_02_02_sol_r()$problems_chr)
  })

  ###### F_03_02: Readouts #####################################################
  # Note: Three numbers the lecture argues about, at the sliders as they are.

  output$tiles <- renderUI({
    sol_lst <- F_02_02_sol_r()
    if (!isTRUE(sol_lst$ok_lgl)) return(NULL)
    # stage 2.5e reads both sides: the two troughs and their model's score
    if (identical(input$stage, "2.5e")) {
      sw_lst <- F_02_05_sw_r()
      if (!isTRUE(sw_lst$ok_lgl)) return(NULL)
      em_df <- C_04_01_irf_fn(sol_lst, "em", 25L)
      sw_df <- C_04_01_irf_fn(sw_lst, "em", 25L)
      sw_now_lgl <- F_02_06_sw_score_r()
      return(T_04_03_row_fn(
        T_04_01_tile_fn(
          "Output Trough, Our Model",
          paste0(T_02_05_num_fn(min(em_df$y), 2), "%"),
          sprintf("At quarter %d after a tightening",
                  which.min(em_df$y) - 1L)),
        T_04_01_tile_fn(
          "Output Trough, Smets and Wouters",
          paste0(T_02_05_num_fn(min(sw_df$y), 2), "%"),
          sprintf("At quarter %d after a tightening",
                  which.min(sw_df$y) - 1L)),
        T_04_01_tile_fn(
          "Whelan's Results Matched, Theirs",
          sprintf("%d of %d", sum(sw_now_lgl %in% TRUE), length(sw_now_lgl)),
          "Their model at your settings",
          class = if (all(sw_now_lgl %in% TRUE)) "good" else "")))
    }
    em_df   <- C_04_01_irf_fn(sol_lst, "em", 25L)
    ea_df   <- C_04_01_irf_fn(sol_lst, "ea", 2L)
    now_lgl <- F_02_04_score_r()
    T_04_03_row_fn(
      T_04_01_tile_fn(
        "Output Trough After a Tightening",
        paste0(T_02_05_num_fn(min(em_df$y), 2), "%"),
        sprintf("At quarter %d; quarter 0 is the impact",
                which.min(em_df$y) - 1L)),
      T_04_01_tile_fn(
        "Hours on Impact of a Technology Shock",
        paste0(T_02_05_num_fn(ea_df$lab[1L], 2), "%"),
        if (ea_df$lab[1L] < 0) "They fall, as Gal&iacute; finds"
        else "They rise, as in the RBC of Part 7",
        class = if (ea_df$lab[1L] < 0) "good" else "bad"),
      T_04_01_tile_fn(
        "Whelan's Results Matched",
        sprintf("%d of %d", sum(now_lgl), length(now_lgl)),
        sprintf("At your settings; %d of %d as estimated",
                sum(F_05_01_post_lgl %in% TRUE), length(F_05_01_post_lgl)),
        class = if (all(now_lgl %in% TRUE)) "good" else ""))
  })

  ###### F_04_01: One Figure, One Builder ######################################
  # Note: The screen and the save buttons call the same function, so a
  #   download always matches what is on screen. pick_chr is passed in
  #   rather than read from input here, so every read of a radio stays
  #   inside a block that req()s it. The hours bridge is drawn at the
  #   estimate, not the sliders, because its three lines are named "As
  #   Estimated", "Real Frictions Off" and "All Frictions Off".

  F_04_00_est_sol <- C_03_02_solve_fn(B_01_02_default_lst)

  F_04_01_build_fn <- function(id_chr, pick_chr = NULL) {
    sol_lst <- F_02_02_sol_r()
    if (identical(id_chr, "plot_ladder")) {
      return(D_02_02_ladder_fn(F_02_01_par_r()))
    }
    if (!isTRUE(sol_lst$ok_lgl)) {
      return(D_09_01_blank_fn("No stable solution at these values."))
    }
    if (grepl("^plot_vs_", id_chr)) {
      sw_lst <- F_02_05_sw_r()
      if (!isTRUE(sw_lst$ok_lgl)) {
        return(D_09_01_blank_fn(paste("Smets and Wouters' model has no",
                                      "stable solution at these values.")))
      }
      return(switch(id_chr,
        plot_vs_irf       = D_06_01_vs_irf_fn(sol_lst, sw_lst, pick_chr,
                                              vars_vec = B_03_08_real_vec),
        plot_vs_irf_nom   = D_06_01_vs_irf_fn(sol_lst, sw_lst, pick_chr,
                                              vars_vec = B_03_08_nominal_vec),
        plot_vs_fevd      = D_06_02_vs_fevd_fn(sol_lst, sw_lst, pick_chr,
                                               c(1L, 10L)),
        plot_vs_fevd_long = D_06_02_vs_fevd_fn(sol_lst, sw_lst, pick_chr,
                                               c(40L, 100L))))
    }
    switch(id_chr,
      plot_friction   = D_02_01_friction_fn(sol_lst, F_02_03_ghost_r()),
      plot_irf        = D_03_01_irf_fn(sol_lst, pick_chr, F_02_03_ghost_r(),
                                       vars_vec = B_03_08_real_vec),
      plot_irf_nom    = D_03_01_irf_fn(sol_lst, pick_chr, F_02_03_ghost_r(),
                                       vars_vec = B_03_08_nominal_vec),
      plot_demand     = D_03_02_demand_fn(sol_lst,
                                          vars_vec = B_03_08_real_vec),
      plot_demand_nom = D_03_02_demand_fn(sol_lst,
                                          vars_vec = B_03_08_nominal_vec),
      plot_fevd       = D_04_01_fevd_fn(sol_lst, pick_chr,
                                        B_03_05b_app_horizon_vec,
                                        label_min_num = 0.06,
                                        text_size_num = 2.8),
      plot_bridge     = D_05_01_bridge_fn(F_04_00_est_sol, B_01_02_default_lst),
      plot_s2_irf     = D_03_01_irf_fn(sol_lst, pick_chr, F_02_03_ghost_r()),
      plot_s2_fevd    = D_04_01_fevd_fn(sol_lst, "y",
                                        B_03_05b_app_horizon_vec,
                                        label_min_num = 0.06,
                                        text_size_num = 2.8),
      D_09_01_blank_fn("No figure by that name."))
  }

  ###### F_04_02: Drawn on Screen ##############################################
  # Note: Through T_02_01c_draw_fn, which lifts the caption under the card;
  #   res = 96 so the type is drawn at the size a browser draws its own
  #   text. No title on screen: the card header names the figure, so where
  #   the title carried something live (the shock or variable chosen) it
  #   moves to the front of the caption (lead_lgl). The req() on each picker
  #   is load-bearing; V_10 in the tests drives the app with the input unset.
  F_04_02_draw_fn <- function(p_obj, lead_lgl = FALSE) {
    if (inherits(p_obj, "ggplot")) {
      if (lead_lgl && !is.null(p_obj$labels$title)) {
        p_obj$labels$caption <- paste0(p_obj$labels$title, ". ",
                                       p_obj$labels$caption)
      }
      p_obj$labels$title <- NULL
    }
    T_02_01c_draw_fn(p_obj)
  }

  output$plot_friction <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_friction")), res = 96)

  output$plot_ladder <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_ladder")), res = 96)

  output$plot_irf <- renderPlot({
    req(input$shock)
    F_04_02_draw_fn(F_04_01_build_fn("plot_irf", input$shock),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_irf_nom <- renderPlot({
    req(input$shock)
    F_04_02_draw_fn(F_04_01_build_fn("plot_irf_nom", input$shock),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_demand <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_demand")), res = 96)

  output$plot_demand_nom <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_demand_nom")), res = 96)

  output$plot_fevd <- renderPlot({
    req(input$target)
    F_04_02_draw_fn(F_04_01_build_fn("plot_fevd", input$target),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_bridge <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_bridge")), res = 96)

  output$plot_s2_irf <- renderPlot({
    req(input$shock2)
    F_04_02_draw_fn(F_04_01_build_fn("plot_s2_irf", input$shock2),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_s2_fevd <- renderPlot(
    F_04_02_draw_fn(F_04_01_build_fn("plot_s2_fevd"), lead_lgl = TRUE),
    res = 96)

  output$plot_vs_irf <- renderPlot({
    req(input$shock5)
    F_04_02_draw_fn(F_04_01_build_fn("plot_vs_irf", input$shock5),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_vs_irf_nom <- renderPlot({
    req(input$shock5)
    F_04_02_draw_fn(F_04_01_build_fn("plot_vs_irf_nom", input$shock5),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_vs_fevd <- renderPlot({
    req(input$target5)
    F_04_02_draw_fn(F_04_01_build_fn("plot_vs_fevd", input$target5),
                    lead_lgl = TRUE)
  }, res = 96)

  output$plot_vs_fevd_long <- renderPlot({
    req(input$target5)
    F_04_02_draw_fn(F_04_01_build_fn("plot_vs_fevd_long", input$target5),
                    lead_lgl = TRUE)
  }, res = 96)

  ###### F_04_03: Save PNG ####################################################
  # Note: The toolkit's handlers, T_07_07h, writing through
  #   T_02_03c_export_fn at the deck's full-width 1600 x 800 px. The file
  #   names are B_03_07's.

  F_04_03_pick_lst <- list(plot_irf = "shock", plot_irf_nom = "shock",
                           plot_fevd = "target", plot_s2_irf = "shock2",
                           plot_vs_irf = "shock5", plot_vs_irf_nom = "shock5",
                           plot_vs_fevd = "target5",
                           plot_vs_fevd_long = "target5")

  for (id_chr in names(B_03_07_figure_lst)) local({
    this_chr <- id_chr
    pick_chr <- F_04_03_pick_lst[[this_chr]]
    T_07_07h_exports_fn(
      output, this_chr,
      plot_fn = function() {
        if (is.null(pick_chr)) return(F_04_01_build_fn(this_chr))
        req(input[[pick_chr]])
        F_04_01_build_fn(this_chr, input[[pick_chr]])
      },
      stem = sub("\\.png$", "", B_03_07_file_fn(this_chr)),
      pair = FALSE)
  })

  ###### F_05_01: What the Model Is Held Against ###############################
  # Note: Two columns, B_05_03. The posterior mode is scored once at
  #   start-up, since it never moves; the sliders whenever they do.

  F_05_01_post_lgl <- B_05_02_score_fn(C_03_02_solve_fn(B_01_02_default_lst))
  F_05_01_sw_post_lgl <- B_05_02_score_fn(C_03_02_solve_fn(
    B_01_02b_sw_fn(B_01_02_default_lst)))

  output$test_ui <- renderUI({
    cap_chr <- B_05_03b_caption_fn(F_05_01_post_lgl, F_05_01_sw_post_lgl,
                                   input$stage)
    if (identical(input$stage, "2.5e")) {
      return(B_05_03_table_fn(list(
        "Ours, Estimated"          = F_05_01_post_lgl,
        "Theirs, Estimated"        = F_05_01_sw_post_lgl,
        "Theirs, At Your Settings" = F_02_06_sw_score_r()), cap_chr))
    }
    B_05_03_table_fn(list("Estimated" = F_05_01_post_lgl,
                          "At Your Settings" = F_02_04_score_r()), cap_chr)
  })

  ###### F_05_03: The Composites, Under the Equations Panel's Names ############
  # Note: HTML, so the letters read as the Equations panel prints them: pi_1
  #   to pi_3 for the solver's p1 to p3, and so on.

  F_05_03_display_vec <- c(
    c1 = "c<sub>1</sub>", c2 = "c<sub>2</sub>", c3 = "c<sub>3</sub>",
    i1 = "i<sub>1</sub>", i2 = "i<sub>2</sub>", q1 = "q<sub>1</sub>",
    z1 = "z<sub>1</sub>", k1 = "k<sub>1</sub>", k2 = "k<sub>2</sub>",
    p1 = "&pi;<sub>1</sub>", p2 = "&pi;<sub>2</sub>", p3 = "&pi;<sub>3</sub>",
    w1 = "w<sub>1</sub>", w2 = "w<sub>2</sub>", w3 = "w<sub>3</sub>",
    w4 = "w<sub>4</sub>")

  ###### F_05_04: Every Composite Is Stored ####################################
  # Note: C_01_03 stores z_1, k_1 and k_2 with the rest, so the panel reads
  #   each from d_lst.

  F_05_04_value_fn <- function(d_lst, nm_chr) d_lst[[nm_chr]]

  output$diagnostic_ui <- renderUI({
    sol_lst <- F_02_02_sol_r()
    share_lst <- C_05_02_shares_fn(sol_lst$d_lst)
    res_lst <- if (isTRUE(sol_lst$ok_lgl)) {
      C_05_01_residual_fn(sol_lst, n_draw_int = 60L)
    } else list(system_num = NA_real_, expect_num = NA_real_)

    tagList(
      div(class = "eq-group-title", "Blanchard and Kahn, Counted"),
      tags$ul(
        tags$li(sprintf("%d stable roots, %d unstable, %d forecast errors.",
                        sol_lst$ns_int %||% 0L, sol_lst$nu_int %||% 0L,
                        sol_lst$neta_int %||% 0L)),
        tags$li(sprintf("The forecast errors span %d of the %d unstable
                        directions.", sol_lst$rank_int %||% 0L,
                        sol_lst$nu_int %||% 0L)),
        tags$li(sprintf("Condition number of the stable basis: %.0f.",
                        sol_lst$cond_num %||% NA_real_))),
      div(class = "eq-group-title", "Does the Solution Satisfy the Model?"),
      tags$ul(
        tags$li(sprintf("Canonical system residual: %.2e.",
                        res_lst$system_num)),
        tags$li(sprintf("Rational expectations, period by period: %.2e.",
                        res_lst$expect_num))),
      div(class = "eq-group-title", "The Steady State"),
      tags$ul(
        tags$li(sprintf("Consumption %.4f, investment %.4f, spending %.4f.",
                        share_lst$ccy_num, share_lst$ciy_num,
                        share_lst$gy_num)),
        tags$li(sprintf("They sum to %.10f, which %s.", share_lst$sum_num,
                        if (share_lst$ok_lgl) "is right"
                        else "is wrong: the steady state does not exist")),
        tags$li(sprintf("Rental rate of capital: %.6f per quarter.",
                        sol_lst$d_lst$crk))),
      # composites under the Equations panel's names
      div(class = "eq-group-title", "Composite Coefficients at These Values"),
      tags$ul(lapply(names(F_05_03_display_vec), function(nm_chr)
        tags$li(HTML(sprintf("%s = %.6f", F_05_03_display_vec[[nm_chr]],
                             F_05_04_value_fn(sol_lst$d_lst, nm_chr)))))),
      # stage 2.5e lists Smets and Wouters' composites beside ours
      if (identical(input$stage, "2.5e") && isTRUE(F_02_05_sw_r()$ok_lgl)) {
        tagList(
          div(class = "eq-group-title",
              "Smets and Wouters' Composites at These Values"),
          tags$ul(lapply(names(F_05_03_display_vec), function(nm_chr)
            tags$li(HTML(sprintf("%s = %.6f", F_05_03_display_vec[[nm_chr]],
                                 F_05_04_value_fn(F_02_05_sw_r()$d_lst,
                                                  nm_chr)))))))
      })
  })

  ###### F_06_01: The Model So Far #############################################
  # Note: The toolkit's three panels. This app's stages are "2.5a" to "2.5e"
  #   rather than numbers, so the items in force are picked here (T_06_03
  #   wants numeric stage keys) and handed to the toolkit's own Equations
  #   and In Words builders, T_06_04 and T_06_06.

  output$eq_title <- renderText({
    trimws(T_05_04_stage_name_fn(B_02_01_stage_vec, input$stage))
  })

  F_06_02_items_r <- reactive({
    rank_int  <- B_02_03_rank_fn(input$stage)
    keep_lst  <- Filter(function(it_lst)
      min(vapply(names(it_lst$versions), B_02_03_rank_fn, integer(1))) <=
        rank_int, B_04_01_equation_lst)
    lapply(keep_lst, function(it_lst) {
      key_chr  <- names(it_lst$versions)
      key_int  <- vapply(key_chr, B_02_03_rank_fn, integer(1))
      now_int  <- max(key_int[key_int <= rank_int])
      now_chr  <- key_chr[key_int == now_int]
      # the assembly stage badges every item new; see B_02_03b
      status_chr <-
        if (identical(input$stage, B_02_03b_assembly_chr)) "new" else
        if (now_int != rank_int) "" else
        if (now_int == min(key_int)) "new" else "changed"
      was_chr <- if (identical(status_chr, "changed")) {
        it_lst$versions[[key_chr[key_int == max(key_int[key_int < now_int])]]]
      }
      list(group = it_lst$group, label = it_lst$label,
           tex = it_lst$versions[[now_chr]], status = status_chr,
           was = was_chr, note = it_lst$notes[[now_chr]])
    })
  })

  # stage 2.5e shows "was" on the Equations tab too, through F_00_02
  output$eq_model <- renderUI({
    if (identical(input$stage, "2.5e")) {
      return(tagList(
        F_00_02_bridge_model_fn(F_06_02_items_r(), B_04_02_group_vec),
        tags$p(class = "stat-caption", HTML(B_04_04b_bridge_note_chr))))
    }
    tagList(
      T_06_04_model_fn(F_06_02_items_r(), B_04_02_group_vec,
                       "Nothing here yet.", cols = 1L),
      tags$p(class = "stat-caption", HTML(B_04_04_equation_note_chr)))
  })

  ###### F_06_05: Panel Two, What Each Symbol Means ############################
  # Note: The toolkit's layout, one column per group, and only what the
  #   student has reached. Built here rather than by T_06_05 for the same
  #   reason as above, and because two descriptions carry HTML.

  output$eq_notation <- renderUI({
    rank_int <- B_02_03_rank_fn(input$stage)
    item_lst <- Filter(function(x_lst)
      B_02_03_rank_fn(x_lst$from) <= rank_int, B_04_03_notation_lst)
    col_fn <- function(grp_chr, title_chr) {
      its_lst <- Filter(function(x_lst) x_lst$grp %in% grp_chr, item_lst)
      div(div(class = "eq-group-title", title_chr),
          tags$table(class = "nota-table", lapply(its_lst, function(x_lst) {
            tags$tr(tags$td(T_06_01_mj_fn(x_lst$sym)),
                    tags$td(HTML(paste0(toupper(substr(x_lst$txt, 1L, 1L)),
                                        substring(x_lst$txt, 2L))),
                            if (identical(B_02_03_rank_fn(x_lst$from),
                                          rank_int) && rank_int > 1L) {
                              tags$span(class = "eq-flag eq-new", "new")
                            }))
          })))
    }
    withMathJax(layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(4, 4, 4)),
      col_fn("var", "Variables"),
      col_fn("par", "Parameters"),
      col_fn(c("tgt", "shk"), "Targets, Thresholds and Shocks")))
  })

  ###### F_06_06: Panel Three, Every Equation With Its Explanation #############
  # Note: The toolkit's In Words builder on the same items.

  output$eq_explain <- renderUI({
    T_06_06_explain_fn(F_06_02_items_r(), B_04_02_group_vec)
  })
}

################################################################################
## G: Run ######################################################################
################################################################################
# Note: Launch.

#### G_01: Launch ##############################################################
# Note: Returns the app object.

###### G_01_01: Build App ######################################################
# Note: UI from E, server from F.

shinyApp(E_02_06_page_ui, F_01_00_server_fn)

#--------------------------------- Script End ---------------------------------#
