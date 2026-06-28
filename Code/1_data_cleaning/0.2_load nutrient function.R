# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 0.2_load nutrient function.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Defines the HPFS nutrient-loading function used to assemble biennial or questionnaire-cycle nutrition variables, including calories, blank FFQ items, alcohol, dietary cholesterol, and saturated fat.
# ==============================================================================

#' Load HPFS Nutrient Data (Updated)
#'
#' This function loads in data from nutrients files.
#'
#' Included in the current version of this function are data on total calories
#' consumed, number of blank items on the food frequency questionnaire (FFQ),
#' and alcohol, dietary cholesterol, and saturated fat consumption.
#'
#' The data here come from derived nutrient files. For example, the alcohol data
#' here reflects total alcohol consumption, in grams/day.
#' Across questionnaire cycles in HPFS, individual types of alcoholic
#' beverages are queried (e.g., beer, wine). Nutritional information in units of
#' grams of alcohol consumed per day has been derived using individual alcoholic
#' beverages. This function reads in just the total alcohol consumed in
#' grams.
#'
#' Individual food and alcohol items queried on the FFQ can be read in from the
#' main questionnaires.
#'
#' Note that nutritional data are queried on HPFS questionnaires every 4 years
#' when the FFQ is asked, rather than biennially like many other variables.
#' Therefore, these data are available in 1986, 1990, 1994, 1998, 2002, 2006,
#' 2010, 2014, and 2018.
#'
#'
#' @return
#' Tibble of >400,000 observations times 7 variables, with one line per person
#' per questionnaire cycle where data on nutrient data are available.
#' The nutrient variables are missing if no (long) questionnaire was returned
#' during that cycle.
#'  * \code{id} HPFS ID, 6 digits.
#'  * \code{cycle} Questionnaire cycle.
#'  * \code{calor} Total calories, kcal/day, numeric.
#'  * \code{nblnk} Number of blank items, numeric.
#'  * \code{sat} Saturated fat, g/day, numeric.
#'  * \code{chol} Cholesterol, mg/day, numeric.
#'  * \code{alco} Alcohol consumed, g/day, numeric.
#'
#' @export
#'
#' @examples
#' nutrients <- get_hpfs_nutrients()
#' nutrients
#'
#' # Number of participants with data per cycle and variable
#' nutrients |>
#'   dplyr::group_by(cycle) |>
#'   dplyr::summarize(
#'     dplyr::across(
#'       .cols = c(-id),
#'       .fns = \(x) sum(!is.na(x))
#'     )
#'   )
get_hpfs_nutrients <- function() {
  # Call in nutrition files where nutrient data is derived.
  # For reference, these are nutrition files, rather than the main questionnaire
  # files. These are available every 4 years and end in ".nts".
  # Use readr::read_fwf to read in, using columns defined in the corresponding
  # input file (/n/hpnh/ReadMacros/HPFS/input/). To find the file version called
  # by the SAS macro, go to /n/hpnh/ReadMacros/HPFS/hpstools/sasautos and
  # find the file name at the top of the program.
  
  h86_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h86_.040925.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor86 = c(7, 10),
      sat86 = c(108, 111),
      chol86 = c(130, 133),
      nblnk86 = c(350, 351),
      alco86 = c(1051, 1055)
    ),
    show_col_types = FALSE
  )
  
  h90_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h90_.040925.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor90 = c(7, 10),
      nblnk90 = c(639, 640),
      sat90 = c(223, 226),
      chol90 = c(374, 377),
      alco90 = c(474, 478)
    ),
    show_col_types = FALSE
  )
  
  h94_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h94_.040825.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor94 = c(7, 10),
      nblnk94 = c(813, 814),
      sat94 = c(194, 198),
      chol94 = c(357, 362),
      alco94 = c(472, 477)
    ),
    show_col_types = FALSE
  )
  
  h98_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h98_.032426.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor98 = c(7, 10),
      nblnk98 = c(786, 787),
      sat98 = c(261, 265),
      chol98 = c(490, 495),
      alco98 = c(595, 600)
    ),
    show_col_types = FALSE
  )
  
  h02_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h02_.031926.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor02 = c(7, 10),
      nblnk02 = c(1125, 1126),
      sat02 = c(345, 350),
      chol02 = c(651, 656),
      alco02 = c(799, 804)
    ),
    show_col_types = FALSE
  )
  
  h06_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h06_.031826.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor06 = c(7, 10),
      nblnk06 = c(11, 12),
      sat06 = c(412, 416),
      chol06 = c(628, 633),
      alco06 = c(843, 848)
    ),
    show_col_types = FALSE
  )
  
  h10_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h10_.032426.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor10 = c(7, 10),
      nblnk10 = c(1527, 1529),
      sat10 = c(555, 560),
      chol10 = c(766, 771),
      alco10 = c(990, 995)
    ),
    show_col_types = FALSE
  )
  
  h14_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h14_.031926.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor14 = c(7, 10),
      nblnk14 = c(1647, 1648),
      sat14 = c(566, 571),
      chol14 = c(877, 882),
      alco14 = c(1105, 1110)
    ),
    show_col_types = FALSE
  )
  
  h18_nts <- readr::read_fwf(
    file = "/n/hpnh/HPFS/Qdata/hp_dat_dtf/h18_.031926.nts",
    col_positions = readr::fwf_cols(
      id = c(1, 6),
      calor18 = c(7, 10),
      nblnk18 = c(1348, 1349),
      sat18 = c(495, 500),
      chol18 = c(598, 603),
      alco18 = c(799, 804)
    ),
    show_col_types = FALSE
  )
  
  # Merge dataframes
  list(
    h86_nts, h90_nts, h94_nts, h98_nts, h02_nts, h06_nts, h10_nts, h14_nts,
    h18_nts
  ) |>
    purrr::reduce(.f = \(x, y) dplyr::left_join(x, y, by = "id")) |>
    # Pivot from wide to long format
    tidyr::pivot_longer(
      cols = !"id",
      names_to = c(".value", "cycle"),
      names_pattern = "(...*)([[:digit:]]{2})$"
    ) |>
    dplyr::mutate(cycle = convert_cycle_to_factor(.data$cycle)) |>
    labelled::set_variable_labels(
      id = "Participant ID",
      cycle = "Questionnaire cycle",
      calor = "Total calories, kcal/d",
      nblnk = "Number of blank items",
      sat = "Saturated fat, g/d",
      chol = "Cholesterol, mg/d",
      alco = "Alcohol, g/d"
    )
}