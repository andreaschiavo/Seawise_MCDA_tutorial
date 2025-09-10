library(jsonlite)

options <- list(
  criteria = list(
    ssb = list(
      weight = 1,
      sub_criteria = list(
        HKE = list(weight = 1, function_ = "logistic2p", ref_points = c(800., 2000., 0.1, 0.99)),
        MUT = list(weight = 1, function_ = "logistic2p", ref_points = c(500., 2000., 0.1, 0.99)),
        SOL = list(weight = 1, function_ = "logistic2p", ref_points = c(200., 2000., 0.1, 0.99))
      ),
      function_ = NULL,
      ref_points = NULL
    ),
    f = list(
      weight = 1,
      sub_criteria = list(
        HKE = list(weight = 1, function_ = "logistic2p", ref_points = c(0.206, 0.412, 0.99, 0.1)),
        MUT = list(weight = 1, function_ = "logistic2p", ref_points = c(0.37, 0.74, 0.99, 0.1)),
        SOL = list(weight = 1, function_ = "logistic2p", ref_points = c(0.937, 1.874, 0.99, 0.1))
      ),
      function_ = NULL,
      ref_points = NULL
    ),
    employment = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(75, 150, 0.1, 0.5)),
    wage = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(750, 1500, 0.1, 0.5)),
    gva = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(9000, 18000, 0.1, 0.5)),
    rsl = list(weight = 1, sub_criteria = NULL, function_ = "beta1p", ref_points = c(0.5, 4)),
    co2 = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(1, 13, 0.99, 0.5)),
    mml = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(30, 45.5, 0.1, 0.5)),
    avpb = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(275, 550, 0.1, 0.5)),
    rbs = list(weight = 1, sub_criteria = NULL, function_ = "logistic2p", ref_points = c(0.35, 0.7, 0.1, 0.5))
  ),
  data = list(
    ssb = "data/MCDA_stock_dummy.csv",
    f = "data/MCDA_stock_dummy.csv",
    employment = "data/MCDA_dummy.csv",
    wage = "data/MCDA_dummy.csv",
    gva = "data/MCDA_dummy.csv",
    rsl = "data/MCDA_dummy.csv",
    co2 = "data/MCDA_dummy.csv",
    mml = "data/MCDA_dummy.csv",
    avpb = "data/MCDA_dummy.csv",
    rbs = "data/MCDA_dummy.csv"
  )
)



write_json(options, "data/MCDA_options_dummy.json", pretty = TRUE)


