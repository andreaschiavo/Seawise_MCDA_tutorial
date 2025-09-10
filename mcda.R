#####
## Author: Andrea Schiavo
## Date: 17.02.2025
## Title: MCDA
## email: andrea.schiavo94@gmail.com
##
## functions to perform a multi-cretiria decision analysis (MCDA) based on Multi-attribute utility theory (MAUT).
##
## version: v.0
#####

library(dplyr)



# Utility functions
logistic_2points <- function(data, x_0, x_1, y_0, y_1) {
  b <- (x_0 - x_1) / (log(1 / y_0 - 1) - log(1 / y_1 - 1))
  a <- b * log(1 / y_0 - 1) - x_0
  return(1 / (1 + exp((data + a) / b)))
}

log_1point <- function(data, ref, u_ref, alpha = 1) {
  b <- log(alpha / (alpha - u_ref)) / ref
  return(1 - exp(-b * data))
}

beta_1point <- function(data, max, k = 4) {
  a <- max * (k - 2) + 1
  b <- (1 - max) * (k - 2) + 1
  dbeta(data, a, b) / dbeta(max, a, b)
}

exp_1point <- function(data, ref, u_ref) {
  a <- log(u_ref + 1) / ref
  return(u_ref * (exp(a * data) - 1) / (exp(a * ref) - 1))
}

negexp_1point <- function(data, ref, u_ref) {
  a <- -log(u_ref) / ref
  return(exp(-a * data))
}

utility_function <- function(data, fun, params) {
  if (fun == "logistic2p") {
    return(logistic_2points(data, params[[1]], params[[2]], params[[3]], params[[4]]))
  } else if (fun == "log1p") {
    return(logistic_1point(data, params[[1]], params[[2]])) # nolint
  } else if (fun == "beta1p") {
    return(beta_1point(data, params[[1]], params[[2]]))
  } else if (fun == "exp1p") {
    return(exp_1point(data, params[[1]], params[[2]]))
  } else if (fun == "negexp1p") {
    return(negexp_1point(data, params[[1]], params[[2]]))
  } else {
    stop("Function not implemented")
  }
}

# Main function to create MCDA object
Mcda <- function(time_period,scenarios, criteria, sub_criteria = c("all"), options_file = "data/MCDA_options.json") {
  # read JSON option file
  options <- jsonlite::fromJSON(options_file)
  
  # Build criteria
  if (length(criteria) == 1) {
    if (criteria[1] == 'all') {
      criteria_list <- names(options$criteria)
    }
  } else if (check_criteria(criteria, options$criteria)) {
    criteria_list <- criteria
  } else {
    stop("One or more criteria are not in the list")
  }
  
  # Build Sub-criteria
  if (is.null(sub_criteria)) {
    sub_criteria_list <- list()
    }
  else if (length(sub_criteria)==1 && sub_criteria[1] == "all") {
    sub_criteria_list <- lapply(criteria_list, function(crt) {
    sub_crit <- options$criteria[[crt]]$sub_criteria
    if (is.null(sub_crit)) {
      return(NULL)  
    } else {
      return(names(sub_crit))
    }
  })
    names(sub_criteria_list) <- criteria_list
    
    sub_criteria_list <- lapply(sub_criteria_list, function(x) {
      if (identical(x, character(0))) return(NULL) else return(x)
    })
    
  } else if (check_subcriteria(sub_criteria, options$criteria)) {
    sub_criteria_list <- lapply(criteria_list, function(crt) {
      if (!is.null(options$criteria[[crt]]$sub_criteria)) {
        return(sub_criteria[[crt]])
      } else {
        return(NULL)
      }
    })
    names(sub_criteria_list) <- criteria_list

  } else {
    stop("Sub-criteria not in the list")
  }

  # Build MCDA object
  mcda <- list(
    options = options,
    criteria = criteria_list,
    sub_criteria = sub_criteria_list,
    data = setNames(vector("list", length(criteria_list)), criteria_list),
    weights = list(),
    alternatives = list(),
    parents = list()
  )
  
  # build alternatives and weights
  mcda$alternatives <- build_alternatives(mcda, time_period, scenarios)
  mcda$weights <- build_weights(mcda)
  
  return(mcda)
}

# Function to create alternatives
build_alternatives <- function(mcda, time_period, scenarios) {
  alternatives <- list()
  all_combinations <- expand.grid(time_period = time_period, scenario = scenarios)
  
  for (i in seq_len(nrow(all_combinations))) { #1:nrow(all_combinations)
    alt <- all_combinations[i, ]
    label <- paste(alt$time_period, alt$scenario, sep = "_")
    alternatives[[label]] <- list(
      time_period = alt$time_period,
      scenario = alt$scenario,
      utilities = setNames(vector("list", length(mcda$criteria)), mcda$criteria),
      total_utility = NULL
    )
  }
  return(alternatives)
}

# Function to create weights
build_weights <- function(mcda) {
  weights <- list()
  for (c in mcda$criteria) {
    weights[[c]] <- mcda$options$criteria[[c]]$weight
    if (!is.null(mcda$options$criteria[[c]]$sub_criteria)) {
      for (s in mcda$sub_criteria[[c]]) {
        weights[[s]] <- mcda$options$criteria[[c]]$sub_criteria[[s]]$weight
      }
    }
  }
  return(weights)
}

# check criteria
check_criteria <- function(criteria, criteria_options) {
  all(sapply(criteria, function(c) c %in% names(criteria_options)))
}

# check sub-criteria
check_subcriteria <- function(subcriteria, criteria_options) {
  for (c in names(subcriteria)) {
    for (s in subcriteria[[c]]) {
      if (!(s %in% names(criteria_options[[c]]$sub_criteria))) {

        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# Load data
load_data <- function(mcda) {
  for (c in mcda$criteria) {
    mcda$data[[c]] <- read.csv(mcda$options$data[[c]], row.names = 1)
  }
  return(mcda)
}

# function to calculate sub-utilities
calculate_sub_utilities <- function(mcda, data, criteria) {
  sub_criteria_list <- names(mcda$options$criteria[[criteria]]$sub_criteria)
  utilities <- list()
  
  for (sub_criteria in mcda$sub_criteria[[criteria]]) {
    data_sel <- data[sub_criteria_list == sub_criteria, ]
    utilities[[sub_criteria]] <- utility_function(
      data_sel[[criteria]],
      mcda$options$criteria[[criteria]]$sub_criteria[[sub_criteria]]$function_,
      mcda$options$criteria[[criteria]]$sub_criteria[[sub_criteria]]$ref_points
    )
  }
  return(utilities)
}

# function to calculate utilities
calculate_utilities <- function(mcda, ... ) {
  
  args <- list(...)
  
  for (alt in names(mcda$alternatives)) {
    for (criteria in mcda$criteria) {
      data <- mcda$data[[criteria]]
      # Basic filter
      alt_data <- data %>% filter(time_period == mcda$alternatives[[alt]]$time_period,
                                  scenario == mcda$alternatives[[alt]]$scenario)
      # Advanced filter using ... arguments
      if (length(args) > 0) {
        for (arg in names(args)) {
          alt_data <- alt_data %>% filter(.data[[arg]] == args[[arg]])
        }
      }
      if (is.null(mcda$sub_criteria[[criteria]])) {
        mcda$alternatives[[alt]]$utilities[[criteria]] <- utility_function(
          alt_data[[criteria]],
          mcda$options$criteria[[criteria]]$function_,
          mcda$options$criteria[[criteria]]$ref_points
        )
      } else {
        mcda$alternatives[[alt]]$utilities[[criteria]] <- calculate_sub_utilities(mcda, alt_data, criteria)
      }
    }
  }
  return(mcda)
}

# function to aggregate utilities
aggregate_utilities <- function(mcda, criteria, label, method = "sum", set_weight = 1.0, sub_criteria = FALSE, weighted = FALSE) {
  
  # this function aggregate the utilities in different ways:
  #   - if only 1 criterion is passed, and it has sub criteria to aggregate along, it is possible to chose the method and perform a weighted operation (weighted sum, or weighted geometric mean)
  #   - if two or more criteria are passed, if they share the same sub criteria, they can be aggregated (setting sub_criteria = TRUE) in one new criterion having the same number and label of sub criteria of the original ones passed
  #       in this case it is possible to chose one of the two methods (sum or gmean of the criteria), however weights are NOT implemented
  #   - if two or more criteria are passed and they do NOT have sub criteria, then it is possible to perform the same operations of the first case
  #   - it's also possible to create new macro criteria by passing the related critierion with sub_criteria set to FALSE. the function will return a new criterion, named after the label passes, witht the same utility of the original criterion
  # Default weight for the new label is set to 1
  # Two methods are implemented: "sum" and "gmean"
  # by setting weighted = TRUE the function performs a weighted arithmetic mean ("sum") or a weighted geometric mean ("gmean"). the result is normalized between 0 and 1
  
  if (!is.character(criteria)) {
    stop("The 'criteria' parameter must be a character vector.")
  }

  for (alt in names(mcda$alternatives)) {
    utilities <- mcda$alternatives[[alt]]$utilities
 
    # Initialize aggregation based on chosen method
    aggregate <- ifelse(method == "gmean", 1, ifelse(method == "sum", 0, stop("Method not implemented")))
    
    if (sub_criteria) {
      # aggregation of sub-criteria
      sub_utilities <- utilities[[criteria[1]]]
      weights <- if (weighted) sapply(names(sub_utilities), function(s) mcda$weights[[s]]) else rep(1, length(sub_utilities))
      weights <- weights / sum(weights)

      new_sub_utilities <- list()

      for (sub in names(sub_utilities)) {
        
        if (length(criteria) == 1) {
          # sub-criteria aggregation for 1 criterion
          value <- sub_utilities[[sub]]
          
          aggregate <- if (method == "gmean") {
            aggregate * (value * weights[[sub]])
          } else if (method == 'sum') {
            aggregate + (value * weights[[sub]])
          }
        } else {
          # sub-criteria aggregation for 2 criteria sharing the same sub-criteria
          value <- if (method == "gmean") {
            prod(sapply(criteria, function(c) utilities[[c]][[sub]]))
          } else if (method == 'sum'){
            sum(sapply(criteria, function(c) utilities[[c]][[sub]]))
          }

          new_sub_utilities[[sub]] <- value
        }
      }

      # add the aggregation result to the utilities
      if (length(criteria) == 1) {
        utilities[[label]] <- aggregate
      } else {
        utilities[[label]] <- new_sub_utilities
      }
    } else {
      # criteria aggregation
      weights <- if (weighted) sapply(criteria, function(c) mcda$weights[[c]]) else rep(1, length(criteria))
      weights <- weights / sum(weights)

      for (c in criteria) {
        value <- utilities[[c]]
        aggregate <- if (method == "gmean") {
          aggregate * (value * weights[criteria == c])
        } else if (method == "sum"){
          aggregate + (value * weights[criteria == c])
        }
      }

      # add the aggregation result to the utilities
      utilities[[label]] <- as.numeric(aggregate)
    }

    # update alternative utilities
    mcda$alternatives[[alt]]$utilities <- utilities
  }

  # add new weight
  if (!label %in% names(mcda$weights)) {
    mcda$weights[[label]] <- set_weight
  }
  mcda$parents[[label]] <- criteria
  return(mcda)
}

# Function to execute maut analysis
run_maut <- function(mcda, criteria) {
  
  # the MAUT method performs a weighted mean of the criteria's utilities that are relevant to the stakeholders. the resultant total utility provide a summary estimate of the alternative value
  
  if (length(criteria) == 1) {
    if (criteria == 'all') {
      criteria_list <- mcda$criteria
  } else {
    criteria_list <- criteria
  }}
    else {
      criteria_list <- criteria
    }
  w <- unlist(lapply(criteria_list, function(c) mcda$weights[[c]]))
  w <- w / sum(w)
  
  names(mcda$alternatives)

  for (alt in names(mcda$alternatives)) {
    U_tot <- 0
    for (i in seq_along(criteria_list)) {
      if (criteria_list[i] %in% names(mcda$alternatives[[alt]]$utilities)){
        U_tot <- U_tot + mcda$alternatives[[alt]]$utilities[[criteria_list[i]]] * w[i]
    }}
    mcda$alternatives[[alt]]$total_utility <- U_tot
  }
  return(mcda)
}

generate_df <- function(mcda, criteria, sub_criteria = NULL, facet_cols = TRUE) {
  
  # funciton to exctract the final scores form the mcda object.
  # the funciton return a dataframe with scenarios on the rows and
  # the selected criteria on the columns
  # it is possible to extract sub_criteria as well
  # if facet_cols argument is set to TRUE (default),
  # 3 additional columns will be extracted from the MCDA object with
  # time period (time), climate change scenario (cc) and management scenario (management)
  
  df <- data.frame()
  U_tot <- c()
  scenario <- c()
  
  for (alt in names(mcda$alternatives)) {
    alternative <- mcda$alternatives[[alt]]
    
    U_tot <- c(U_tot, alternative$total_utility)
    utilities <- alternative$utilities
    
    values <- sapply(criteria, function(c) utilities[[c]])
    keys <- criteria
    
    if (!is.null(sub_criteria)) {
      sub_values <- unlist(lapply(sub_criteria, function(s) unlist(utilities[[s]])))
      sub_labels <- unlist(lapply(sub_criteria, function(s) paste0(substr(s, 1, 3), "_", names(utilities[[s]]))))
      
      values <- c(values, sub_values)
      keys <- c(keys, sub_labels)
    }
    
    row_ <- as.data.frame(t(values))
    colnames(row_) <- keys
    df <- rbind(df, row_)
    scenario <- c(scenario, alt)
  }
  
  df$scenario <- scenario
  df$total_utility <- U_tot
  df <- df %>% select(scenario, everything())

  if (facet_cols){

    time <- list()
    cc <- list()
    mgmt <- list()

    for (alt in mcda$alternatives) {
      time <- c(time, as.character(alt$time_period))
      sce <- as.character(alt$scenario[1])
      manage <- strsplit(sce, "_")[[1]][1]
      climate <- strsplit(sce, "_")[[1]][2]
      cc <- c(cc, climate)
      mgmt <- c(mgmt, manage)
    }
    df$time_period <- as.character(time)
    df$cc <- as.character(cc)
    df$management <- as.character(mgmt)
  }

  return(df)
}

top_scenarios <- function(dataframe, criterion, num = 5, ascending = FALSE) {
  
  dataframe %>%
    arrange(if (ascending) !!sym(criterion) else desc(!!sym(criterion))) %>%
    head(num)
}

set_weights <- function(mcda, criteria, new_weights) { 
  i<-1
  for (name in criteria) { 
    mcda$weights[[name]] <- new_weights[i]
    i <- i+1
  }
  return(mcda)
} 

#### VISUALIZATON #####

plot_matrix <- function(df, criteria = "all", scenarios = "all", color = "Blues", style = "minimal") {

  # Filter on criteria (cols) & scenarios (rows)
  if (length(criteria) > 1) {
    df <- df[, c(criteria, "scenario"), drop = FALSE]
  }
  if (length(scenarios) > 1) {
    df <- df[df$scenario %in% scenarios, , drop = FALSE]
  }
  
  df_long <- reshape2::melt(df, id.vars = "scenario", variable.name = "Criteria", value.name = "Value")
  df_long$scenario <- factor(df_long$scenario, levels = rev(unique(df_long$scenario)))
  df_long$Value <- as.numeric(df_long$Value)

  if (style == "dark") {
    p <- ggplot2::theme_set(ggplot2::theme_dark())
  } else if (style == "minimal") {
    p <- ggplot2::theme_set(ggplot2::theme_minimal())
  } else {
    p <- ggplot2::theme_set(ggplot2::theme_gray())
  }

## heatmap
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Criteria, y = scenario, fill = Value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colors = scales::brewer_pal(palette = color)(9), limits = c(0, 1)) +
    ggplot2::labs(x = "Criteria", y = "Scenarios", fill = "Utility\nvalue") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

 print(p)
}

## Parallel graph
plot_parallel_graph <- function(df, criteria, highlight = "top", top_bottom_criterion = "total_utility", save = FALSE) {
  
  ascending <- ifelse(highlight == "top", FALSE, TRUE)

  top_df <- top_scenarios(df, criterion = top_bottom_criterion, num = 5, ascending = ascending)

  colors <- RColorBrewer::brewer.pal(5, "Set1")

  df_long <- reshape2::melt(df, id.vars = "scenario", measure.vars = criteria, variable.name = "Criteria", value.name = "Utility")
  
  top_df_long <- reshape2::melt(top_df, id.vars = "scenario", measure.vars = criteria, variable.name = "Criteria", value.name = "Utility")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Criteria, y = Utility, group = scenario)) +
    ggplot2::geom_line(color = "grey", alpha = 0.5, linewidth = 0.5) +
    ggplot2::geom_line(data = top_df_long, ggplot2::aes(color = scenario), linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = seq_along(criteria), color = "black", linetype = "solid", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), color = "#a7a2a2", linetype = "dashed", linewidth = 0.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(y = "Utility", color = "Highlighted Scenarios") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  print(p)
}  

plot_histogram <- function(df, criteria, weights, scenarios = NULL, subcriteria = FALSE, parents = NULL, color = "Set3", facet = TRUE, save = FALSE) {

  # Filtrare i dati per gli scenari selezionati
  if (!is.null(scenarios)) {
    data_filtered <- df[df$scenario %in% scenarios, ]
  } else {
    data_filtered <- df
  }

  # Normalizzare i pesi dei macro-criteri
  w <- sapply(criteria, function(c) weights[[c]])
  w <- w / sum(w)
  y <- data_filtered[, criteria, drop = FALSE]
  y <- sweep(y, 2, w, "*") # Calcolare le utility dei macro-criteri
  y <- cbind(scenario = df$scenario, time = df$time_period, cc = df$cc, mgmt = df$management, y)

  # Palette colori per i macro-criteri
  macro_colors <- if (length(color) == 1) {
    setNames(RColorBrewer::brewer.pal(n = length(criteria), name = color), criteria)
  } else {
    setNames(color, criteria)
  }

  if (subcriteria) {
    # Estrarre i sottocriteri e calcolare la loro utility
    sub_criteria <- unlist(lapply(criteria, function(c) parents[[c]]))
    y_sub <- data_filtered[, sub_criteria, drop = FALSE]

    for (macro in criteria) {
      sub_set <- parents[[macro]]
      y_sub[, sub_set] <- sweep(y_sub[, sub_set, drop = FALSE], 1, rowSums(y_sub[, sub_set, drop = FALSE]), FUN = "/") * y[, macro]
    }

    df_long <- reshape2::melt(
      cbind(scenario = df$scenario, time = df$time_period, cc = df$cc, mgmt = df$management, y_sub),
      id.vars = c("scenario", "time", "cc", "mgmt"))

    # Creare le sfumature di colore per i sottocriteri
    sub_colors <- unlist(lapply(1:length(criteria), function(i) {
      shades <- colorRampPalette(c("white", macro_colors[i]))(length(parents[[criteria[i]]]) + 1)[-1]
      names(shades) <- parents[[criteria[i]]]
      return(shades)
    }))

    df_long$MacroCriterion <- unlist(lapply(df_long$variable, function(x) {
      macro <- names(parents)[sapply(parents, function(subs) x %in% subs)]
      return(ifelse(length(macro) > 0, macro, NA))
    }))

    all_colors <- c(sub_colors)

    if (facet) {
      p <- ggplot2::ggplot(df_long, ggplot2::aes(x = mgmt, y = value, fill = variable))
    } else {
      p <- ggplot2::ggplot(df_long, ggplot2::aes(x = scenario, y = value, fill = variable))
    }

    p <- p +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::scale_fill_manual(values = all_colors, name = "Subcriteria") +
      ggplot2::labs(x = "Scenarios", y = "Utility") +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            legend.position = "right",
            legend.title = ggplot2::element_text(face = "bold"),
            panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.5),
            panel.grid.minor = ggplot2::element_line(color = "gray90", linewidth = 0.3),
            panel.background = ggplot2::element_blank()) +
      ggplot2::ylim(0,1)
    if (facet){
      p <- p + ggplot2::facet_grid(rows = cc, cols = time)
    }

  } else {
    df_long <- reshape2::melt(
      cbind(scenario = df$scenario, time = df$time_period, cc = df$cc, mgmt = df$management, y),
      id.vars = c("scenario", "time", "cc", "mgmt"))

    if (facet) {
      p <- ggplot2::ggplot(df_long, ggplot2::aes(x = mgmt, y = value, fill = variable))
    } else {
      p <- ggplot2::ggplot(df_long, ggplot2::aes(x = scenario, y = value, fill = variable))
    }
    p <- p +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::scale_fill_manual(values = macro_colors, name = "Criteria") +
      ggplot2::labs(x = "Scenarios", y = "Utility") +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            legend.position = "right",
            legend.title = ggplot2::element_text(face = "bold"),
            panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.5),
            panel.grid.minor = ggplot2::element_line(color = "gray90", linewidth = 0.3),
            panel.background = ggplot2::element_blank()) +
      ggplot2::ylim(0,1)

      if (facet) {
        p <- p + ggplot2::facet_grid(rows = vars(cc), cols = vars(time))
      }
  }
  
  return(p)
}
