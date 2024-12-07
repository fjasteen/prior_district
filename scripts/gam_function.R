apply_gam
function(df,
         y_var,
         eval_years,
         year = "year",
         taxonKey = "taxonKey",
         type_indicator = "observations",
         baseline_var = NULL,
         p_max = 0.1,
         taxon_key = NULL,
         name = NULL,
         df_title = NULL,
         x_label = "jaar",
         y_label = "aantal km²-hokken",
         saveplot = FALSE,
         dir_name = NULL,
         verbose = FALSE) {
  if (is.numeric(taxon_key)) {
    taxon_key <- as.character(taxon_key)
  }
  
  # Check right type of inputs
  assert_that(is.data.frame(df),
              msg = paste(
                paste(as.character(df), collapse = ","),
                "is not a data frame.",
                "Check value of argument df."
              )
  )
  map2(
    list(y_var, year, taxonKey, type_indicator, x_label, y_label),
    c("y_var", "year", "taxonKey", "type_indicator", "x_label", "y_label"),
    function(x, y) {
      # Check right type of inputs
      assert_that(is.character(x),
                  msg = paste0(
                    paste(as.character(x), collapse = ","),
                    " is not a character vector.",
                    " Check value of argument ", y, "."
                  )
      )
      # Check y_var, year, taxonKey, type_indicator have length 1
      assert_that(length(x) == 1,
                  msg = paste0(
                    "Multiple values for argument ",
                    paste0(y, collapse = ","),
                    " provided."
                  )
      )
    }
  )
  assert_that(is.numeric(eval_years),
              msg = paste(
                paste(as.character(eval_years), collapse = ","),
                "is not a numeric or integer vector.",
                "Check value of argument eval_years."
              )
  )
  
  map2(
    list(baseline_var, taxon_key, name, df_title, dir_name),
    c("baseline_var", "taxon_key", "name", "df_title", "dir_name"),
    function(x, y) {
      # check argument type
      assert_that(is.null(x) | is.character(x),
                  msg = paste0(
                    paste(as.character(x), collapse = ","),
                    " is not a character vector.",
                    " Check value of argument ", y, "."
                  )
      )
      # check length
      assert_that(length(x) < 2,
                  msg = paste(
                    "Multiple values for argument",
                    y, "provided."
                  )
      )
    }
  )
  
  map2(
    list(saveplot, verbose),
    c("saveplot", "verbose"),
    function(x, y) {
      assert_that(is.logical(x),
                  msg = paste(
                    paste(as.character(x), collapse = ","),
                    "is not a logical vector.",
                    "Check value of argument saveplot.",
                    "Did you maybe use quotation marks?"
                  )
      )
      assert_that(length(x) == 1,
                  msg = paste("Multiple values for argument", y, "provided.")
      )
    }
  )
  
  map2(
    list(y_var, year, taxonKey),
    c("y_var", "year", "taxonKey"),
    function(x, y) {
      # Check y_var, year, taxonKey are present in df
      assert_that(x %in% names(df),
                  msg = paste0(
                    "The column ", x,
                    " is not present in df. Check value of",
                    " argument ", y, "."
                  )
      )
    }
  )
  
  if (!is.null(baseline_var)) {
    # Check baseline_var is present in df
    assert_that(
      baseline_var %in% names(df),
      msg = paste0(
        "The column ", baseline_var,
        " is not present in df. Check value of argument baseline_var."
      )
    )
    method_em <- "correct_baseline"
  } else {
    method_em <- "basic"
  }
  
  if (isFALSE(saveplot)) {
    if (!is.null(dir_name)) {
      warning(paste(
        "saveplot is FALSE: plots are not saved.",
        "Argument dir_name ignored."
      ))
    }
  } else {
    if (!is.null(dir_name)) {
      dir.create(dir_name, showWarnings = FALSE)
    } else {
      # current directory
      dir_name <- "./"
    }
  }
  
  year <- vars_pull(names(df), !!enquo(year))
  taxonKey <- vars_pull(names(df), !!enquo(taxonKey))
  
  # Check eval_year is present in column year
  assert_that(all(eval_years %in% df[[year]]),
              msg = paste(
                "One or more evaluation years",
                "not present in df.",
                "Check value of argument eval_years."
              )
  )
  
  assert_that(is.numeric(p_max) && p_max >= 0 && p_max <= 1,
              msg = paste(
                "p_max is a p-value: it has to be a",
                "number between 0 and 1."
              )
  )
  
  # Check type_indicator is one of the two allowed values
  assert_that(type_indicator %in% c("aantal km²-hokken", "observations"),
              msg = paste(
                "Invalid type_indicator.",
                "type_indicator has to be one of:",
                "aantal km²-hokken, observations."
              )
  )
  
  if (verbose == TRUE) {
    print(paste0("Analyzing: ", name, "(", taxon_key, ")"))
  }
  
  if (nrow(df) > 0) {
    # Maximum minimum time series (year)
    fyear <- min(df[[year]], na.rm = TRUE) # first year
    lyear <- max(df[[year]], na.rm = TRUE) # last year
    
    # Define model to use for GAM
    maxk <- max(round((lyear - fyear) / 10, digits = 0), 5) # max number of knots
  }
  if (method_em == "correct_baseline") {
    fm <- paste0(
      y_var,
      " ~ s(",
      year,
      ", k = maxk, m = 3, bs = \"tp\") + s(",
      baseline_var,
      ")"
    )
    fm <- formula(fm)
  } else {
    method_em <- "basic"
    fm <- paste0(
      y_var,
      " ~ s(",
      year,
      ", k = maxk, m = 3, bs = \"tp\")"
    )
    fm <- formula(fm)
  }
  
  # Initialization
  output_model <- as_tibble(df)
  output_model <-
    output_model %>%
    mutate(
      fit = NA_real_,
      ucl = NA_real_,
      lcl = NA_real_,
      em1 = NA_real_,
      em2 = NA_real_,
      em = NA_real_,
      em_status = NA_real_,
      growth = NA_real_,
      method = method_em
    )
  model <- deriv1 <- deriv2 <- plot_gam <- summary_pv <- p_ok <- NULL
  emerging_status_output <-
    output_model %>%
    filter(!!sym(year) %in% eval_years) %>%
    select(
      !!sym(taxonKey),
      .data$year,
      .data$em_status,
      .data$growth,
      .data$method
    )
  
  if (nrow(df) > 3 & sum(df[[y_var]][2:nrow(df)]) != 0) {
    result <- tryCatch(expr = {
      model <- gam(
        formula = fm,
        family = nb(),
        data = df,
        method = "REML"
      )
      # Check that p-value of at least one smoother < 0.1
      summary_pv <- summary.gam(model)$s.pv
      p_ok <- ifelse(any(summary_pv < p_max), TRUE, FALSE)
    }, error = function(e) e, warning = function(w) w)
    
    if (class(result)[1] %in% c("simpleWarning", "simpleError")) {
      if (verbose) {
        warning(paste0(
          "GAM (",
          method_em,
          ") cannot be performed or cannot converge.\n"
        ))
      }
    } else {
      if (isFALSE(p_ok)) {
        if (verbose) {
          warning(paste0(
            "GAM output cannot be used: ",
            "p-values of all GAM smoothers are above ",
            p_max, ".\n"
          ))
        }
      } else {
        output_model <- df
        # Add method
        output_model <-
          output_model %>%
          mutate(method = method_em)
        # Predict to new data (5 values per year)
        temp <- predict(
          object = model,
          newdata = output_model,
          type = "iterms",
          interval = "prediction",
          se.fit = TRUE
        )
        
        # Calculate confidence intervals & backtransform to real scale
        intercept <- unname(model$coefficients[1])
        output_model$fit <- model$family$linkinv(temp$fit[, 1] + intercept)
        output_model$ucl <- model$family$linkinv(temp$fit[, 1] + intercept + temp$se.fit[, 1] * 1.96)
        output_model$lcl <- model$family$linkinv(temp$fit[, 1] + intercept - temp$se.fit[, 1] * 1.96)
        
        # Check that fit ucl and lcl are all above zero
        output_model <-
          output_model %>%
          mutate(
            fit = ifelse(.data$fit < 0, 0, .data$fit),
            ucl = ifelse(.data$ucl < 0, 0, .data$ucl),
            lcl = ifelse(.data$lcl < 0, 0, .data$lcl)
          )
        
        # Calculate first and second derivative + conf. interval
        deriv1 <- derivatives(model,
                              type = "central", order = 1, level = 0.8,
                              n = nrow(output_model), eps = 1e-4
        )
        deriv2 <- derivatives(model,
                              type = "central", order = 2, level = 0.8,
                              n = nrow(output_model), eps = 1e-4
        )
        
        # Emerging status based on first and second derivative
        em1 <-
          deriv1 %>%
          as_tibble() %>%
          filter(.data$var == year) %>%
          mutate(em1 = case_when(
            .data$lower < 0 & .data$upper <= 0 ~ -1,
            .data$lower < 0 & .data$upper > 0 ~ 0,
            .data$lower >= 0 & .data$upper > 0 ~ 1
          )) %>%
          select(!!sym(year) := .data$data, .data$em1) %>%
          mutate(!!sym(year) := round(!!sym(year)))
        
        em2 <- deriv2 %>%
          as_tibble() %>%
          filter(.data$var == year) %>%
          mutate(em2 = case_when(
            .data$lower < 0 & .data$upper <= 0 ~ -1,
            .data$lower < 0 & .data$upper > 0 ~ 0,
            .data$lower >= 0 & .data$upper > 0 ~ 1
          )) %>%
          select(!!sym(year) := .data$data, .data$em2) %>%
          mutate(!!sym(year) := round(!!sym(year)))
        
        em_level_gam <- full_join(em1, em2, by = year) %>%
          mutate(em = case_when(
            .data$em1 == 1 & .data$em2 == 1 ~ 4,
            .data$em1 == 1 & .data$em2 == 0 ~ 3,
            .data$em1 == 1 & .data$em2 == -1 ~ 2,
            .data$em1 == 0 & .data$em2 == 1 ~ 1,
            .data$em1 == 0 & .data$em2 == 0 ~ 0,
            .data$em1 == 0 & .data$em2 == -1 ~ -1,
            .data$em1 == -1 & .data$em2 == 1 ~ -2,
            .data$em1 == -1 & .data$em2 == 0 ~ -3,
            .data$em1 == -1 & .data$em2 == -1 ~ -4
          ))
        
        # Emerging status
        em_levels <-
          em_level_gam %>%
          mutate(em_status = case_when(
            .data$em < 0 ~ 0, # not emerging
            .data$em == 0 ~ 1, # unclear
            .data$em < 3 ~ 2, # potentially emerging
            .data$em >= 3 ~ 3 # emerging
          ))
        
        output_model <- left_join(output_model, em_levels, by = year)
        
        # Lower value of first dedrivative (minimal guaranted growth) if positive
        lower_deriv1 <-
          deriv1 %>%
          filter(.data$var == year) %>%
          rename(!!sym(year) := .data$data) %>%
          mutate(!!sym(year) := round(!!sym(year), digits = 0)) %>%
          mutate(growth = model$family$linkinv(.data$lower)) %>%
          select(!!sym(year), .data$growth)
        
        # Add lower value of first derivative
        output_model <- left_join(output_model, lower_deriv1, by = "year")
        
        # Get emergin status summary for output
        emerging_status_output <-
          output_model %>%
          filter(!!sym(year) %in% eval_years) %>%
          select(
            !!sym(taxonKey),
            .data$year,
            .data$em_status,
            .data$growth,
            .data$method
          )
        
        # Create plot with conf. interval + colour for status
        if (!is.null(name)) {
          ptitle <- name
        } else {
          ptitle <- "No Species Name Provided"
        }
        plot_gam <- plot_ribbon_em(
          df_plot = output_model,
          x_axis = year,
          y_axis = y_var,
          x_label = x_label,
          y_label = y_label,
          ptitle = ptitle
        )
        if (saveplot == TRUE) {
          if (str_ends(dir_name, pattern = "/")) {
            # remove "/" at the end
            dir_name <- str_sub(dir_name, end = -2)
          }
          file_name <- paste0(dir_name, "/", ptitle, ".png")
          if (isTRUE(verbose)) {
            print(paste("Output plot:", file_name))
          }
          ggsave(filename = file_name, plot_gam)
        }
      }
    }
  } else {
    if (verbose) {
      if (!is.null(name) & !is.null(taxon_key)) {
        warning(paste0(
          "Too few data for applying GAM (",
          method_em,
          ") to ", name, " (", taxon_key, ").\n"
        ))
      } else {
        if (!is.null(name)) {
          warning(paste0(
            "Too few data for applying GAM (",
            method_em, ") to ", name, ".\n"
          ))
        } else {
          if (!is.null(taxon_key)) {
            warning(paste0(
              "Too few data for applying GAM (",
              method_em,
              ") to taxon key: ", taxon_key, ".\n"
            ))
          } else {
            warning(paste0(
              "Too few data for applying GAM (",
              method_em, ").\n"
            ))
          }
        }
      }
    }
  }
  
  return(list(
    em_summary = emerging_status_output,
    model = model,
    output = output_model,
    first_derivative = deriv1,
    second_derivative = deriv2,
    plot = plot_gam
  ))
}