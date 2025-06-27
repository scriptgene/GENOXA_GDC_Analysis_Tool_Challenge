# Source additional R scripts
source("01_library.R")
source("02_plot_utils.R")
source("03_utils.R")
source("04_GSEA.R")

# Split data frame --
# Split the data frame to create training and test data
split_data <- function(model_df,
                       train_size = 0.80,
                       seed = NULL) {
  set.seed(seed)
  train_index <- createDataPartition(model_df$group,
                                     p = train_size,
                                     list = FALSE
  )
  # Use the train_index to subset the data frame
  train_df <- model_df[train_index, ]
  test_df <- model_df[-train_index, ]
  
  # Remove rownames
  rownames(train_df) <- NULL
  rownames(test_df) <- NULL
  
  # Create a list with test and training data frames
  split_dataframes <- list()
  split_dataframes[[1]] <- train_df
  split_dataframes[[2]] <- test_df
  
  # Rename the items of the list
  names(split_dataframes) <- c("training", "test")
  return(split_dataframes)
}

# Train models -------------------------------------------------------------
# Train machine learning models on training data
train_models <- function(split_df,
                         resample_method = "repeatedcv",
                         resample_iterations = 10,
                         num_repeats = 3,
                         algorithm_list,
                         seed = NULL,
                         ...) {
  
  # If algorithm_list is not provided, use the default list of algorithms.
  if (missing(algorithm_list)) {
    algorithm_list <- c("rf", "naive_bayes", "glmboost", "xgbLinear")
  }
  
  # Set trainControl parameters for resampling
  set.seed(seed)
  fit_control <- trainControl(
    method = resample_method,
    number = resample_iterations,
    repeats = num_repeats,
    classProbs = TRUE
  )
  
  # Extract the training data set from the split_df object
  training_data <- split_df$training
  
  # Train models using ML algorithms from the algorithm_list.
  model_list <- lapply(
    setNames(algorithm_list, algorithm_list),
    function(x) {
      tryCatch(
        {
          set.seed(seed)
          message(paste0("\n", "Running ", x, "...", "\n"))
          train(group ~ .,
                data = training_data,
                trControl = fit_control,
                method = x,
                ...
          )
        },
        error = function(e) {
          message(paste0(
            x,
            " failed."
          ))
        }
      )
    }
  )
  
  message(paste0("Done!"))
  
  # Drop models that failed to build from the list
  model_list <- Filter(Negate(is.null), model_list)
  return(model_list)
}

# Get model results matrix
library(caret)

extract_and_sort_metrics <- function(model_list, metrics = c("Accuracy", "Kappa"), output_file = NULL)
  {
  # Define a function to safely extract metrics
  extract_model_metrics <- function(model, metric_names) {
    results <- model$results
    if (is.null(results)) {
      return(NULL)
    }
    
    # Check which metrics are present in the results
    available_metrics <- intersect(metric_names, colnames(results))
    
    # If no metrics are found, return NULL
    if (length(available_metrics) == 0) {
      return(NULL)
    }
    
    # Extract metrics, defaulting to NA if the metric is not available
    metrics_values <- sapply(available_metrics, function(metric) {
      if (metric %in% colnames(results)) {
        # Get the best result for each metric
        best_index <- which.max(results$Accuracy)
        results[best_index, metric]
      } else {
        NA
      }
    })
    
    # Return as a data frame
    data.frame(
      Model = model$method,
      Accuracy = metrics_values["Accuracy"],
      Kappa = metrics_values["Kappa"],
      stringsAsFactors = FALSE
    )
  }
  
  # Extract metrics from each model and create a data frame
  metrics_list <- lapply(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    metrics_df <- extract_model_metrics(model, metrics)
    if (!is.null(metrics_df)) {
      metrics_df$Model <- model_name
      return(metrics_df)
    } else {
      return(NULL)
    }
  })
  
  # Combine list into a data frame
  metrics_df <- do.call(rbind, metrics_list)
  
  # Ensure that metrics_df is not empty before sorting
  if (nrow(metrics_df) > 0) {
    # Sort the data frame based on Accuracy in descending order
    sorted_metrics_df <- metrics_df[order(-metrics_df$Accuracy, na.last = TRUE), ]
    
    # Print the sorted metrics
    print(sorted_metrics_df)
    
    # Optionally, save the results to a file
    if (!is.null(output_file)) {
      write.csv(sorted_metrics_df, file = output_file, row.names = FALSE)
      message(paste("Results saved to", output_file))
    }
    
    return(sorted_metrics_df)
  } else {
    message("No metrics were extracted. Check your model results.")
    return(NULL)
  }
}


# Model Performance plots - DOT PLOTS
perf_plot <- function(model_list,
                      text_size = 10,
                      palette = "plasma",
                      ggtitle = NULL){
  
  #binding for global variable
  method <- value <- NULL
  
  # Extract resample data from the model_list object
  resample_data <- resamples(model_list)
  
  # Convert the performance data into long-form format for plotting.
  plot_data <- reshape2::melt(resample_data$values)
  
  # Extract method and metric information from the variable column and add them
  # to new columns.
  plot_data$method <- sapply(
    strsplit(
      as.character(plot_data$variable),
      "~"
    ),
    "[", 1
  )
  plot_data$metric <- sapply(
    strsplit(
      as.character(plot_data$variable),
      "~"
    ),
    "[", 2
  )
  
  # Remove the variable column as it is no longer needed.
  plot_data$variable <- NULL
  
  # Assign palette color
  pal_col <- set_col(palette, 1)
  
  # Make plots
  text_size = text_size
  # Make dot plots
  perform_plot <- ggplot2::ggplot(
    plot_data,
    aes(
      x = method,
      y = value
    )
  ) +
    ggplot2::stat_summary(
      fun.data = "mean_se",
      lwd = text_size * 0.2,
      color = pal_col
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "pointrange",
      size = text_size * 0.2,
      pch = 16,
      position = "identity",
      color = pal_col
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~metric,
                        scales = "free"
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    extra_facet_theme() +
    ggtitle(ggtitle) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = text_size, face = "bold", hjust = 0.5),
      plot.title = ggplot2::element_text(color = "black", size = text_size * 1.2, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.text.y = element_text(size = text_size * 0.8, face = "bold"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  return(perform_plot)
}

# prob list
# --------------------Test the models ------------------------------------
## Test machine learning models on test data
test_models <- function(model_list,
                        split_df,
                        type = "raw",
                        save_confusionmatrix = FALSE,
                        file_path = NULL,
                        ...) {
  
  # Extract test data from the split_df object
  test_data <- split_df$test
  
  # Predict test data
  pred_list <- lapply(
    model_list,
    function(x) {
      message(paste0(
        "\n",
        "Testing ",
        x$method,
        "...",
        "\n"
      ))
      predict(x,
              test_data,
              type = type
      )
    }
  )
  message(paste0("\n", "Done!"))
  
  # Get confusion matrices and associated statistics
  if (type == "raw") {
    cm_list <- lapply(
      pred_list,
      function(x) {
        confusionMatrix(
          x,
          test_data$group
        )
      }
    )
    
    # Convert c.matrices to long-form data frames
    cm_df <- lapply(
      cm_list,
      function(x) reshape2::melt(as.table(x))
    )
    
    #Get the list of methods
    method_list <- names(cm_df)
    
    #Add method names to the data frames
    cm_dfm <- lapply(
      seq_along(method_list),
      function(x) {
        cm_df[[x]]["method"] <- method_list[x]
        cm_df[[x]]
      }
    )
    
    #Combine all data frames into one
    cm_dfm_long <- do.call(
      "rbind",
      cm_dfm
    )
    
    #Add column names
    colnames(cm_dfm_long) <- c("Prediction", "Reference", "Value", "Method")
    
    # Save data if required
    if (save_confusionmatrix) {
      if (is.null(file_path)) {
        file_path <- tempdir()
      }
      write.table(cm_dfm_long,
                  file = paste0(file_path, "/Confusion_matrices.txt"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE
      )
    }
    
    return(cm_dfm_long)
  }
  
  return(pred_list)
}

# ROC plot
roc_plot <- function(probability_list,
                     split_df,
                     ...,
                     multiple_plots = TRUE,
                     text_size = 10,
                     palette = "viridis",
                     title = NULL,
                     save = FALSE,
                     file_path = NULL,
                     file_name = "ROC_plot",
                     file_type = "pdf",
                     plot_width = 7,
                     plot_height = 7,
                     dpi = 80) {
  
  #binding for global variable
  method <- auc <- NULL
  
  # Generate ROC objects and add them to a list
  roc_list <- lapply(
    probability_list,
    function(x) {
      pROC::roc(x[[1]],
                response = split_df$test$group,
                percent = TRUE,
                ...
      )
    }
  )
  
  # Extract auc values and create a dataframe for plot labels
  auc_list <- lapply(
    roc_list,
    function(x) {
      round(
        x$auc,
        1
      )
    }
  )
  auc_labels <- as.data.frame(do.call(
    "rbind",
    auc_list
  ))
  auc_labels$method <- rownames(auc_labels)
  rownames(auc_labels) <- NULL
  names(auc_labels) <- c(
    "auc",
    "method"
  )
  
  # Extract sensitivities and specificities from ROC object and make a data
  # frame for plotting
  sens_list <- lapply(
    roc_list,
    function(x) x$sensitivities
  )
  sens <- do.call(
    "rbind",
    sens_list
  )
  sens_melted <- reshape2::melt(sens)
  names(sens_melted) <- c(
    "method",
    "number",
    "sensitivity"
  )
  
  spec_list <- lapply(
    roc_list,
    function(x) x$specificities
  )
  spec <- do.call(
    "rbind",
    spec_list
  )
  spec_melted <- reshape2::melt(spec)
  names(spec_melted) <- c(
    "method",
    "number",
    "specificity"
  )
  
  # Merge to create the final data frame for plotting
  roc_plotdata <- merge(
    sens_melted,
    spec_melted
  )
  
  # Order data frame by sensitivity. Important for geom_step
  roc_plotdata <- roc_plotdata[order(
    roc_plotdata$method,
    roc_plotdata$sensitivity
  ), ]
  
  # Make ROC curves
  rocplots <- ggplot2::ggplot(roc_plotdata, aes(
    x = specificity,
    y = sensitivity,
    colour = method
  )) +
    ggplot2::geom_step(
      direction = "hv",
      lwd = text_size * 0.1
    ) +
    ggplot2::scale_x_reverse(lim = c(100, 0), ) +
    viridis::scale_color_viridis(
      discrete = TRUE,
      option = palette,
      begin = 0.3,
      end = 0.7
    ) +
    ggplot2::xlab("Specificity (%)") +
    ggplot2::ylab("Sensitivity (%)") +
    ggplot2::geom_abline(
      intercept = 100,
      slope = 1,
      color = "grey60",
      linetype = 2,
      show.legend = FALSE
    )
  
  if (multiple_plots == FALSE) {
    rocplots1 <- rocplots +
      extra_theme() +
      ggplot2::ggtitle(title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = text_size * 1.2, face = "bold", hjust = 0.5),
        strip.text = ggplot2::element_text(color = "black", size = text_size, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(
          size = text_size * 0.9,
          face = "bold"
        ),
        axis.title = element_text(size = text_size)
      )
  } else {
    rocplots1 <- rocplots +
      ggplot2::facet_wrap(~method) +
      extra_facet_theme() +
      ggplot2::ggtitle(title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(color = "black", size = text_size * 1.2, face = "bold", hjust = 0.5),
        strip.text = ggplot2::element_text(color = "black", size = text_size, face = "bold", hjust = 0.5),
        text = element_text(size = text_size),
        legend.position = "none",
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(
          size = text_size,
          angle = 90
        )
      ) +
      ggplot2::geom_label(
        data = auc_labels,
        aes(label = paste0(
          "AUC: ",
          auc,
          "%"
        )),
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = -0.5,
        inherit.aes = FALSE,
        label.size = NA
      )
  }
  
  #Set temporary file_path if not specified
  if(is.null(file_path)){
    file_path <- tempdir()
  }
  
  if (save == TRUE) {
    ggplot2::ggsave(paste0(
      file_path,
      "/",
      file_name,
      ".",
      file_type
    ),
    rocplots1,
    dpi = dpi,
    width = plot_width,
    height = plot_height
    ) +
      ggplot2::ggtitle(title)
  }
  
  return(rocplots1)
}