# ML model options for UI:

ML_model_options <- selectInput("gmgmt5", "Select ML algorithms:",
                  c(
                    "Naive Bayes" = "naive_bayes",
                    "Flexible Discriminant Analysis" = "fda",
                    "Support Vector Machines" = "svmRadial",
                    "Boosted Generalized Linear Model" = "glmboost",
                    "Random Forest" = "rf",
                    "C5.0" = "C5.0",
                    "glmnet" = "glmnet",
                    "Partial Least Squares" = "kernelpls"
                    ),
                    multiple = TRUE, selected = c("naive_bayes", "fda", "svmRadial","glmboost" )
                  )