# ML model options for UI:

ML_model_options <- selectInput("gmgmt5", "Select ML algorithms for prediction:",
                  c(
"Bagged CART" = "treebag",
"Bagged FDA using gCV Pruning" = "bagFDAGCV",
"Bagged Flexible Discriminant Analysis" = "bagFDA",
"Boosted Generalized Linear Model" = "glmboost",
"C5.0" = "C5.0",
"CART" = "rpart",
"Conditional Inference Random Forest" = "cforest",
"eXtreme Gradient Boosting_xgbLinear" = "xgbLinear",
"eXtreme Gradient Boosting_xgbTree" = "xgbTree",
"Flexible Discriminant Analysis" = "fda",
"Generalized Linear Model" = "glm",
"glmnet" = "glmnet",
"Parallel Random Forest" = "parRF",
"Partial Least Squares" = "kernelpls",
"Random Forest" = "rf",
"Regularized Random Forest" = "RRF",
"Regularized Random Forest_RRFglobal" = "RRFglobal",
"Rotation Forest" = "rotationForest",
"Rotation Forest_rotationForestCp" = "rotationForestCp",
"Single C5.0 Ruleset" = "C5.0Rules",
"Single C5.0 Tree" = "C5.0Tree",
"Sparse Distance Weighted Discrimination" = "sdwd"
), 
multiple = TRUE, selected = c("glmnet", "rf"))

