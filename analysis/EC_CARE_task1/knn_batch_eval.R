library(class)

setwd("//wsl$/Ubuntu/home/rashedul/project/ProSEC/analysis/EC_CARE_task1/")

# List all model output files
dist_files <- list.files("./fastas", pattern = "^sample_30-50_2000_.*\\.csv$", full.names = TRUE)

# Load training metadata
training_metadata <- read.csv("protein_train.csv")

# Define test metadata files
test_metadata_files <- c("30-50_protein_test.csv", "30_protein_test.csv")

# Prepare results storage
results <- data.frame(
  Model = character(),
  TestSet = character(),
  TotalPredictions = integer(),
  CorrectPredictions = integer(),
  Accuracy = numeric(),
  stringsAsFactors = FALSE
)

for (dist_file in dist_files) {
  # Extract model name from file
  model_name <- sub(".*sample_30-50_2000_(.*)\\.csv$", "\\1", dist_file)
  
  # Read distance matrix
  print(paste("Processing model:", model_name))
  dist_mat <- read.csv(dist_file, row.names = 1)
  dist_mat_proteins <- rownames(dist_mat)
  training_proteins <- training_metadata$Entry
  
  # Find common proteins for training
  common_train_proteins <- intersect(dist_mat_proteins, training_proteins)
  train_mat <- dist_mat[common_train_proteins, common_train_proteins]
  train_labels <- training_metadata$EC1[match(rownames(train_mat), training_metadata$Entry)]
  
  for (test_metadata_file in test_metadata_files) {
    test_metadata <- read.csv(test_metadata_file)
    test_proteins <- test_metadata$Entry
    common_test_proteins <- intersect(dist_mat_proteins, test_proteins)
    test_mat <- dist_mat[common_test_proteins, common_train_proteins]
    test_labels <- test_metadata$EC1[match(rownames(test_mat), test_metadata$Entry)]
    
    # Run k-NN 5 times and calculate average accuracy
    accuracies <- numeric(5)
    total_predictions <- 0
    correct_predictions <- 0
    for (i in 1:5) {
      knn_pred <- knn(train = train_mat, test = test_mat, cl = train_labels, k = 5)
      conf_matrix <- table(Predicted = knn_pred, Actual = test_labels)
      total_predictions <- sum(conf_matrix)
      correct_predictions <- sum(diag(conf_matrix))
      accuracies[i] <- (correct_predictions / total_predictions) * 100
    }
    avg_accuracy <- mean(accuracies)
    
    # Store results
    results <- rbind(
      results,
      data.frame(
        Model = model_name,
        TestSet = test_metadata_file,
        TotalPredictions = total_predictions,
        CorrectPredictions = correct_predictions,
        Accuracy = round(avg_accuracy, 2),
        stringsAsFactors = FALSE
      )
    )
  }
}

# Save results to CSV
write.csv(results, "knn_accuracy_summary_sample_30-50_2000.csv", row.names = FALSE)