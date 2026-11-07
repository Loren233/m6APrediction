#' Encode DNA Sequences
#'
#' @description Convert DNA sequences into nucleotide-position factors.
#' @param dna_strings A character vector of DNA sequences.
#' @return A data frame containing nucleotide-position factors.
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "GGACT"))
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Sites for Multiple Observations
#'
#' @description Predict m6A probabilities and statuses using a Random Forest model.
#' @param ml_fit A trained Random Forest model.
#' @param feature_df A data frame containing the seven input features.
#' @param positive_threshold The probability threshold for positive classification.
#' @return A data frame with predicted probabilities and statuses.
#' @examples
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction"
#' ))
#'
#' example_data <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction"
#' ))
#'
#' prediction_multiple(rf_fit, example_data,
#'                     positive_threshold = 0.6)
#'
#' @export
#' @import randomForest
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {

  stopifnot(all(c("gc_content", "RNA_type", "RNA_region",
                  "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer")
                %in% colnames(feature_df)))
  pred_df <- feature_df
  pred_df$RNA_type <- factor(
    pred_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  pred_df$RNA_region <- factor(
    pred_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )
  pred_df <- cbind(pred_df, dna_encoding(pred_df$DNA_5mer))
  prob <- predict(ml_fit, newdata = pred_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- prob
  feature_df$predicted_m6A_status <- ifelse(
    prob > positive_threshold,
    "Positive",
    "Negative"
  )
  return(feature_df)
}


#' Predict m6A for a Single Observation
#'
#' @description Predict m6A probability and status from one set of input features.
#' @param ml_fit A trained Random Forest model.
#' @param gc_content GC content of the site.
#' @param RNA_type RNA type.
#' @param RNA_region RNA region.
#' @param exon_length Exon length.
#' @param distance_to_junction Distance to exon junction.
#' @param evolutionary_conservation Evolutionary conservation score.
#' @param DNA_5mer A five-nucleotide DNA sequence.
#' @param positive_threshold The probability threshold for positive classification.
#' @return A named vector containing the predicted probability and status.
#' @examples
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction"
#' ))
#'
#' prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type,
                      levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region,
                        levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  pred <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  returned_vector <- c(
    pred$predicted_m6A_prob[1],
    pred$predicted_m6A_status[1]
  )
  names(returned_vector) <- c(
    "predicted_m6A_prob",
    "predicted_m6A_status"
  )
  return(returned_vector)
}
