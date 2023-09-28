library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data_subset.csv')
read_expression_table <- function(filename) {
  expression_csv <- read_csv("data/example_intensity_data_subset.csv") %>%
    pivot_longer(cols = c(starts_with("GSM")), #or col=-1
                 names_to = "subject_id") %>%
    pivot_wider(names_from = c("probe"))
  return (expression_csv)
}
#solution2
read_expression_table <- function(filename) {
  expr_mat <- read.table(
    filename,
    header = TRUE,
    sep = ",",
    row.names = 1
  )
  expr_mat <-  tibble::as_tibble(
    t(expr_mat),
    rownames = "subject_id"
  )
  return(expr_mat)
}


#' Load Metadata from Specified CSV File
#'
#' This function reads the provided CSV file into a dataframe.
#'
#' @param filepath (character) The path to the CSV file.(data/proj_metadata.csv)
#'
#' @return A dataframe containing the loaded metadata.
#'

load_metadata <- function(filepath) {
  
  # TODO: Use appropriate function to read in the CSV file
  #metadata <-   
  
  # Return the loaded metadata
  #return(metadata)
  metadata <- read_csv(filepath) %>% as_tibble()
  return(metadata)
}



#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {
  output_str <- gsub("\\.", "_", str)
  #or stringr::str_replace_all(str, "\\.", "_") or [.] instead of \\.
  return (output_str)
}


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  #rename the column names
  #col_meta <- colnames(metadata)
  #col_meta <- period_to_underscore(col_meta)
  #colnames(metadata) <- col_meta
  colnames(data) <-
    lapply(colnames(data), period_to_underscore)
  #select columns
  data <- select(data, Sex, Age=Age_at_diagnosis, TNM_Stage, Tumor_Location, geo_accession, KRAS_Mutation,
                 Subtype=SixSubtypesClassification, Batch=normalizationcombatbatch)
  #or 
  #out <- dplyr::rename(data, "age" = Age_at_diagnosis,
  #                     "Subtype" = SixSubtypesClassification,
  #                     "Batch" = normalizationcombatbatch) %>% dplyr::select(...)
  return (data)
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
  #data$Stage <-factor(paste("stage", data$TNM_Stage, sep = ""), levels = unique(paste("stage", data$TNM_Stage, sep = "")))
  data <- data %>%
    mutate(Stage = factor(paste("stage", TNM_Stage, sep = ""), levels = unique(paste("stage", TNM_Stage, sep = "")))) #or paste0() ignore the space between 
  #solution2
  out <- dplyr::mutate(data, 
                       Stage = as.factor(
                         stringr::str_c("stage", TNM_Stage) #similar to paste0()
                       ))
  return (data)
}


#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  mean_age <- data %>% filter(Sex == sex) %>% summarise(mean(Age)) %>% pull()
  #solution2 - baseR
  data_sex <- data[data$Sex == sex,]
  avg <- mean(data_sex$Age)
  #solution3 -- stringR
  s <- stringr::str_which(data$Sex, sex)
  avg <- mean(data$Age[s])
  return (mean_age)
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage. Name the newly created column containing the average, 'mean_avg'
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  mean_avg <- data %>% group_by(Stage) %>% summarise(mean(Age))
  return (mean_avg)
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  cross_tab <- data %>% group_by(Stage) %>% 
    summarize(C3 = sum(Subtype == "C3"), C4 = sum(Subtype == "C4"))
  #table() function
  #count() function
  #solution2
  out <- dplyr::group_by(data, Subtype, Stage) %>%
    dplyr::summarise("n" = dplyr::n()) %>%
    tidyr::pivot_wider(
      names_from = Subtype,
      values_from = n,
      values_fill = 0
    )
  return (cross_tab)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
summarize_expression <- function(exprs) {
  main_exp <- exprs %>% select(2:1001) %>% apply( 2,mean) # select(!c("subject_id))/select(-subject_id)
  variance <- exprs %>% select(2:1001) %>% apply( 2,var)
  sum_tb <- data.frame(main_exp = main_exp, variance=variance, 
                       probe=colnames(exprs %>% select(ends_with("at")))) %>% as_tibble()
  #solution2- tidyvirse
  summarized <- pivot_longer(
    exprs,
    cols = -1,
    names_to = "probe",
    values_to = "expression"
  ) %>%
    group_by(probe) %>%
    summarise(mean = mean(expression), var = var(expression))
  #solution3 - baser
  datal <- exprs %>% dplyr::select(-subject_id)
  mean_expr <- colMeans(datal)
  var_expr <- sapply(datal, var)
  tbs <- tibble(
    mean_exp = mean_expr,
    variance = var_expr,
    probe = colnames(exprs)[-1] #remove the first element
  )
  return (sum_tb)
}
