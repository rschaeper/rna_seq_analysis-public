# script to clean up metadata
# outputs metadata.xlsx, metadata.csv

# load libraries
library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)
library(openxlsx)

# column descriptions -----------------------------------------------------
# description dictionary
description_dict <- list(
  CEGAT_ID = "Identifier used by CeGaT for bulk RNAseq",
  INTERNAL_ID = "Identifier used by us for bulk RNAseq",
  STUDY_ID = "Identifier for each subject",
  SAMPLE_ID = "Identifier for each sample extraction",
  STUDY = "Identifier for study season",
  VISIT = "Visit number",
  DAY = "Day within vaccination kinetic",
  NOMINAL_DAY = "Day on which the sample was idealy collected according to study design (V1 = 0, V2 = 1, V3 = 7, V4 = 30, V5 = 90, V6 = 180)",
  STUDY_GROUP = "Identifier for disease state (A = HC, B = SMM/MM, C = NSCLC)",
  SUBJECT_GROUP_NUMBER = "Identifier for subject within study group",
  DIAGNOSIS = "Diagnosis of the subject (HC = Healthy, SMM = Smouldering Multiple Myeloma, MM = Multiple Myeloma, NSCLC = Lung Cancer)",
  SEX = "Biological sex of the subject, as stated by anamnesis questionnaire",
  YEAR_OF_BIRTH = "Year of birth, as stated by anamnesis questionnaire",
  AGE = "Age of the subject in years (2022 - YEAR_OF_BIRTH for SYS01 and 2023 - YEAR_OF_BIRTH for SYS03)",
  AGE_BELOW_60 = "Whether the subject has an AGE below 60 and thus receives a lower vaccine dose than those with AGE above 60",
  HEIGHT = "Height of the subject in cm",
  WEIGHT = "Weight of the subject in kg",
  TRI_ELISA = "Titer Response Index (TRI). Collapses multiple antibody‐titer readouts (HAI and neutralization, across several vaccine antigens) into a single scalar per person. Only available for V3, V4 and V5.",
  TRI_ELISA_V3 = "Titer Response Index (TRI) at V3",
  TRI_ELISA_V4 = "Titer Response Index (TRI) at V4",
  TRI_ELISA_V5 = "Titer Response Index (TRI) at V5",
  IgG = "IgG antibody levels as Area Under The Curve (AUC) value calculated from raw ECL-ELISA counts (measured as blood plasma IgG capable of binding a plate coated with vaccine antigen; MSD platform)",
  IgG_V1 = "IgG antibody levels at V1 as Area Under The Curve (AUC) value calculated from raw ECL-ELISA counts (measured as blood plasma IgG capable of binding a plate coated with vaccine antigen; MSD platform)",
  IgG_V3 = "IgG antibody levels at V3 as Area Under The Curve (AUC) value calculated from raw ECL-ELISA counts (measured as blood plasma IgG capable of binding a plate coated with vaccine antigen; MSD platform)",
  IgG_V1_to_V3_l2fc = "IgG antibody level increase from V1 baseline to V3 (log2(V3_auc/V1_auc))",
  IgG_V4 = "IgG antibody levels at V4 as Area Under The Curve (AUC) value calculated from raw ECL-ELISA counts (measured as blood plasma IgG capable of binding a plate coated with vaccine antigen; MSD platform)",
  IgG_V1_to_V4_l2fc = "IgG antibody level increase from V1 baseline to V4 (log2(V4_auc/V1_auc))",
  VISIT_DATE = "Date of sample extraction, as stated by sample documentation",
  RNA_SAMPLE_SEQUENCING_DATE = "Date of bulk RNAseq",
  RNA_SAMPLE_MONTHS_SINCE_VISIT = "Months since bulk RNAseq (RNA_SAMPLE_SEQUENCING_DATE - VISIT_DATE)",
  SUBJECT_LOCATION = "Geographical location of subject, determined by SUBJECT_GROUP_NUMBER",
  HSCT_DATE = "Date of hematopoietic stem cell transplantation, as stated by clinical information",
  LENALIDOMIDE_START_DATE = "Date on which Lenalidomide was first given to subject, as stated by clinical information",
  LENALIDOMIDE_END_DATE = "Date on which Lenalidomide was last given to subject, as stated by clinical information",
  MONTHS_SINCE_HSCT = "Months since hematopoietic stem cell transplantation (VISIT_DATE - HSCT_DATE)",
  MONTHS_ON_LENALIDOMIDE = "Months on Lenalidomide treatment before sample extraction (VISIT_DATE - LENALIDOMIDE_START_DATE); value of 0 means that Lenalidomide had been given, but treatment finished before our systems immunology trial",
  LC_TYPE = "Detailed lung cancer type for C subjects, be aware that Adeno and Pleca are not exclusive to each other",
  LC_STAGE = "Stage of lung cancer for C subjects",
  LC_THERAPY_START_DATE = "Date on which lung cancer therapy was started, as stated by clinical information",
  HEMOGLOBIN_LEVEL = "Hemoglobin level of the participant in g/dL",
  CANCER_THERAPY_AT_V1 = "The treatment given to MM or NSCLC subjects at V1 (w&w = 'watch and wait', NSCLC treatment before current treatment has been discarded)",
  CRIZOTINIB_AT_V1 = "Whether subject receives lung cancer drug Crizotinib at V1",
  ISATUXIMAB_AT_V1 = "Whether subject receives multiple myeloma drug Isatuximab at V1",
  LENALIDOMIDE_AT_V1 = "Whether subject receives multiple myeloma drug Lenalidomide at V1",
  PEMBROLIZUMAB_AT_V1 = "Whether subject receives lung cancer drug Pembrolizumab at V1",
  PEMETREXED_AT_V1 = "Whether subject receives lung cancer drug Pemetrexed at V1",
  DEXAMETHASONE_AT_V1 = "Whether subject receives multiple myeloma drug Dexamethasone at V1",
  CARBOPLATIN_AT_V1 = "Whether subject receives lung cancer drug Carboplatin at V1",
  PACLITAXEL_AT_V1 = "Whether subject receives chemotherapy drug Paclitaxel at V1",
  CEMIPLIMAB_AT_V1 = "Whether subject receives lung cancer drug Cemiplimab at V1",
  DATO_DX_AT_V1 = "Whether subject receives lung cancer drug Dato-DXd at V1",
  MONTHS_ON_LC_THERAPY = "Months on lung cancer therapy before sample extraction (VISIT_DATE - LC_THERAPY_START_DATE)",
  CHEMOTHERAPY_AT_V1 = "Whether the subject receives any kind of chemotherapy drug (Pemetrexed, Paclitaxel or Carboplatin)",
  PD1_INHIBITOR_AT_V1 = "Whether the subject receives a PD-1 inhibiting therapeutic antibody (Pembrolizumab or Cemiplimab)",
  INFLUENZA_VACCINATIONS_LAST_5_YEARS = "Number of influenza vaccinations received by the subject within 5 years prior to our systems immunology trial",
  COVID19_VACCINATIONS = "Total number of COVID19 vaccinations received by the subject",
  HISTORY_OF_CANCER = "Whether subject had a cancer diagnosis in the past or not",
  REGULAR_MEDICATION = "Whether subject receives regular medication apart from cancer specific treatments",
  RNA_SAMPLE_VISITS_COMPLETE = "Whether a full kinetic (V1 to V6) was obtained for subject during bulk RNAseq",
  RNA_SAMPLE_INPUT = "Amount of RNA that was available after blood RNA extraction for bulk RNAseq at CeGaT",
  RNA_SAMPLE_DNA_DIGESTION = "Whether an additional DNA digestion step was performed using DNAse at CeGaT",
  RNA_SAMPLE_SUSPECTED_CONTAMINATION_WITH = "Other samples of which traces could be found within given sample",
  RNA_SAMPLE_RESEQUENCING = "Whether resequencing step was performed at CeGaT",
  AGE_CENTERED_SCALED = "AGE column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  HEIGHT_CENTERED_SCALED = "HEIGHT column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  WEIGHT_CENTERED_SCALED = "WEIGHT column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  IgG_CENTERED_SCALED = "IgG column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  IgG_V1_CENTERED_SCALED = "IgG_V1 column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  IgG_V3_CENTERED_SCALED = "IgG_V3 column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  IgG_V4_CENTERED_SCALED = "IgG_V4 column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  HEMOGLOBIN_LEVEL_CENTERED_SCALED = "HEMOGLOBIN_LEVEL column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  MONTHS_SINCE_HSCT_CENTERED_SCALED = "MONTHS_SINCE_HSCT column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  MONTHS_ON_LENALIDOMIDE_CENTERED_SCALED = "MONTHS_ON_LENALIDOMIDE column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  MONTHS_ON_LC_THERAPY_CENTERED_SCALED = "MONTHS_ON_LC_THERAPY column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)",
  RNA_SAMPLE_MONTHS_SINCE_VISIT_CENTERED_SCALED = "RNA_SAMPLE_MONTHS_SINCE_VISIT column centered around zero (subtract mean from value) and scaled to unit variance (divide by standard deviation)"
)

# note on date format:
# in some cases information was available down to the specific day, sometimes just to the month.
# in order to preserve data, a format of YYYY-MM-DD was chosen, with the option of DD being NA

# centered and scaled versions of covariates are appended as additional columns with suffix _CENTERED_SCALED to column name

# helper functions --------------------------------------------------------

# function to clean the dates (use with sapply to transform column)
clean_date <- function(date_string) {
  if (is.na(date_string) || date_string == "" || is.null(date_string)) {
    return(NA_character_) # handle missing or empty values
  }
  
  # try different formats, including "Sep/ 23"
  parsed_date <- parse_date_time(date_string, c("b %d, %Y", "b-%y", "d/m/Y", "b/ %y"), quiet = TRUE)
  
  if (!is.na(parsed_date)) {
    if (day(parsed_date) == 1 && grepl("Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|Jan|Feb|Mär|Apr|Mai|Jun|Jul|Aug|Sep|Okt|Nov|Dez", date_string)) {
      return(format(parsed_date, "%Y-%m-NA")) # no day, format as YYYY-MM-NA
    } else {
      return(format(parsed_date, "%Y-%m-%d")) # format as YYYY-MM-DD
    }
  } else {
    return(NA_character_) # if it doesn't match any format, return NA
  }
}

# function to calculate month differences (from format YYYY-MM-DD), blind for the day as it is sometimes NA (use with mapply)
calculate_month_diff <- function(date1, date2) {
  if (is.na(date1) || is.na(date2) || date1 == "None" || date2 == "None") {
    return(NA_character_) # Return NA if either date is missing or "None"
  }
  
  # extract year and month
  year1 <- as.numeric(substr(date1, 1, 4))
  month1 <- as.numeric(substr(date1, 6, 7))
  year2 <- as.numeric(substr(date2, 1, 4))
  month2 <- as.numeric(substr(date2, 6, 7))
  
  # calculate the difference in months
  month_diff <- (year2 - year1) * 12 + (month2 - month1)
  return(month_diff)
}

# assemble metadata -------------------------------------------------------

# construct from excel file mapping CeGaT IDs to external IDs
metadata_sys01 <- read.csv("data/meta/S11664_Liste_CeGaT_ID_External_ID.csv", sep = ";")
metadata_sys03 <- read.csv("data/meta/SYS03_S11947_Liste_CeGaT ID_External ID.csv", sep = ";")

# make sys03's column names match sys01's
colnames(metadata_sys03) <- colnames(metadata_sys01)

# bind the rows together
metadata <- rbind(metadata_sys01, metadata_sys03)

# rename columns CeGat_ID and External_ID for consistency (style: COLUMN_CASE)
# also, rename EXTERNAL_ID to INTERNAL_ID, as it is now what we would refer to
metadata <- dplyr::rename(metadata, CEGAT_ID = CeGaT_ID, INTERNAL_ID = External_ID)

# discard day information in SAMPLE_ID
metadata$SAMPLE_ID <- stringr::str_extract(metadata$INTERNAL_ID, "^[^_]+(?:_[^_]+){3}")
# remove spaces within SAMPLE_ID string (unfortunately it does occur!)
metadata$SAMPLE_ID <- gsub("\\s+", "", metadata$SAMPLE_ID)

# split SAMPLE_ID to create STUDY, STUDY_ID and VISIT columns
metadata <- tidyr::separate(data = metadata, col = SAMPLE_ID, into = c("STUDY", "STUDY_ID", "VISIT"), sep = "_(?=[A-Z])", extra = "merge", remove = FALSE)

# STUDY_ID should be of the form SYSXX_Y_ZZZ, thus, we concatenate STUDY and STUDY_ID and keep those columns
metadata <- metadata %>% tidyr::unite(col = STUDY_ID, STUDY, STUDY_ID, sep = "_", remove = FALSE)

# add STUDY_GROUP column
metadata$STUDY_GROUP <- stringr::str_extract(metadata$STUDY_ID, "A|B|C")

# add DIAGNOSIS column
# add diagnosis for group A subjects
metadata$DIAGNOSIS[metadata$STUDY_GROUP == "A"] <- "HC"
# note: consider differentiating lung cancer diagnoses more precisely and
# add diagnosis for group C subjects
metadata$DIAGNOSIS[metadata$STUDY_GROUP == "C"] <- "NSCLC"

# add column DAY
metadata$DAY <- stringr::str_extract(metadata$INTERNAL_ID, "\\d+$")

# make sure to remove leading and trailing whitespace from every character column
metadata <- metadata %>%
  mutate(across(where(is.character), trimws))


# add visit dates, rna sample sequencing dates and subject location ---------------------------------------------------------

# initialize list to store sample documentations
sample_docs <- list()

sample_docs[["SYS01"]] <- read.csv("data/meta/SYS01_sample documentation.csv", sep = ";")
sample_docs[["SYS03"]] <- read.csv("data/meta/SYS03_sample_documentation MODIFIED.csv", sep = ";")

# the code section below can be used to check which kind of changes were applied by hand
# library(diffdf)
# 
# df_old <- read.csv("data/meta/SYS03_sample_documentation.csv", sep = ";")
# df_new <- read.csv("data/meta/SYS03_sample_documentation MODIFIED.csv", sep = ";")
# 
# # this yields data.frame of all differences
# diffs <- diffdf(df_old, df_new)
# 
# # see summary
# print(diffs)

# loop over sample docs
for (name in names(sample_docs)) {
  sample_docs[[name]] <- sample_docs[[name]] %>%
    # remove rows where everything is empty or NA
    filter(!if_all(everything(), ~ .x %in% c("", NA))) %>%
    # trim whitespace from character columns
    mutate(across(where(is.character), trimws)) %>%
    mutate(
      Visit = str_extract(Visit.day, "V\\d+"),
      SAMPLE_ID = paste0(Participant.ID, "_", Visit),
      VISIT_DATE = if_else(
        is.na(dmy(Collection.date)), #try day month year
        if_else(is.na(ymd(Collection.date)), # if day month year fails, try year month day.
                as.character(Collection.date), # if all fails, return the value as is.
                as.character(ymd(Collection.date))),
        as.character(dmy(Collection.date))
      ),
      VISIT_DATE = ymd(VISIT_DATE), #convert to date object.
      VISIT_DATE = format(VISIT_DATE, "%Y-%m-%d") #format to year month day
    ) %>%
    dplyr::select(-Participant.ID, -Visit.day)
}

sample_doc <- bind_rows(sample_docs)

metadata <- dplyr::left_join(metadata, sample_doc[c("SAMPLE_ID", "VISIT_DATE")], by = "SAMPLE_ID")

# create NOMINAL_DAY column
# define nominal days for each visit
nominal <- tibble(
  VISIT = c("V1","V2","V3","V4","V5","V6"),
  nominal_day = c(0, 1, 7, 30, 90, 180)
)

metadata <- metadata %>%
  left_join(nominal, by = "VISIT") %>%
  rename(NOMINAL_DAY = nominal_day)

# create RNA_SAMPLE_SEQUENCING_DATE column
# note: precision is currently just the moth
metadata <- metadata %>%
  mutate(
    RNA_SAMPLE_SEQUENCING_DATE = case_when(
      STUDY == "SYS01" ~ "2024-10-NA",
      STUDY == "SYS03" ~ "2025-01-NA",
      TRUE ~ NA_character_  # for everything else, leave it blank
    )
  )

# create RNA_SAMPLE_MONTHS_SINCE_VISIT
metadata$RNA_SAMPLE_MONTHS_SINCE_VISIT <- mapply(calculate_month_diff, metadata$VISIT_DATE, metadata$RNA_SAMPLE_SEQUENCING_DATE)

# create SAMPLE_LOCATION column based on STUDY_ID
metadata <- metadata %>%
  mutate(
    # extract the SUBJECT_GROUP_NUMBER part using regular expressions
    SUBJECT_GROUP_NUMBER = as.numeric(gsub(".*_(\\d{3})$", "\\1", STUDY_ID)),
    # create SUBJECT_LOCATION based on the SUBJECT_GROUP_NUMBER
    SUBJECT_LOCATION = case_when(
      (SUBJECT_GROUP_NUMBER >= 0 & SUBJECT_GROUP_NUMBER < 200) ~ "Tübingen",
      SUBJECT_GROUP_NUMBER >= 300 & SUBJECT_GROUP_NUMBER < 400 ~ "Stuttgart",
      SUBJECT_GROUP_NUMBER >= 500 & SUBJECT_GROUP_NUMBER < 600 ~ "Esslingen",
      TRUE ~ NA_character_ # if the SUBJECT_GROUP_NUMBER doesn't match, put NA
    )
  )

# add B subject info ---------------------------------------------------

# columns created or added to: DIAGNOSIS, HSCT_DATE, "LENALIDOMIDE_START_DATE", "LENALIDOMIDE_END_DATE", "MONTHS_ON_LENALIDOMIDE", "HEMOGLOBIN_LEVEL"

# initialize list of b_subject infos
b_subjects_infos <- list()

# --- SYS01 ---
b_subjects_infos[["SYS01"]] <- read.csv("data/meta/sys01_gruppeB_treatmentInfo.csv", sep = ";")
# --- SYS03 ---
b_subjects_infos[["SYS03_tuebingen"]] <- read.csv("data/meta/SYS03_clinical-Information_GroupB_(S)MM_Tübingen MODIFIED.csv", sep = ";")
b_subjects_infos[["SYS03_stuttgart"]] <- read.csv("data/meta/SYS03_clinical-Information_GroupB_(S)MM_Stuttgart.csv", sep = ";")

for (name in names(b_subjects_infos)) {
  # remove leading and trailing whitespace from every character column
  b_subjects_infos[[name]] <- b_subjects_infos[[name]] %>%
    mutate(across(where(is.character), trimws))
  
  # format Study.ID
  if (name == "SYS01") {
    b_subjects_infos[[name]] <- b_subjects_infos[[name]] %>%
      mutate(Study.ID = gsub("^B", "SYS01_B_", Study.ID))
  } else {
    b_subjects_infos[[name]] <- b_subjects_infos[[name]] %>%
      mutate(Study.ID = gsub("^B", "SYS03_B_", Study.ID))
  }
  
  # rename Entität to DIAGNOSIS
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Entität"] <- "DIAGNOSIS"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Entität..MM.SMM."] <- "DIAGNOSIS"
  
  # rename stem cell transplant date columns to HSCT_DATE
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Zeit.seit.HSCT..M."] <- "HSCT_DATE"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Datum.HSCT"] <- "HSCT_DATE"
  
  # rename first Lenalidomide date columns to LENALIDOMIDE_START_DATE
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Zeit.seit.Len..M."] <- "LENALIDOMIDE_START_DATE"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Datum.erste.Len.Gabe"] <- "LENALIDOMIDE_START_DATE"
  
  # rename columns containing information about therapy at vaccination to CANCER_THERAPY_AT_V1
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Therapie.1..Impfung"] <- "CANCER_THERAPY_AT_V1"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Therapie.Impfung"] <- "CANCER_THERAPY_AT_V1"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Erhält.Len.zum..Zeitpunkt.der.Impfung"] <- "CANCER_THERAPY_AT_V1"
  
  # rename columns containing information about IVIG dates to IVIG_DATES
  # names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "Datum.IVIG.30.Tage.vor.und.100.Tage.nach.Impfung"] <- "IVIG_DATE"
  # names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "IVIG.Daten"] <- "IVIG_DATE"
  # names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "IVIG"] <- "IVIG_DATE"
  
  # rename columns containing information about hemoglobin to HEMOGLOBIN_LEVEL
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "HB"] <- "HEMOGLOBIN_LEVEL"
  names(b_subjects_infos[[name]])[names(b_subjects_infos[[name]]) == "HB.Wert"] <- "HEMOGLOBIN_LEVEL"
  
}

b_subjects_info <- bind_rows(b_subjects_infos)

# remove rows where Study.ID is empty or NA
b_subjects_info <- b_subjects_info %>%
  filter(!is.na(Study.ID) & Study.ID != "")

# rename Study.ID to STUDY_ID
b_subjects_info <- dplyr::rename(b_subjects_info, "STUDY_ID" = "Study.ID")

# format HSCT_DATE column
b_subjects_info$HSCT_DATE <- sapply(b_subjects_info$HSCT_DATE, clean_date)
b_subjects_info$LENALIDOMIDE_START_DATE <- sapply(b_subjects_info$LENALIDOMIDE_START_DATE, clean_date)

# format CANCER_THERAPY_AT_V1 column
b_subjects_info$CANCER_THERAPY_AT_V1[b_subjects_info$CANCER_THERAPY_AT_V1 == "Z.n. Lena-Erhaltung" | b_subjects_info$CANCER_THERAPY_AT_V1 == "nein"] <- "w&w"
# identify synonyms for "Lenalidomide"
b_subjects_info$CANCER_THERAPY_AT_V1[b_subjects_info$CANCER_THERAPY_AT_V1 == "Lena Erhaltung" |
                                       b_subjects_info$CANCER_THERAPY_AT_V1 == "Lena-Erhaltung" |
                                       b_subjects_info$CANCER_THERAPY_AT_V1 == "ja"] <- "Lenalidomide"
# treat remaining two cases of combination therapy
b_subjects_info$CANCER_THERAPY_AT_V1[b_subjects_info$CANCER_THERAPY_AT_V1 == "Lena/Dexa"] <- "Lenalidomide,Dexamethasone"
b_subjects_info$CANCER_THERAPY_AT_V1[b_subjects_info$CANCER_THERAPY_AT_V1 == "Isatuximab/Lena"] <- "Isatuximab,Lenalidomide"

# create LENALIDOMIDE_START_DATE and LENALIDOMIDE_END_DATE for the three MM subjects that have finished their Lenalidomide treatment in the past
b_subjects_info <- b_subjects_info %>%
  mutate(
    LENALIDOMIDE_START_DATE = ifelse(
      STUDY_ID %in% c("SYS03_B_114", "SYS03_B_116", "SYS03_B_120"),
      case_when(
        STUDY_ID == "SYS03_B_114" ~ "2014-07-NA",
        STUDY_ID == "SYS03_B_116" ~ "2018-07-NA",
        STUDY_ID == "SYS03_B_120" ~ "2018-03-NA",
        TRUE ~ LENALIDOMIDE_START_DATE # this line should not be reached, it keeps the data already in place
      ),
      LENALIDOMIDE_START_DATE # keep the existing value for other STUDY_IDs
    )
  )

b_subjects_info <- b_subjects_info %>%
  mutate(
    LENALIDOMIDE_END_DATE = case_when(
      STUDY_ID == "SYS03_B_114" ~ "2020-09-NA",
      STUDY_ID == "SYS03_B_116" ~ "2020-07-NA",
      STUDY_ID == "SYS03_B_120" ~ "2019-09-NA",
      TRUE ~ NA_character_  # for everyone else, leave it blank
    )
  )


# b_subjects_info to metadata
metadata <- dplyr::left_join(metadata, b_subjects_info[c("STUDY_ID", "DIAGNOSIS", "HSCT_DATE", "HEMOGLOBIN_LEVEL", "CANCER_THERAPY_AT_V1", "LENALIDOMIDE_START_DATE", "LENALIDOMIDE_END_DATE")], by = "STUDY_ID") %>%
  # combine DIAGNOSIS columns
  mutate(DIAGNOSIS = coalesce(DIAGNOSIS.x, DIAGNOSIS.y)) %>%
  dplyr::select(-DIAGNOSIS.x, -DIAGNOSIS.y)  # remove the old DIAGNOSIS columns

# create MONTHS_SINCE_HSCT covariate from HSCT_DATE and VISIT_DATE
metadata$MONTHS_SINCE_HSCT <- mapply(calculate_month_diff, metadata$HSCT_DATE, metadata$VISIT_DATE)

# create MONTHS_ON_LENALIDOMIDE covariate from LENALIDOMIDE_START_DATE and VISIT_DATE
metadata$MONTHS_ON_LENALIDOMIDE <- mapply(calculate_month_diff, metadata$LENALIDOMIDE_START_DATE, metadata$VISIT_DATE)

# for four more subjects we have a LENALIDOMIDE_START_DATE but no LENALIDOMIDE_END_DATE
# all these are MM and have CANCER_THERAPY_AT_V1 "w&w". Set MONTHS_ON_LENALIDOMIDE to 0 for them
metadata <- metadata %>%
  mutate(
    MONTHS_ON_LENALIDOMIDE = ifelse(
      !is.na(LENALIDOMIDE_START_DATE) & is.na(LENALIDOMIDE_END_DATE) & CANCER_THERAPY_AT_V1 == "w&w" & DIAGNOSIS == "MM",
      0, # Set MONTHS_ON_LENALIDOMIDE to 0
      MONTHS_ON_LENALIDOMIDE # Keep existing values for others
    )
  )


# make sure to set the months on lenalidomide to NA or 0 for the three MM subjects that have finished Lena treatment before vaccination
# create LENALIDOMIDE_START_DATE and LENALIDOMIDE_END_DATE for the three MM subjects that have finished their Lenalidomide treatment in the past
metadata <- metadata %>%
  mutate(
    MONTHS_ON_LENALIDOMIDE = ifelse(
      STUDY_ID %in% c("SYS03_B_114", "SYS03_B_116", "SYS03_B_120"),
      case_when(
        STUDY_ID == "SYS03_B_114" ~ "0",
        STUDY_ID == "SYS03_B_116" ~ "0",
        STUDY_ID == "SYS03_B_120" ~ "0",
        TRUE ~ MONTHS_ON_LENALIDOMIDE # this line should not be reached, it keeps the data already in place
      ),
      MONTHS_ON_LENALIDOMIDE # keep the existing value for other STUDY_IDs
    )
  )

# add C subject info ---------------------------------------------------
# initialize list to store C subjects infos
c_subjects_infos <- list()

c_subjects_infos[["SYS01"]] <- read.csv("data/meta/sys01_clinical_data_C MODIFIED.csv", sep = ";")
c_subjects_infos[["SYS03_tuebingen"]] <- read.csv("data/meta/SYS03_clinical-Information_GroupC_NSCLC_Tübingen.csv", sep = ";")
c_subjects_infos[["SYS03_stuttgart"]] <- read.csv("data/meta/SYS03_clinical-Information_GroupC_NSCLC_Stuttgart MODIFIED.csv", sep = ";")
c_subjects_infos[["SYS03_esslingen"]] <- read.csv("data/meta/SYS03_clinical-Information_GroupC_NSCLC_Esslingen.csv", sep = ";")

for (name in names(c_subjects_infos)) {
  # remove leading and trailing whitespace from every character column
  c_subjects_infos[[name]] <- c_subjects_infos[[name]] %>%
    mutate(across(where(is.character), trimws))
  
  
  # rename Participant.ID column to Study.ID
  if (name == "SYS01") {
    c_subjects_infos[[name]] <- c_subjects_infos[[name]] %>%
      mutate(STUDY_ID = Participant.ID,
             CANCER_THERAPY_AT_V1 = Therapie.bei.Impfung,
             LC_THERAPY_START_DATE = Datum.Beginn.der.Therapie) %>%
      # split Krankheitsstadium.bei.Impfung into LC_TYPE and LC_STAGE
      separate(Krankheitsstadium.bei.Impfung, into = c("LC_TYPE", "LC_STAGE"), sep = " Stadium ")
    
    # # clean up LC_TYPE
    # c_subjects_infos[[name]]$LC_TYPE <- str_replace_all(c_subjects_infos[[name]]$LC_TYPE, "NSCLC \\(|\\)", "")
    c_subjects_infos[[name]]$LC_TYPE <- str_replace(c_subjects_infos[[name]]$LC_TYPE, "PLECA", "Pleca")
    # 
    # #clean up LC_STAGE
    # c_subjects_infos[[name]]$LC_STAGE <- str_replace_all(c_subjects_infos[[name]]$LC_STAGE, "Stadium ", "")
  }
  
  if (name == "SYS03_tuebingen") {
    c_subjects_infos[[name]] <- c_subjects_infos[[name]] %>%
      mutate(
        STUDY_ID = Study.ID,
        CANCER_THERAPY_AT_V1 = Therapie.an.Tag.0,
        LC_THERAPY_START_DATE = Datum.Beginn.der.Therapie,
        HEMOGLOBIN_LEVEL = HB.Wert..g.dL.,
        LC_TYPE = Entität,
        LC_STAGE = Stadium
      )
  }
  
  if (name == "SYS03_stuttgart") {
    c_subjects_infos[[name]] <- c_subjects_infos[[name]] %>%
      mutate(
        STUDY_ID = paste0("SYS03_C_", gsub("C", "", Study.ID)),
        CANCER_THERAPY_AT_V1 = Therapie.an.Tag.0,
        LC_THERAPY_START_DATE = Datum.Beginn..der.Therapie,
        HEMOGLOBIN_LEVEL = HB.Wert..g.dL.,
        LC_TYPE = Entität,
        LC_STAGE = Stadium
      )
    
    # Clean up the LC_STAGE column
    c_subjects_infos[[name]]$LC_STAGE <- str_replace_all(c_subjects_infos[[name]]$LC_STAGE, "4", "IV")
    c_subjects_infos[[name]]$LC_STAGE <- str_replace_all(c_subjects_infos[[name]]$LC_STAGE, "3", "III")
    c_subjects_infos[[name]]$LC_STAGE <- str_replace_all(c_subjects_infos[[name]]$LC_STAGE, "a", "A")
    c_subjects_infos[[name]]$LC_STAGE <- str_replace_all(c_subjects_infos[[name]]$LC_STAGE, "b", "B")
  }
  
  if (name == "SYS03_esslingen") {
    c_subjects_infos[[name]] <- c_subjects_infos[[name]] %>%
      mutate(STUDY_ID = paste0("SYS03_C_", gsub("C", "", Study.ID)),
             CANCER_THERAPY_AT_V1 = Therapie.an.Tag.0,
             LC_THERAPY_START_DATE = Datum.Beginn.der.Therapie,
             HEMOGLOBIN_LEVEL = HB.Wert..g.dL.,
             LC_TYPE = Entität,
             LC_STAGE = Stadium
      )
    
    # clean up columns
    c_subjects_infos[[name]]$LC_TYPE <- str_replace(c_subjects_infos[[name]]$LC_TYPE, "Plattenepithel", "Pleca")
  }
  
}

c_subjects_info <- bind_rows(c_subjects_infos)

# format CANCER_THERAPY_AT_V1 column
c_subjects_info$CANCER_THERAPY_AT_V1 <- c_subjects_info$CANCER_THERAPY_AT_V1 %>%
  # 1. remove everything after "nach" or "-Erhaltung"
  str_replace_all("-Erhaltung.*", "") %>%
  str_replace_all(" - Erhaltung.*", "") %>%
  str_replace_all("nach.*", "") %>%
  # 2. standardize drug names and capitalize
  str_replace_all("Pemprolizumab", "Pembrolizumab") %>%
  str_replace_all("Premetrexed", "Pemetrexed") %>%
  str_replace_all("nab-Paclitaxel", "Paclitaxel") %>%
  str_replace_all("  /", "/") %>% #remove double space
  str_replace_all("Radiochemotherapie","Radiotherapy") %>%
  str_replace_all("- ","/") %>%
  str_replace_all(" / ","/") %>%
  str_replace_all(" , ", ",") %>%
  str_replace_all(", ", ",") %>%
  str_replace_all(" ,", ",") %>%
  str_replace_all(" /", "/") %>%
  str_replace_all("/ ", "/") %>%
  str_replace_all(" ","") %>%
  # 3. replace comma and whitespace with only comma.
  str_replace_all(",\\s+",",") %>%
  # 4. split combinations, sort, and re-combine
  str_split("/") %>%
  lapply(sort) %>%
  lapply(paste, collapse = ",") %>%
  unlist() %>%
  str_split(",") %>%
  lapply(sort) %>%
  lapply(paste, collapse = ",") %>%
  unlist() %>%
  # 5. replace empty strings with "NA"
  str_replace_all("^$", NA_character_) %>%
  #6. capitalize the first letter of each word.
  str_to_title()


# trim whitespace from LC_STAGE
c_subjects_info$LC_STAGE <- str_trim(c_subjects_info$LC_STAGE)

# Remove ALL internal spaces
c_subjects_info$LC_STAGE <- str_replace_all(c_subjects_info$LC_STAGE, " ", "")

# remove rows where STUDY_ID is empty or NA
c_subjects_info <- c_subjects_info %>%
  filter(!is.na(STUDY_ID) & STUDY_ID != "")

# unify date formats
c_subjects_info$LC_THERAPY_START_DATE <- sapply(c_subjects_info$LC_THERAPY_START_DATE, clean_date)

# add c_subjects_info to metadata
metadata <- dplyr::left_join(metadata, c_subjects_info[c("STUDY_ID", "CANCER_THERAPY_AT_V1", "LC_TYPE", "LC_STAGE", "LC_THERAPY_START_DATE", "HEMOGLOBIN_LEVEL")], by = "STUDY_ID") %>%
  # combine HEMOGLOBIN_LEVEL columns
  mutate(HEMOGLOBIN_LEVEL = coalesce(HEMOGLOBIN_LEVEL.x, HEMOGLOBIN_LEVEL.y)) %>%
  # combine CANCER_THERAPY_AT_V1 columns
  mutate(CANCER_THERAPY_AT_V1 = coalesce(CANCER_THERAPY_AT_V1.x, CANCER_THERAPY_AT_V1.y)) %>%
  dplyr::select(-CANCER_THERAPY_AT_V1.x, -CANCER_THERAPY_AT_V1.y, -HEMOGLOBIN_LEVEL.x, -HEMOGLOBIN_LEVEL.y)  # remove the old HEMOGLOBIN_LEVEL and CANCER_THERAPY_AT_V1 columns

# with the CANCER_THERAPY_AT_V1 column being completed with B and C subject info, split the treatments into binary indicator variables ("Y" or "N")
# 1. Split the drug strings
split_drugs <- strsplit(tolower(as.character(metadata$CANCER_THERAPY_AT_V1)), ",")

# 2. Find unique drugs
all_drugs <- unlist(split_drugs)
unique_drugs <- unique(trimws(all_drugs)) # trimws removes whitespace

# Remove "w&w" and NA from unique_drugs
unique_drugs <- unique_drugs[!unique_drugs %in% c("w&w", NA)]

# 3. Create binary columns
for (drug in unique_drugs) {
  col_name <- paste0(toupper(drug), "_AT_V1") # Create column name, replace spaces with underscores
  metadata[, col_name] <- ifelse(grepl(drug, tolower(metadata$CANCER_THERAPY_AT_V1), fixed = TRUE), "Y", "N")
}

# make sure to name "DATO-DX_AT_V1" to "DATO_DX_AT_V1" (using base R)
colnames(metadata)[colnames(metadata) == "DATO-DX_AT_V1"] <- "DATO_DX_AT_V1"

# create binary column CHEMOTHERAPY_AT_V1 from cancer treatment info just generated
metadata <- metadata %>%
  mutate(
    CHEMOTHERAPY_AT_V1 = case_when(
      PEMETREXED_AT_V1 == "Y" ~ "Y",
      PACLITAXEL_AT_V1 == "Y" ~ "Y",
      CARBOPLATIN_AT_V1 == "Y" ~ "Y",
      TRUE ~ "N"  # for everyone else, leave it at N
    )
  )

# create binary column PD1_INHIBITOR_AT_V1 from cancer treatment info just generated
metadata <- metadata %>%
  mutate(
    PD1_INHIBITOR_AT_V1 = case_when(
      PEMBROLIZUMAB_AT_V1 == "Y" ~ "Y",
      CEMIPLIMAB_AT_V1 == "Y" ~ "Y",
      TRUE ~ "N"  # for everyone else, leave it at N
    )
  )


# calculate MONTHS_ON_LC_THERAPY from LC_THERAPY_START_DATE and VISIT_DATE
metadata$MONTHS_ON_LC_THERAPY <- mapply(calculate_month_diff, metadata$LC_THERAPY_START_DATE, metadata$VISIT_DATE)

# clean up hemoglobin levels column and make it in english decimal format
metadata$HEMOGLOBIN_LEVEL <- sapply(metadata$HEMOGLOBIN_LEVEL, function(x) {
  x <- gsub(",", ".", x)  # Replace commas with dots
  x <- gsub("\\s.*", "", x) # Remove everything after the first space
  return(x)
})

# add anamnesis questionnaire info ----------------------------------------
anamnesis_info_sys01 <- read.csv("data/meta/Anamnesis Questionnaire SYS01 - Season 2223 MODIFIED.csv", sep = ";", fileEncoding = "ISO-8859-1") # the encoding is necessary to handle the "µ" char
anamnesis_info_sys03 <- read.csv("data/meta/SYS03_AnamnesisQuestionnaire-Season2324.csv", sep = ";")

# make sure to remove leading and trailing whitespace from every character column and NA rows
anamnesis_info_sys01 <- anamnesis_info_sys01 %>%
  mutate(across(where(is.character), trimws)) %>%
  filter(!is.na(Study.ID) & Study.ID != "")

# also remove columns containing just NA values in SYS03 table
anamnesis_info_sys03 <- anamnesis_info_sys03 %>%
  mutate(across(where(is.character), trimws)) %>%
  filter(!is.na(Study.ID) & Study.ID != "") %>%
  dplyr::select(where(~ !all(is.na(.))))

# rename the columns to match preferred format COLUMN_CASE
rename_dict_sys01 <- c(
  "YEAR_OF_BIRTH" = "Date_of_Birth",
  "SEX" = "Gender",
  "STUDY_ID" = "Subject",
  "REGULAR_MEDICATION" = "Regular_medication",
  #"MEDICATION" = "Drug_classes_used_for_regular_medication",
  "INFLUENZA_VACCINATIONS_LAST_5_YEARS" = "No_of_Influenca_vaccinations_since_17.18",
  "COVID19_VACCINATIONS" = "No_of_COVID.19_vaccinations",
  "HISTORY_OF_CANCER" = "History_of_cancer"
)
rename_dict_sys03 <- c(
  "YEAR_OF_BIRTH" = "Year_of_Birth",
  "SEX" = "Gender",
  "STUDY_ID" = "Subject",
  "REGULAR_MEDICATION" = "Regular_medication",
  #"MEDICATION" = "Drug_classes_used_for_regular_medication",
  "INFLUENZA_VACCINATIONS_LAST_5_YEARS" = "How_many_since_2018.19",
  "COVID19_VACCINATIONS" = "No_of_COVID19_vaccinations",
  "HISTORY_OF_CANCER" = "History_of_cancer",
  "HEIGHT" = "Height..cm.",
  "WEIGHT" = "Weight.kg."
)

anamnesis_info_sys01 <- anamnesis_info_sys01 %>% dplyr::rename(all_of(rename_dict_sys01))
anamnesis_info_sys03 <- anamnesis_info_sys03 %>% dplyr::rename(all_of(rename_dict_sys03))

# convert all occurences of "/" to value 0 in INFLUENZA_VACCINATIONS_LAST_5_YEARS
# treat INFLUENZA_VACCINATIONS_LAST_5_YEARS column as factor for now (beware of how many levels are introduced)
anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS <- gsub("/", "0", as.character(anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS))
anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS <- as.factor(anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS)

anamnesis_info_sys03$INFLUENZA_VACCINATIONS_LAST_5_YEARS <- gsub("/", "0", as.character(anamnesis_info_sys03$INFLUENZA_VACCINATIONS_LAST_5_YEARS))
anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS <- as.factor(anamnesis_info_sys01$INFLUENZA_VACCINATIONS_LAST_5_YEARS)

# bind rows together
anamnesis_info <- bind_rows(anamnesis_info_sys01, anamnesis_info_sys03)

# join with metadata
metadata <- dplyr::left_join(metadata,
                             anamnesis_info[c("STUDY_ID", "YEAR_OF_BIRTH", "SEX", "HEIGHT", "WEIGHT", "INFLUENZA_VACCINATIONS_LAST_5_YEARS", "COVID19_VACCINATIONS", "HISTORY_OF_CANCER", "REGULAR_MEDICATION")],
                             by = "STUDY_ID")

# set the one subject where INFLUENZA_VACCINATIONS_LAST_5_YEARS was stated as 3,4 to 3 for now
# perhaps exclude it in analyses involving this factor or put it to NA
metadata$INFLUENZA_VACCINATIONS_LAST_5_YEARS <- gsub("3,4", "3", as.character(metadata$INFLUENZA_VACCINATIONS_LAST_5_YEARS))

# set the one subject where COVID19_VACCINATIONS was stated as >4 to 5 for now
# because not a single other subject has received 6 or more vaccinations, we can assume that this subject received 5 for now
# perhaps exclude it in analyses involving this factor or put it to NA
metadata$COVID19_VACCINATIONS <- gsub(">4", "5", as.character(metadata$COVID19_VACCINATIONS))


# convert YEAR_OF_BIRTH character vector into AGE covariat,
# format German decimal to English decimal in WEIGHT column,
# substitute "/" with NA in HISTORY_OF_CANCER column
metadata <- metadata %>%
  mutate(AGE = if_else(STUDY == "SYS01", 
                       2022 - as.numeric(YEAR_OF_BIRTH), # the study SYS01 took place in 22/23
                       2023 - as.numeric(YEAR_OF_BIRTH)), # the study SYS03 took place in 23/24
         WEIGHT = gsub(",", ".", WEIGHT),
         SEX = gsub("w", "f", as.character(SEX)),
         HISTORY_OF_CANCER = gsub("/", NA_character_, as.character(HISTORY_OF_CANCER)))

# create binary column AGE_BELOW_60
metadata <- metadata %>%
  mutate(AGE_BELOW_60 = ifelse(AGE < 60, "Y", "N"))


# add results from other experiments --------------------------------------
# add antibody data
antibody_info_sys01 <- read.csv("data/antibody/SYS01_MSD-ELISA-230509.csv", sep = ";") %>% filter(!row_number() == 1) # remove second row containing the days
antibody_info_sys03 <- read.csv("data/antibody/AUC of SYS03_all.csv", sep = ",") #%>% filter(!row_number() == 1) # remove second row containing the days
antibody_info_sys03$X <- trimws(antibody_info_sys03$X)

# SYS01
# convert to long format (one value for each study_id and visit combination)
antibody_info_sys01 <- antibody_info_sys01 %>%
  pivot_longer(cols = -Study.ID, names_to = "VISIT", values_to = "IgG") %>%
  mutate(
    VISIT = gsub("Visite\\.", "V", VISIT), # convert Visite.1 to V1
    SAMPLE_ID = paste0(Study.ID, "_", VISIT) # match metadata format
  )

# SYS03

# extract the row where X is "Total Area"
area_data <- antibody_info_sys03 %>% filter(X == "Total Area")

# remove the "X" column as it's no longer needed
area_data <- area_data %>% dplyr::select(-X)

# get the original column names (excluding X)
original_names <- colnames(area_data)

# function to transform column names, necessary to match SAMPLE_ID in metadata later on
convert_name <- function(name) {
  parts <- unlist(strsplit(name, "[.]")) # split at "."
  if (length(parts) == 2) {
    sample_id <- substr(parts[1], 2, nchar(parts[1])) # remove first letter
    paste0("SYS03_", substr(parts[1], 1, 1), "_", sample_id, "_", parts[2])
  } else {
    name # if it doesn't match expected pattern, leave it unchanged
  }
}

# rename the columns
new_names <- sapply(original_names, convert_name)
colnames(area_data) <- new_names

# convert to long format
antibody_info_sys03 <- area_data %>%
  pivot_longer(cols = everything(), names_to = "SAMPLE_ID", values_to = "IgG")

# convert datatype from character and double to integer
antibody_info_sys01$IgG <- as.integer(antibody_info_sys01$IgG)
antibody_info_sys03$IgG <- as.integer(antibody_info_sys03$IgG)

# bind rows
antibody_info <- bind_rows(antibody_info_sys01, antibody_info_sys03)

# join with metadata
metadata <- dplyr::left_join(metadata, antibody_info[c("SAMPLE_ID", "IgG")], by = join_by(SAMPLE_ID))

# create columns IgG_V3, IgG_V1_to_V3_l2fc, IgG_V4 and IgG_V1_to_V4_l2fc
metadata <- metadata %>%
  group_by(STUDY_ID) %>%
  dplyr::mutate(
    IgG_V1 = IgG[VISIT == "V1"][1], # get IgG value where VISIT is "V1"
    IgG_V3 = IgG[VISIT == "V3"][1], # get IgG value where VISIT is "V3"
    IgG_V1_to_V3_l2fc = log2(IgG[VISIT == "V3"][1]/IgG[VISIT == "V1"][1]), # compute a l2fc for V3
    IgG_V4 = IgG[VISIT == "V4"][1], # get IgG value where VISIT is "V4",
    IgG_V1_to_V4_l2fc = log2(IgG[VISIT == "V4"][1]/IgG[VISIT == "V1"][1]), # compute a l2fc for V4
  ) %>%
  ungroup()


# add titer response index TRI_ELISA
tri_info <- read.csv("data/antibody/MSD_SYS01+SYS03_TRI-ELISA.csv", sep = ";")

# create SAMPLE_ID column and replace comma with decimal dot in tri_info
tri_info <- tri_info %>%
  mutate("SAMPLE_ID" = paste0(Study.ID, "_", Visit),
         "TRI_ELISA" = as.numeric( gsub(",", ".", TRI_ELISA)))

# join with metadata
metadata <- dplyr::left_join(metadata, tri_info[c("SAMPLE_ID", "TRI_ELISA")], by = join_by(SAMPLE_ID))

# create columns TRI_ELISA_V3, TRI_ELISA_V4 and TRI_ELISA_V5
metadata <- metadata %>%
  group_by(STUDY_ID) %>%
  dplyr::mutate(
    TRI_ELISA_V3 = TRI_ELISA[VISIT == "V3"][1], # get TRI_ELISA value where VISIT is "V3"
    TRI_ELISA_V4 = TRI_ELISA[VISIT == "V4"][1], # get TRI_ELISA value where VISIT is "V4"
    TRI_ELISA_V5 = TRI_ELISA[VISIT == "V4"][1], # get TRI_ELISA value where VISIT is "V5"
  ) %>%
  ungroup()

# add technical details ---------------------------------------------------

# lost due to processing errors
lost_list <- c(
  # SYS01
  "S11664Nr184",
  "S11664Nr205",
  "S11664Nr192",
  "S11664Nr203",
  "S11664Nr188",
  "S11664Nr197",
  "S11664Nr172",
  "S11664Nr178",
  "S11664Nr199",
  "S11664Nr183",
  "S11664Nr207",
  "S11664Nr180",
  "S11664Nr201",
  "S11664Nr194",
  "S11664Nr189",
  "S11664Nr182",
  "S11664Nr198",
  "S11664Nr186",
  "S11664Nr190",
  "S11664Nr176",
  "S11664Nr195",
  "S11664Nr187",
  "S11664Nr173",
  "S11664Nr181",
  "S11664Nr175",
  "S11664Nr200",
  "S11664Nr193",
  "S11664Nr185",
  "S11664Nr206",
  "S11664Nr174",
  "S11664Nr202",
  "S11664Nr177",
  "S11664Nr179",
  "S11664Nr204",
  "S11664Nr196",
  "S11664Nr191",
  # SYS03
  "S11947Nr45",
  "S11947Nr231",
  "S11947Nr239",
  "S11947Nr320",
  "S11947Nr509"
)

# V1 missing - could also be done in code via comparison, or check which have a complete sample row
v1_missing_list <- c(
  # SYS03
  "SYS03_B_104"
)

# create RNA_SAMPLE_VISITS_COMPLETE column
metadata <- metadata %>%
  group_by(STUDY_ID) %>% # Group the data by STUDY_ID
  mutate(
    RNA_SAMPLE_VISITS_COMPLETE = ifelse(
      all(VISITS %in% VISIT), # check if all visits are present
      "Y",
      "N"
    )
  ) %>%
  ungroup()

# RNA input amount record for RNA_SAMPLE_INPUT column
# by default, 100 ng of RNA is used as input for sequencing
default_rna_input = 100.0
# cases for which less material was available at maximum input volume are documented below
rna_input_dict <- c(
  # SYS01
  "S11664Nr101" = 34.0,
  "S11664Nr103" = 93.0,
  # SYS03
  "S11947Nr52" = 64.0,
  "S11947Nr136" = 86.0,
  "S11947Nr235" = 82.0,
  "S11947Nr370" = 82.0,
  "S11947Nr323" = 68.0,
  "S11947Nr410" = 61.0,
  "S11947Nr418" = 69.0
)

# create RNA_SAMPLE_INPUT column
metadata$RNA_SAMPLE_INPUT <- ifelse(metadata$CEGAT_ID %in% names(rna_input_dict),
                                    rna_input_dict[metadata$CEGAT_ID], 
                                    default_rna_input)


# DNAse digestion performed record for RNA_SAMPLE_DNA_DIGESTION column
dna_digestion_list <- c(
  # SYS01
  "S11664Nr49",
  "S11664Nr50",
  "S11664Nr52",
  "S11664Nr54",
  "S11664Nr58",
  "S11664Nr61",
  "S11664Nr62",
  "S11664Nr67",
  "S11664Nr73",
  "S11664Nr77",
  "S11664Nr85",
  # SYS03
  "S11947Nr212",
  "S11947Nr300",
  "S11947Nr417",
  "S11947Nr433",
  "S11947Nr448"
)

# create RNA_SAMPLE_DNA_DIGESTION column
metadata$RNA_SAMPLE_DNA_DIGESTION <- ifelse(metadata$CEGAT_ID %in% dna_digestion_list,
                                            "Y", 
                                            "N")

# contamination record for RNA_SAMPLE_SUSPECTED_CONTAMINATION_WITH column
contamination_dict <- c(
  "S11947Nr54" = "S11947Nr55",
  "S11947Nr55" = "S11947Nr54"
)

# create RNA_SAMPLE_SUSPECTED_CONTAMINATION_WITH column
metadata$RNA_SAMPLE_SUSPECTED_CONTAMINATION_WITH <- ifelse(metadata$CEGAT_ID %in% names(contamination_dict),
                                                           contamination_dict[metadata$CEGAT_ID], 
                                                           NA_character_)

# create RNA_SAMPLE_RESEQUENCING
resequencing_failed_list <- c(
  # SYS03
  "S11947Nr322"
)

# create RNA_SAMPLE_RESEQUENCING column
metadata$RNA_SAMPLE_RESEQUENCING <- ifelse(metadata$CEGAT_ID %in% resequencing_failed_list,
                                           "N", 
                                           "Y")

# order columns
metadata <- metadata %>%
  dplyr::select(CEGAT_ID,
                INTERNAL_ID,
                STUDY_ID,
                SAMPLE_ID,
                STUDY,
                VISIT,
                DAY,
                NOMINAL_DAY,
                STUDY_GROUP,
                SUBJECT_GROUP_NUMBER,
                DIAGNOSIS,
                SEX,
                YEAR_OF_BIRTH,
                AGE,
                WEIGHT,
                HEIGHT,
                TRI_ELISA,
                TRI_ELISA_V3,
                TRI_ELISA_V4,
                TRI_ELISA_V5,
                IgG,
                IgG_V1,
                IgG_V3,
                IgG_V4, 
                everything())  # reorder specific columns, keep others

# create covariates -------------------------------

# convert everything into a factor except for those which will be treated as covariates
# note: RNA_SAMPLE_INPUT could be treated as a covriate, but because of the few levels it makes more sense to see it as a factor in plotting

covariates <- c("AGE", "HEIGHT", "DAY", "NOMINAL_DAY", "WEIGHT", "TRI_ELISA", "TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5", "IgG", "IgG_V1", "IgG_V3", "IgG_V1_to_V3_l2fc", "IgG_V4", "IgG_V1_to_V4_l2fc", "HEMOGLOBIN_LEVEL", "MONTHS_SINCE_HSCT", "MONTHS_ON_LENALIDOMIDE", "MONTHS_ON_LC_THERAPY", "RNA_SAMPLE_MONTHS_SINCE_VISIT")

metadata <- metadata %>%
  dplyr::mutate(across(-all_of(covariates), as.factor))

# center and scale covariates
for (covar in covariates) {
  # Convert to numeric (if possible)
  metadata[[covar]] <- as.numeric(metadata[[covar]])
  
  if (covar == "DAY") next # centering an scaling makes not much sense
  if (covar == "NOMINAL_DAY") next # centering an scaling makes not much sense
  if (covar == "TRI_ELISA") next # already centered and scaled
  if (covar == "TRI_ELISA_V3") next # already centered and scaled
  if (covar == "TRI_ELISA_V4") next # already centered and scaled
  if (covar == "TRI_ELISA_V5") next # already centered and scaled
  
  # Center and scale
  scaled_name <- paste0(covar, "_CENTERED_SCALED")
  metadata[[scaled_name]] <- scale(metadata[[covar]], center = TRUE, scale = TRUE)
}

# note: continue assembling MEDICATION once SYS03 medication is clear
# # split MEDICATION multi-label categorical variable into binary indicator variables
# # step 0: make sure whitespace is trimmed
# metadata$MEDICATION <- trimws(metadata$MEDICATION)
# # step 1: identify unique medication types
# medication_list <- unique(unlist(strsplit(paste(metadata$MEDICATION, collapse = ", "), ", "))) # creates a vector of all unique medication types
# # remove "none" and "NA" value, as it does not make much sense for a column
# medication_list <- setdiff(medication_list, c("none", "NA"))
# print(medication_list)
# # step 2: generate binary indicator variables
# # loop through each medication and generate binary column
# for (med in medication_list) {
#   col_name <- toupper(gsub(" ", "_", med)) # replaces spaces for underscores and capitalizes letters
#   metadata[[col_name]] <- as.factor(ifelse(grepl(med, metadata$MEDICATION), "Y", "N"))
# }
# remove MEDICATION column and instead grab the regular medication binary column
# metadata <- metadata %>% dplyr::select(-MEDICATION)

# finish metadata by setting CEGAT_IDs as rownames
rownames(metadata) <- metadata$CEGAT_ID

head(metadata)


# save metadata as .csv and .xlsx -----------------------------------------

# 1. saving sample based data

# loop through and add descriptions as attributes
for (col_name in names(description_dict)) {
  if (col_name %in% names(metadata)) {
    attr(metadata[[col_name]], "description") <- description_dict[[col_name]]
  }
}

# create a data frame for the Excel sheet of descriptions
column_descriptions <- data.frame(
  COLUMN = names(metadata),
  DESCRIPTION = sapply(names(metadata), function(col) {
    desc <- attr(metadata[[col]], "description")
    if (is.null(desc)) {
      return("") # return an empty string if there's no description
    } else {
      return(desc) # return the description if it exists
    }
  })
)

# create a new workbook
wb <- createWorkbook()

# add the metadata sheet
addWorksheet(wb, "sample_metadata")
# add the actual data to the sheet
writeData(wb, "sample_metadata", metadata, rowNames = FALSE, keepNA = TRUE, na.string = "NA")

# add the column descriptions sheet
addWorksheet(wb, "column_descriptions")
writeData(wb, "column_descriptions", column_descriptions, rowNames = FALSE)

# save the workbook
saveWorkbook(wb, "cleaned_output/sample_metadata.xlsx", overwrite = TRUE)

#save csv
write.csv2(metadata, file = "cleaned_output/sample_metadata.csv", row.names = FALSE, na = "NA")

# save R object for script 03
# note: is this still necessary, as we are now retrieving everything from dds?
saveRDS(metadata, file = "output/metadata.rds")

# 2. saving subject specific data

# columns to keep (subject-specific)
subject_cols <- c("STUDY_ID",
                  "STUDY",
                  "STUDY_GROUP",
                  "SUBJECT_GROUP_NUMBER",
                  "DIAGNOSIS",
                  "VISIT",
                  "VISIT_DATE",
                  "SEX",
                  "YEAR_OF_BIRTH",
                  "AGE",
                  "AGE_BELOW_60",
                  "WEIGHT",
                  "HEIGHT",
                  "TRI_ELISA_V3",
                  "TRI_ELISA_V4",
                  "TRI_ELISA_V5",
                  "IgG_V1",
                  "IgG_V3",
                  "IgG_V1_to_V3_l2fc",
                  "IgG_V4",
                  "IgG_V1_to_V4_l2fc",
                  "CANCER_THERAPY_AT_V1",
                  "CRIZOTINIB_AT_V1",
                  "ISATUXIMAB_AT_V1",
                  "LENALIDOMIDE_AT_V1",
                  "PEMBROLIZUMAB_AT_V1",
                  "PEMETREXED_AT_V1",
                  "DEXAMETHASONE_AT_V1",
                  "CARBOPLATIN_AT_V1",
                  "PACLITAXEL_AT_V1",
                  "CEMIPLIMAB_AT_V1",
                  "DATO_DX_AT_V1",
                  "HEMOGLOBIN_LEVEL",
                  "HSCT_DATE",
                  "MONTHS_SINCE_HSCT",
                  "LENALIDOMIDE_START_DATE",
                  "LENALIDOMIDE_END_DATE",
                  "MONTHS_ON_LENALIDOMIDE",
                  "LC_TYPE",
                  "LC_STAGE",
                  "LC_THERAPY_START_DATE",
                  "MONTHS_ON_LC_THERAPY",
                  "CHEMOTHERAPY_AT_V1",
                  "PD1_INHIBITOR_AT_V1",
                  "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
                  "COVID19_VACCINATIONS",
                  "HISTORY_OF_CANCER",
                  "REGULAR_MEDICATION",
                  "SUBJECT_LOCATION")

# group by STUDY_ID and keep subject-specific columns
subject_metadata <- metadata %>%
  group_by(STUDY_ID) %>%
  filter(VISIT == "V1") %>% # keep only rows where VISIT is "V1" to avoid duplicates
  slice(1) %>% # keep only the first row for each subject
  dplyr::select(all_of(subject_cols)) %>%
  ungroup()

# filter description_dict to only include subject_cols
subject_description_dict <- description_dict[names(description_dict) %in% subject_cols]

# loop through and add descriptions as attributes to subject_metadata
for (col_name in names(subject_description_dict)) {
  attr(subject_metadata[[col_name]], "description") <- subject_description_dict[[col_name]]
}

# create a data frame for the Excel sheet of descriptions
column_descriptions <- data.frame(
  COLUMN = names(subject_metadata),
  DESCRIPTION = sapply(names(subject_metadata), function(col) {
    desc <- attr(subject_metadata[[col]], "description")
    if (is.null(desc)) {
      return("") # return an empty string if there's no description
    } else {
      return(desc) # return the description if it exists
    }
  })
)

# create a new workbook
wb <- createWorkbook()

# add the subject metadata sheet
addWorksheet(wb, "subject_metadata")
writeData(wb, "subject_metadata", subject_metadata, rowNames = FALSE, keepNA = TRUE, na.string = "NA")

# add the column descriptions sheet
addWorksheet(wb, "column_descriptions")
writeData(wb, "column_descriptions", column_descriptions, rowNames = FALSE)

# ave the workbook
saveWorkbook(wb, "cleaned_output/subject_metadata.xlsx", overwrite = TRUE)

#save csv
write.csv2(subject_metadata, file = "cleaned_output/subject_metadata.csv", row.names = FALSE, na = "NA")

# save R object for script 03
saveRDS(subject_metadata, file = "output/subject_metadata.rds")

# example to see a description in R:
attr(metadata$AGE, "description")

message("Metadata cleaned up.")
