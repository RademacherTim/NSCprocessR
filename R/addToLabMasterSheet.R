#' Create or append lab master sheet with the processed samples
#'
#' Function to create or add to the lab master sheet.
#' @param data Tibble with the data that is to be added to the master sheet. The data was created using the processNSC function.
#' @return
#' @import openxlsx
#' @import readxl
#' @import tidyverse
#' @export
addToLabMasterSheet <- function (data, fileDir, fileName, IDs = c ('SampleID')) {

  # get the number of ID coumns
  #--------------------------------------------------------------------------------------
  nIDs <- length (IDs)

  # Paste directory name and file name to get a full path
  #--------------------------------------------------------------------------------------
  fullPath <- paste (fileDir,'/',fileName, sep = '')

  # Check whether master sheet exists for lab and create it
  #--------------------------------------------------------------------------------------
  if (!file.exists (fullPath)) { # If master sheet does NOT exist, create it with the data
    openxlsx::write.xlsx (x         = data,
                          file      = fullPath,
                          sheetName = "Master",
                          col.names = TRUE,
                          row.names = FALSE,
                          append    = FALSE,
                          keepNA    = FALSE)
  } else { # If it does exist

    # Read the master sheet
    #------------------------------------------------------------------------------------
    types <- c (rep ('text', (nIDs+1)), rep ('date', 3), rep ('numeric', 19), 'text',
                rep ('numeric', 2), 'text', 'numeric','text', rep ('numeric', 6), 'text',
                rep ('numeric', 6), 'text')
    master <- readxl::read_excel (fullPath, col_types = types)

    # Join the two tibble
    #------------------------------------------------------------------------------------
    temp <- dplyr::full_join (x = master,
                              y = data,
                              by = c ("RCLabNumber", "SampleID", "Tissue", "BatchID",
                                      "SampleLocation", "DateOfSampleCollection",
                                      "DateOfSugarAnalysis", "DateOfStarchAnalysis",
                                      "MassOfEmptyTube", "MassOfTubeAndSample",
                                      "Absorbance490_1", "Absorbance490_2",
                                      "Absorbance490_Blank", "Absorbance525_1",
                                      "Absorbance525_2", "DilutionFactorSugar",
                                      "VolumeSugar", "DilutionFactorStarch",
                                      "VolumeStarch", "MassSample", "MeanAbsorbance490",
                                      "MeanAbsorbance525", "CorrectedMeanAbsorbance490",
                                      "SDAbsorbance490", "SDAbsorbance525",
                                      "CVAbsorbance490", "CVAbsorbance525", "HighCV",
                                      "SlopeSugar", "TBHigh",
                                      "CorrectedMeanAbsorbance525", "LowAbsorbance525",
                                      "SlopeStarch","ConcentrationSugar", "MassSugar",
                                      "ConcentrationStarch", "MassStarch","SRFHigh",
                                      "MeanStarchRecovery", "CorrectedMassStarch",
                                      "ConcentrationSugarMgG", "ConcentrationStarchMgG",
                                      "ConcentrationSugarPerDW",
                                      "ConcentrationStarchPerDW","LCSOakDeviation"))

    # Delete duplicate rows in temp, in case the data had already been added
    #------------------------------------------------------------------------------------
    temp <- distinct (temp, .keep_all = TRUE)

    # Overwrite the old master sheet with the new version
    #------------------------------------------------------------------------------------
    openxlsx::write.xlsx (x         = temp,
                          file      = fullPath,
                          sheetName = "Master",
                          col.names = TRUE,
                          row.names = FALSE,
                          append    = FALSE, # Do not append the already existing sheet
                          keepNA    = FALSE) # Empty cells for NAs
  }

  # Return zero if the function ran smoothly and data has been successfully appended to
  # the master sheet
  #--------------------------------------------------------------------------------------
  return (0)
}
# To-do list:
