#' Process sugar data
#'
#' Function to process raw nonstructural carbohydrate data as read in using NSCprocessR::readRawNSCData ()
#' @param rawData Tibble containing the raw data, as read in by the readRawNSCData function.
#' @param cvLimitSample Limit for the coefficient of variation within a sample at which the sample if flagged for re-measurement.
#' @param cvLimitTube Limit for the coefficient of variation within a tube at which the sample is flagged for re-measurement.
#' @param forceIntercept Logical variable describing whether to force the calibration curve intercept through 0 or not.
#' @return The processed data, results of the anyalsis and a summary.
#' @export
processNSCs <- function (rawData, cvLimitSample = 0.25, cvLimitTube = 0.05, forceIntercept = FALSE) {

  # Calculate the sample weight
  rawData [['MassSample']] <- rawData [['MassOfTubeAndSample']] - rawData [['MassOfEmptyTube']]
  # What to do with negative mass of sample??? Maybe jjst set them to 0?

  # Create two new colums with the average absorbance at 490 and 525 nm
  absorbances490 <- cbind (rawData [['Absorbance490_1']], rawData [['Absorbance490_2']])
  absorbances525 <- cbind (rawData [['Absorbance525_1']], rawData [['Absorbance525_2']])
  rawData [['MeanAbsorbance490']] <- rowMeans (absorbances490)
  rawData [['MeanAbsorbance525']] <- rowMeans (absorbances525)

  # Calculate the within-sample coefficient of variation
  rawData [['SDAbsorbance490']] <- apply (absorbances490, 1, sd)
  rawData [['SDAbsorbance525']] <- apply (absorbances525, 1, sd)
  rawData [['CVAbsorbance490']] <- rawData [['SDAbsorbance490']] / rawData [['MeanAbsorbance490']]
  rawData [['CVAbsorbance525']] <- rawData [['SDAbsorbance525']] / rawData [['MeanAbsorbance525']]

  # Flag samples with high CV for re-measuring
  rawData [['HighCV']] <- 'N'
  rawData [['HighCV']] [rawData [["CVAbsorbance490"]] > cvLimitTube] <- 'Sugar'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube & rawData [['HighCV']] == 'N'] <- 'Starch'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube & rawData [['HighCV']] == 'Sugar'] <- 'Sugar and starch'

  # Extract blanks with with solution
  blanks <- rawData [rawData [['SampleID']] == 'B', ]

  # Extract blanks with only empty tubes
  tubeBlanks <- rawData [rawData [['SampleID']] == 'TB', ]

  # Extract Lab Control Standards (LCS)
  labControlStandards <- rawData [substr (rawData [['SampleID']], 1, 3) == 'LCS', ]

  # Extract reference values for calibration curve for each batch
  referenceValues <- rawData [substr (rawData [['SampleID']], 1, 3) == 'REF', ]

  # Return the processed data
  return (processedData)
}
