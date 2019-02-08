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

  # Get indices for samples with negative weight
  indices <- which (rawData [['MassSample']] < 0.0)

  # Set 'MassSample' to zero for blanks (B) and tube blanks (TB) that have negative weight
  rawData [['MassSample']] [rawData [['SampleID']] [indices] == 'B' | rawData [['SampleID']] [indices] == 'TB'] <- 0.0

  # Error for non-blanks with negativ weight.
  if (length (which (rawData [['MassSample']] < 0.0)) > 0) {
    if (which (rawData [['MassSample']] < 0.0)) stop ('Error: There is a non-blank sample with a negative weight.')
  }

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

  # Extract blanks with solution
  blanks <- rawData [rawData [['SampleID']] == 'B', ]

  # Extract blanks with only empty tubes
  tubeBlanks <- rawData [rawData [['SampleID']] == 'TB', ]

  # Extract Lab Control Standards (LCS)
  labControlStandards <- rawData [substr (rawData [['SampleID']], 1, 3) == 'LCS', ]

  # Compile a list of batches
  batches <- unique (rawData [['BatchID']])

  # Make and save a calibration curve for each batch
  for (batch in batches) {
    for (NSC in c ('sugar','starch')) {
      condition <- substr (rawData [['SampleID']], 1, 3) == 'REF' & rawData [['BatchID']] == batch
      referenceValues <- rawData [condition, ]
      concentrations <- substr (rawData [['SampleID']] [condition], 4, nchar (rawData [['SampleID']] [condition]))
      if (NSC == 'sugar') {
        concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      } else if (NSC == 'starch'){
        concentrations [concentrations == '100/200'] <- 200.0 # 200 for starrch
      }
      concentrations <- as.numeric (concentrations)

      # Get the slope and intercept
      if (forceIntercept) { # Get slope for intercepts forced through zero
        if (NSC == 'sugar') {
          fitSugar  <- lm (referenceValues [['MeanAbsorbance490']] ~ 0 + concentrations)
        } else if (NSC == 'starch') {
          fitStarch <- lm (referenceValues [['MeanAbsorbance525']] ~ 0 + concentrations)
        }
      } else { # Get intercept and slope
        if (NSC == 'sugar') {
          fitSugar  <- lm (referenceValues [['MeanAbsorbance490']] ~ concentrations)
        } else if (NSC == 'starch') {
          fitStarch <- lm (referenceValues [['MeanAbsorbance525']] ~ concentrations)
        }
      }

      # Add the slopes and intercepts to respective rows in the tibble
      if (forceIntercept) { # intercept forced through 0.0
        if (NSC == 'sugar') {
          rawData [['InterceptSugar']]  [rawData [['BatchID']] == batch] <- 0.0
          rawData [['SlopeSugar']]      [rawData [['BatchID']] == batch] <- fitSugar$coefficients
        } else if (NSC == 'starch') {
          rawData [['InterceptStarch']] [rawData [['BatchID']] == batch] <- 0.0
          rawData [['SlopeStarch']]     [rawData [['BatchID']] == batch] <- fitStarch$coefficients
        }
      } else { # variable intercept
        if (NSC == 'sugar') {
          rawData [['InterceptSugar']]  [rawData [['BatchID']] == batch] <- fitSugar$coefficients [1]
          rawData [['SlopeSugar']]      [rawData [['BatchID']] == batch] <- fitSugar$coefficients [2]
        } else if (NSC == 'starch') {
          rawData [['InterceptStarch']] [rawData [['BatchID']] == batch] <- fitStarch$coefficients [1]
          rawData [['SlopeStarch']]     [rawData [['BatchID']] == batch] <- fitStarch$coefficients [2]
        }
      }
    }
  }

  # Determine concentrations from absorbance values for sugar
  rawData [['ConcentrationSugarM']] <- (rawData [['MeanAbsorbance490']] - rawData [['InterceptSugar']]) /
                                        rawData [['SlopeSugar']]
  rawData [['ConcentrationSugarB']] <- (rawData [['Absorbance490_Blank']] - rawData [['InterceptSugar']]) /
                                        rawData [['SlopeSugar']]


  # Correct for background signal using no phenol
  rawData [['ConcentrationSugar']] <- (rawData [['ConcentrationSugarM']] - rawData [['ConcentrationSugarB']]) *
                                       rawData [['DilutionFactorSugar']]

  # Determine concentrations from absorbance values for starch
  rawData [['ConcentrationStarch']] <- (rawData [['MeanAbsorbance525']] - rawData [['InterceptStarch']]) /
                                        rawData [['SlopeStarch']]
  rawData [['ConcentrationStarch']] <- rawData [['ConcentrationStarch']] * rawData [['DilutionFactorStarch']]

  # Once converted to concentration what do I need to do???

  # Declare data as processed
  processedData <- rawData

  # Return the processed data
  return (processedData)
}
