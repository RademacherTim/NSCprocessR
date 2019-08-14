#' Process sugar data
#'
#' Function to process raw nonstructural carbohydrate data as read in using NSCprocessR::readRawNSCData ()
#' @param rawData Tibble containing the raw data, as read in by the readRawNSCData function.
#' @param cvLimitSample Limit for the coefficient of variation within a sample at which the sample if flagged for re-measurement.
#' @param cvLimitTube Limit for the coefficient of variation within a tube at which the sample is flagged for re-measurement.
#' @param forceIntercept Logical variable describing whether to force the calibration curve intercept through 0 or not.
#' @param minimumSampleMass Threshold weight in milligram below which samples are dropped because they are unreliable. The threshold is set to 5 mg by default.
#' @return The processed data, results of the anyalsis and a summary.
#' @import tibble
#' @export
processNSCs <- function (rawData,
                         cvLimitSample = 0.25,
                         cvLimitTube = 0.05, # Should this be 0.10 or 10% error, as suggested by Jim?
                         forceIntercept = FALSE,
                         maxStarchRecoveryFraction = 0.9,
                         LCS = 'Oak',
                         minimumSampleMass = 5.0) {

  # Load dependencies
  #--------------------------------------------------------------------------------------
  if (!existsFunction ('tibble')) library ('tidyverse')

  # Calculate the sample weight [g] and convert to [mg]
  #--------------------------------------------------------------------------------------
  rawData [['MassSample']] <- (rawData [['MassOfTubeAndSample']] -
                               rawData [['MassOfEmptyTube']]) * 1e3

  # check whether a sample weight is below 10 mg and drop it if it is.
  #--------------------------------------------------------------------------------------
  indicesToDrop <- which (rawData [['MassSample']]                 <  minimumSampleMass &
                          rawData [['SampleID']]                   != 'B'               &
                          rawData [['SampleID']]                   != 'TB'              &
                          substring (rawData [['SampleID']], 1, 3) != 'REF'             &
                          substring (rawData [['SampleID']], 1, 3) != 'LCS')
  rawData <- rawData [-indicesToDrop, ]
  if (length (indicesToDrop) >= 1) {
    warning (cat (sprintf ('Warning: %s samples are excluded from the analysis,',
                           length (indicesToDrop)),
                  sprintf ('because their sample mass is below the threshold of %s mg.',
                           minimumSampleMass)))
  }

  # Get indices for samples with negative weight
  #--------------------------------------------------------------------------------------
  indices <- which (rawData [['MassSample']] < 0.0)

  # Set 'MassSample' to zero for blanks (B) and tube blanks (TB) that have negative weight
  #--------------------------------------------------------------------------------------
  for (i in indices) {
    if (rawData [['SampleID']] [i] == 'B' |
        rawData [['SampleID']] [i] == 'TB') rawData [['MassSample']] [i] <- 0.0
  }

  # Error for non-blanks with negativ weight.
  #--------------------------------------------------------------------------------------
  if (length (which (rawData [['MassSample']] < 0.0)) > 0) {
    if (which (rawData [['MassSample']] < 0.0)) {
      stop ('Error: There is a non-blank sample with a negative weight.')
    }
  }

  # Create two new colums with the average absorbance at 490 and 525 nm
  #--------------------------------------------------------------------------------------
  absorbances490 <- cbind (rawData [['Absorbance490_1']], rawData [['Absorbance490_2']])
  absorbances525 <- cbind (rawData [['Absorbance525_1']], rawData [['Absorbance525_2']])
  rawData [['MeanAbsorbance490']] <- rowMeans (absorbances490, na.rm = T)
  rawData [['MeanAbsorbance525']] <- rowMeans (absorbances525, na.rm = T)
  rawData [['CorrectedMeanAbsorbance490']] <- rawData [['MeanAbsorbance490']] -
                                              rawData [['Absorbance490_Blank']] # TTR Some corrected mean absorbances are negative.

  # Calculate the within-sample coefficient of variation
  #--------------------------------------------------------------------------------------
  rawData [['SDAbsorbance490']] <- apply (absorbances490, 1, sd)
  rawData [['SDAbsorbance525']] <- apply (absorbances525, 1, sd)
  rawData [['CVAbsorbance490']] <- rawData [['SDAbsorbance490']] /
                                   rawData [['MeanAbsorbance490']]
  rawData [['CVAbsorbance490']] [rawData [['CVAbsorbance490']] < 0.0   |
                                 is.na  (rawData [['CVAborbance490']]) |
                                 is.nan (rawData [['CVABsorbance490']])] <- NA
  rawData [['CVAbsorbance525']] <- rawData [['SDAbsorbance525']] /
                                   rawData [['MeanAbsorbance525']]
  rawData [['CVAbsorbance525']] [rawData [['CVAbsorbance525']] < 0.0   |
                                   is.na  (rawData [['CVAborbance525']]) |
                                   is.nan (rawData [['CVABsorbance525']])] <- NA

  # Flag samples with high CV for re-measuring
  #--------------------------------------------------------------------------------------
  rawData [['HighCV']] <- 'N'
  rawData [['HighCV']] [rawData [["CVAbsorbance490"]] > cvLimitTube] <- 'Sugar'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube &
                        rawData [['HighCV']] == 'N'] <- 'Starch'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube &
                        rawData [['HighCV']] == 'Sugar'] <- 'Sugar and starch'

  # Compile a list of unique batches
  #--------------------------------------------------------------------------------------
  batches <- unique (rawData [['BatchID']])
  for (batch in batches) {
    dates <- unique (rawData [['DateOfSugarAnalysis']] [rawData [['BatchID']] == batch])
    if (batch == batches [1]) {
      extractionsSugar <- tibble::tibble (batch = batch, date = dates)
    } else {
      extractionsSugar <- tibble::add_row (extractionsSugar, batch = batch, date = dates)
    }
  }
  # Delete rows that do not have calibration curves
  #--------------------------------------------------------------------------------------
  if (sum (is.na (extractionsSugar [['date']])) > 0) {
    extractionsSugar <- extractionsSugar [-which (is.na (extractionsSugar [['date']])), ]
  }

  # Compile a list of unique batches and dates for sugar extractions
  #--------------------------------------------------------------------------------------
  for (batch in batches) {
    dates <- unique (rawData [['DateOfStarchAnalysis']] [rawData [['BatchID']] == batch])
    if (batch == batches [1]) {
      extractionsStarch <- tibble (batch = batch, date = dates)
    } else {
      extractionsStarch <- add_row (extractionsStarch, batch = batch, date = dates)
    }
  }

  # Delete rows that do not have calibration curves
  #--------------------------------------------------------------------------------------
  if (sum (is.na (extractionsStarch [['date']])) > 0) {
    extractionsStarch <- extractionsStarch [-which (is.na (extractionsStarch [['date']])), ]
  }

  # Make and save a sugar calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  for (extraction in 1:(dim (extractionsSugar) [1])) {

    # Get date of analysis and batch number
    #--------------------------------------------------------------------------------------
    analysisDate <- extractionsSugar [['date']]  [extraction]
    batch        <- extractionsSugar [['batch']] [extraction]

    # Get absorbances for reference values to create a calibration curve for each batch
    #----------------------------------------------------------------------------------
    refCondition <- substr (rawData [['SampleID']], 1, 3)     == 'REF' &
                            rawData [['BatchID']]             == batch &
                            rawData [['DateOfSugarAnalysis']] == analysisDate
    batchCondition <- rawData [['BatchID']]             == batch &
                      rawData [['DateOfSugarAnalysis']] == analysisDate
    referenceValues <- rawData [refCondition, ]

    # Get reference solution concentrations
    #----------------------------------------------------------------------------------
    concentrations <- substr (rawData [['SampleID']] [refCondition],
                              4,
                              nchar (rawData [['SampleID']] [refCondition]))

    # Set reference solution concentrations to 100 for sugars
    #----------------------------------------------------------------------------------
    concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
    concentrations <- as.numeric (concentrations)

    # Get the slope and intercept
    #----------------------------------------------------------------------------------
    if (forceIntercept) { # Get slope for intercepts forced through zero
      # sugar <- slope * absorbance + intercept
      fitSugar  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance490']])
    } else { # Get intercept and slope
      fitSugar  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance490']])
    }

    # Add the slopes and intercepts to respective rows in the tibble
    #----------------------------------------------------------------------------------
    if (forceIntercept) { # intercept forced through 0.0
      rawData [['InterceptSugar']]  [batchCondition] <- 0.0
      rawData [['SlopeSugar']]      [batchCondition] <- fitSugar$coefficients
    } else { # variable intercept
      rawData [['InterceptSugar']]  [batchCondition] <- fitSugar$coefficients [1]
      rawData [['SlopeSugar']]      [batchCondition] <- fitSugar$coefficients [2]
    }
  } # End of extraction loop

  # Make and save a starch calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  for (extraction in 1:(dim (extractionsStarch) [1])) {

    # Get date of analysis and batch number
    #--------------------------------------------------------------------------------------
    analysisDate <- extractionsStarch [['date']]  [extraction]
    batch        <- extractionsStarch [['batch']] [extraction]

    # Get absorbances for reference values to create a calibration curve for each batch
    #----------------------------------------------------------------------------------
    refCondition <- substr (rawData [['SampleID']], 1, 3) == 'REF'        &
                    rawData [['BatchID']]                 == batch        &
                    rawData [['DateOfStarchAnalysis']]    == analysisDate &
                    !is.na (rawData [['SampleID']])
    batchCondition <- rawData [['BatchID']]              == batch &
                      rawData [['DateOfStarchAnalysis']] == analysisDate
    referenceValues <- rawData [refCondition, ]

    # Get reference solution concentrations
    #----------------------------------------------------------------------------------
    concentrations <- substr (rawData [['SampleID']] [refCondition],
                              4,
                              nchar (rawData [['SampleID']] [refCondition]))

    # Set reference solution concentrations to 200 for starch
    #----------------------------------------------------------------------------------
    concentrations [concentrations == '100/200'] <- 200.0 # 200 for starch
    concentrations <- as.numeric (concentrations)

    # Get the batch's mean absorbance at 525nm for tube blanks
    #--------------------------------------------------------------------------------
    if (sum (rawData [['SampleID']] == 'TB' & batchCondition, na.rm = T) > 0.0) {
      batchTBAbsorbance <- mean (rawData [['MeanAbsorbance525']] [rawData [['SampleID']] == 'TB' &
                                                                  batchCondition],
                                 na.rm = T)
    } else { # In case there are no tube blanks we just assume no interference form tubes
      batchTBAbsorbance <- 0.0
    }


    # Determine correction factor from TB, unless they are larger than the sample
    # absorbances at 525nm. If the TB is larger than sample absorbance at 525nm
    #--------------------------------------------------------------------------------
    batchCorrection <- min (batchTBAbsorbance,
                            rawData [['MeanAbsorbance525']] [batchCondition])

    # Flag for higher TB than MeanAbsorbance525
    #--------------------------------------------------------------------------------
    if (batchTBAbsorbance != batchTBAbsorbance) {
      rawData [['TBHigh']] [batchCondition] <- 'Y'
    } else {
      rawData [['TBHigh']] [batchCondition] <- 'N'
    }

    # Correct mean absorbance values at 525nm
    #--------------------------------------------------------------------------------
    rawData [['CorrectedMeanAbsorbance525']] [batchCondition] <-
      rawData [['MeanAbsorbance525']] [batchCondition] - batchCorrection

    # Set flag for low absorbance in samples
    #--------------------------------------------------------------------------------
    rawData [['LowAbsorbance525']] [batchCondition] <- 'N'
    rawData [['LowAbsorbance525']] [rawData [['MeanAbsorbance525']] < batchTBAbsorbance &
                                    batchCondition] <- 'Y'

    # Correct reference values for tube blank absorbance at 525nm
    #--------------------------------------------------------------------------------
    referenceValues [['CorrectedMeanAbsorbance525']] <-
      referenceValues [['MeanAbsorbance525']] - batchCorrection

    # Get the slope and intercept (startch = slope * absorbance + intercept)
    #----------------------------------------------------------------------------------
    if (forceIntercept) { # Get slope for intercepts forced through zero
      fitStarch  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance525']])
    } else { # Get intercept and slope
      fitStarch  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance525']])
    }

    # Drop 250 from calibration curve if R2 is below 0.9
    #----------------------------------------------------------------------------------
    if (summary (fitStarch)$r.squared < 0.9) {
      indexToDrop <- which (concentrations == 250.0)
      concentrations <- concentrations [-indexToDrop]
      referenceValues <- referenceValues [['CorrectedMeanAbsorbance525']] [-indexToDrop]
      if (forceIntercept) { # Get slope for intercepts forced through zero
        fitStarch  <- lm (concentrations ~ 0 + referenceValues)
      } else { # Get intercept and slope
        fitStarch  <- lm (concentrations ~ referenceValues)
      }
    }

    # check whether calibration curve is sufficiently precise
    #----------------------------------------------------------------------------------
    if (summary (fitStarch)$r.squared < 0.9) {
      warning (paste ('Warning: the calibration curve for batch ',batch,' on the ',
                      analysisDate,' has an R2 lower than 0.9.', sep = ''))
    }

    # Add the slopes and intercepts to respective rows in the tibble
    #----------------------------------------------------------------------------------
    if (forceIntercept) { # intercept forced through 0.0
      rawData [['InterceptStarch']] [batchCondition] <- 0.0
      rawData [['SlopeStarch']]     [batchCondition] <- fitStarch$coefficients
    } else { # variable intercept
      rawData [['InterceptStarch']] [batchCondition] <- fitStarch$coefficients [1]
      rawData [['SlopeStarch']]     [batchCondition] <- fitStarch$coefficients [2]
    }

  } # End of starch extractions loop

  # Determine concentrations from absorbance values for sugar # Call it a concentrations [mg ml-1]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugar']] <- rawData [['CorrectedMeanAbsorbance490']] *
                                      rawData [['SlopeSugar']] +
                                      rawData [['InterceptSugar']]
  rawData [['ConcentrationSugar']] [rawData [['ConcentrationSugar']] < 0.0   |
                                    is.na  (rawData [['ConcentrationSugar']])|
                                    is.nan (rawData [['ConcentrationSugar']])] <- NA

  # Correct sugar amount [mg ml-1] for background signal using no phenol
  #--------------------------------------------------------------------------------------
  rawData [['MassSugar']] <- rawData [['ConcentrationSugar']] *
                             rawData [['DilutionFactorSugar']] *
                             rawData [['VolumeSugar']]

  # Determine concentration of strach [mg ml-1] from absorbance
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationStarch']] <- rawData [['CorrectedMeanAbsorbance525']] *
                                       rawData [['SlopeStarch']] +
                                       rawData [['InterceptStarch']]
  rawData [['ConcentrationStarch']] [rawData [['ConcentrationStarch']] < 0.0   |
                                      is.na  (rawData [['ConcentrationStarch']])|
                                      is.nan (rawData [['ConcentrationStarch']])] <- NA
  rawData [['MassStarch']] <- rawData [['ConcentrationStarch']] *
                              rawData [['DilutionFactorStarch']] *
                              rawData [['VolumeStarch']]

  # Calculate starch recovery for each combination of batch and date
  #--------------------------------------------------------------------------------------
  for (extraction in 1:(dim (extractionsStarch) [1])) {

    # Get date of analysis and batch number
    #--------------------------------------------------------------------------------------
    analysisDate <- extractionsStarch [['date']]  [extraction]
    batch        <- extractionsStarch [['batch']] [extraction]

    # Get absorbances for potato starch for each combination of batch and analysisDate
    #----------------------------------------------------------------------------------
    refCondition <- substr (rawData [['SampleID']], 1, 10) == 'LCS Potato' &
                    rawData [['BatchID']]                  == batch        &
                    rawData [['DateOfStarchAnalysis']]     == analysisDate &
                    !is.na (rawData [['SampleID']])
    batchCondition <- rawData [['BatchID']]              == batch        &
                      rawData [['DateOfStarchAnalysis']] == analysisDate

    # Set flag for unrealistically high starch recovery fraction
    #----------------------------------------------------------------------------------
    rawData [['SRFHigh']] [batchCondition] <- 'N'

    # Check whether potato standard was run at all, otherwise use maxStarchRecoveryFraction
    #----------------------------------------------------------------------------------
    if (sum (refCondition, na.rm = T) == 0) {
      rawData [['MeanStarchRecovery']] [batchCondition] <- maxStarchRecoveryFraction
    } else {
      LCSPotato <- rawData [refCondition, ]
      meanPotatoMassRecovered <- mean (LCSPotato [['MassStarch']])
      if (meanPotatoMassRecovered >
          mean (LCSPotato [['MassSample']]) * maxStarchRecoveryFraction * 1000.0) { # Maybe I should just drop to one high value?
        meanPotatoMass <- mean (LCSPotato [['MassSample']])
        rawData [['SRFHigh']] [batchCondition] <- 'Y'
      } else {
        meanPotatoMass <- mean (LCSPotato [['MassSample']]) * maxStarchRecoveryFraction
      }
      meanPotatoMass <- meanPotatoMass * 1000.0 # Convert from g to mg
      meanRecoveryPer <- meanPotatoMassRecovered / meanPotatoMass * 100.0
      if (meanRecoveryPer > 100.0) rawData [['SRFHigh']] [batchCondition] <- 'Y'
      meanRecoveryPer <- min (meanRecoveryPer, 100.0) # Hack to avoid recovery rate above 100%.

      # Set batch's mean starch recovery rate
      #----------------------------------------------------------------------------------
      rawData [['MeanStarchRecovery']] [batchCondition] <- meanRecoveryPer
    }
  }

  # Correct all samples for the mean starch recovery percentage of each batch
  #------------------------------------------------------------------------------------
  rawData [['CorrectedMassStarch']] <- rawData [['MassStarch']] /
                                      (rawData [['MeanStarchRecovery']] / 100.0)

  # Convert mass to concentrations [mg g-1]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugarMgG']]  <- rawData [['MassSugar']] /
                                          rawData [['MassSample']]
  rawData [['ConcentrationSugarMgG']] [rawData [['ConcentrationSugarMgG']] < 0.0   |
                                       is.na (rawData [['ConcentrationSugarMgG']]) |
                                       is.nan (rawData [['ConcentrationSugarMgG']])] <- NA
  rawData [['ConcentrationStarchMgG']] <- rawData [['CorrectedMassStarch']] /
                                          rawData [['MassSample']]

  # Convert concentrations from [mg g-1] to [% dry weight]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugarPerDW']]  <- rawData [['ConcentrationSugarMgG']]  / 10.0
  rawData [['ConcentrationStarchPerDW']] <- rawData [['ConcentrationStarchMgG']] / 10.0

  # Check whether Lab Control Standard for Oak wood is high
  #----------------------------------------------------------------------------------
  if (LCS == 'Oak') {
    for (extraction in 1:(dim (extractionsSugar) [1])) {

      # Get date of analysis and batch number
      #--------------------------------------------------------------------------------------
      analysisDate <- extractionsStarch [['date']]  [extraction]
      batch        <- extractionsStarch [['batch']] [extraction]

      # Get absorbances for potato starch for each combination of batch and analysisDate
      #----------------------------------------------------------------------------------
      refCondition <- substr (rawData [['SampleID']], 1, 7) == 'LCS Oak'    &
                      rawData [['BatchID']]                 == batch        &
                      rawData [['DateOfStarchAnalysis']]    == analysisDate &
                      !is.na (rawData [['SampleID']])

      # check whether thebatch has a LCS Oak
      #----------------------------------------------------------------------------------
      if (sum (refCondition) == 0) {
        warning (paste ('Warning: There is no LCS Oak for batch ',batch,
                        ' analysed on the ',analysisDate,'.', sep = ''))
      } else { # Check the LCS Oak is low or high
        # get LCS Oak standard and compare it against threshold
        #--------------------------------------------------------------------------------
        LCSOakSugar  <- rawData [['ConcentrationSugarPerDW']]  [refCondition]
        LCSOakStarch <- rawData [['ConcentrationStarchPerDW']] [refCondition]

        # Flag for high or low oak lab standard
        #--------------------------------------------------------------------------------
        rawData [['LCSOakDeviation']] <- 'N'

        # oak sugar was 3.59+-0.38 % DW at Harvard and 3.28+-0.51 at NAU thus far
        #--------------------------------------------------------------------------------
        if (mean (LCSOakSugar, na.rm = T) < 3.59-0.38 &
            mean (LCSOakSugar, na.rm = T) > 3.59+0.38) {
          rawData [['LCSOakDeviation']] <- 'sugar'
        }
        # oak starch was 2.45+-0.21 % DW at Harvard and 2.91+-0.33 at NAU thus far
        #--------------------------------------------------------------------------------
        if (mean  (LCSOakStarch, na.rm = T) < 2.45-0.21 &
            mean (LCSOakStarch, na.rm = T) > 2.45+0.21) {
          if (mean (LCSOakSugar, na.rm = T) < 3.59-0.38 &
              mean (LCSOakSugar, na.rm = T) > 3.59+0.38) {
            rawData [['LCSOakDeviation']] <- 'sugar & starch'
          } else {
            rawData [['LCSOakDeviation']] <- 'starch'
          }
        }
      } # End sum (refCondition) == 0 condition
    } # End extraction loop
  } # End LCS == Oak condition

  # Declare data as processed
  #--------------------------------------------------------------------------------------
  processedData <- rawData

  # Return the processed data
  #--------------------------------------------------------------------------------------
  return (processedData)
}
# To-do-list:
# TTR Test that TBHigh flag works
# - meanStarchRecovery can be more than 100 percent, which does not make sense. As a hack I limited it to 100% as max!!!
