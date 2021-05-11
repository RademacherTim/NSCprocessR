#' Process sugar data
#'
#' Function to process raw nonstructural carbohydrate data as read in using NSCprocessR::readRawNSCData ()
#' @param rawData Tibble containing the raw data, as read in by the readRawNSCData function.
#' @param cvLimitSample Limit for the coefficient of variation within a sample at which the sample if flagged for re-measurement.
#' @param cvLimitTube Limit for the coefficient of variation within a tube at which the sample is flagged for re-measurement.
#' @param prescribedStarchRecoveryFraction This parameter allows to prescribe the starch recovery fraction.
#' @param maxStarchRecoveryFraction Parameter setting the maximum starch recovery fraction, normally set to 90% as this is realistically recoverable.
#' @param LCS Name of the Laboratory control standards, which would be 'Oak' for the RC lab.
#' @param minimumSampleMass Threshold weight in milligram below which samples are dropped because they are unreliable. The threshold is set to 5 mg by default.
#' @param correctAbsorbance490 Boolean indicating whether absobance at 490 nm for soluble sugars is corrected using the sample blanks of the absorbance490_blank column.
#' @return The processed data, results of the anyalsis and a summary.
#' @import tibble
#' @export
processNSCs <- function (rawData,
                         cvLimitSample = 0.25,
                         cvLimitTube = 0.05, # Should this be 0.10 or 10% error, as suggested by Jim?
                         prescribedStarchRecoveryFraction = NA,
                         maxStarchRecoveryFraction = 0.9,
                         LCS = 'Oak',
                         minimumSampleMass = 5.0,
                         correctAbsorbance490 = TRUE) {

  # Load dependencies
  #--------------------------------------------------------------------------------------
  if (!existsFunction ('tibble')) library ('tidyverse')

  # Calculate the sample weight [g] and convert to [mg]
  #--------------------------------------------------------------------------------------
  rawData [['MassSample']] <- (rawData [['MassOfTubeAndSample']] -
                               rawData [['MassOfEmptyTube']]) * 1e3

  # Check whether a sample weight is below 10 mg and drop it if it is.
  #--------------------------------------------------------------------------------------
  indicesToDrop <- which (rawData [['MassSample']]              <  minimumSampleMass &
                          rawData [['SampleID']]                != 'B'               &
                          rawData [['SampleID']]                != 'TB'              &
                          substr (rawData [['SampleID']], 1, 3) != 'REF'             &
                          substr (rawData [['SampleID']], 1, 3) != 'Ref'             &
                          substr (rawData [['SampleID']], 1, 3) != 'LCS')
  if (length (indicesToDrop) >= 1) {
    rawData <- rawData [-indicesToDrop, ]
    warning (paste ('Warning: ',length (indicesToDrop),' samples are excluded because ',
                    'their sample mass is below the threshold of ',minimumSampleMass,
                    ' mg.', sep = ''))
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
  if (sum (!is.na (absorbances490)) > 0) {
    rawData [['MeanAbsorbance490']] <- rowMeans (absorbances490, na.rm = TRUE)
    if (correctAbsorbance490) {
      rawData [['CorrectedMeanAbsorbance490']] <- rawData [['MeanAbsorbance490']] -
        min (rawData [['Absorbance490_Blank']], rawData [['MeanAbsorbance490']], na.rm = TRUE)
    } else {
      rawData [['CorrectedMeanAbsorbance490']] <- NA
    }
  }
  absorbances525 <- cbind (rawData [['Absorbance525_1']], rawData [['Absorbance525_2']])
  if (sum (!is.na (absorbances525)) > 0) {
    rawData [['MeanAbsorbance525']] <- rowMeans (absorbances525, na.rm = TRUE)
  }

  # Calculate the within-sample coefficient of variation
  #--------------------------------------------------------------------------------------
  if (sum (!is.na (absorbances490)) > 0) {
    rawData [['SDAbsorbance490']] <- apply (absorbances490, 1, sd)
    rawData [['CVAbsorbance490']] <- rawData [['SDAbsorbance490']] /
      rawData [['MeanAbsorbance490']]
    rawData [['CVAbsorbance490']] [rawData [['CVAbsorbance490']] < 0.0   |
                                   is.na  (rawData [['CVAborbance490']]) |
                                   is.nan (rawData [['CVABsorbance490']])] <- NA
  } else {
    rawData [['CVAbsorbance490']] <- NA
  }
  if (sum (!is.na (absorbances525)) > 0) {
    rawData [['SDAbsorbance525']] <- apply (absorbances525, 1, sd)
    rawData [['CVAbsorbance525']] <- rawData [['SDAbsorbance525']] /
      rawData [['MeanAbsorbance525']]
    rawData [['CVAbsorbance525']] [rawData [['CVAbsorbance525']] < 0.0   |
                                   is.na  (rawData [['CVAborbance525']]) |
                                   is.nan (rawData [['CVABsorbance525']])] <- NA
  } else {
    rawData [['CVAbsorbance525']] <- NA
  }

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
    if (batch == batches [1] & !is.na (dates [1])) {
      extractionsSugar <- tibble::tibble (batch = batch, date = dates)
    } else if (!is.na (dates [1])) {
      extractionsSugar <- tibble::add_row (extractionsSugar, batch = batch, date = dates)
    }
  }

  # Delete rows that do not have calibration curves
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsSugar')) {
    if (sum (is.na (extractionsSugar [['date']])) > 0) {
      extractionsSugar <- extractionsSugar [-which (is.na (extractionsSugar [['date']])), ]
    }

    # Check that there is still one or more extractions or delete the variable
    #--------------------------------------------------------------------------------------
    if (dim (extractionsSugar) [1] == 0) rm (extractionsSugar)
  }

  # Compile a list of unique batches and dates for starch extractions
  #--------------------------------------------------------------------------------------
  for (batch in batches) {
    dates <- unique (rawData [['DateOfStarchAnalysis']] [rawData [['BatchID']] == batch])
    if (batch == batches [1] & !is.na (dates [1])) {
      extractionsStarch <- tibble (batch = batch, date = dates)
    } else if (!is.na (dates [1])) {
      extractionsStarch <- add_row (extractionsStarch, batch = batch, date = dates)
    }
  }

  # Delete rows that do not have calibration curves
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsStarch')) {
    if (sum (is.na (extractionsStarch [['date']])) > 0) {
      extractionsStarch <- extractionsStarch [-which (is.na (extractionsStarch [['date']])), ]
    }

    # Check that there is still one or more extractions or delete the variable
    #--------------------------------------------------------------------------------------
    if (dim (extractionsStarch) [1] == 0) rm (extractionsStarch)

  }

  # Make and save a sugar calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsSugar')) {
    for (extraction in 1:(dim (extractionsSugar) [1])) {

      # Get date of analysis and batch number
      #----------------------------------------------------------------------------------
      analysisDate <- extractionsSugar [['date']]  [extraction]
      batch        <- extractionsSugar [['batch']] [extraction]

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      refCondition <- (substr (rawData [['SampleID']], 1, 3)    == 'REF' |
                       substr (rawData [['SampleID']], 1, 3)    == 'Ref') &
                              rawData [['BatchID']]             == batch  &
                              rawData [['DateOfSugarAnalysis']] == analysisDate
      if (correctAbsorbance490) {
        referenceValues <- rawData  [['CorrectedMeanAbsorbance490']] [refCondition]
      } else {
        referenceValues <- rawData  [['MeanAbsorbance490']] [refCondition]
      }

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (rawData [['SampleID']] [refCondition], 4,
                              nchar (rawData [['SampleID']] [refCondition]))

      # Set reference solution concentrations to 100 for sugars
      #----------------------------------------------------------------------------------
      concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      concentrations [concentrations == '100/100'] <- 100.0 # 100 for sugar
      concentrations <- as.numeric (concentrations)

      # Sort out reference values with an absorbance above 1.0
      #----------------------------------------------------------------------------------
      indicesToKeep <- which (referenceValues >= -0.1 &
                              referenceValues <=  1.1)
      referenceValues <- referenceValues [indicesToKeep]
      concentrations  <- concentrations  [indicesToKeep]

      # Get the slope of a linear curve forced through 0 to make sure no values are
      # estimated to be negative values
      #----------------------------------------------------------------------------------
      # sugar <- slope * absorbance
      fitSugarAll <- lm (concentrations ~ 0 + referenceValues)

      # Check for outliers (residual > 2.0 * sigma) and take them out
      #----------------------------------------------------------------------------------
      indicesToDrop <- which (fitSugarAll$residuals > 2.0 * sigma (fitSugarAll))
      if (length (indicesToDrop) > 0) {
        concentrations  <- concentrations  [-indicesToDrop]
        referenceValues <- referenceValues [-indicesToDrop]
      }

      # Get the slope (startch = slope * absorbance)
      #----------------------------------------------------------------------------------
      fitSugar <- lm (concentrations ~ 0 + referenceValues)

      # Add the slopes to the tibble
      #----------------------------------------------------------------------------------
      batchCondition <- rawData [['BatchID']]              == batch &
                        rawData [['DateOfSugarAnalysis']] == analysisDate
      rawData [['SlopeSugar']] [batchCondition] <- fitSugar$coefficients

    } # End of extraction loop
  } # End exists ('extractionsSugar')

  # Make and save a starch calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsStarch')) {
    for (extraction in 1:(dim (extractionsStarch) [1])) {

      # Get date of analysis and batch number
      #--------------------------------------------------------------------------------------
      analysisDate <- extractionsStarch [['date']]  [extraction]
      batch        <- extractionsStarch [['batch']] [extraction]

      # Get the batch's mean absorbance at 525nm for tube blanks
      #--------------------------------------------------------------------------------
      batchCondition <- rawData [['BatchID']]              == batch &
                        rawData [['DateOfStarchAnalysis']] == analysisDate
      if (sum (rawData [['SampleID']] == 'TB' & batchCondition, na.rm = T) > 0.0) {
        batchTBAbsorbance <- mean (rawData [['MeanAbsorbance525']] [rawData [['SampleID']] == 'TB' &
                                                                    batchCondition],
                                   na.rm = T)
      } else { # In case there are no tube blanks we just assume no interference form tubes
        batchTBAbsorbance <- 0.0
      }

      # Determine starch correction factor from TB, unless they are larger than the sample
      # absorbances at 525nm. If the TB is larger than sample absorbance at 525nm, use
      # the smallest mean absorbance to avoid negetive numbers due to the correction.
      #--------------------------------------------------------------------------------
      condition <- batchCondition                                 &
                   substr (rawData [['SampleID']], 1, 3) != 'REF' &
                   substr (rawData [['SampleID']], 1, 3) != 'Ref' &
                   substr (rawData [['SampleID']], 1, 2) != 'TB'  &
                   substr (rawData [['SampleID']], 1, 1) != 'B'
      batchCorrection <- min (batchTBAbsorbance,
                              rawData [['MeanAbsorbance525']] [condition],
                              na.rm = TRUE)

      # Flag for higher TB than MeanAbsorbance525
      #--------------------------------------------------------------------------------
      if (batchCorrection != batchTBAbsorbance) {
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

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      refCondition <- (substr (rawData [['SampleID']], 1, 3) == 'REF' |
                       substr (rawData [['SampleID']], 1, 3) == 'Ref')      &
                      rawData [['BatchID']]                 == batch        &
                      rawData [['DateOfStarchAnalysis']]    == analysisDate &
                      !is.na (rawData [['SampleID']])
      referenceValues <- rawData [['CorrectedMeanAbsorbance525']] [refCondition]

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (rawData [['SampleID']] [refCondition], 4,
                                nchar (rawData [['SampleID']] [refCondition]))

      # Set reference solution concentrations to 200 for starch
      #----------------------------------------------------------------------------------
      concentrations [concentrations == '100/200'] <- 200.0 # 200 for starch
      concentrations [concentrations == '100/100'] <- 100.0 # 100 for starch # TR Waiting to hear from Jim to confirm this.
      concentrations <- as.numeric (concentrations)

      # Sort out reference values with an absorbance above 1.0
      #----------------------------------------------------------------------------------
      indicesToKeep <- which (referenceValues >= -0.1 &
                              referenceValues <= 1.1)

      # Get the slope (startch = slope * absorbance)
      #------------------------------------------------------------------------------------
      fitStarchAll <- lm (concentrations ~ 0 + referenceValues)

      # Check for outlier (residual > 2.0 * sigma) and take them out
      #------------------------------------------------------------------------------------
      indicesToDrop <- which (fitStarchAll$residuals > 2.0 * sigma (fitStarchAll))
      if (length (indicesToDrop) > 0) {
        concentrations  <- concentrations  [-indicesToDrop]
        referenceValues <- referenceValues [-indicesToDrop]
      }

      # Get the slope (startch = slope * absorbance)
      #------------------------------------------------------------------------------------
      fitStarch <- lm (concentrations ~ 0 + referenceValues)

      # Add the slopes to the tibble
      #----------------------------------------------------------------------------------
      rawData [['SlopeStarch']] [batchCondition] <- fitStarch$coefficients

    } # End of starch extractions loop
  } # End exists (extractionsStarch) condition

  # Determine concentrations from absorbance values for sugar # Call it a concentrations [mg ml-1]
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsSugar')) {
    if (correctAbsorbance490) {
      rawData [['ConcentrationSugar']] <- rawData [['CorrectedMeanAbsorbance490']] *
        rawData [['SlopeSugar']]
    } else {
      rawData [['ConcentrationSugar']] <- rawData [['MeanAbsorbance490']] *
        rawData [['SlopeSugar']]
    }
    rawData [['ConcentrationSugar']] [rawData [['ConcentrationSugar']] < 0.0   |
                                      is.na  (rawData [['ConcentrationSugar']])|
                                      is.nan (rawData [['ConcentrationSugar']])] <- NA

    # Correct sugar amount [mg ml-1] for background signal using no phenol
    #------------------------------------------------------------------------------------
    rawData [['MassSugar']] <- rawData [['ConcentrationSugar']] *
                               rawData [['DilutionFactorSugar']] *
                               rawData [['VolumeSugar']]

    # Convert mass to concentrations [mg g-1]
    #------------------------------------------------------------------------------------
    rawData [['ConcentrationSugarMgG']] <- rawData [['MassSugar']] /
      rawData [['MassSample']]
    rawData [['ConcentrationSugarMgG']] [rawData [['ConcentrationSugarMgG']] < 0.0   |
                                           is.na (rawData [['ConcentrationSugarMgG']]) |
                                           is.nan (rawData [['ConcentrationSugarMgG']])] <- NA

    # Convert concentrations from [mg g-1] to [% dry weight]
    #--------------------------------------------------------------------------------------
    rawData [['ConcentrationSugarPerDW']]  <- rawData [['ConcentrationSugarMgG']]  / 10.0
  }

  # Determine concentration of starch [mg ml-1] from absorbance
  #--------------------------------------------------------------------------------------
  if (exists ('extractionsStarch')) {
    rawData [['ConcentrationStarch']] <- rawData [['CorrectedMeanAbsorbance525']] *
                                         rawData [['SlopeStarch']]
    rawData [['ConcentrationStarch']] [rawData [['ConcentrationStarch']] < 0.0   |
                                       is.na  (rawData [['ConcentrationStarch']])|
                                       is.nan (rawData [['ConcentrationStarch']])] <- NA
    rawData [['MassStarch']] <- rawData [['ConcentrationStarch']] *
                                rawData [['DilutionFactorStarch']] *
                                rawData [['VolumeStarch']]

    # Calculate starch recovery for each combination of batch and date, unless it is prescribed
    #------------------------------------------------------------------------------------
    if (is.na (prescribedStarchRecoveryFraction)) {
      for (extraction in 1:(dim (extractionsStarch) [1])) {

        # Get date of analysis and batch number
        #------------------------------------------------------------------------------------
        analysisDate <- extractionsStarch [['date']]  [extraction]
        batch        <- extractionsStarch [['batch']] [extraction]

        # Get absorbances for potato starch for each combination of batch and analysisDate
        #--------------------------------------------------------------------------------
        refCondition <- substr (rawData [['SampleID']], 1, 10) == 'LCS Potato' &
                        rawData [['BatchID']]                  == batch        &
                        rawData [['DateOfStarchAnalysis']]     == analysisDate &
                        !is.na (rawData [['SampleID']])
        batchCondition <- rawData [['BatchID']]              == batch        &
                          rawData [['DateOfStarchAnalysis']] == analysisDate

        # Set flag for unrealistically high starch recovery fraction
        #--------------------------------------------------------------------------------
        rawData [['SRFHigh']] [batchCondition] <- 'N'

        # Check whether potato standard was run at all, otherwise use maxStarchRecoveryFraction
        #--------------------------------------------------------------------------------
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
          #--------------------------------------------------------------------------------
          rawData [['MeanStarchRecovery']] [batchCondition] <- meanRecoveryPer
        }
      }
    } else {
      # Set batch's mean starch recovery rate
      #--------------------------------------------------------------------------------
      rawData [['MeanStarchRecovery']] <- prescribedStarchRecoveryFraction
    }

    # Correct all samples for the mean starch recovery percentage of each batch
    #----------------------------------------------------------------------------------
    rawData [['CorrectedMassStarch']] <- rawData [['MassStarch']] /
                                        (rawData [['MeanStarchRecovery']] / 100.0)

    # Convert mass to concentrations [mg g-1]
    #------------------------------------------------------------------------------------
    rawData [['ConcentrationStarchMgG']] <- rawData [['CorrectedMassStarch']] /
                                            rawData [['MassSample']]

    # Convert concentrations from [mg g-1] to [% dry weight]
    #------------------------------------------------------------------------------------
    rawData [['ConcentrationStarchPerDW']] <- rawData [['ConcentrationStarchMgG']] / 10.0
  }

  # Check whether Lab Control Standard for Oak wood is high
  #----------------------------------------------------------------------------------
  if (LCS == 'Oak') {
    if (exists ('extractionsSugar') & exists ('extractionsStarch')) {
      for (extraction in 1:(dim (extractionsSugar) [1])) {

        # Get date of analysis and batch number
        #------------------------------------------------------------------------------------
        analysisDate <- extractionsStarch [['date']]  [extraction]
        batch        <- extractionsStarch [['batch']] [extraction]

        # Get absorbances for potato starch for each combination of batch and analysisDate
        #--------------------------------------------------------------------------------
        refCondition <- substr (rawData [['SampleID']], 1, 7) == 'LCS Oak'    &
                        rawData [['BatchID']]                 == batch        &
                        rawData [['DateOfStarchAnalysis']]    == analysisDate &
                        !is.na (rawData [['SampleID']])

        # check whether the batch has a LCS Oak
        #--------------------------------------------------------------------------------
        if (sum (refCondition) == 0) {
          warning (paste ('Warning: There is no LCS Oak for batch ',batch,
                          ' analysed on the ',analysisDate,'.', sep = ''))

          # Flag for high or low oak lab standard
          #------------------------------------------------------------------------------
          rawData [['LCSOakDeviation']] [rawData [['BatchID']] == batch] <- NA

        } else { # Check the LCS Oak is low or high

          # get LCS Oak standard and compare it against threshold
          #--------------------------------------------------------------------------------
          LCSOakSugar  <- rawData [['ConcentrationSugarPerDW']]  [refCondition]
          LCSOakStarch <- rawData [['ConcentrationStarchPerDW']] [refCondition]

          # Flag for high or low oak lab standard
          #--------------------------------------------------------------------------------
          rawData [['LCSOakDeviation']] [rawData [['BatchID']] == batch] <- 'N'

          # oak sugar was 3.59+-0.38 % DW at Harvard and 3.28+-0.51 at NAU thus far
          #--------------------------------------------------------------------------------
          if (mean (LCSOakSugar, na.rm = T) < 3.59-0.38 &
              mean (LCSOakSugar, na.rm = T) > 3.59+0.38) {
            rawData [['LCSOakDeviation']] [rawData [['BatchID']] == batch] <- 'sugar'
          }
          # oak starch was 2.45+-0.21 % DW at Harvard and 2.91+-0.33 at NAU thus far
          #--------------------------------------------------------------------------------
          if (mean  (LCSOakStarch, na.rm = T) < 2.45-0.21 &
              mean (LCSOakStarch, na.rm = T) > 2.45+0.21) {
            if (mean (LCSOakSugar, na.rm = T) < 3.59-0.38 &
                mean (LCSOakSugar, na.rm = T) > 3.59+0.38) {
              rawData [['LCSOakDeviation']] [rawData [['BatchID']] == batch] <- 'sugar & starch'
            } else {
              rawData [['LCSOakDeviation']] [rawData [['BatchID']] == batch] <- 'starch'
            }
          }
        } # End sum (refCondition) == 0 condition
      } # End extraction loop
    } # End exists ('extractionsSugar') & exists ('extractionsStarch')
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
