#' Process sugar data
#'
#' Function to process raw nonstructural carbohydrate data as read in using NSCprocessR::readRawNSCData ()
#' @param rawData Tibble containing the raw data, as read in by the readRawNSCData function.
#' @param cvLimitSample Limit for the coefficient of variation within a sample at which the sample if flagged for re-measurement.
#' @param cvLimitTube Limit for the coefficient of variation within a tube at which the sample is flagged for re-measurement.
#' @param forceIntercept Logical variable describing whether to force the calibration curve intercept through 0 or not.
#' @return The processed data, results of the anyalsis and a summary.
#' @export
processNSCs <- function (rawData,
                         cvLimitSample = 0.25,
                         cvLimitTube = 0.10,
                         forceIntercept = FALSE,
                         maxStarchRecoveryFraction = 0.9) {

  # Calculate the sample weight [g] and convert to [mg]
  #--------------------------------------------------------------------------------------
  rawData [['MassSample']] <- (rawData [['MassOfTubeAndSample']] -
                               rawData [['MassOfEmptyTube']]) * 1e3

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
  rawData [['MeanAbsorbance490']] <- rowMeans (absorbances490)
  rawData [['MeanAbsorbance525']] <- rowMeans (absorbances525)
  rawData [['CorrectedMeanAbsorbance490']] <- rawData [['MeanAbsorbance490']] -
                                              rawData [['Absorbance490_Blank']]

  # Calculate the within-sample coefficient of variation
  #--------------------------------------------------------------------------------------
  rawData [['SDAbsorbance490']] <- apply (absorbances490, 1, sd)
  rawData [['SDAbsorbance525']] <- apply (absorbances525, 1, sd)
  rawData [['CVAbsorbance490']] <- rawData [['SDAbsorbance490']] /
                                   rawData [['MeanAbsorbance490']]
  rawData [['CVAbsorbance525']] <- rawData [['SDAbsorbance525']] /
                                   rawData [['MeanAbsorbance525']]

  # Flag samples with high CV for re-measuring
  #--------------------------------------------------------------------------------------
  rawData [['HighCV']] <- 'N'
  rawData [['HighCV']] [rawData [["CVAbsorbance490"]] > cvLimitTube] <- 'Sugar'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube &
                        rawData [['HighCV']] == 'N'] <- 'Starch'
  rawData [['HighCV']] [rawData [["CVAbsorbance525"]] > cvLimitTube &
                        rawData [['HighCV']] == 'Sugar'] <- 'Sugar and starch'

  # Extract blanks with solution
  #--------------------------------------------------------------------------------------
  #blanks <- rawData [rawData [['SampleID']] == 'B', ]

  # Extract blanks with only empty tubes
  #--------------------------------------------------------------------------------------
  #tubeBlanks <- rawData [rawData [['SampleID']] == 'TB', ]

  # Extract Lab Control Standards (LCS)
  #--------------------------------------------------------------------------------------
  #labControlStandards <- rawData [substr (rawData [['SampleID']], 1, 3) == 'LCS', ]

  # Compile a list of batches
  #--------------------------------------------------------------------------------------
  batches <- unique (rawData [['BatchID']])

  # Make and save a calibration curve for each batch
  #--------------------------------------------------------------------------------------
  for (batch in batches) {
    for (NSC in c ('sugar','starch')) {

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      condition <- substr (rawData [['SampleID']], 1, 3) == 'REF' &
                           rawData [['BatchID']] == batch
      referenceValues <- rawData [condition, ]

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (rawData [['SampleID']] [condition],
                                4,
                                nchar (rawData [['SampleID']] [condition]))

      # Set reference solution concentrations depending on whether we are calibrating
      # starch or sugar
      #----------------------------------------------------------------------------------
      if (NSC == 'sugar') {
        concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      } else if (NSC == 'starch'){
        concentrations [concentrations == '100/200'] <- 200.0 # 200 for starrch

        # Get the batch's mean absorbance at 525nm for tube blanks
        #--------------------------------------------------------------------------------
        batchTBAbsorbance <- rawData [['MeanAbsorbance525']] [rawData [['SampleID']] == 'TB']

        # Check that TB absorbance is not larger than the sample, otherwise set flag
        #--------------------------------------------------------------------------------
        #if (batchTBAbsorbance > ) {}

        # Correct reference value for tube blank absorbance
        #--------------------------------------------------------------------------------
        rawData [['CorrectedMeanAbsorbance525']] [rawData [['BatchID']] == batch] <- rawData [['MeanAbsorbance525']] -
                                                                                     batchTBAbsorbance
      }
      concentrations <- as.numeric (concentrations)

      # Get the slope and intercept
      #----------------------------------------------------------------------------------
      if (forceIntercept) { # Get slope for intercepts forced through zero
        if (NSC == 'sugar') { # sugar <- slope * absorbance + intercept
          fitSugar  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance490']])
        } else if (NSC == 'starch') {
          fitStarch <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance525']])
        }
      } else { # Get intercept and slope
        if (NSC == 'sugar') {
          fitSugar  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance490']])
        } else if (NSC == 'starch') {
          fitStarch <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance525']])
        }
      }

      # Add the slopes and intercepts to respective rows in the tibble
      #----------------------------------------------------------------------------------
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

  # Determine concentrations from absorbance values for sugar # Call it a concentrations [mg ml-1]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugar']] <- rawData [['CorrectedMeanAbsorbance490']] *
                                      rawData [['SlopeSugar']] +
                                      rawData [['InterceptSugar']]

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
  rawData [['MassStarch']] <- rawData [['ConcentrationStarch']] *
                              rawData [['DilutionFactorStarch']] *
                              rawData [['VolumeStarch']]

  # Calcualte starch recovery
  #--------------------------------------------------------------------------------------
  for (batch in batches) {
    # Get absorbances for potato starch for each batch
    #----------------------------------------------------------------------------------
    condition <- substr (rawData [['SampleID']], 1, 3) == 'LCS potato' &
                         rawData [['BatchID']] == batch
    LCSPotato <- rawData [condition, ]
    meanPotatoMassRecovered <- mean (LCSPotato [['MassStarch']])
    meanPotatoMass <- mean (LCSPotato [['MassSample']])
    meanRecoveryPer <- meanPotatoMassRecovered /
                      (meanPotatoMass * maxStarchRecoveryFraction) * 100.0

    # Set batch's mean starch recovery rate
    #----------------------------------------------------------------------------------
    rawData [['MeanStarchRecovery']] [rawData [['BatchID']] == batch] <- meanRecoveryPer
  }

  # Correct all samples for the mean starch recovery percentage of each batch
  #------------------------------------------------------------------------------------
  rawData [['CorrectedMassStarch']] <- rawData [['MassStarch']] /
                                      (rawData [['MeanStarchRecovery']] / 100.0)

  # Convert mass to concentrations [mg g-1]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugarMgG']]  <- rawData [['MassSugar']] /
                                          rawData [['MassSample']]
  rawData [['ConcentrationStarchMgG']] <- rawData [['CorrectedMassStarch']] /
                                          rawData [['MassSample']]

  # Convert concentrations from [mg g-1] to [% dry weight]
  #--------------------------------------------------------------------------------------
  rawData [['ConcentrationSugarPerDW']]  <- rawData [['ConcentrationSugarMgG']]  / 10.0
  rawData [['ConcentrationStarchPerDW']] <- rawData [['ConcentrationStarchMgG']] / 10.0


  # Declare data as processed
  #--------------------------------------------------------------------------------------
  processedData <- rawData

  # Return the processed data
  #--------------------------------------------------------------------------------------
  return (processedData)
}
# To-do-list:
# TTR Correct starch concentration using the TB and percentage of extracted potato starch
