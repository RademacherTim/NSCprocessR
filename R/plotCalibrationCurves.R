#' Plot calibration curves
#'
#' Function to plot the calibration curves returned by the NSCprocessR::processNSC () function.
#' @param data Tibble with the processed data using the NSCprocessR::processNSC function.
#' @param epsilon This is a very small number (default = 0.01). When the sample absorbance range is within an order of magnitude of this number REF0 is included in the calibration curve.
#' @return pdf file with a graph of each batch's calibration curve.
#' @import tidyverse
#' @export
plotCalibrationCurves <- function (data,
                                   forceIntercept = FALSE,
                                   epsilon = 0.01) {

  # Compile a list of unique batches and dates for sugar extractions
  #--------------------------------------------------------------------------------------
  batches <- unique (data [['BatchID']])
  for (batch in batches) {
    dates <- unique (data [['DateOfSugarAnalysis']] [data [['BatchID']] == batch])
    if (batch == batches [1]) {
      extractionsSugar <- tibble (batch = batch, date = dates)
    } else {
      extractionsSugar <- add_row (extractionsSugar, batch = batch, date = dates)
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
    dates <- unique (data [['DateOfStarchAnalysis']] [data [['BatchID']] == batch])
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
      analysisDate <- extractionsSugar [['date']] [extraction]
      batch <- extractionsSugar [['batch']] [extraction]

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      condition <- substr (data [['SampleID']], 1, 3)     == 'REF' &
                           data [['BatchID']]             == batch &
                           data [['DateOfSugarAnalysis']] == analysisDate
      referenceValues <- data [condition, ]

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (data [['SampleID']] [condition],
                                4,
                                nchar (data [['SampleID']] [condition]))

      # Set reference solution concentrations to 100 for sugars
      #----------------------------------------------------------------------------------
      concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      concentrations <- as.numeric (concentrations)

      # Get the range of absorbance values for the samples
      #----------------------------------------------------------------------------------
      condition <- data [['BatchID']]                    == batch        &
                   substring (data [['SampleID']], 1, 3) != 'REF'        &
                   substring (data [['SampleID']], 1, 3) != 'LCS'        &
                   substring (data [['SampleID']], 1, 1) != 'B'          &
                   substring (data [['SampleID']], 1, 2) != 'TB'         &
                   data [['DateOfSugarAnalysis']]        == analysisDate &
                   !is.na (data [['CorrectedMeanAbsorbance490']])
      if (sum (condition, na.rm = TRUE) > 0) {
        absorbanceRange <- range (data [['CorrectedMeanAbsorbance490']] [condition],
                                  na.rm = TRUE)

        # Only use refence values that are within an order of magnitude of the sample
        # absorbances
        #--------------------------------------------------------------------------------
        fitRange <- c (absorbanceRange [1] / 10.0, absorbanceRange [2] * 10.0)
        if (fitRange [1] < epsilon) { # If the lower range is a very small number include REF0
          indicesToDrop <- which ((referenceValues [['CorrectedMeanAbsorbance490']] < fitRange [1] |
                                   referenceValues [['CorrectedMeanAbsorbance490']] > fitRange [2]) &
                                  referenceValues [['SampleID']] != 'REF0')
        } else {
          indicesToDrop <- which (referenceValues [['CorrectedMeanAbsorbance490']] < fitRange [1] |
                                  referenceValues [['CorrectedMeanAbsorbance490']] > fitRange [2])
        }
        if (length (indicesToDrop) > 0) {
          referenceValues <- referenceValues [-indicesToDrop, ]
          concentrations  <- concentrations  [-indicesToDrop]
        }
      } else {
        warning (paste ('Warning: There are no sugar samples in batch ',batch,
                        ' analysed on the ',analysisDate, sep = ''))
      }

      # Get the slope and intercept
      #----------------------------------------------------------------------------------
      if (forceIntercept) { # Get slope for intercepts forced through zero
      # sugar <- slope * absorbance + intercept
        fitSugar  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance490']])
      } else { # Get intercept and slope
        fitSugar  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance490']])
      }

      # Create fileName for the pdf
      #----------------------------------------------------------------------------------
      fileName <- paste ("calibrationCurve_",format (analysisDate, "%Y-%m-%d"),
                         "_batch",batch,"_sugar.pdf", sep = "")

      # Open the pdf
      #----------------------------------------------------------------------------------
      pdf (file = fileName)

      # Plot the sugar calibration curve
      #----------------------------------------------------------------------------------
      plot (x = referenceValues [['CorrectedMeanAbsorbance490']],
            y = concentrations,
            main = paste ('calibration curve for sugar (batch ',batch,'; ',analysisDate,')', sep = ''),
            las = 1,
            xlab = 'absorbance at 490 nm',
            ylab = 'sugar (mg / ml)') # TTR Should this be per ml
      points (x = referenceValues [['CorrectedMeanAbsorbance490']],
              y = concentrations,
              col = '#91b9a499',
              pch = 19)
      # Error bars are not really meaningful, because they are based on colour determinations for most of the values but on repetitions for REF100
      #if (length (unique (concentrations)) > length (concentrations)) {
      #  aggregate (referenceValues [['CorrectedMeanAbsorbance490']], list (concentrations), sd)
      #  arrows ()
      #}
      abline (fitSugar,
              col = 'grey',
              lwd = 2, lty = 2)
      text (x = mean (referenceValues [['CorrectedMeanAbsorbance490']]),
            y = 20,
            labels = expression (paste (R^2,' = ', sep = '')),
            pos = 4)
      text (x = mean (referenceValues [['CorrectedMeanAbsorbance490']]) * 1.2,
            y = 20,
            labels = round (summary (fitSugar)$r.squared, 3),
            pos = 4)

      # close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()
  }

  # Make and save a starch calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  for (extraction in 1:(dim (extractionsStarch) [1])) {

    # Get date of analysis and batch number
    #--------------------------------------------------------------------------------------
    analysisDate <- extractionsStarch [['date']] [extraction]
    batch <- extractionsStarch [['batch']] [extraction]

    # Get absorbances for reference values to create a calibration curve for each batch
    #------------------------------------------------------------------------------------
    condition <- substr (data [['SampleID']], 1, 3)      == 'REF' &
                         data [['BatchID']]              == batch &
                         data [['DateOfStarchAnalysis']] == analysisDate
    referenceValues <- data [condition, ]

    # Get reference solution concentrations
    #------------------------------------------------------------------------------------
    concentrations <- substr (data [['SampleID']] [condition],
                              4,
                              nchar (data [['SampleID']] [condition]))

    # Set reference solution concentrations to 200 for starch
    #------------------------------------------------------------------------------------
    concentrations [concentrations == '100/200'] <- 200.0 # 200 for starch
    concentrations <- as.numeric (concentrations)

    # Get the range of absorbance values for the samples
    #----------------------------------------------------------------------------------
    condition <- data [['BatchID']]                    == batch  &
                 substring (data [['SampleID']], 1, 3) != 'REF'  &
                 substring (data [['SampleID']], 1, 3) != 'LCS'  &
                 substring (data [['SampleID']], 1, 1) != 'B'    &
                 substring (data [['SampleID']], 1, 2) != 'TB'   &
                 data [['DateOfStarchAnalysis']] == analysisDate &
                 !is.na (data [['CorrectedMeanAbsorbance525']])
    if (sum (condition, na.rm = T) > 0) {
      absorbanceRange <- range (data [['CorrectedMeanAbsorbance525']] [condition],
                                na.rm = TRUE)

      # Only use refence values that are within an order of magnitude of the sample
      # absorbances
      #--------------------------------------------------------------------------------
      fitRange <- c (absorbanceRange [1] / 10.0, absorbanceRange [2] * 10.0)
      if (fitRange [1] < epsilon) { # If the lower range is a very small number include REF0
        indicesToDrop <- which ((referenceValues [['CorrectedMeanAbsorbance525']] < fitRange [1] |
                                 referenceValues [['CorrectedMeanAbsorbance525']] > fitRange [2]) &
                                referenceValues [['SampleID']] != 'REF0')
      } else {
        indicesToDrop <- which (referenceValues [['CorrectedMeanAbsorbance525']] < fitRange [1] |
                                referenceValues [['CorrectedMeanAbsorbance525']] > fitRange [2])
      }
      if (length (indicesToDrop) > 0) {
        referenceValues <- referenceValues [-indicesToDrop, ]
        concentrations  <- concentrations  [-indicesToDrop]
      }
    } else {
      warning (paste ('Warning: There are no sugar samples in batch ',batch,
                      ' analysed on the ',analysisDate, sep = ''))
    }

    # Get the slope and intercept (startch = slope * absorbance + intercept)
    #------------------------------------------------------------------------------------
    if (forceIntercept) { # Get slope for intercepts forced through zero
      fitStarch  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance525']])
    } else { # Get intercept and slope
      fitStarch  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance525']])
    }

    # Drop 250 from calibration curve if R2 is below 0.9
    #------------------------------------------------------------------------------------
    if (summary (fitStarch)$r.squared < 0.9) {
      indexToDrop <- which (concentrations == 250.0)
      concentrations <- concentrations [-indexToDrop]
      referenceValues <- referenceValues [-indexToDrop, ]
      if (forceIntercept) { # Get slope for intercepts forced through zero
        fitStarch  <- lm (concentrations ~ 0 + referenceValues [['CorrectedMeanAbsorbance525']])
      } else { # Get intercept and slope
        fitStarch  <- lm (concentrations ~ referenceValues [['CorrectedMeanAbsorbance525']])
      }
    }

    # check whether calibration curve is sufficiently precise
    #------------------------------------------------------------------------------------
    if (summary (fitStarch)$r.squared < 0.9) {
      warning (paste ('Warning: the calibration curve for batch ',batch,' on the ',
                   analysisDate,' has an R2 lower than 0.9.', sep = ''))
    }
    # Create fileName for the pdf
    #------------------------------------------------------------------------------------
    fileName <- paste ("calibrationCurve_",format (analysisDate, "%Y-%m-%d"),
                       "_batch",batch,"_starch.pdf", sep = "")

    # Open the pdf
    #------------------------------------------------------------------------------------
    pdf (file = fileName)

    # Plot the starch calibration curve
    #------------------------------------------------------------------------------------
    plot (x = referenceValues [['CorrectedMeanAbsorbance525']],
          y = concentrations,
          main = paste ('calibration curve for starch (batch ',batch,'; ',analysisDate,')', sep = ''),
          las = 1,
          xlab = 'absorbance at 525 nm',
          ylab = 'glucose equivalent (mg / ml)')
    points (x = referenceValues [['CorrectedMeanAbsorbance525']],
            y = concentrations,
            col = '#91b9a499',
            pch = 19)
    abline (fitStarch,
            col = 'grey',
            lwd = 2, lty = 2)
    text (x = mean (referenceValues [['CorrectedMeanAbsorbance525']]),
          y = 20,
          labels = expression (paste (R^2,' = ', sep = '')),
          pos = 4)
    text (x = mean (referenceValues [['CorrectedMeanAbsorbance525']]) * 1.1,
          y = 20,
          labels = round (summary (fitStarch)$r.squared, 3),
          pos = 4)

    # close graphics device
    #------------------------------------------------------------------------------------
    dev.off ()
  }

  # Return zero exit status, if it ran smoothly
  #--------------------------------------------------------------------------------------
  return (0)
}
# To do list:
# TTR discuss calibrations curves with Jim
# TTR add error bars as standard deviation
