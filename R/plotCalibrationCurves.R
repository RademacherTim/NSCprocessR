#' Plot calibration curves
#'
#' Function to plot the calibration curves returned by the NSCprocessR::processNSC () function.
#' @param data Tibble with the processed data using the NSCprocessR::processNSC function.
#' @return pdf file with a graph of each batch's calibration curve.
#' @import tidyverse
#' @export
plotCalibrationCurves <- function (data) {

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

  # Compile a list of unique batches and dates for starch extractions
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

  # Create two new colums with the average absorbance at 490 and 525 nm
  #--------------------------------------------------------------------------------------
  absorbances490 <- cbind (data [['Absorbance490_1']], data [['Absorbance490_2']])
  if (sum (!is.na (absorbances490)) > 0) {
    data [['MeanAbsorbance490']] <- rowMeans (absorbances490, na.rm = TRUE)
    data [['CorrectedMeanAbsorbance490']] <- data [['MeanAbsorbance490']] -
      min (data [['Absorbance490_Blank']], data [['MeanAbsorbance490']], na.rm = TRUE)
  }
  absorbances525 <- cbind (data [['Absorbance525_1']], data [['Absorbance525_2']])
  if (sum (!is.na (absorbance525)) > 0) {
    data [['MeanAbsorbance525']] <- rowMeans (absorbances525, na.rm = TRUE)
  }

  # Make and save a sugar calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  if (dim (extractionsSugar) [1] != 0) {
    for (extraction in 1:(dim (extractionsSugar) [1])) {

      # Get date of analysis and batch number
      #--------------------------------------------------------------------------------------
      analysisDate <- extractionsSugar [['date']] [extraction]
      batch <- extractionsSugar [['batch']] [extraction]

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      refCondition <- (substr (data [['SampleID']], 1, 3)    == 'REF' |
                       substr (data [['SampleID']], 1, 3)    == 'Ref' ) &
                              data [['BatchID']]             == batch   &
                              data [['DateOfSugarAnalysis']] == analysisDate
      referenceValues <- data  [['CorrectedMeanAbsorbance490']] [refCondition]

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (data [['SampleID']] [refCondition], 4,
                                nchar (data [['SampleID']] [refCondition]))

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

      # Get the slope and intercept of a linear curve forced through 0 to make sure no
      # are estimated to be negative values
      #----------------------------------------------------------------------------------
      # sugar <- slope * absorbance + intercept
      fitSugarAll <- lm (concentrations ~ 0 + referenceValues)

      # Check for outliers (residual > 2.0 * sigma) and take them out
      #------------------------------------------------------------------------------------
      allReferenceValues <- referenceValues
      allConcentrations  <- concentrations
      indicesToDrop <- which (fitSugarAll$residuals > 2.0 * sigma (fitSugarAll))
      if (length (indicesToDrop) > 0) {
        concentrations  <- concentrations  [-indicesToDrop]
        referenceValues <- referenceValues [-indicesToDrop]
      }

      # Get the slope and intercept (startch = slope * absorbance + intercept)
      #------------------------------------------------------------------------------------
      fitSugar <- lm (concentrations ~ 0 + referenceValues)

      # Create fileName for the pdf
      #----------------------------------------------------------------------------------
      fileName <- paste ("calibrationCurve_",format (analysisDate, "%Y-%m-%d"),
                         "_batch",batch,"_sugar.pdf", sep = "")

      # Open the pdf
      #----------------------------------------------------------------------------------
      pdf (file = fileName)

      # Plot the sugar calibration curve
      #----------------------------------------------------------------------------------
      plot (x = allReferenceValues,
            y = allConcentrations,
            main = paste ('calibration curve for sugar (batch ',batch,'; ',analysisDate,')', sep = ''),
            las = 1,
            xlab = 'absorbance at 490 nm',
            ylab = 'sugar (mg / ml)',
            xlim = c (0, 1))
      points (x = referenceValues,
              y = concentrations,
              col = '#91b9a499',
              pch = 19)
      abline (fitSugarAll,
              col = 'gray',
              lwd = 1, lty = 2)
      abline (fitSugar,
              col = '#feb24c',
              lwd = 2, lty = 2)
      text (x = 0.2,
            y = max (allConcentrations) * 0.9,
            labels = expression (paste (R^2,' = ', sep = '')),
            pos = 2)
      text (x = 0.2,
            y = max (allConcentrations) * 0.9,
            labels = round (summary (fitSugar)$r.squared, 4),
            pos = 4)

      # close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()
    }
  }

  # Make and save a starch calibration curve for each combination of batch and date
  #--------------------------------------------------------------------------------------
  for (extraction in 1:(dim (extractionsStarch) [1])) {

    # Get date of analysis and batch number
    #--------------------------------------------------------------------------------------
    analysisDate <- extractionsStarch [['date']] [extraction]
    batch <- extractionsStarch [['batch']] [extraction]

    # Get the batch's mean absorbance at 525nm for tube blanks
    #--------------------------------------------------------------------------------
    batchCondition <- data [['BatchID']]              == batch &
                      data [['DateOfStarchAnalysis']] == analysisDate
    if (sum (data [['SampleID']] == 'TB' & batchCondition, na.rm = T) > 0.0) {
      batchTBAbsorbance <- mean (data [['MeanAbsorbance525']] [data [['SampleID']] == 'TB' &
                                                                    batchCondition],
                                 na.rm = T)
    } else { # In case there are no tube blanks we just assume no interference form tubes
      batchTBAbsorbance <- 0.0
    }

    # Determine correction factor from TB, unless they are larger than the sample
    # absorbances at 525nm. If the TB is larger than sample absorbance at 525nm, use
    # the smallest mean absorbance to avoid negetive numbers due to the correction.
    #--------------------------------------------------------------------------------
    condition <- batchCondition                              &
                 substr (data [['SampleID']], 1, 3) != 'REF' &
                 substr (data [['SampleID']], 1, 3) != 'Ref' &
                 substr (data [['SampleID']], 1, 2) != 'TB'  &
                 substr (data [['SampleID']], 1, 1) != 'B'
    batchCorrection <- min (batchTBAbsorbance,
                            data [['MeanAbsorbance525']] [condition],
                            na.rm = TRUE)

    # Correct mean absorbance values at 525nm
    #--------------------------------------------------------------------------------
    data [['CorrectedMeanAbsorbance525']] [batchCondition] <-
      data [['MeanAbsorbance525']] [batchCondition] - batchCorrection

    # Get reference solution concentrations
    #------------------------------------------------------------------------------------
    refCondition <- (substr (data [['SampleID']], 1, 3) == 'REF' |
                     substr (data [['SampleID']], 1, 3) == 'Ref' ) &
                    data [['BatchID']]                 == batch    &
                    data [['DateOfStarchAnalysis']]    == analysisDate
    concentrations <- substr (data [['SampleID']] [refCondition], 4,
                              nchar (data [['SampleID']] [refCondition]))
    #------------------------------------------------------------------------------------
    referenceValues <- data [['CorrectedMeanAbsorbance525']] [refCondition]

    # Set reference solution concentrations to 200 for starch
    #------------------------------------------------------------------------------------
    concentrations [concentrations == '100/200'] <- 200.0 # 200 for starch
    concentrations [concentrations == '100/100'] <- 100.0 # 100 for starch
    concentrations <- as.numeric (concentrations)

    # Sort out reference values with an absorbance above 1.0
    #----------------------------------------------------------------------------------
    indicesToKeep <- which (referenceValues >= -0.1 &
                            referenceValues <= 1.1)

    # Get the slope and intercept (startch = slope * absorbance + intercept)
    #------------------------------------------------------------------------------------
    fitStarchAll <- lm (concentrations ~ 0 + referenceValues)

    # Check for outlier (residual > 2.0 * sigma) and take them out
    #------------------------------------------------------------------------------------
    allReferenceValues <- referenceValues
    allConcentrations  <- concentrations
    indicesToDrop <- which (fitStarchAll$residuals > 2.0 * sigma (fitStarchAll))
    if (length (indicesToDrop) > 0) {
      concentrations  <- concentrations  [-indicesToDrop]
      referenceValues <- referenceValues [-indicesToDrop]
    }

    # Get the slope and intercept (startch = slope * absorbance + intercept)
    #------------------------------------------------------------------------------------
    fitStarch <- lm (concentrations ~ 0 + referenceValues)

    # Create fileName for the pdf
    #------------------------------------------------------------------------------------
    fileName <- paste ("calibrationCurve_",format (analysisDate, "%Y-%m-%d"),
                       "_batch",batch,"_starch.pdf", sep = "")

    # Open the pdf
    #------------------------------------------------------------------------------------
    pdf (file = fileName)

    # Plot the starch calibration curve
    #------------------------------------------------------------------------------------
    plot (x = allReferenceValues,
          y = allConcentrations,
          main = paste ('calibration curve for starch (batch ',batch,'; ',analysisDate,')', sep = ''),
          las = 1,
          xlab = 'absorbance at 525 nm',
          ylab = 'glucose equivalent (mg / ml)',
          xlim = c (0, 1))
    points (x = referenceValues,
            y = concentrations,
            col = '#91b9a499',
            pch = 19)
    abline (fitStarchAll,
            col = 'gray',
            lwd = 1, lty = 2)
    abline (fitStarch,
            col = '#feb24c',
            lwd = 2, lty = 2)
    text (x = 0.2,
          y = max (allConcentrations) * 0.9,
          labels = expression (paste (R^2,' = ', sep = '')),
          pos = 2)
    text (x = 0.2,
          y = max (allConcentrations) * 0.9,
          labels = round (summary (fitStarch)$r.squared, 4),
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
