#' Plot calibration curves
#'
#' Function to plot the calibration curves returned by the NSCprocessR::processNSC () function.
#' @param data Tibble with the processed data using the NSCprocessR::processNSC function.
#' @return pdf file with a graph of each batch's calibration curve.
#' @export
plotCalibrationCurves <- function (data, forceIntercept = FALSE) {

  # Compile a list of unique batches # TTR Might ahve to include analysis dates in the future
  #--------------------------------------------------------------------------------------
  batches <- unique (data [['BatchID']])

  # Make and save a calibration curve for each batch # TTR Should be dependent on BatchID and Date, because BatchID might be repeated.
  #--------------------------------------------------------------------------------------
  for (batch in batches) {

    for (NSC in c ('sugar','starch')) {

      # Get date of analysis
      #--------------------------------------------------------------------------------------
      if (NSC == 'sugar') {
        analysisDate <- unique (data [['DateOfSugarAnalysis']] [data [['BatchID']] == batch])
      } else if (NSC == 'starch') {
        analysisDate <- unique (data [['DateOfStarchAnalysis']] [data [['BatchID']] == batch])
      }
      if (length (analysisDate) > 1) stop ('Error: The batch has multiple associated analysis dates!')

      # Get absorbances for reference values to create a calibration curve for each batch
      #----------------------------------------------------------------------------------
      condition <- substr (data [['SampleID']], 1, 3) == 'REF' &
                           data [['BatchID']] == batch
      referenceValues <- data [condition, ]

      # Get reference solution concentrations
      #----------------------------------------------------------------------------------
      concentrations <- substr (data [['SampleID']] [condition],
                                4,
                                nchar (data [['SampleID']] [condition]))

      # Set reference solution concentrations depending on whether we are calibrating
      # starch or sugar
      #----------------------------------------------------------------------------------
      if (NSC == 'sugar') {
        concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      } else if (NSC == 'starch'){
        concentrations [concentrations == '100/200'] <- 200.0 # 200 for starrch
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

      # Create fileName for the pdf
      #----------------------------------------------------------------------------------
      fileName <- paste ("calibrationCurveBatch",batch,"_",format (analysisDate, "%Y-%m-%d"),
                         "_",NSC,".pdf", sep = "")

      # Open the pdf
      #----------------------------------------------------------------------------------
      pdf (file = fileName)

      # Plot the sugar calibration curve
      #----------------------------------------------------------------------------------
      if (NSC == 'sugar') {
        plot (x = referenceValues [['CorrectedMeanAbsorbance490']],
              y = concentrations,
              main = paste ('calibration curve for ',NSC,' (batch ',batch,'; ',analysisDate,')', sep = ''),
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
      } else if (NSC == 'starch') { # Plot starch calibration curve
        plot (x = referenceValues [['CorrectedMeanAbsorbance525']],
              y = concentrations,
              main = paste ('calibration curve for ',NSC,' (batch ',batch,'; ',analysisDate,')', sep = ''),
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
      }

      # close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()

    } # End NSC (sugar versus starch) loop

  } # End of batch loop

  # Return zero exit status, if it ran smoothly
  return (0)
}
# To do list: Make a separate calibration curve for each extraction date. At the moment, it only depends on batch number, which will eventually be repeated.
# TTR test function, produce pdf and discuss calibrations curves with Jim
# TTR add error bars as standard deviation
# TTR create plot of REF100 values over time with mean and standard deviation
# TTR create plot of LCS over time (Oak wood and potato)
