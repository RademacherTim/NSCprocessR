#' Plot calibration curves
#'
#' Function to plot the calibration curves returned by the NSCprocessR::processNSC () function.
#' @param data Tibble with the processed data.
#' @return pdf file with a graph of each batch's calibration curve.
#' @export
plotCalibrationCurve <- function (data) {

  # Compile a list of batches
  #--------------------------------------------------------------------------------------
  batches <- unique (data [['BatchID']])

  # Make and save a calibration curve for each batch TTR I need to adjust the absorbance of the calibration curve for the no phenol blanks
  #--------------------------------------------------------------------------------------
  for (batch in batches) {
    for (NSC in c ('sugar','starch')) {
      condition <- substr (data [['SampleID']], 1, 3) == 'REF' & data [['BatchID']] == batch
      referenceValues <- data [condition, ]
      concentrations <- substr (data [['SampleID']] [condition], 4, nchar (data [['SampleID']] [condition]))
      if (NSC == 'sugar') {
        concentrations [concentrations == '100/200'] <- 100.0 # 100 for sugar
      } else if (NSC == 'starch'){
        concentrations [concentrations == '100/200'] <- 200.0 # 200 for starrch
      }
      concentrations <- as.numeric (concentrations)

      # Get the slope and intercept
      #----------------------------------------------------------------------------------
      if (forceIntercept) { # Get slope for intercepts forced through zero
        if (NSC == 'sugar') {
          fitSugar  <- lm (referenceValues [['MeanAbsorbance490']] -
                           referenceValues [['Absorbance490_Blank']] ~ 0 + concentrations)
        } else if (NSC == 'starch') {
          fitStarch <- lm (referenceValues [['MeanAbsorbance525']] ~ 0 + concentrations)
        }
      } else { # Get intercept and slope
        if (NSC == 'sugar') {
          fitSugar  <- lm (referenceValues [['MeanAbsorbance490']] -
                           referenceValues [['Absorbance490_Blank']] ~ concentrations)
        } else if (NSC == 'starch') {
          fitStarch <- lm (referenceValues [['MeanAbsorbance525']] ~ concentrations)
        }
      }

      # Add the slopes and intercepts to respective rows in the tibble
      #----------------------------------------------------------------------------------
      if (forceIntercept) { # intercept forced through 0.0
        if (NSC == 'sugar') {
          data [['InterceptSugar']]  [data [['BatchID']] == batch] <- 0.0
          data [['SlopeSugar']]      [data [['BatchID']] == batch] <- fitSugar$coefficients
        } else if (NSC == 'starch') {
          data [['InterceptStarch']] [data [['BatchID']] == batch] <- 0.0
          data [['SlopeStarch']]     [data [['BatchID']] == batch] <- fitStarch$coefficients
        }
      } else { # variable intercept
        if (NSC == 'sugar') {
          data [['InterceptSugar']]  [data [['BatchID']] == batch] <- fitSugar$coefficients [1]
          data [['SlopeSugar']]      [data [['BatchID']] == batch] <- fitSugar$coefficients [2]
        } else if (NSC == 'starch') {
          data [['InterceptStarch']] [data [['BatchID']] == batch] <- fitStarch$coefficients [1]
          data [['SlopeStarch']]     [data [['BatchID']] == batch] <- fitStarch$coefficients [2]
        }
      }
    }

    # Create fileName for the pdf
    fileName <- paste ("calibrationCurvesBatch",batch,".pdf")

    # Open the pdf
    pdf (file = fileName)

    # Plot the sugar calibration curve
    if (NSC == 'sugar') {
      plot (x = concentrations,
            y = referenceValues [['MeanAbsorbance490']] - referenceValues [['Absorbance490_Blank']],
            las = 1,
            xlab = 'sugar (mg)',
            ylab = 'absorbance at 490 nm')
      points (x = concentrations,
              y = referenceValues [['MeanAbsorbance490']] - referenceValues [['Absorbance490_Blank']],
              col = '#91b9a499',
              pch = 19)
      abline (fitSugar,
              col = 'grey',
              lwd = 2, lty = 2, )
      text (x = 30,
            y = max (referenceValues [['MeanAbsorbance490']] - referenceValues [['Absorbance490_Blank']]) * 0.9,
            labels = expression (paste (R^2,' = ', sep = '')),
            pos = 4)
      text (x = 50,
            y = max (referenceValues [['MeanAbsorbance490']] - referenceValues [['Absorbance490_Blank']]) * 0.9,
            labels = round (summary (fitSugar)$r.squared, 3),
            pos = 4)
    } else if (NSC == 'starch') { # Plot starch calibration curve
      plot (x = concentrations,
            y = referenceValues [['MeanAbsorbance525']],
            las = 1,
            xlab = 'glucose equivalent (mg)',
            ylab = 'absorbance at 525 nm')
      points (x = concentrations,
              y = referenceValues [['MeanAbsorbance525']],
              col = '#91b9a499',
              pch = 19)
      abline (fitStarch,
              col = 'grey',
              lwd = 2, lty = 2, )
      text (x = 30,
            y = max (referenceValues [['MeanAbsorbance525']]) * 0.9,
            labels = expression (paste (R^2,' = ', sep = '')),
            pos = 4)
      text (x = 50,
            y = max (referenceValues [['MeanAbsorbance525']]) * 0.9,
            labels = round (summary (fitStarch)$r.squared, 3),
            pos = 4)
    }

    # close graphics device
    dev.off ()

  } # End of batch loop

  # Return zero exit status, if it ran smoothly
  return (0)
}
# To do list: Make a separate calibration curve for each extraction date. At the moment, it only depends on batch number, which will eventually be repeated.
