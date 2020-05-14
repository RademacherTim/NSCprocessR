#' Plot laboratory control standards over time
#'
#' Function to plot the laboratory control standards over time.
#' @param data Tibble with the processed data using the NSCprocessR::processNSC function. It can also be from the already processed data in the master file read in using NSCprocessR::readLabMasterSheet () function.
#' @return pdf file with a graph of each lab control standard over time.
#' @export
plotLabControlStandards <- function (data) {

  # Compile a list of unique LCS, whose SampleID must start with LCS
  #--------------------------------------------------------------------------------------
  labControlStandards <- unique (data [['SampleID']] [substr (data [['SampleID']], 1, 3) == 'LCS'])

  # Make graph of
  #--------------------------------------------------------------------------------------
  for (standard in labControlStandards) {

    # Check whether it is a starch or sugar standard
    #----------------------------------------------------------------------------------
    if (standard == 'LCS Oak') {
      datesSugar  <- data [['DateOfSugarAnalysis']] [data [['SampleID']] == standard]
      datesStarch <- data [['DateOfStarchAnalysis']] [data [['SampleID']] == standard]
      concentrationsOfSugarStds <- data [['ConcentrationSugarPerDW']] [data [['SampleID']] == standard]
      concentrationsOfStarchStds <- data [['ConcentrationStarchPerDW']] [data [['SampleID']] == standard]
    } else if (standard == 'LCS Potato') {
      datesStarch <- data [['DateOfStarchAnalysis']] [data [['SampleID']] == standard]
      concentrationsOfStarchStds <- data [['ConcentrationStarchPerDW']] [data [['SampleID']] == standard]
    }

    # Check whether this standard is for sugar
    #----------------------------------------------------------------------------------
    if (exists ('datesSugar')) {

      # Create fileName for the pdf
      #----------------------------------------------------------------------------------
      fileName <- paste (standard,"_SolubleSugar.png", sep = "")

      # Open the png
      #----------------------------------------------------------------------------------
      png (file = fileName)

      # Plot the sugar standard over time
      #----------------------------------------------------------------------------------
      plot (x = datesSugar,
            y = concentrationsOfSugarStds,
            main = paste (standard,'over time'),
            xlab = 'time',
            ylab = 'soluble sugar concentration (% dry weight)',
            las = 1,
            col = '#91b9a499',
            pch = 19,
            ylim  = c (0, 1.5*max (concentrationsOfSugarStds)))
      if (standard == 'LCS Oak') {

        # Add morgan's interval of variation
        #--------------------------------------------------------------------------------
        rect (xleft   = min (datesSugar) - 86400*365,
              xright  = max (datesSugar) + 86400*365,
              ybottom = 3.59-0.38,
              ytop    = 3.59+0.38,
              col = "#33333333",
              lty = 0)
        abline (h = 3.59,
                col = '#66666666')

        # Add jim's interval of variation
        #--------------------------------------------------------------------------------
        rect (xleft   = min (datesSugar) - 86400*365,
              xright  = max (datesSugar) + 86400*365,
              ybottom = 2.45-0.21,
              ytop    = 2.45+0.21,
              col = "#91b9a444",
              lty = 0)
        abline (h = 2.45,
                col = '#10647066',
                lty = 2)
      }

      # Add the mean and standard deviation
      #----------------------------------------------------------------------------------
      means <- aggregate (concentrationsOfSugarStds, by = list (datesSugar), FUN = mean)
      sds <- aggregate (concentrationsOfSugarStds, by = list (datesSugar), FUN = sd)
      arrows (x0 = means [, 1],
              y0 = means [, 2] - sds [, 2],
              x1 = means [, 1],
              y1 = means [, 2] + sds [, 2],
              col = '#DC143C99',
              length = 0.05,
              angle = 90,
              code = 3)
      points (x = means [, 1],
              y = means [, 2],
              pch = 15,
              col = '#DC143C99')

      # Add legend
      #----------------------------------------------------------------------------------
      legend (x = min (dates),
              y = 1.5*max (concentrationsOfSugarStds),
              legend = c ('mean ± standard deviation','individual measurement'),
              pch = c (15, 19),
              col = c ('#DC143C99','#91b9a499'),
              box.lty = 0,
              bg = 'transparent')

      # Close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()

    } # End if sugar exists condition

    # Check whether a starch standard exists
    #----------------------------------------------------------------------------------
    if (exists ('datesStarch')) {

      # Create fileName for the png file
      #----------------------------------------------------------------------------------
      fileName <- paste (standard,"_Starch.png", sep = "")

      # Open the png file
      #----------------------------------------------------------------------------------
      png (file = fileName)

      # Plot the starch standard over time
      #----------------------------------------------------------------------------------
      plot (x = datesStarch,
            y = concentrationsOfStarchStds,
            main = paste (standard,'over time'),
            xlab = 'time',
            ylab = 'soluble sugar concentration (% dry weight)',
            las = 1,
            col = '#91b9a499',
            pch = 19,
            ylim  = c (0, 1.5*max (concentrationsOfStarchStds)))
      if (standard == 'LCS Potato') {
        abline (h = 100,
                col = '#66666666',
                lty = 1)
      }

      # Add the mean and standard deviation
      #----------------------------------------------------------------------------------
      means <- aggregate (concentrationsOfStarchStds, by = list (datesStarch), FUN = mean)
      sds <- aggregate (concentrationsOfStarchStds, by = list (datesStarch), FUN = sd)
      arrows (x0 = means [, 1],
              y0 = means [, 2] - sds [, 2],
              x1 = means [, 1],
              y1 = means [, 2] + sds [, 2],
              col = '#DC143C99',
              length = 0.05,
              angle = 90,
              code = 3)
      points (x = means [, 1],
              y = means [, 2],
              pch = 15,
              col = '#DC143C99')

      # Add legend
      #----------------------------------------------------------------------------------
      legend (x = min (datesStarch),
              y = 1.5*max (concentrationsOfStarchStds),
              legend = c ('mean ± standard deviation','individual measurement'),
              pch = c (15, 19),
              col = c ('#DC143C99','#91b9a499'),
              box.lty = 0,
              bg = 'transparent')

      # Close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()

    } # End if starch exists condition

  } # End the standards loop

  # Return zero exit status, if it ran smoothly
  #--------------------------------------------------------------------------------------
  return (0)
}
