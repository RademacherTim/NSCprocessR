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

      # Create fileName for the pdf
      #----------------------------------------------------------------------------------
      fileName <- paste (standard,".png", sep = "")

      # Open the png
      #----------------------------------------------------------------------------------
      png (file = fileName)

      # Check whether it is a starch or sugar standard
      #----------------------------------------------------------------------------------
      if (standard == 'LCS Oak') {
        dates <- data [['DateOfSugarAnalysis']] [data [['SampleID']] == standard]
        concentrationsOfStds <- data [['ConcentrationSugarPerDW']] [data [['SampleID']] == standard]
      } else if (standard == 'LCS Potato') {
        dates <- data [['DateOfStarchAnalysis']] [data [['SampleID']] == standard]
        concentrationsOfStds <- data [['ConcentrationStarchPerDW']] [data [['SampleID']] == standard]
      }

      # plot the standard over time
      #----------------------------------------------------------------------------------
      plot (x = dates,
            y = concentrationsOfStds,
            main = paste (standard,'over time'),
            xlab = 'time',
            ylab = 'concentration (% dry weight)',
            las = 1,
            col = '#91b9a499',
            pch = 19,
            ylim  = c (0, 1.5*max (concentrationsOfStds)))
      # Add the mean and standard deviation
      #----------------------------------------------------------------------------------
      means <- aggregate (concentrationsOfStds, by = list (dates), FUN = mean)
      sds <- aggregate (concentrationsOfStds, by = list (dates), FUN = sd)
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

      # close graphics device
      #----------------------------------------------------------------------------------
      dev.off ()

    } # End the standards loop

  # Return zero exit status, if it ran smoothly
  #--------------------------------------------------------------------------------------
  return (0)
}