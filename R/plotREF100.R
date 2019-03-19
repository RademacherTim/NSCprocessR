#' Plot REF100 values over time
#'
#' Function to plot the REF100 values over time.
#' @param data Tibble with the processed data using the NSCprocessR::processNSC function. It can also be from the already processed data in the master file read in using NSCprocessR::readLabMasterSheet () function.
#' @return pdf file with a graph of the REF100 value over time.
#' @export
plotREF100 <- function (data) {

  # Create fileName for the pdf
  #------------------------------------------------------------------------------------
  fileName <- paste ("REF100.png", sep = "")

  # Open the png
  #------------------------------------------------------------------------------------
  png (file = fileName)

    # Check whether it is a strach or sugar standard
    #----------------------------------------------------------------------------------
    dates <- data [['DateOfSugarAnalysis']] [data [['SampleID']] == 'REF100' |
                                             data [['SampleID']] == 'REF100/200']
    concentrationsOfStds <- data [['ConcentrationSugar']] [
                                  data [['SampleID']] == 'REF100' |
                                  data [['SampleID']] == 'REF100/200']

    # plot the REF100 values over time
    #----------------------------------------------------------------------------------
    plot (x = dates,
          y = concentrationsOfStds,
          main = 'REF100 over time',
          xlab = 'time',
          ylab = 'REF 100 concentration (% dry weight)',
          las = 1,
          col = '#91b9a499',
          pch = 19,
          ylim  = c (70, 130))
    # Add 100 line and 5 and 10% error shades
    #----------------------------------------------------------------------------------
    rect (xleft = min (dates) - 0.1 * (max (dates) - min (dates)),
          xright = max (dates) + 0.1 * (max (dates) - min (dates)),
          ybottom = 95,
          ytop = 105,
          col = '#bbbbbb66',
          lty = 0)
    rect (xleft = min (dates) - 0.1 * (max (dates) - min (dates)),
          xright = max (dates) + 0.1 * (max (dates) - min (dates)),
          ybottom = 90,
          ytop = 110,
          col = '#bbbbbb66',
          lty = 0)
    abline (h = 100,
            col = '#666666')
    points (x = dates,
            y = concentrationsOfStds,
            col = '#91b9a499',
            pch = 19)

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
  #--------------------------------------------------------------------------------------
  dev.off ()

  # Return zero exit status, if it ran smoothly
  #--------------------------------------------------------------------------------------
  return (0)
}
# To-do list:
# - Need to test this function.
