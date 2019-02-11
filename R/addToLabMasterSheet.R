#' Create or append lab master sheet with the processed samples
#'
#' Function to create or add to the lab master sheet.
#' @param data Tibble with the data that is to be added to the master sheet. The data was created using the processNSC function.
#' @return
#' @import openxlsx
#' @import readxl
#' @export
addToLabMasterSheet <- function (data, fileDir, fileName) {

  # Paste directory name and file name to get a full path
  fullPath <- paste (fileDir,'/',fileName, sep = '')

  # Check whether master sheet exists for lab and create it
  if (!file.exists (fullPath)) {
    openxlsx::write.xlsx (x = data, # TTR Need to go over these to check them! Is all still old from java dependent package.
                          file = fullPath,
                          sheetName = "Master",
                          col.names = TRUE,
                          row.names = TRUE,
                          append = FALSE) # If not create it with the data
  } else { # If it does exist

    # Read the master sheet
    #------------------------------------------------------------------------------------
    types <- c (rep ('text', (nIDs+1)), rep ('date', 3), rep ('numeric', 11)) # TTR Need to adjust col_types to include processed columns!!!
    data <- readxl::read_excel (filePath, col_types = types)

    # Check whether data is already in the data sheet
    #------------------------------------------------------------------------------------

  }



  #
  return ()
}
