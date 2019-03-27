#' Read raw sugar and starch data
#'
#' Function to read excel file with raw data of nonstructural carbohydrates
#' @param workDir Path to the directory containing the excel spreadsheets with raw data for the nonstructural carbohydrate assay.
#' @param fileName Character string with the file name of the excel spreadsheet containing the raw data.
#' @param IDs Vector of character strings with sample ID specifications, such as c ('SampleID','BatchID','Tissue','Treatment'). It needs to contain at least one column called 'SampleID'. The identification columns are expected to preceed any other columns.
#' @return The raw data from the nonstructural carbohydrate assay.
#' @import readxl
#' @export
readLabMasterSheet <- function (fileDir, fileName, IDs = c ('SampleID')){

  # Get the filePath to the relevant excel spreadsheet
  #---------------------------------------------------------------------------------------
  fullPath <- paste (fileDir, fileName, sep = "/")

  # Check that at least one ID column has been provide
  #---------------------------------------------------------------------------------------
  nIDs <- length (IDs) # Get number of ID columns
  if (nIDs < 1) stop ('Error: Not a single column name has been provided to identify the samples.')

  # Read the master sheet
  #------------------------------------------------------------------------------------
  types <- c (rep ('text', (nIDs+1)), rep ('date', 3), rep ('numeric', 19), 'text',
              rep ('numeric', 2), 'text', 'numeric','text', rep ('numeric', 6), 'text',
              rep ('numeric', 6), 'text')
  master <- readxl::read_excel (fullPath, col_types = types)  [, 1:(41 + nIDs)]

  # Return the tibble with all necessary information
  #---------------------------------------------------------------------------------------
  return (master)
}
