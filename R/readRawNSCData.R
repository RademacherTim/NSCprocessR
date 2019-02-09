#' Read raw sugar and starch data
#'
#' Function to read excel file with raw data of nonstructural carbohydrates
#' @param workDir Path to the directory containing the excel spreadsheets with raw data for the nonstructural carbohydrate assay.
#' @param fileName Character string with the file name of the excel spreadsheet containing the raw data.
#' @param IDs Vector of character strings with sample ID specifications, such as c ('SampleID','BatchID','Tissue','Treatment'). It needs to contain at least one column called 'SampleID'. The identification columns are expected to preceed any other columns.
#' @return The raw data from the nonstructural carbohydrate assay.
#' @import readxl
#' @export
readRawNSCData <- function (workDir, fileName, IDs = c ('SampleID')){

 # Get the filePath to the relevant excel spreadsheet
 #---------------------------------------------------------------------------------------
 filePath <- paste (workDir, fileName, sep = "/")

 # Check that at least one ID column has been provide
 #---------------------------------------------------------------------------------------
 nIDs <- length (IDs) # Get number of ID columns
 if (nIDs < 1) stop ('Error: Not a single column name has been provided to identify the samples.')

 # Read the relevant columns in the excel spreadsheet and organise the tibble
 #---------------------------------------------------------------------------------------
 if (file.exists (filePath)){
   types <- c (rep ('text', (nIDs+1)), rep ('date', 3), rep ('numeric', 11))
   temp <- readxl::read_excel (filePath, col_types = types) [, 1:(15 + nIDs)] # Only read the relevant columns
 } else {
   stop ("Error: The file your are trying to read does not exist.")
 }

 # Check that ID columns exist and are properly named
 #---------------------------------------------------------------------------------------
 for (i in 1:nIDs) {
   if (length (which (names (temp) == IDs [i])) > 1) {
      stop (paste ('There are multiple ID column with the name: ',IDs [i],'.', sep = ''))
   } else if (length (which (names (temp) == IDs [i])) < 1) {
      stop (paste ('The ID column ',IDs [i], 'could not be found. Maybe it was miss-spelled?', sep = ''))
   }
 }

 # Check that sample location is provided for the sample collection and analysis
 #---------------------------------------------------------------------------------------
 if (length (which (names (temp) == 'SampleLocation')) > 1) {
   stop ('There are multiple columns with the sample location.')
 } else if (length (which (names (temp) == 'SampleLocation')) < 1) {
   warning ('There is no column with the samnple location.')
 }

 # Check that dates are provided for the sample collection and analysis
 #---------------------------------------------------------------------------------------
 for (date in c ('DateOfSampleCollection','DateOfSugarAnalysis','DateOfStarchAnalysis')) {
   if (length (which (names (temp) == date)) > 1) {
     stop (paste ('There are multiple columns with the ',date, sep = ''))
   } else if (length (which (names (temp) == date)) < 1) {
     warning (paste ('There is no column with the ',date, sep = ''))
   }
 }

 # Check that mass columns are provided
 #---------------------------------------------------------------------------------------
 for (mass in c ('MassOfEmptyTube','MassOfTubeAndSample')) {
    if (length (which (names (temp) == mass)) > 1) {
       stop (paste ('There are multiple columns with data for the ', mass, sep = ''))
    } else if (length (which (names (temp) == mass)) < 1) {
       warning (paste ('There is no column with data for the ', mass, sep = ''))
    }
 }

 # Check that absorbance columns are provided
 #---------------------------------------------------------------------------------------
 for (absorbance in c ('Absorbance490_1','Absorbance490_2','Absorbance490_Blank','Absorbance525_1','Absorbance525_2')) {
    if (length (which (names (temp) == absorbance)) > 1) {
       stop (paste ('There are multiple columns with data for the ', absorbance, sep = ''))
    } else if (length (which (names (temp) == absorbance)) < 1) {
       warning (paste ('There is no column with data for the ', absorbance, sep = ''))
    }
 }

 # Check that dillution columns are provided
 #---------------------------------------------------------------------------------------
 for (dilution in c ('DilutionFactorSugar','VolumeSugar','DilutionFactorStarch','VolumeStarch')) {
    if (length (which (names (temp) == dilution)) > 1) {
       stop (paste ('There are multiple columns with data for the ', dilution, sep = ''))
    } else if (length (which (names (temp) == dilution)) < 1) {
       warning (paste ('There is no column with data for the ', dilution, sep = ''))
    }
 }

 # Return the tibble with all necessary information
 #---------------------------------------------------------------------------------------
 return (temp)
}
