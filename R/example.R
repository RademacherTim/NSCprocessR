# Example script


# Read data from exampleData.xlsx spreadsheet
rawData <- readRawNSCData (workDir  = getwd (),
                           fileName = 'exampleData.xlsx',
                           IDs      = c ('RCLabNumber','SampleID','Tissue','BatchID'))

