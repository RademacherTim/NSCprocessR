#======================================================================================
# plot starch recovery fraction
#--------------------------------------------------------------------------------------

# read the lab master sheet and plot standards over time
#--------------------------------------------------------------------------------------
master <- readLabMasterSheet (fileDir = getwd (), fileName = 'RCLabNSCMasterSheet.xlsx', IDs = c ('SampleID', 'Tissue', 'BatchID', 'SampleLocation'))

# Compile a list of unique batches
#--------------------------------------------------------------------------------------
batches <- unique (master [['BatchID']])

# Compile a list of unique batches and dates for sugar extractions
#--------------------------------------------------------------------------------------
for (batch in batches) {
  dates <- unique (master [['DateOfStarchAnalysis']] [master [['BatchID']] == batch])
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

# plot potato starch recovery over time
#--------------------------------------------------------------------------------------
plot (x = master$DateOfStarchAnalysis [1],
      y = master$MeanStarchRecovery   [1],
      ylab = 'mean batch starch recovery fraction (%)',
      xlab = 'date',
      col = 'white',
      xlim = c (as.POSIXct ("2019-01-01", format = "%Y-%m-%d"),
                as.POSIXct ("2019-04-01", format = "%Y-%m-%d")),
      ylim = c (60, 110),
      main = 'Potato starch recovery fraction over time')

for (extraction in 1:(dim (extractionsStarch) [1])) {

  analysisDate <- extractionsStarch$date [extraction]
  batch <- extractionsStarch$batch       [extraction]

  condition <- master$DateOfStarchAnalysis == analysisDate &
               master$BatchID              == batch

  points (x = unique (master$DateOfStarchAnalysis [condition]),
          y = unique (master$MeanStarchRecovery   [condition]),
          col = '#91b9a499',
          pch = 19)
}

abline (h = 100,
        col = '#66666666')
#======================================================================================
