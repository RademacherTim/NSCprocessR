# read the lab master sheet and plot standards over time
master <- readLabMasterSheet (fileDir = getwd (), fileName = 'RCLabNSCMasterSheet.xlsx', IDs = c ('SampleID', 'Tissue', 'BatchID', 'SampleLocation'))

# plot standards over time
plotLabControlStandards (master)

# plot REF100 values over time
plotREF100 (master)

# plot sugar calibration curves over time
plot (1,
      main = 'calibration curves for soluble sugars over time',
      las = 1,
      col = 'white',
      xlim = c (0, 2.5),
      ylim = c (0, 250),
      xlab = 'absorbance at 490 nm',
      ylab = 'glucose equivalent (mg / ml)')
mapply (function (x, y, z) {
           abline (x, y, col = z)},
           unique (master [['InterceptSugar']]),
           unique (master [['SlopeSugar']]),
           c ('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))

#plot strach calibration curves over time
plot (1,
      main = 'calibration curves for starch over time',
      las = 1,
      col = 'white',
      xlim = c (0, 1.5),
      ylim = c (0, 250),
      xlab = 'absorbance at 525 nm',
      ylab = 'glucose equivalent (mg / ml)')
mapply (function (x, y, z) {
  abline (x, y, col = z)},
  unique (master [['InterceptStarch']]),
  unique (master [['SlopeStarch']]),
  c ('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))

