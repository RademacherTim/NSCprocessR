# Example script

# load dependencies
library ('lubridate')

# Read data from exampleData.xlsx spreadsheet
rawData  <- readRawNSCData (fileDir  = '../data/NSCs/',
                            fileName = 'RCLab_2019_02_09_HF_Leaves_Exp2017.xlsx',
                            IDs      = c ('RCLabNumber','SampleID','Tissue','BatchID'))
rawData1 <- readRawNSCData (fileDir  = '../data/NSCs/',
                            fileName = 'RCLab_2019_01_26_HF_Leaves_Exp2017.xlsx',
                            IDs      = c ('RCLabNumber','SampleID','Tissue','BatchID'))
rawData2 <- readRawNSCData (fileDir  = '../data/NSCs/',
                            fileName = 'RCLab_2019_03_22_HF_Roots_Exp2017.xlsx',
                            IDs      = c ('RCLabNumber','SampleID','Tissue','BatchID'))

# Process the raw data
processedData  <- processNSCs (rawData = rawData,
                               cvLimitSample = 0.25,
                               cvLimitTube = 0.05,
                               forceIntercept = FALSE,
                               LCS = 'Oak')
processedData1 <- processNSCs (rawData = rawData1,
                               cvLimitSample = 0.25,
                               cvLimitTube = 0.05,
                               forceIntercept = FALSE,
                               LCS = 'Oak')
processedData2 <- processNSCs (rawData = rawData2,
                               cvLimitSample = 0.25,
                               cvLimitTube = 0.05,
                               forceIntercept = FALSE,
                               LCS = 'Oak')

# Produce pdf file with calibration curves
res <- plotCalibrationCurves (data = processedData)
if (res != 0) print ('Error: plotCalibrationCurves () did not work!')
res <- plotCalibrationCurves (data = processedData1)
if (res != 0) print ('Error: plotCalibrationCurves () did not work!')
res <- plotCalibrationCurves (data = processedData2)
if (res != 0) print ('Error: plotCalibrationCurves () did not work!')

# Add both datasets to the lab master sheet
res <- addToLabMasterSheet (data = processedData,
                            fileDir = getwd (),
                            fileName = 'RCLabNSCMasterSheet.xlsx',
                            IDs = c ('RCLabNumber','SampleID','Tissue','BatchID'))
if (res != 0) print ('Error: addToLabMasterSheet () did not work!')
res <- addToLabMasterSheet (processedData1,
                            fileDir = getwd (),
                            fileName = 'RCLabNSCMasterSheet.xlsx',
                            IDs = c ('RCLabNumber','SampleID','Tissue','BatchID'))
if (res != 0) print ('Error: addToLabMasterSheet () did not work!')
res <- addToLabMasterSheet (processedData2,
                            fileDir = getwd (),
                            fileName = 'RCLabNSCMasterSheet.xlsx',
                            IDs = c ('RCLabNumber','SampleID','Tissue','BatchID'))
if (res != 0) print ('Error: addToLabMasterSheet () did not work!')

# Add columns for study, treeID, treatment to tibble
processedData [['Study']] <- NA
processedData [['Tree']] <- NA
processedData [['Treatment']] <- NA
processedData1 [['Study']] <- NA
processedData1 [['Tree']] <- NA
processedData1 [['Treatment']] <- NA
processedData2 [['Study']] <- NA
processedData2 [['Tree']] <- NA
processedData2 [['Treatment']] <- NA

# Function to extract study, treeID and treatment
for (i in 1:length (processedData [['SampleID']])) {
  processedData [['Study']] [i] <- unlist (strsplit (processedData [['SampleID']] [i], '-')) [1]
  string <- unlist (strsplit (processedData [['SampleID']] [i], '-')) [2]
  processedData [['Tree']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [1])
  processedData [['Treatment']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [2])
}
for (i in 1:length (processedData1 [['SampleID']])) {
  processedData1 [['Study']] [i] <- unlist (strsplit (processedData1 [['SampleID']] [i], '-')) [1]
  string <- unlist (strsplit (processedData1 [['SampleID']] [i], '-')) [2]
  processedData1 [['Tree']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [1])
  processedData1 [['Treatment']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [2])
}
for (i in 1:length (processedData2 [['SampleID']])) {
  processedData2 [['Study']] [i] <- unlist (strsplit (processedData2 [['SampleID']] [i], '-')) [1]
  string <- unlist (strsplit (processedData2 [['SampleID']] [i], '-')) [2]
  processedData2 [['Tree']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [1])
  processedData2 [['Treatment']] [i] <- as.numeric (strsplit (x = string, split = '.', fixed = T) [[1]] [2])
}

# Make scatterplot of all sugar concentrations
par (mar = c (5, 5, 1, 1))
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 1
plot (y = processedData [['ConcentrationSugarPerDW']] [condition],
      x = processedData [['DateOfSampleCollection']] [condition] - 3*86400,
      las = 1,
      pch = 19,
      col = '#8dd3c799',
      xlim = c (as.POSIXct ('2017-06-15', format = '%Y-%m-%d'),
                as.POSIXct ('2017-11-15', format = '%Y-%m-%d')),
      xlab = 'date',
      ylab = 'sugar concentration (% dry weight)',
      ylim = c (0, 12))
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 2
points (y = processedData [['ConcentrationSugarPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition],
        pch = 19,
        col = '#ffffb399')
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 3
points (y = processedData [['ConcentrationSugarPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition] + 3*86400,
        pch = 19,
        col = '#bebada99')
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 4
points (y = processedData [['ConcentrationSugarPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition] + 6*86400,
        pch = 19,
        col = '#fb807299')
legend (x = as.POSIXct ('2017-06-15', format = '%Y-%m-%d'),
        y = 12.4,
        box.lty = 0,
        pch = 19,
        col = c ('#8dd3c799', '#ffffb399', '#bebada99', '#fb807299'),
        legend = c ('control','girdled','compression','double compressed'),
        bg = 'transparent',
        cex = 0.8)


# Make scatterplot of all Starch concentrations
par (mar = c (5, 5, 1, 1))
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 1
plot (y = processedData [['ConcentrationStarchPerDW']] [condition],
      x = processedData [['DateOfSampleCollection']] [condition] - 3*86400,
      las = 1,
      pch = 19,
      col = '#8dd3c799',
      xlim = c (as.POSIXct ('2017-06-15', format = '%Y-%m-%d'),
                as.POSIXct ('2017-11-15', format = '%Y-%m-%d')),
      xlab = 'date',
      ylab = 'starch concentration (% dry weight)',
      ylim = c (0.0, 8.0))
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 2
points (y = processedData [['ConcentrationStarchPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition],
        pch = 19,
        col = '#ffffb399')
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 3
points (y = processedData [['ConcentrationStarchPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition] + 3*86400,
        pch = 19,
        col = '#bebada99')
condition <- processedData [['Study']] == 'Exp' & processedData [['Treatment']] == 4
points (y = processedData [['ConcentrationStarchPerDW']] [condition],
        x = processedData [['DateOfSampleCollection']] [condition] + 6*86400,
        pch = 19,
        col = '#fb807299')
legend (x = as.POSIXct ('2017-10-01', format = '%Y-%m-%d'),
        y = 8,
        box.lty = 0,
        pch = 19,
        col = c ('#8dd3c799', '#ffffb399', '#bebada99', '#fb807299'),
        legend = c ('control','girdled','compression','double compressed'),
        bg = 'transparent',
        cex = 0.8)

# Scatterplot of July and Augsut repeat measruements
condition  <- processedData  [['Study']] == 'Exp' &
             (month (processedData  [['DateOfSampleCollection']]) == 7 |
              month (processedData  [['DateOfSampleCollection']]) == 8)
plot (x = processedData1 [['ConcentrationSugarPerDW']] [condition],
      y = processedData  [['ConcentrationSugarPerDW']] [condition],
      col = 'white',
      xlim = c (0, 11),
      ylim = c (0, 11),
      las = 1,
      xlab = 'sugar concentration [% dry weight] in January',
      ylab = 'sugar concentration [% dry weight] in February')
for (mon in 7:8) {
  if (mon == 7) {
    colour = '#91b9a4aa'
  } else {
    colour = '#DC143Caa'
  }
  for (i in 1:41) {
    condition  <- processedData  [['Study']] == 'Exp' &
                  month (processedData  [['DateOfSampleCollection']]) == mon &
                  processedData [['Tree']] == i
    condition1 <- processedData1  [['Study']] == 'Exp' &
                  month (processedData1  [['DateOfSampleCollection']]) == mon &
                  processedData1 [['Tree']] == i
    if (sum (condition) != sum (condition1) | sum (condition) == 0) next
    points (x = processedData [['ConcentrationSugarPerDW']] [condition],
            y = processedData1 [['ConcentrationSugarPerDW']] [condition1])
    points (x = processedData [['ConcentrationSugarPerDW']] [condition],
            y = processedData1 [['ConcentrationSugarPerDW']] [condition1],
            pch = 19,
            col = colour)
    if (i == 1) {
      adata <- tibble (jan = processedData  [['ConcentrationSugarPerDW']] [condition],
                       feb = processedData1 [['ConcentrationSugarPerDW']] [condition1])
    } else {
      adata <- add_row (adata, jan = processedData  [['ConcentrationSugarPerDW']] [condition],
                               feb = processedData1 [['ConcentrationSugarPerDW']] [condition1])
    }
  }
}
abline (0, 1, col = 'grey', lty = 2, lwd = 2)
legend (x = 0,
        y = 10,
        col = c ('#91b9a4aa','#DC143Caa'),
        box.lty = 0,
        bg = 'transparent',
        legend = c ('july', 'august'),
        pch = 19)

fit <- lm (adata$jan ~ adata$feb)
summary (fit)
round (summary (fit)$r.squared, 2)
text (x = 7,
      y = 1,
      labels = expression (paste (R^2,'= 0.62')))
abline (fit,
        col = '#DC143Caa',
        lwd = 2)

# plot root sugar concentrations
par (mar = c (5, 5, 1, 1))
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 1
plot (y = processedData2 [['ConcentrationSugarPerDW']] [condition],
      x = processedData2 [['DateOfSampleCollection']] [condition] - 3*86400,
      las = 1,
      pch = 19,
      col = '#8dd3c799',
      xlim = c (as.POSIXct ('2017-06-15', format = '%Y-%m-%d'),
                as.POSIXct ('2017-11-15', format = '%Y-%m-%d')),
      xlab = 'date',
      ylab = 'root sugar concentration (% dry weight)',
      ylim = c (0.0, 3.2))
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 2
points (y = processedData2 [['ConcentrationSugarPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition],
        pch = 19,
        col = '#ffffb399')
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 3
points (y = processedData2 [['ConcentrationSugarPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition] + 3*86400,
        pch = 19,
        col = '#bebada99')
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 4
points (y = processedData2 [['ConcentrationSugarPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition] + 6*86400,
        pch = 19,
        col = '#fb807299')
legend (x = as.POSIXct ('2017-07-01', format = '%Y-%m-%d'),
        y = 3.0,
        box.lty = 0,
        pch = 19,
        col = c ('#8dd3c799', '#ffffb399', '#bebada99', '#fb807299'),
        legend = c ('control','girdled','compression','double compressed'),
        bg = 'transparent')

# TTR There seems to be an issue with index c (196, 375) which are the samples
# Exp-10.3-R-Oct_2017 and LCS Potato in batch 4 and 6 respectively.

# plot root starch concentrations
par (mar = c (5, 5, 1, 1))
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 1
plot (y = processedData2 [['ConcentrationStarchPerDW']] [condition],
      x = processedData2 [['DateOfSampleCollection']] [condition] - 3*86400,
      las = 1,
      pch = 19,
      col = '#8dd3c799',
      xlim = c (as.POSIXct ('2017-06-15', format = '%Y-%m-%d'),
                as.POSIXct ('2017-11-15', format = '%Y-%m-%d')),
      xlab = 'date',
      ylab = 'root starch concentration (% dry weight)',
      ylim = c (0, 3.0))
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 2
points (y = processedData2 [['ConcentrationStarchPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition],
        pch = 19,
        col = '#ffffb399')
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 3
points (y = processedData2 [['ConcentrationStarchPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition] + 3*86400,
        pch = 19,
        col = '#bebada99')
condition <- processedData2 [['Study']] == 'Exp' & processedData2 [['Treatment']] == 4
points (y = processedData2 [['ConcentrationStarchPerDW']] [condition],
        x = processedData2 [['DateOfSampleCollection']] [condition] + 6*86400,
        pch = 19,
        col = '#fb807299')
legend (x = as.POSIXct ('2017-07-01', format = '%Y-%m-%d'),
        y = 3.0,
        box.lty = 0,
        pch = 19,
        col = c ('#8dd3c799', '#ffffb399', '#bebada99', '#fb807299'),
        legend = c ('control','girdled','compression','double compressed'),
        bg = 'transparent')

# There seems to be an issue with the indices c (375, 376, 471, 473:479),
# as they give ridiculously high values. They are from LCS Potato from batch 6
# and all samples except for one duplicate in the re-run batch 12.
