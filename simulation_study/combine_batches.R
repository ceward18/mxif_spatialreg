################################################################################
# combine batch runs
################################################################################

outputFolder <- 'output'
resultsFolder <- 'results'

outputFiles <- sort(list.files(paste0('./', outputFolder)))

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles[1]))

for (i in 2:length(outputFiles)) {
    if (i %% 100 == 0) print(i)
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}

saveRDS(res_all, paste0('./', resultsFolder, '/res_all.rds'))





