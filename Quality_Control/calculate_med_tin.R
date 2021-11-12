##### Calculating median TIN for the GBE samples ########

#Data lives here:
#/Users/amandan/Desktop/Dworkin/background_effects/GBE_TIN


setwd("/Users/amandan/Desktop/Dworkin/background_effects/GBE_TIN")

my_files <- list.files(pattern = ".xls")
my_files

tin_sum <- lapply(my_files, function(i){
  x = read.table(i, header = TRUE )
  x$file = i
  x
})
tin_sum[[1]]
tin_sum = do.call("rbind.data.frame", tin_sum)

#Okay... first step would be to extract the line number, make that its own column
#Then, group by strain and calculate median TIN value per strain

strain <- rep(NA, times = nrow(tin_sum))

for (i in 1:length(tin_sum$file)){
  file.name = tin_sum$file[i]
  strain[i] <- paste(strsplit(file.name, split = "_")[[1]][1],
                     strsplit(file.name, split = "_")[[1]][2],
                     strsplit(file.name, split = "_")[[1]][3],
                     sep = "_")
}

tin_sum$strain <- strain #Now I added strain as a variable

median_tin <- aggregate(tin_sum$TIN, by = list(tin_sum$strain), FUN = median)
colnames(median_tin) <- c("strain", "median_TIN")

#Median TIN of all the data (including dgrp not used in my project)
hist(median_tin$median_TIN, main = "Hist. of median TIN values in sample")
