args <- commandArgs(T)

setwd(args)

file <- dir()

for(i in file){
  
  load(i)
  
}

rm(i)

rm(file)

save.image("input_data.Rdata")