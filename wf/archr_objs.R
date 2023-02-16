library(ArchR)

args = commandArgs(trailingOnly=TRUE)
print('hi')

project_name <- args[1]
genome <- args[2]
threads <- as.integer(args[3])
tile_size <- as.integer(args[4])
min_TSS <- as.numeric(args[5])
min_frags <- as.integer(args[6])

runs = strsplit(args[7:length(args)], ',')
inputs = c()
for (run in runs) {inputs[run[1]] = run[2]}

addArchRGenome(genome)
addArchRThreads(threads=threads)

ArrowFiles <- createArrowFiles(
   inputFiles = inputs,
   sampleNames = names(inputs),
   minTSS = min_TSS,
   minFrags = min_frags,
   maxFrags = 1e+07,
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   offsetPlus = 0,
   offsetMinus = 0,
   TileMatParams = list(tileSize=tile_size)
)

proj <- ArchRProject(
  ArrowFiles=ArrowFiles, 
  outputDirectory= paste0(project_name, "_ArchRProject")
)

# Add an additional Conditions column
for (run in runs) {proj$Condition[proj$Sample==run[1]] <- run[3]}

saveArchRProject(ArchRProj = proj)

