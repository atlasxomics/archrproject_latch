library(ArchR)

args = commandArgs(trailingOnly=TRUE)

inputFile <- args[1]
run_id <- args[2]
genome <- args[3]
threads <- as.integer(args[4])
tile_size <- as.integer(args[5])
min_TSS <- as.numeric(args[6])
min_frags <- as.integer(args[7])

addArchRGenome(genome)
addArchRThreads(threads=threads)

ArrowFiles <- createArrowFiles(
   inputFiles = inputFile,
   sampleNames = run_id,
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
  outputDirectory= paste0(run_id, "_ArchRProject")
)

