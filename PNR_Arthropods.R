# ---- Powdermill Nature Preserve Community Data ----


pnr2015 <- read.csv("Data/Processed/PNR2015_InvertebrateCommunity.csv")

pnr2022 <- read.csv("Data/Processed/Powdermill_Community_2022.csv")

traploc <- read.csv("Data/Processed/PNR_PitfallTrapLocations.csv")

pnr15 <- pnr2015[pnr2015$Quadrat >= 41, ]
