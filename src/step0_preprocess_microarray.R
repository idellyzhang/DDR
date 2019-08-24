MR <- read.csv("../data/microarray_clean.csv", header = T,check.names = F,row.names = 1)


########MAPPING#######################
G4 <- c(1,2,4,6,7,8,13,14,15,18,19,20,21,22,25,27,28,33,36, 38,39,40,45,46,49,50, 53,54,55,56,57,58,68, 69, 70, 71, 73, 75, 76)
G3 <- c(9, 10, 11, 12, 17, 23, 24,34, 35, 37, 47, 52, 59, 60, 62, 63)
WNT <- c(3, 16, 26, 48, 61, 64, 65, 66)
SHH <- c(5, 29, 30, 31, 32, 42, 43, 44, 51, 72, 74)

G4_m <-MR[, G4]
G3_m <- MR[, G3]
WNT_m <- MR[, WNT]
SHH_m <- MR[, SHH]

complete <- cbind(WNT_m, SHH_m, G3_m, G4_m)

write.csv(complete, "processedTable.csv")