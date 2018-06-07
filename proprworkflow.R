library(propr)
library(ALDEx2)

#read in your table and check to see format, taxonomy columns are not supported in the clr function
m <- read.table("ebimsseq.tsv")

#you may want to eliminate sequence variants with low counts using the below line of code (set as 10 for this example), but it is not required.
#You also want the counts table to not contain any of the taxonomy information for calculating the clr, so select the columns without them
m.n0 <- m[,2:ncol(m)]

#now take the sequence table (full or reduced version), and enter it into the clr function--selecting the rows with read counts. for this function "conds" doesn't affect the calculations
m.clr <- aldex.clr(m.n0, conds = rep("x", ncol(m.n0)))

#now that you have your clr values, you can use the aldex2propr function to generate the propr object, which containsL: the counts for the features across the samples[samples,features], the clr values for the features across the samples[samples,features], and the matrix of the rho values [features,features]
m.propr <- aldex2propr(m.clr)

#this function is to choose your rho cutoff of interest and subset the pairs of features satisfying the cutoff into the @pairs portion of the propr object
m.propr <- m.propr[">", 0.7]

#to reduce the size of the object to contain only the features (aka sequence variants) which result in pairs, use the simplify function
m.propr.simp <- simplify(m.propr)

#since the @pairs is populated with indices of the matrix, instead of the coordinates (which are needed to get the identity of the features, I wrote a function to generate the coordinates)
rhopairs <- data.frame()
coords <- for(i in m.propr.simp@pairs) {
     if(i %% ncol(m.propr.simp@logratio) == 0) {
         rhopairs <- rbind(rhopairs, c(ncol(m.propr.simp@logratio), i %/% ncol(m.propr.simp@logratio)))
     }
     else {
         rhopairs <- rbind(rhopairs, c(i %% ncol(m.propr.simp@logratio) , i %/% ncol(m.propr.simp@logratio) +1))
     }
}

#calculate the mean clr values for the features
clrmeans <- colMeans(m.propr.simp@logratio)
clrmeans <- as.vector(clrmeans)

#calculate the differences in the mean clr values for the pairs
clrdiffer <- function(x) {
  clrmeans[rhopairs[,1]] - clrmeans[rhopairs[,2]]
}
clrdiff <- sapply(1, clrdiffer)

#generate a numeric vector of the row names/numbers of the features that generated the pairs
# if using the dada2 count table derived from the workflow, then as.numeric() can be removed
taxpairs <- as.numeric(colnames(m.propr.simp@logratio))

#use the row names as a call for the taxonomy column in the original original count table
taxa <- m[taxpairs,1]
taxa <- as.character(taxa)
final1 <- cbind(colnames(m.propr.simp@logratio[,rhopairs[,1]]), colnames(m.propr.simp@logratio[,rhopairs[,2]]), format(round(m.propr.simp@matrix[m.propr.simp@pairs], 4), nsamll=4), format(round(clrdiff, 4), nsmall=4), taxa[rhopairs[,1]], taxa[rhopairs[,2]])
colnames(final1) <- c("OTU/ASV_1", "OTU/ASV_2", "Rho_Value", "CLR_Difference", "Taxonomy_of_OTU/ASV_1", "Taxonomy_of_OTU/ASV 2")
View(final1)
#inspect table, looking for clr differences around 5 to 6

### WORK IN PROGRESS--removing the OTUs/ASVs which are technical variants from the original count table
# this chunk is designed for the EBI pipeline which has row numbers, not names, in its count table
#obtain rows with clr difference greater than 5 (current arbitrary cutoff)
final5 <- final1[which(abs(as.numeric(final1[,4])) > 5),]
#Obtaining row numbers for the technical variants
m.neg <- as.numeric(final5[which(as.numeric(final5[,4]) < 0), 1])
m.pos <- as.numeric(final5[which(as.numeric(final5[,4]) > 0), 2])
#concatenating and removing redundant row numbers
m.remove <- unique(sort(c(m.pos, m.neg)))
#remove rows from the count table which are technical variants
m.reduced <- m[-m.remove,]
