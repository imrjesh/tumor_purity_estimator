install.packages("estimate")
install.packages("estimate", repos="http://R-Forge.R-project.org")
###load the library
library(estimate)

#####set path to the file, my case its in download folder
OvarianCancerExpr <- "/Users/kumarr9/Desktop/rajesh_projects/tumor_purity/DSP_Q3_stroma_PD2.txt"

filterCommonGenes(input.f=OvarianCancerExpr, output.f="OV_10412genes.gct", id="GeneSymbol")

estimateScore("OV_10412genes.gct", "OV_estimate_score.gct", platform="affymetrix")

#in.file <- system.file("OV_estimate_score.gct", package="estimate")

in.file <- "OV_estimate_score.gct"
plotPurity(in.file)

