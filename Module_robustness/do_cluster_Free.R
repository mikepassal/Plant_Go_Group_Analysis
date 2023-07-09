h5ls('/data/CoCoCoNet/networks/arabidopsis_prioAggNet.hdf5')
bulk_network =  h5read("/data/CoCoCoNet/networks/arabidopsis_prioAggNet.hdf5","agg")
cols = h5read("/data/CoCoCoNet/networks/arabidopsis_prioAggNet.hdf5","col")
rows = h5read("/data/CoCoCoNet/networks/arabidopsis_prioAggNet.hdf5","row")
library(WGCNA)
options(stringsAsFactors = FALSE);
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(bulk_network, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-o
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")