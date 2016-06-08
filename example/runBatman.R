library(batman)

# set to the directory containing this script
setwd('/home/rstudio/example')

# Configure the parameters of Batman from the batmanOptions.txt file inside
# the runBATMAN/BatmanInput directory

# Multiple spectra data
bm <- batman(BrukerDataDir='brucker_input', figBatmanFit=FALSE)

# Plot spectral fit and inferred relative concentration
plotBatmanFit(bm, showPlot=FALSE)
plotRelCon(bm, showPlot=FALSE)
plotMetaFit(bm, showPlot=FALSE)

# plot invididual metabolite fit
plotMetaFit(bm, metaName='Alanine', showPlot=FALSE, from=1.46, to=3.76)
plotMetaFit(bm, metaName='Ethanol', showPlot=FALSE)
plotMetaFit(bm, metaName='Glucose', showPlot=FALSE)
plotMetaFit(bm, metaName='Glycine', showPlot=FALSE)
plotMetaFit(bm, metaName='Histidine', showPlot=FALSE)
plotMetaFit(bm, metaName='Isoleucine', showPlot=FALSE)
plotMetaFit(bm, metaName='Lactate', showPlot=FALSE)
plotMetaFit(bm, metaName='Leucine', showPlot=FALSE)
plotMetaFit(bm, metaName='Phenylalanine', showPlot=FALSE)
plotMetaFit(bm, metaName='Tyrosine', showPlot=FALSE)
plotMetaFit(bm, metaName='Valine', showPlot=FALSE)
