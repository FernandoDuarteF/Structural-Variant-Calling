install.packages("circlize")
library("circlize")
library(rtracklayer)
library(StructuralVariantAnnotation)

palette("Dark 2")

data = generateRandomBed(nc=2, nr = 50)



mants_vcf <- VariantAnnotation::readVcf(system.file("extdata", "representations.vcf", package = "StructuralVariantAnnotation"))

manta_bp <- breakpointRanges(mants_vcf)
manta_be <- breakendRanges(mants_vcf)

manta_bp$svtype

#Why there are two coordinates for every breakend? From Lumpy github: 
#"Lumpy assigns a probability distribution to the breakpoint, reflecting
#the relative uncertainty in the definition of the start and end breakpoints"

#First for breakpoints

manta_bp_bedpe= breakpointgr2bedpe(manta_bp)

manta_bp_bedpe$name <- gsub(":.*", "",
                       gsub("Manta", "",
                            manta_bp_bedpe[,7]))

manta_bp_bedpe_notrans = manta_bp_bedpe[manta_bp_bedpe$name != "BND",]

manta_bp_bedpe_trans = manta_bp_bedpe[manta_bp_bedpe$name == "BND",]

manta_bp_bedpe_notrans$end1 = manta_bp_bedpe_notrans$end2

#Now for breakends

manta_be$score = manta_be$QUAL

manta_be = as.data.frame(manta_be)

manta_be = manta_be[,c(1:3,12)]

#manta_be$seqnames <- "chr2"

circos.initializeWithIdeogram(species = "hg19")

par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))

#Deletions, insertions, duplications and inversions

typeE=c("DEL","INS", "DUP","INV")

colE=c("blue","yellow","green", "red")

for (i in 1:4) {
  bed_list=manta_bp_bedpe_notrans[manta_bp_bedpe_notrans[,7]==typeE[i],]
  bed_list=bed_list[,1:3]
  circos.genomicTrackPlotRegion(bed_list, bg.border = NA, bg.col = add_transparency("grey", 0.8), stack=TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, cex = 0.5, pch = 16, col = add_transparency(i, 0.5), border = i, lwd = "2",...)
  })
}

#Samples in the same track. Might be useful for

#same_track = manta_bp_bedpe_notrans[,c(1:3,7)]

#circos.genomicTrack(same_track, stack = TRUE, 
#                    panel.fun = function(region, value, ...) {
#                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = factor(value$name),...)
#                    })

#Now for breakends (if exists)

circos.genomicTrackPlotRegion(manta_be, bg.border = NA, bg.col = add_transparency("grey", 0.8), stack=TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, cex = 0.5, pch = 16, col = 5, border = 5, lwd = "2",...)
})

#Now for translocation breakpoints

bed_1 = manta_bp_bedpe_trans[,1:3]

bed_1[,2:3] <- lapply(bed_1[,2:3], as.integer)

bed_2 = manta_bp_bedpe_trans[,4:6]

circos.genomicLink(bed_1, bed_2, col = 6, 
                   border = 6, lwd = "2")

circos.clear()

#legend(0.7,1,legend=c("DELETION","INSERTION", "DUPLICATION","INVERSION"),col=unique(manta_bp_bedpe_notrans$name),pch=c(15,15,15,15),cex=0.75,bty='n')

legend("topleft", right,legend=c("DELETION","INSERTION", "DUPLICATION", "INVERSION", "SINGLE BND", "TRANSLOCATION"),col=1:6,pch=c(15,15,15,15,15,NA),cex=0.75,bty='n',lty=c(NA,NA,NA,NA,NA,1),lwd=c(NA,NA,NA,NA,NA,1))

palette("default")
