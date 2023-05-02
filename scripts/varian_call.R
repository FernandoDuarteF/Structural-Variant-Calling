BiocManager::install("StructuralVariantAnnotation")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("ggbio")

library(ggbio)
library(StructuralVariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


browseVignettes("StructuralVariantAnnotation")

vcf = VariantAnnotation::readVcf("SV.vcf", "hg19")

vcf[, 1:3]

gt = geno(vcf)

header(vcf)

remove(vcf1)


gr = breakendRanges(vcf)

head(gt$AF)

head(gr)


help(geno)

################

vcf = VariantAnnotation::readVcf("diploid_sed.vcf")

gr = breakpointRanges(vcf)

partner(gr)


q=GRanges(seqnames="8",
          ranges=IRanges(start = 1, end = 146364022))

chr8 = subsetByOverlaps(biovizBase::hg19sub, q)

seqlevels(chr8) <- "8"

#chr8 <- chr8[seqnames(chr8) %in% seqlevels(gr)]

chr8

mcols(gr)$to.gr <- granges(partner(gr))

p <- ggbio() +
  circle(gr, geom="link", linked.to="to.gr", aes(color = svtype)) +
  circle(chr8, geom='ideo', fill='gray70') +
  circle(chr8, geom='scale', size=2) +
  circle(chr8, geom='text', aes(label=seqnames), vjust=0, size=3)
p
