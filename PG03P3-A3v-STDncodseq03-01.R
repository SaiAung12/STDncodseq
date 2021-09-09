# commands to load libraries and datasets
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(vegan)
library(tidyr)
library(dplyr)
library(compositions)
library(ggplot2)
library(DESeq2)
library(dendextend)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(ggplot2)
library(DESeq2)
library(dendextend)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(ggsci)
library(scales)


##### STDnCoda_pipeline

### READING IN OUR DATA  ## need to set this
setwd("/home/sylvester/test/PG03P3-A3v-STDncodseq03-01")

#### load Rdata
load("/home/sylvester/test/PG03P3-A3v-STDncodseq02-01/dada2f230/202106-pg03P3-A3v-16S-rmchloMTunc-Q1f230-STDnCoda01-minprop-001minOc03OTUs1417_allphyloseq.Rdata")

##### set params
fpond = "pg03P3-A3v"
ftiss = "16S-rmchloMTunc"
prop = "Q1f230-STDnCoda01-minprop"
prop1 = "Q1f230-STDnCoda01-cutoffs"

##output directory
outdir = "/home/sylvester/test/PG03P3-A3v-STDncodseq03-01/dada2f230/"

if (file.exists(outdir)){
    setwd(file.path(outdir))
} else {
    dir.create(file.path(outdir))
    setwd(file.path(outdir))
}

testcolk = pal_simpsons("springfield", alpha = 0.75)(16)

wkcolpal = testcolk[c(1, 5, 6,
                     15, 2, 3, 9, 8,
                     11, 10, 16, 13, 12,
                     7, 4, 14)]

names(wkcolpal) = c("n0", "n4", "n6", 0:12)
show_col(wkcolpal)

fminprop = "001"
fminoc =  "03"

##nbrtotal OTU
nbrOTUs = dim(fcount_tab)[1]
nbrOTUs

#[1] 1417

otu_tabx <- abundances(f2.PG03p3Q1_NonChl_physeq)

pdf(paste(fileprefx, "-cntchloMTunc-rarecurve.pdf", sep=""), width=11, height=8)

vegan::rarecurve(t(otu_tabx), step=20, col=filt_sample_info_tab$wk, 
                 lwd=2, ylab="# of ASVs", xlab="# of Sequences",
                 main = paste(fileprefx, " rare curve", sep = "\n"),
                 label = FALSE)
								 
vegan::ordilabel(cbind(rowSums(t(otu_tabx)), vegan::specnumber(t(otu_tabx))), 
                 labels = rownames(t(otu_tabx)), 
                 fill = "NA")
								 
abline(v=(min(rowSums(t(otu_tabx)))))


################################## no chloroplast and MTunc here #######


### remove chlo
fcount_tab.chlo.rowid = which(row.names(fcount_tab) %in% tax_tab.chlo.asvid)
filt_fcount_tab = fcount_tab[-fcount_tab.chlo.rowid , ]
chlor_fcount_tab = fcount_tab[fcount_tab.chlo.rowid , ]

fcount_tab.chloMTunc.rowid = which(row.names(filt_fcount_tab) %in% tax_tab.MTunc.asvid)
MTunc_fcount_tab = filt_fcount_tab[fcount_tab.chloMTunc.rowid , ]
if(length(fcount_tab.chloMTunc.rowid) > 0){
  filt_fcount_tab = filt_fcount_tab[-fcount_tab.chloMTunc.rowid , ]
}else{
  filt_fcount_tab = filt_fcount_tab
}

filt_tax_tab = tax_tab[row.names(filt_fcount_tab), ]

chlor_tax_tab = tax_tab[row.names(chlor_fcount_tab), ]

MTunc_tax_tab = tax_tab[row.names(MTunc_fcount_tab), ]


#reset the number of OTUs after removed
### this is the number of OTUs left in filt_fcount_tab
nbrOTUs = dim(filt_tax_tab)[1]
nbrOTUs

#1327


###########  
## file prefix
fileprefx = paste("202106-", fpond, "-", 
                  ftiss, "-", 
                  prop, "-", fminprop, "minOc", fminoc, "OTUs", nbrOTUs, sep = "")

### removed the chloroplast  and MTunc
fcount_tab = filt_fcount_tab

# replace 0 values with an estimate  ######### f.n0.count_tab
f.n0.count_tab <- cmultRepl(t(fcount_tab), method="CZM", label=0)

# generate the CLR values for plotting later
f.n0.clr <- codaSeq.clr(f.n0.count_tab)   

####Returns a matrix of center log-ratio transformed data with samples by row. 
##Each value is equivalent to log(x/gx) where gx is the geometic mean of the row vector X.
head(t(f.n0.clr), n=2)


###### reconstruct phyloseq again with the raw count and CZM method (replace 0 with an estimated).
#Transform into matrixes otu and tax tables
otu_mat <- as.matrix(fcount_tab)
tax_mat <- tax_tab[row.names(fcount_tab),]


OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(filt_sample_info_tab)  ### change to filter for ponds/tissues

f2.PG03p3Q1_noChloMTunc_physeq <- phyloseq(OTU, TAX, samples)
f2.PG03p3Q1_noChloMTunc_physeq   ### phyloseq with the counts.


f2relCZM.PG03p3Q1_noChloMTunc_physeq = f2.PG03p3Q1_noChloMTunc_physeq

otu_mat <- as.matrix(t(f.n0.count_tab))  ### transposed f.n0

otu_table(f2relCZM.PG03p3Q1_noChloMTunc_physeq) = otu_table(otu_mat, taxa_are_rows = TRUE)

f2relCZM.PG03p3Q1_noChloMTunc_physeq ## phyloseq with estimated 0 count from CZM method

### codaSeq.clr results
f2relCodaCRL.PG03p3Q1_noChloMTunc_physeq = f2.PG03p3Q1_noChloMTunc_physeq

otu_mat <- as.matrix(t(f.n0.clr))  ### transposed f.n0.clr

otu_table(f2relCodaCRL.PG03p3Q1_noChloMTunc_physeq) = otu_table(otu_mat, taxa_are_rows = TRUE)

f2relCodaCRL.PG03p3Q1_noChloMTunc_physeq ## phyloseq with CodaCRL


### microbiome transform cleer
f2rel.microCLR.PG03p3Q1_noChloMTunc_physeq = microbiome::transform(f2.PG03p3Q1_noChloMTunc_physeq, 'clr')
sample_sums(f2rel.microCLR.PG03p3Q1_noChloMTunc_physeq)
sample_sums(f2relCodaCRL.PG03p3Q1_noChloMTunc_physeq)


### save

save(
  fcount_tab,
  f2.PG03p3Q1physeq, 
  f2.PG03p3Q1_noChloMTunc_physeq,
  f2relCZM.PG03p3Q1_noChloMTunc_physeq,
  f2relCodaCRL.PG03p3Q1_noChloMTunc_physeq,
  #f2rel.microCLR.PG03p3Q1_noChloMTunc_physeq,
  filt_count_tab,
  filt_sample_info_tab,
  phyq1PG03A3v,
  tax_tab,
  tax_tab.chlo,
  tax_tab.nochlo,
  tax_tab.MT,
  #tax_tab.rick.unc,
  tax_tab.chlo.asvid, 
  tax_tab.MTunc.asvid,
  tax_tab.MT.asvid,
  #tax_tab.rick.unc.asvid,
  
  file = paste(fileprefx, "_allphyloseq.Rdata", sep = ""))

####
######  test aldex
samples = sample_data(filt_sample_info_tab)

samples$srtest = samples$pdfx  ##should be changed everytime by P Ay

sample_data(f2.PG03p3Q1_noChloMTunc_physeq) = samples

f2.conds = sample_data(f2.PG03p3Q1_noChloMTunc_physeq)[colnames(otu_table(f2.PG03p3Q1_noChloMTunc_physeq))]$srtest

# estimate the distribution of CLR values
f2.x <- aldex.clr(as(otu_table(f2.PG03p3Q1_noChloMTunc_physeq), "matrix"), f2.conds, 
                  mc.samples=1000, verbose=FALSE)
									

# generate the expected CLR value for each OTU
# along with expected value of effect sizes
f2.e <- aldex.effect(f2.x, f2.conds,
                     include.sample.summary=TRUE,
                     verbose=FALSE, useMC = FALSE)

###### ordination
# get a matrix of E(CLR) values
# samples by row
# use for ordination and dendrogram
E.E.clr <- t(f2.e[,grep("rab.sample", colnames(f2.e))])
rownames(E.E.clr) <- gsub("rab.sample.", "", rownames(E.E.clr))
exp <- apply(E.E.clr, 1, function(x) 2^x)
E.clr <- t(apply(exp, 2, function(x) log2(x) - mean(log2(x)) ))


save(
  f2.conds,
  f2.x,
  f2.e,
  E.E.clr,
  exp,
  file = paste(fileprefx, "_aldextest.Rdata", sep = ""))


####
# perform a singular value decomposition
pcx <- prcomp(E.clr)

# plot a PCA biplot
# calculate percent variance explained for the axis labels
pc1 <- round(pcx$sdev[1]^2/sum(pcx$sdev^2),2)
pc2 <- round(pcx$sdev[2]^2/sum(pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")

# colobuf = as(filt_sample_info_tab[names(pcx$x[,1]), "colorbuf"], "vector")
# colobuf = colobuf[[1]]
# namecolbuf = as(filt_sample_info_tab[names(pcx$x[,1]), "extbuf"], "vector")
# names(colobuf) = namecolbuf[[1]]

## samples

muldifsample = data.frame(pc1=pcx$x[,1], pc2=pcx$x[,2])

muldifsample = merge(muldifsample, filt_sample_info_tab, by=0, all.x=TRUE)

row.names(muldifsample) = muldifsample$Row.names

muldifsample = muldifsample[-1]

muldifsample$samid = row.names(muldifsample)

muldifsample$idno = paste(muldifsample$shrimpid, muldifsample$site, muldifsample$filter,  sep="-")


### samples plot
psamp = ggscatter(data = muldifsample, x="pc1", y="pc2", 
                  xlab=xlab, ylab = ylab, size = 3,
                  main = paste(fileprefx, "Nochlo MTunc PCA samples", sep = "\n"),
                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk", 
                  shape = "pdfx",
                  label = "idno", 
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg03), alpha = 0.1, geom = "polygon") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))

ggsave(filename=paste(fileprefx, "-cntchlo-PCA", "_sample.pdf", sep = ""), 
       width=210, height=230, units="mm")


psamp = ggscatter(data = muldifsample, x="pc1", y="pc2", 
                  xlab=xlab, ylab = ylab, size = 3,
                  main = paste(fileprefx, "Nochlo MTunc PCA samples", sep = "\n"),
                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk", 
                  shape = "pdfx",
                  #label = "idno", 
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg03), alpha = 0.1, geom = "polygon") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))

ggsave(filename=paste(fileprefx, "-cntchlo-PCA", "_xlab_sample.pdf", sep = ""), 
       width=210, height=230, units="mm")


### NMDS

sigcol2 = col2rgb(c("grey", "#E85285"))

unique(as(tax_table(f2.PG03p3Q1_noChloMTunc_physeq)[,"Species"], "vector"))

colorsiggenus = c(rgb(sigcol2[1,1], sigcol2[2,1], sigcol2[3,1], max = 255, alpha = 200), 
                  rgb(sigcol2[1,2], sigcol2[2,2], sigcol2[3,2], max = 255, alpha = 155))
#names(colorsiggenus) = unique(as(tax_table(f2.PG03P3Q1_NonChl_physeq)[,"Species"], "vector"))

f2.p01physeq.nmdsord <- ordinate(f2.PG03p3Q1_noChloMTunc_physeq, "NMDS", "bray", trymax = 100)

p2 = plot_ordination(f2.PG03p3Q1_noChloMTunc_physeq, f2.p01physeq.nmdsord,
                     type = "samples",
                     shape = "pd",
                     color = "wk",
                     title=paste(fileprefx, "cntchlo-NonChlMTunc-NMDS-bray", "_sample", sep = ""))
										 
p2 + scale_colour_manual(values =  c(wkcolpal), aesthetics = "colour") + 
  geom_point(size=3) + 
  geom_text(mapping = aes(label = site), size = 5, color=colorsiggenus[1], hjust = 0, vjust = -0.25)  +
  #  scale_colour_manual(palette = wkcolpal) +
  theme_bw()

ggsave(filename=paste(fileprefx, "-cntchlo-NonChl-MTunc-NMDS-bray", "_sample.pdf", sep = ""), 
       width=210, height=230, units="mm")

p2 = plot_ordination(f2.PG03p3Q1_noChloMTunc_physeq, f2.p01physeq.nmdsord,
                     type = "species", #shape = "Genus",
                     color = "Phylum",
                     title=paste(fileprefx, "-cnt-NonChl-NMDS-bray", "_ASVPhylum", sep = ""))
										 
p2 + #scale_colour_manual(values =  c(colowfsg, colorsiggenus), aesthetics = "colour") +
  geom_point(size=3) +
  #geom_text(mapping = aes(label = shrimpno), size = 8, color="red", hjust = 0, vjust = -0.25)  +
  theme_bw()
	
ggsave(filename=paste(fileprefx, "-cntchlo-NonChl-NMDS-bray", "_ASVPhylum.pdf", sep = ""),
       width=840, height=230, units="mm")


p2 = plot_ordination(f2.PG03p3Q1_noChloMTunc_physeq, f2.p01physeq.nmdsord,
                     type = "species", #shape = "Genus",
                     color = "Class",
                     title=paste(fileprefx, "-cnt-NonChl-NMDS-bray", "_ASVClass", sep = ""))
										 
p2 + #scale_colour_manual(values =  c(colowfsg, colorsiggenus), aesthetics = "colour") +
  geom_point(size=3) +
  #geom_text(mapping = aes(label = shrimpno), size = 8, color="red", hjust = 0, vjust = -0.25)  +
  theme_bw()
	
ggsave(filename=paste(fileprefx, "-cntchlo-NonChl-NMDS-bray", "_ASVClass.pdf", sep = ""),
       width=840, height=230, units="mm")


