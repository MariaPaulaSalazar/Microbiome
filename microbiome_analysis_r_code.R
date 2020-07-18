library(dada2)

#cargar las secuencias fastaq

path <- "~/Universidad/Tesis/Secuencias_fasta/Heliconius_desc/"
list.files(path)

# Separar las secuencias forward y reverse

fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Subdirectorio de filtros
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
head(out)


#Calcular errores
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

#Inferencia
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Combinar F R y los filtros

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspeccionar distribución
table(nchar(getSequences(seqtab)))

#Quitar quimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Trackear las lecturas
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
## If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Taxonomia
taxa <- assignTaxonomy(seqtab.nochim, "~/Universidad/Tesis/Secuencias_fasta/Heliconius_desc/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#PHYLOSEQ
## Instalar el programa
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}

BiocManager::install("phyloseq")

## Cargar paquetes para phyloseq y gráficas

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(mergers)

##Gráficas

#Paleta de colores (20 colores)
library(RColorBrewer)

palette = c("#cec88c", "#0164e3", "#e3be00", "#9f008a", "#27e273", "#ff65be", "#487300", "#f890ff","#7d7600", "#faabfc", "#3d5b12", "#ff637b", "#017c64", "#d05800", "#005d93", "#ffaa8f", "#9cc5ff", "#93322c", "#dab6ea", "#7f396e")
palette2 = c('cadetblue1', 'chartreuse1', 'chocolate1', 'tomato3', 'cadetblue', 'blue3', 'bisque3', 'darkgoldenrod1', 'chartreuse4', 'deepskyblue3', 'darkolivegreen2', 'darkorchid1', 'darkred', 'gold4', 'gray1', 'deeppink', 'deeppink4', 'gray52', 'turquoise4', 'yellowgreen', 'snow2', 'slateblue4')

# Variables especie y Lugar
Especie <- c('Heliconius clysonimus', 'Heliconius cydno', 'Heliconius cydno', 'Heliconius cydno', 'Heliconius cydno', 'Heliconius cydno', 'Heliconius clysonimus', 'Heliconius cydno', 'Heliconius cydno', 'Heliconius clysonimus', 'Heliconius cydno', 'Heliconius cydno')
Lugar <- c('BB', 'BB', 'BB','BB', 'BB', 'BB', 'EAB', 'EAB', 'EAB', 'EAB', 'EAB', 'EAB')

#Agregar variables

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "Especie"), `[`, 1)
especies <- substr(subject,1,1)
especies <- as.integer(sapply(strsplit(samples.out, "Especie"), `[`, 1))
samdf <- data.frame(Subject=subject, Especies=especies)
samdf$Sitio <- Lugar
samdf$Especies <- Especie
samdf$Sitio[samdf$Lugar!="BB"] <- "EAB"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

top20_2 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20_2 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20_2 <- prune_taxa(top20_2, ps.top20_2)

#Graficar abundancia

plot_bar(ps.top20_2, x = "sample_Sample", fill = "Genus") + facet_grid(Especies~Sitio) + geom_bar(stat="identity", position="stack")

ab1 = plot_bar(ps.top20_2, "Genus", fill="Genus", facet_grid=Sitio~Especies)
ab1 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

ab2 = plot_bar(ps, fill="Genus", x="Especies")
ab2 + scale_fill_manual(values=palette) + geom_bar(stat="identity", position="stack")

ab3 = plot_heatmap(ps, method=NULL, sample.label="sample_Sample", taxa.label="Genus")

plot_tree(ex1, color = "SampleType", label.tips = "Phylum", ladderize = "left", justify = "left" , size = "Abundance")


#Graficar índices de Shannon y Simpson
plot_richness(ps, measures=c("Shannon", "Simpson"))

###Separado por especie
Alpha = plot_richness(ps, measures=c("Shannon", "Simpson"), color = "Especies", shape = "Sitio") 
Alpha + geom_point(size=4, alpha=0.75)

###Graficar todos los índices alpha
plot_richness(ps)

# Gráfica separada de especies (CORRER PARA ORDENACIONES)

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Genus") + scale_fill_manual(values=palette)

#############################################################################################################################
#ORDENACIONES

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
theme_set(theme_bw())


OD = ps
wh0 = genefilter_sample(OD, filterfun_sample(function(x) x > 5), A=0.5*nsamples(OD))
OD1 = prune_taxa(wh0, OD)

OD1 = transform_sample_counts(OD1, function(x) 1E6 * x/sum(x))

phylum.sum = tapply(taxa_sums(OD1), tax_table(OD1)[, "Genus"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
OD1 = prune_taxa((tax_table(OD1)[, "Genus"] %in% top7phyla), OD1)


OD.ord <- ordinate(OD1, "PCoA", "bray")

O1 = plot_ordination(OD1, OD.ord, type="taxa", color="Genus", title="taxa")
O1 + facet_wrap(~Phylum, 3)

O2 = plot_ordination(OD1, OD.ord, type="split", color="Genus", shape="Especies", label="Sitio", title="biplot")
O2 + geom_point(size=4, alpha=0.75)

#Grafica por especies

O3 = plot_ordination(OD1, OD.ord, type="samples", color="Sitio")
O3 + geom_polygon(aes(fill=Sitio)) + geom_point(size=5) + ggtitle("samples")


O4 = plot_ordination(OD1, OD.ord, type="biplot", color="Genus", shape="Especies", label="Sitio", title="Biplot PCoA")
O4 + geom_point(size=3, alpha=0.75) 

O5 = plot_ordination(OD1, OD.ord, color="Genus", shape="Especies", label="Sitio", title="Bray PCoA")
O5 + geom_point(size=3, alpha=0.75)

############################################################################################################################################################################




