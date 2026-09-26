library(limma)
library(affy)
library(annotate)
library(xlsx)
library(statmod)
library(magrittr)

microarrayDegAnalysis <- function(fileNames, type, outputFileName)
{
  fileNames <- fileNames
  
  type <- type
  
  targets <- data.frame(fileNames, type)
  
  #targets <- readTargets("targets_de_rybel.txt", sep="", row.names="filename") ## Reading from files
  my_data <- ReadAffy(filenames = targets$fileName)
  eset <- rma(my_data) ## Preprocessing of raw data
  expr <- exprs(eset)
  exp <- as.factor(targets$type) %>% relevel(., ref = '0h')
  design = model.matrix(~exp)
  
  design <- cbind(control = 1, treat = targets$type == "72h") ## Creating of design  matrix
  Annotation <- read.delim("GPL198-17390.txt", 
                           comment.char="#", quote="\"", sep = '\t',
                           check.names=FALSE, stringsAsFactors=FALSE)
  annot <- Annotation[, c(1, 9)] 
  ID <- featureNames(eset)
  Symbol <- cbind(ID, AGI = Annotation$`Representative Public ID`)
  fData(eset) <- data.frame(Symbol = Symbol) ## Association of Affy IDs with AG codes
  fit <- lmFit(eset, design)
  fit <- eBayes(fit, trend = T, robust = T) ## Creating of liner model and DEG finding
  All_BH <- topTable(fit, coef=2, number = Inf, adjust = 'BH')
  expr2 <- cbind(Annotation$`Representative Public ID`, expr)
  # write.csv(All_BH, "24h_acc.csv")
  
  degs_name <- paste0(outputFileName, '_DEGs.txt')
  counts_name <- paste0(outputFileName, '_counts.txt')
  
  write.table(All_BH, file=degs_name, sep="\t")
  write.table(expr2, file=counts_name, sep="\t")
}

fileNames <- c('24D_0h_rep1.CEL','24D_0h_rep2.CEL','24D_0h_rep3.CEL',
               '24D_72h_rep1.CEL','24D_72h_rep2.CEL','24D_72h_rep3.CEL')

type <- c('0h','0h', '0h','72h','72h', '72h')

outputFileName <- 'GSM5413995_24D_72h'

microarrayDegAnalysis(fileNames = fileNames, type = type, outputFileName = outputFileName)




