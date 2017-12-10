##################################################################################################################################################################
######################################################################### FONCTIONS ####################################################################
##################################################################################################################################################################


GTF <- function(id){
  #### By Pierre-Louis Stenger (Pierrelouis.stenger@gmail.com) ####
  
  library(splitstackshape)
  
  pdf("GTF.pdf", height=10,width=20) # create the PDF
  id <- list.files(pattern = "*assembled_transcripts].gtf", recursive = TRUE) # catch all the images call "adapter_content.png" in all folders
  nb <- length(list.files(pattern = "*assembled_transcripts].gtf", recursive = TRUE))
  
  for(i in id){ 

    gtf <- read.table(i, header = FALSE, sep = '\t')
    
    dat <- data.frame(gtf$V1, gtf$V2, gtf$V3, gtf$V4, gtf$V5, gtf$V6, gtf$V7, gtf$V8)
    dat2 <- data.frame(gtf$V9)
    dat3 <- cSplit(cSplit(dat2, 'gtf.V9', sep = ';', direction = 'long'),
                   'gtf.V9', sep = ' ')[, dcast(.SD, cumsum(gtf.V9_1 == 'gene_id') ~ gtf.V9_1)]
    dat4 <- data.frame(gtf$V1, gtf$V2, gtf$V3, gtf$V4, gtf$V5, gtf$V6, gtf$V7, gtf$V8, dat3$conf_hi, dat3$conf_lo, dat3$cov, dat3$exon_number, dat3$FPKM, dat3$frac, dat3$gene_id, dat3$transcript_id)
    colnames(dat4) <- c("Seqname",	"Source",	"Feature",	"Start",	"End",	"Score",	"Strand",	"Frame", "conf_hi",
                        "conf_lo", "cov", "exon_number", "FPKM", "frac", "gene_id", "transcript_id")
    
    # Taille des transcrits et des exon
    dat5 <- data.frame(dat4$Seqname, dat4$Feature, dat4$End-dat4$Start)
    
    # Exon
    exon <- subset(dat5, dat4.Feature=="exon")
    colnames(exon) <- c("Seqname", "Feature", "Size")
    
    # Transcri
    transcri <- subset(dat5, dat4.Feature=="transcript")
    colnames(transcri) <- c("Seqname", "Feature", "Size")
    
    ## Position des exons et transcrits
    dat6 <- data.frame(dat4$Seqname, dat4$Feature, ((dat4$End+dat4$Start)/2))
    colnames(dat6) <- c("Seqname", "Feature", "Position")
    
    # Exon
    exon2 <- subset(dat6, Feature=="exon")
    colnames(exon2) <- c("Seqname", "Feature", "Position")
    
    # Transcri
    transcri2 <- subset(dat6, Feature=="transcript")
    colnames(transcri2) <- c("Seqname", "Feature", "Position")
    
    
    ## Quatres graphiques assemblés
    par(mfrow=c(2,2))
    print(plot(density(exon$Size), main ="Distribution of exon lenghts", sub=i))
    abline(v = c(exon$Size[which.max(exon$Size)]),
           col = 'red', lty = 2)
    text(x = exon$Size[which.max(exon$Size)],
         y = 0,
         labels = round(exon$Size[which.max(exon$Size)], 2))
    plot(density(transcri$Size), main ="Distribution of transcripts lenghts", sub=i)
    abline(v = c(transcri$Size[which.max(transcri$Size)]),
           col = 'red', lty = 2)
    text(x = transcri$Size[which.max(transcri$Size)],
         y = 0,
         labels = round(transcri$Size[which.max(transcri$Size)], 2))
    plot(exon2$Position, main ="Exon position", sub=i)
    plot(transcri2$Position, main ="Transcrits position", sub=i)
    
  }
  dev.off() # close the pdf
}




# GTF()




