#@Author=Arnaud NG
#This script plots the results from the high-confidence epitopes analyses

#import libraries
library("ggplot2")
library("seqinr")
library("grid")
library("RColorBrewer")
library("randomcoloR")
library("gplots")
library("lmPerm")
library("ggpubr")
library("gridExtra")
library("RColorBrewer")
library("indicspecies")
library("tidyr")
library("Cairo")
library("parallel")
library("foreach")
library("doParallel")
library("infotheo")
library("VennDiagram")
library("Biostrings")
library("FD")
library("vegan")
library("lme4")
library("lmerTest")
library("MuMIn")
library("EnvStats")
library("session")
#import script arguments
output_workspace <- as.character(commandArgs(TRUE)[1])
nb_cpus <- as.integer(commandArgs(TRUE)[2])
depth_data_wp <- output_workspace

#name of the reference genome fasta file
fasta_refseq_filename <- "MN908947_3.fasta"
#import reference fasta
genome_refseq <- seqinr::getSequence(object = toupper(read.fasta(paste0(output_workspace,fasta_refseq_filename),seqtype = "DNA",as.string = TRUE,forceDNAtolower = FALSE)),as.string = TRUE)[[1]]
v_orfs_of_interest <- c("orf1a","orf1b","S","E","M","N")
df_epitopes <- read.csv2(file = paste0(output_workspace,"Epitopes_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes$Genomic_start <- as.integer(df_epitopes$Genomic_start)
df_epitopes$Genomic_End <- as.integer(df_epitopes$Genomic_End)

v_lst_id_peptide_seq_in_order <- 1:length(sort(unique(df_epitopes$Peptide),decreasing = FALSE))
names(v_lst_id_peptide_seq_in_order) <- sort(unique(df_epitopes$Peptide),decreasing = FALSE)

#position Inter (prevalence >=0.1)
v_positions_inter <- c(241,   313,  1059,  1163,  2416,  2480,  2558,  3037,  4002,
                       4346,  6312,  7540,  8782,  9286,  9477, 10097, 10319, 10376,
                       11083, 11916, 12525, 13536, 13730, 14408, 14708, 14805, 15324,
                       16647, 17247, 17747, 17858, 18060, 18555, 18877, 19839, 20268,
                       22992, 23401, 23403, 23731, 23929, 25429, 25563, 25979, 26144,
                       26735, 27964, 28144, 28311, 28580, 28657, 28725, 28854, 28863,
                       28881, 28882, 28883, 29540, 29692, 29742, 29870)

#Get list of genomic region and positions
v_orfs <- c("5'UTR", "orf1a", "orf1b", "S","ORF3a","ORF3b","ORF3c","E","M","ORF6","ORF7a", "ORF7b","ORF8", "N", "ORF9c","ORF10","3'UTR")
v_start_orfs <- c(1, 266, 13468, 21563, 25393, 25814, 25524,26245, 26523, 27202, 27394, 27756,27894, 28274, 28734, 29558, 29675)
names(v_start_orfs) <- v_orfs
v_end_orfs <- c(265, 13468, 21555, 25384, 26220, 25882, 25697, 26472, 27191, 27387, 27759, 27887,28259, 29533, 28955, 29674, 29903)
names(v_end_orfs) <- v_orfs
find_ORF_of_mutation <- function(the_site_position){
  indx <- which((v_start_orfs<=the_site_position)&(v_end_orfs>=the_site_position))[1]
  if (length(indx)==0){
    return(NA)
  }else{
    return(v_orfs[indx])
  }
}
v_orfs_length <- v_end_orfs - v_start_orfs + 1
palette_orfs_epitopes <- c("orf1a"="red","orf1b"="blue","S"="green3","E"="orange","M"="grey","N"="purple")

v_genes_with_unique_product <- c(paste0("NSP",1:10),paste0("NSP",12:16), "S","ORF3a","ORF3b","ORF3c","E","M","ORF6","ORF7a", "ORF7b","ORF8", "N", "ORF9c", "ORF10")
v_start_genes <- c(265+1,265+541,265+2455,265+8290,265+9790,265+10708,265+11578,265+11827,265+12421,265+12760,265+13176,265+15972,265+17775,265+19356,265+20394,21563, 25393, 25814, 25524, 26245, 26523, 27202, 27394, 27756, 27894, 28274,28734, 29558)
names(v_start_genes) <- v_genes_with_unique_product
v_end_genes <- c(265+540,265+2454,265+8289,265+9789,265+10707,265+11577,265+11826,265+12420,265+12759,265+13176,265+15971,265+17774,265+19355,265+20393,265+21287,25384, 26220,25882, 25697, 26472, 27191, 27387, 27759, 27887, 28259,29533, 28955, 29674)
names(v_end_genes) <- v_genes_with_unique_product
find_gene_of_mutation <- function(the_site_position){
  indx <- which((v_start_genes<=the_site_position)&(v_end_genes>=the_site_position))[1]
  if (length(indx)==0){
    return(NA)
  }else{
    return(v_genes_with_unique_product[indx])
  }
}
v_genes_length <- v_end_genes - v_start_genes + 1


#find protein site from mutation name
find_prot_site_from_mut_name <- function(the_mut){
  the_mut <- gsub(pattern = "Stop",replacement = "*",x = the_mut,fixed = T)
  v_positions_split <- as.vector(gregexpr(pattern = ";",text = the_mut,fixed = T)[[1]])
  return(as.integer(substr(the_mut,v_positions_split[1]+2,v_positions_split[2]-2)))
}
#find protein from mutation name
find_prot_from_mut_name <- function(the_mut){
  the_mut <- gsub(pattern = "Stop",replacement = "*",x = the_mut,fixed = T)
  v_positions_split <- as.vector(gregexpr(pattern = ";",text = the_mut,fixed = T)[[1]])
  return(substr(the_mut,v_positions_split[2]+1,v_positions_split[3]-1))
}

#Get list of S protein domains positions (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7266584/)
v_S_protein_domains <- c("NTD", "RBD", "SD1", "SD2","CR","HR1","CH-BH","SD3","HR2-TM-CT")
v_start_S_protein_domains <- c(18,331,528,589,846,912,985,1072,1163)
names(v_start_S_protein_domains) <- v_S_protein_domains
v_end_S_protein_domains <- c(306,528,589,677,912,985,1072,1163,1273)
names(v_end_S_protein_domains) <- v_S_protein_domains
find_S_protein_domain_of_mutation <- function(the_site_position){
  indx <- which((v_start_S_protein_domains<=the_site_position)&(v_end_S_protein_domains>=the_site_position))[1]
  if (length(indx)==0){
    return(NA)
  }else{
    return(v_S_protein_domains[indx])
  }
}
v_S_protein_domains_length <- v_end_S_protein_domains - v_start_S_protein_domains + 1

df_epitopes$Mapping_region <- vapply(X = df_epitopes$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes$Mapping_region <- factor(as.character(df_epitopes$Mapping_region),intersect(v_orfs,df_epitopes$Mapping_region))
df_epitopes$peptide_id <- paste0(df_epitopes$Mapping_region,"_",unname(v_lst_id_peptide_seq_in_order[df_epitopes$Peptide]))
df_epitopes <- df_epitopes[,c("peptide_id",names(read.csv2(file = paste0(output_workspace,"Epitopes_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)))]
#exclude possible annotation mistakes
df_epitopes <- subset(df_epitopes,vapply(X = 1:nrow(df_epitopes),FUN = function(i) (return(grepl(pattern = df_epitopes$Mapping_region[i],x = df_epitopes$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes$Group <- as.character(df_epitopes$Group)
df_epitopes$RFU <- as.numeric(df_epitopes$RFU)
df_epitopes <- subset(df_epitopes, RFU>=1000)
df_epitopes$protein_start <- unname(ceiling((df_epitopes$Genomic_start - v_start_orfs[as.character(df_epitopes$Mapping_region)] + 1)/3))
df_epitopes$protein_end <- unname(ceiling((df_epitopes$Genomic_End - v_start_orfs[as.character(df_epitopes$Mapping_region)] + 1)/3))
#save table
#write.table(x=df_epitopes,file = paste0(output_workspace,"Mapped_Epitopes.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

df_variants_NCBI_SRA_amplicon <- rbind(readRDS(paste0(output_workspace,"df_variants_SRA_amplicon_first_wave.rds")),readRDS(paste0(output_workspace,"df_variants_SRA_amplicon_second_wave.rds")))
df_variants_NCBI_SRA_amplicon <- subset(df_variants_NCBI_SRA_amplicon,ORF%in%c("orf1a","orf1b","S","E","M","N"))
df_variants_NCBI_SRA_amplicon$is_fixed <- ifelse(test = df_variants_NCBI_SRA_amplicon$VarFreq>0.75,yes = "Yes",no = "No")

df_variants_NCBI_SRA_amplicon$ORF <- vapply(X = df_variants_NCBI_SRA_amplicon$Position,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_variants_NCBI_SRA_amplicon$gene <- vapply(X = df_variants_NCBI_SRA_amplicon$Position,FUN = function(x) return(find_gene_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
#duplicate variants of orf3a if they occur also in orf3b or orf3c
df_subset_orf3b_orf3c <- subset(df_variants_NCBI_SRA_amplicon,subset=(Position>=v_start_orfs["ORF3a"])&(Position<=v_end_orfs["ORF3a"]))
print(paste0("Example positions:",head(df_variants_NCBI_SRA_amplicon$Position)))
if (nrow(df_subset_orf3b_orf3c)>0){
  if (sum((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"]))>0){
    df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"ORF"] <- "ORF3b"
    df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3b"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3b"])),"gene"] <- "ORF3b"
  }
  if (sum((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"]))>0){
    df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"ORF"] <- "ORF3c"
    df_subset_orf3b_orf3c[which((df_subset_orf3b_orf3c$Position>=v_start_orfs["ORF3c"])&(df_subset_orf3b_orf3c$Position<=v_end_orfs["ORF3c"])),"gene"] <- "ORF3c"
  }
  if (nrow(df_subset_orf3b_orf3c)>0){
    df_subset_orf3b_orf3c <- subset(df_subset_orf3b_orf3c,ORF!="ORF3a")
  }
  df_variants_NCBI_SRA_amplicon <- rbind(df_variants_NCBI_SRA_amplicon,df_subset_orf3b_orf3c)
}
#duplicate variants of ORF7a if they occur also in ORF7b
df_subset_orf7b <- subset(df_variants_NCBI_SRA_amplicon,subset=(Position>=v_start_orfs["ORF7a"])&(Position<=v_end_orfs["ORF7a"]))
if (nrow(df_subset_orf7b)>0){
  if (sum((df_subset_orf7b$Position>=v_start_orfs["ORF7b"])&(df_subset_orf7b$Position<=v_end_orfs["ORF7b"]))>0){
    df_subset_orf7b[which((df_subset_orf7b$Position>=v_start_orfs["ORF7b"])&(df_subset_orf7b$Position<=v_end_orfs["ORF7b"])),"ORF"] <- "ORF7b"
  }
  if (nrow(df_subset_orf7b)>0){
    df_subset_orf7b <- subset(df_subset_orf7b,ORF!="ORF7a")
  }
  df_variants_NCBI_SRA_amplicon <- rbind(df_variants_NCBI_SRA_amplicon,df_subset_orf7b)
}
#duplicate variants of N if they occur also in orf9c
df_subset_ORF9c <- subset(df_variants_NCBI_SRA_amplicon,subset=(Position>=v_start_orfs["N"])&(Position<=v_end_orfs["N"]))
if (nrow(df_subset_ORF9c)>0){
  if (sum((df_subset_ORF9c$Position>=v_start_orfs["ORF9c"])&(df_subset_ORF9c$Position<=v_end_orfs["ORF9c"]))>0){
    df_subset_ORF9c[which((df_subset_ORF9c$Position>=v_start_orfs["ORF9c"])&(df_subset_ORF9c$Position<=v_end_orfs["ORF9c"])),"ORF"] <- "ORF9c"
  }
  if (nrow(df_subset_ORF9c)>0){
    df_subset_ORF9c <- subset(df_subset_ORF9c,ORF!="N")
  }
  df_variants_NCBI_SRA_amplicon <- rbind(df_variants_NCBI_SRA_amplicon,df_subset_ORF9c)
}

#function that find the original and mutated codons of a variant
get_ref_and_mutated_codon <- function(the_position,ref_nucl,new_nucl){
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    the_ref_codon <- NA
    the_mut_codon <- NA
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      the_ref_codon <- paste0(ref_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2),sep="")
      the_mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2),sep="")
    }else if (pos_in_codon==2){
      the_ref_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),ref_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1),sep="")
      the_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1),sep="")
    }else if (pos_in_codon==3){
      the_ref_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),ref_nucl,sep="")
      the_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl,sep="")
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  return(list(ref_codon=the_ref_codon,mutated_codon=the_mut_codon))
}
#build function that determines whether a mutation is synonymous or not
is_mutation_synonymous <- function(the_reference_codon,the_mutated_codon){
  if (the_reference_codon %in% c("TAA","TAG","TGA")){
    return(NA)
  }else{
    return(seqinr::translate(seq = unlist(strsplit(the_reference_codon,"")))==seqinr::translate(seq = unlist(strsplit(the_mutated_codon,""))))
  }
}
#build function that determines whether a mutation is synonymous or not
translate_seq <- function(the_codon){
  if (is.na(the_codon)){
    return(NA)
  }else if (the_codon %in% c("TAA","TAG","TGA")){
    return("Stop")
  }else{
    return(seqinr::translate(seq = unlist(strsplit(the_codon,""))))
  }
}
#Original codon and mutated codon
df_variants_NCBI_SRA_amplicon$ref_codon <- NA
df_variants_NCBI_SRA_amplicon$mut_codon <- NA
df_variants_NCBI_SRA_amplicon$pos_in_ORF <- NA
df_variants_NCBI_SRA_amplicon$pos_in_gene <- NA
df_variants_NCBI_SRA_amplicon$pos_in_protein <- NA
for (i in 1:nrow(df_variants_NCBI_SRA_amplicon)){
  df_variants_NCBI_SRA_amplicon$ref_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants_NCBI_SRA_amplicon$Position[i],ref_nucl = df_variants_NCBI_SRA_amplicon$Ref[i],new_nucl = df_variants_NCBI_SRA_amplicon$VarAllele[i]))$ref_codon
  df_variants_NCBI_SRA_amplicon$mut_codon[i] <-(get_ref_and_mutated_codon(the_position = df_variants_NCBI_SRA_amplicon$Position[i],ref_nucl = df_variants_NCBI_SRA_amplicon$Ref[i],new_nucl = df_variants_NCBI_SRA_amplicon$VarAllele[i]))$mutated_codon
  df_variants_NCBI_SRA_amplicon$old_aa[i] <- translate_seq(the_codon = df_variants_NCBI_SRA_amplicon$ref_codon[i])
  df_variants_NCBI_SRA_amplicon$new_aa[i] <- translate_seq(the_codon = df_variants_NCBI_SRA_amplicon$mut_codon[i] )
  df_variants_NCBI_SRA_amplicon$pos_in_ORF[i] <- df_variants_NCBI_SRA_amplicon$Position[i] - v_start_orfs[df_variants_NCBI_SRA_amplicon$ORF[i]] + 1
  df_variants_NCBI_SRA_amplicon$pos_in_gene[i] <- df_variants_NCBI_SRA_amplicon$Position[i] - v_start_genes[df_variants_NCBI_SRA_amplicon$gene[i]] + 1
  df_variants_NCBI_SRA_amplicon$pos_in_protein[i] <- ceiling(df_variants_NCBI_SRA_amplicon$pos_in_gene[i]/3)
  #print(paste0("Iterations ",i," out of ",nrow(df_variants_NCBI_SRA_amplicon)))
}
df_variants_NCBI_SRA_amplicon$mutation_name <- paste0(paste0(df_variants_NCBI_SRA_amplicon$Ref,df_variants_NCBI_SRA_amplicon$Position,df_variants_NCBI_SRA_amplicon$VarAllele,""),";",paste0(df_variants_NCBI_SRA_amplicon$old_aa,df_variants_NCBI_SRA_amplicon$pos_in_protein,df_variants_NCBI_SRA_amplicon$new_aa),";",df_variants_NCBI_SRA_amplicon$ORF,";",df_variants_NCBI_SRA_amplicon$gene)
#Define Nonsense and non-coding mutations
df_variants_NCBI_SRA_amplicon$is_nonsense <- (df_variants_NCBI_SRA_amplicon$new_aa=="Stop")
df_variants_NCBI_SRA_amplicon$is_UTR <- (is.na(df_variants_NCBI_SRA_amplicon$new_aa))
df_variants_NCBI_SRA_amplicon$is_synonymous <- ifelse((is.na(df_variants_NCBI_SRA_amplicon$old_aa)|(df_variants_NCBI_SRA_amplicon$new_aa=="Stop")),yes = NA,no = df_variants_NCBI_SRA_amplicon$old_aa==df_variants_NCBI_SRA_amplicon$new_aa)

df_variants_NCBI_SRA_amplicon$mutation_type <- ifelse(test = df_variants_NCBI_SRA_amplicon$is_UTR,yes = "UTR",no = ifelse(test = df_variants_NCBI_SRA_amplicon$is_nonsense,yes = "Nonsense",no = ifelse(test = df_variants_NCBI_SRA_amplicon$is_synonymous,yes = "Synonymous",no = "Non-Synonymous")))
df_variants_NCBI_SRA_amplicon$S_protein_domain <- ifelse(test=df_variants_NCBI_SRA_amplicon$ORF=="S",yes = vapply(X = df_variants_NCBI_SRA_amplicon$pos_in_protein,FUN = function(x) return(find_S_protein_domain_of_mutation(the_site_position = x)),FUN.VALUE = c("")),no=NA)
v_recurrence_mut_NCBI_SRA_amplicon <- as.vector(table(df_variants_NCBI_SRA_amplicon$mutation_name))
names(v_recurrence_mut_NCBI_SRA_amplicon) <- names(table(df_variants_NCBI_SRA_amplicon$mutation_name))
df_variants_NCBI_SRA_amplicon$is_prevalence_above_transmission_threshold <- df_variants_NCBI_SRA_amplicon$mutation_name%in%(names(v_recurrence_mut_NCBI_SRA_amplicon)[v_recurrence_mut_NCBI_SRA_amplicon>=3])
v_nb_samples_NCBI_SRA_amplicon <- length(unique(df_variants_NCBI_SRA_amplicon$Sample))
#name of the reference genome fasta file
fasta_refseq_filename <- "MN908947_3.fasta"
#import reference fasta
genome_refseq <- seqinr::getSequence(object = toupper(read.fasta(paste0(output_workspace,fasta_refseq_filename),seqtype = "DNA",as.string = TRUE,forceDNAtolower = FALSE)),as.string = TRUE)[[1]]

#Presence of mutations of interest in the S protein (define by Emma B. Hodcroft as of 2020-12-23 and https://virological.org/t/mutations-arising-in-sars-cov-2-spike-on-sustained-human-to-human-transmission-and-human-to-animal-passage/578)
v_S_region_mutations_of_interest <- sort(unique(c("A222V","S477N","S98F","D80Y","N439K","Y453F","N501S","N501T","N501Y","A626S","V1122L","H69","D80A", "D215G", "P681H", "A701V", "T716I", "D1118H", "D614G","A570D", "K417N", "E484K", "N501Y", "S983A")))

#patients group
v_patients_to_group <- apply(as.matrix(table(df_epitopes$patient_ID,df_epitopes$Group)),MARGIN = 1,FUN = function(x) colnames(as.matrix(table(df_epitopes$patient_ID,df_epitopes$Group)))[which(x>0)])
v_patients_to_binary_group <- ifelse(test=apply(as.matrix(table(df_epitopes$patient_ID,df_epitopes$Group)),MARGIN = 1,FUN = function(x) colnames(as.matrix(table(df_epitopes$patient_ID,df_epitopes$Group)))[which(x>0)])%in%c("1","2"),yes="Positive",no="Negative")
names(v_patients_to_binary_group) <- rownames(as.matrix(table(df_epitopes$patient_ID,df_epitopes$Group)))
v_patients_to_binary_group <- v_patients_to_binary_group[c(which(v_patients_to_group=="1"),which(v_patients_to_group=="2"),which(v_patients_to_group=="3"))]
v_covid_pos_patients <- names(v_patients_to_binary_group[v_patients_to_binary_group=="Positive"])
v_covid_neg_patients <- names(v_patients_to_binary_group[v_patients_to_binary_group=="Negative"])
#peptide per group
v_peptide_group1 <- subset(df_epitopes,Group==1)$Peptide
v_peptide_group2 <- subset(df_epitopes,Group==2)$Peptide
v_peptide_group3 <- subset(df_epitopes,Group==3)$Peptide
v_position_peptide_group1 <- NULL
for (i in 1:nrow(subset(df_epitopes,Group==1))){
  v_position_peptide_group1 <- c(v_position_peptide_group1, ((subset(df_epitopes,Group==1)$Genomic_start[i]):(subset(df_epitopes,Group==1)$Genomic_End[i])))
}
v_position_peptide_group1 <- sort(v_position_peptide_group1)
v_position_peptide_group2 <- NULL
for (i in 1:nrow(subset(df_epitopes,Group==2))){
  v_position_peptide_group2 <- c(v_position_peptide_group2, ((subset(df_epitopes,Group==2)$Genomic_start[i]):(subset(df_epitopes,Group==2)$Genomic_End[i])))
}
v_position_peptide_group2 <- sort(v_position_peptide_group2)
v_position_peptide_group3 <- NULL
for (i in 1:nrow(subset(df_epitopes,Group==3))){
  v_position_peptide_group3 <- c(v_position_peptide_group3, ((subset(df_epitopes,Group==3)$Genomic_start[i]):(subset(df_epitopes,Group==3)$Genomic_End[i])))
}
v_position_peptide_group3 <- sort(v_position_peptide_group3)
v_position_with_highest_antibody_response <- NULL
for (i in 1:nrow(subset(df_epitopes,RFU>=quantile(df_epitopes$RFU,probs = 0.5)))){
  v_position_with_highest_antibody_response <- c(v_position_with_highest_antibody_response, subset(df_epitopes,RFU>=quantile(df_epitopes$RFU,probs = 0.5))$Genomic_start[i]:subset(df_epitopes,RFU>=quantile(df_epitopes$RFU,probs = 0.5))$Genomic_End[i])
}
v_position_with_highest_antibody_response <- sort(v_position_with_highest_antibody_response)
palette_patient_groups <- c("1"=alpha("tomato",0.6), "2"=alpha('green3',0.6), "3"=alpha('royalblue',0.6))

#Venn Diagram
venn.diagram(
  x = list(v_peptide_group1,v_peptide_group2,v_peptide_group3),
  category.names = c("Group 1" , "Group 2","Group 3") ,
  filename =  paste0(output_workspace,"Venn_diagram_epitopes_peptides.png"),
  output = TRUE ,
  imagetype="png" ,
  fill = c(alpha("tomato",0.6), alpha('green3',0.6), alpha('deepskyblue',0.6)),
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

#RFU across ORFs and groups
p <- ggboxplot(data = df_epitopes, x = "Mapping_region", y="RFU",color = "Group",add = "jitter")+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Genomic region") +scale_y_continuous(limits = c(0,max(df_epitopes$RFU)+10000),breaks=seq(0,max(df_epitopes$RFU)+10000,10000))
facet(p +  stat_compare_means(), facet.by = "Group", ncol = 1)
ggsave(filename = "RFU_by_ORF_and_Groups.png", path=output_workspace, width = 20, height = 20, units = "cm",dpi = 1200)

#RFU across groups
ggboxplot(data = df_epitopes, x = "Group", y="RFU",color = "Group",add = "jitter") +  stat_compare_means() + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Group")
ggsave(filename = "RFU_by_Groups.png", path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#RFU across Antibody and groups
p <- ggboxplot(data = df_epitopes, x = "Antibody", y="RFU",color = "Group",add = "jitter")+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Antibody") +scale_y_continuous(limits = c(0,max(df_epitopes$RFU)+10000),breaks=seq(0,max(df_epitopes$RFU)+10000,10000))
facet(p +  stat_compare_means(), facet.by = "Group", ncol = 1)
ggsave(filename = "RFU_by_Antibody_and_Groups.png", path=output_workspace, width = 20, height = 20, units = "cm",dpi = 1200)

df_regions_per_group <- as.data.frame(table(df_epitopes$Mapping_region,df_epitopes$Group),stringAsFactors=FALSE)
names(df_regions_per_group) <- c("Mapped_region","Group","Count")
df_regions_per_group$Group <- paste0("Group",as.integer(df_regions_per_group$Group))
df_regions_per_group$Mapped_region <- factor(df_regions_per_group$Mapped_region,intersect(v_orfs,df_regions_per_group$Mapped_region))

p <- ggbarplot(data = df_regions_per_group, x = "Mapped_region", y="Count", fill="Group",ggtheme = theme_light())+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("Number of epitope occurrences") + xlab("Genomic region")
facet(p, facet.by = "Group", ncol = 1)
ggsave(filename = "Number_of_mapped_epitope_peptides_per_genomic_region.png", path=output_workspace, width = 15, height = 18, units = "cm",dpi = 1200)

df_regions_per_group$Density <- df_regions_per_group$Count/v_orfs_length[df_regions_per_group$Mapped_region]
p <- ggbarplot(data = df_regions_per_group, x = "Mapped_region", y="Density", fill="Group",ggtheme = theme_light())+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("Density of mapped peptides (count/ORF length)") + xlab("Genomic region")
facet(p, facet.by = "Group", ncol = 1)
ggsave(filename = "Density_of_mapped_epitope_peptides_per_genomic_region.png", path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#top 100 epitope hotspots
df_epitopes_top100_hotspots <- (unique(df_nb_occurence_epitope_position_by_group[,c("Position","Count")])[order(unique(df_nb_occurence_epitope_position_by_group[,c("Position","Count")])$Count,decreasing = TRUE),])[1:100,]

#function for plotting linear model
ggplotRegression <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit)$coefficients[,4][2]), format = "e", digits = 3)))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P =",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm")
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}
ggplotRegression_export_eps <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit)$coefficients[,4][2]), format = "e", digits = 3)))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P =",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm", device = cairo_ps)
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}
#build function that determines whether a mutation is synonymous or not 
is_mutation_synonymous <- function(the_reference_codon,the_mutated_codon){
  if (the_reference_codon %in% c("TAA","TAG","TGA")){
    return(NA)
  }else{
    return(translate(seq = unlist(strsplit(the_reference_codon,"")))==translate(seq = unlist(strsplit(the_mutated_codon,""))))
  }
}
#build function that determines whether a mutation is synonymous or not 
translate_seq <- function(the_codon){
  if (is.na(the_codon)){
    return(NA)
  }else if (the_codon %in% c("TAA","TAG","TGA")){
    return("Stop")
  }else{
    return(translate(seq = unlist(strsplit(the_codon,""))))
  }
}
# #function that determines if a mutation is in a mutation hotspot, as identified on nextrain April 24, 2020
# is_in_mutation_hotspot <- function(the_variant){
#   return(any(sapply(X = c(4049:4051,11081:11083,13400:13402,14407:14409,21575:21577),FUN = function(x) return(grepl(pattern = as.character(x),x = the_variant,fixed = TRUE)))))
# }
# #Vectorial version of the function "is_in_mutation_hotspot"
# is_in_mutation_hotspot_vec <- function(x){
#   return(vapply(X = x,FUN = function(y) return(is_in_mutation_hotspot(y)),FUN.VALUE = c(FALSE)))
# }

#create a function that returns number of synonymous sites for a single position in the genome
calculate_nb_ss_position_in_genome <- function(the_position){
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      the_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
    }else if (pos_in_codon==2){
      the_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
    }else if (pos_in_codon==3){
      the_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(the_codon)!=3){
    stop("codon length should be 3!")
  }
  possible_single_site_mutated_codons <- rep("",3)
  num_mut_codon <-1
  for (pos_codon in pos_in_codon){
    if (substr(the_codon,start = pos_codon,stop=pos_codon)=="A"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
      
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="T"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="C"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else{#G
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }
  }
  #count the number of synonymous mutations based on the genetic code
  nb_unique_syn_mut_codons <-0 #default initialization
  if (the_codon == "TTT") {
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTC"])
    
  } else if (the_codon == "TTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTT"])
    
  } else if (the_codon == "TTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTG","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCC","TCA","TCG","AGT","AGC")])
  } else if (the_codon == "TCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCC","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCC","AGT","AGC")])
    
  } else if (the_codon == "TAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAC")])
    
  } else if (the_codon == "TAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAT")])
    
  } else if (the_codon == "TGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGC")])
    
  } else if (the_codon == "TGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGT")])
    
  } else if (the_codon == "TGG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "CTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTC","CTA","CTG")])
    
  } else if (the_codon == "CTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTA","CTG")])
    
  } else if (the_codon == "CTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTG")])
    
  } else if (the_codon == "CTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTA")])
    
  } else if (the_codon == "CCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCC","CCA","CCG")])
    
  } else if (the_codon == "CCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCA","CCG")])
    
    
  } else if (the_codon == "CCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCG")])
    
  } else if (the_codon == "CCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCA")])
    
  } else if (the_codon == "CAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAC")])
    
  } else if (the_codon == "CAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAT")])
    
  } else if (the_codon == "CAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAG")])
    
  } else if (the_codon == "CAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAA")])
    
  } else if (the_codon == "CGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGC","CGA","CGG")])
    
  } else if (the_codon == "CGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGG")])
    
  } else if (the_codon == "CGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGC","CGG")])
    
  } else if (the_codon == "CGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGC")])
    
  } else if (the_codon == "ATT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATA")])
    
  } else if (the_codon == "ATC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATT","ATA")])
    
  } else if (the_codon == "ATA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATT")])
    
  } else if (the_codon == "ATG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "ACT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACC","ACA","ACG")])
    
    
  } else if (the_codon == "ACC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACA","ACG")])
    
  } else if (the_codon == "ACA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACG")])
    
    
  } else if (the_codon == "ACG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACA")])
    
  } else if (the_codon == "AAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAC")])
    
  } else if (the_codon == "AAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAT")])
    
  } else if (the_codon == "AAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAG")])
    
  } else if (the_codon == "AAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAA")])
    
  } else if (the_codon == "AGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGC","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGT","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGG")])
    
  } else if (the_codon == "AGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGA")])
    
  } else if (the_codon == "GTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTG")])
    
  } else if (the_codon == "GTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTT","GTA","GTG")])
    
  } else if (the_codon == "GTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTT","GTG")])
    
  } else if (the_codon == "GTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTT")])
    
  } else if (the_codon == "GCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCG")])
    
  } else if (the_codon == "GCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCT","GCA","GCG")])
    
  } else if (the_codon == "GCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCT","GCG")])
    
  } else if (the_codon == "GCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCT")])
    
  } else if (the_codon == "GAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAC")])
    
  } else if (the_codon == "GAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAT")])
    
  } else if (the_codon == "GAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAG")])
    
  } else if (the_codon == "GAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAA")])
    
  } else if (the_codon == "GGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGG")])
    
  } else if (the_codon == "GGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGT","GGA","GGG")])
    
  } else if (the_codon == "GGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGT","GGG")])
    
    
  } else if (the_codon == "GGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGT")])
  }
  return((nb_unique_syn_mut_codons/3))
}

#create a function that returns number of possible SINGLE-SITE synonymous mutations divided by 3 for a CODON
calculate_third_of_possible_ns_codon <- function(the_codon){
  the_codon <- toupper(the_codon)
  if (nchar(the_codon)!=3){
    stop("codon length should be 3!")
  }
  possible_single_site_mutated_codons <- rep("",9)
  num_mut_codon <-1
  for (pos_codon in 1:3){
    if (substr(the_codon,start = pos_codon,stop=pos_codon)=="A"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
      
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="T"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="C"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else{#G
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }
  }
  #count the number of synonymous mutations based on the genetic code
  nb_unique_syn_mut_codons <-0 #default initialization
  if (the_codon == "TTT") {
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTC"])
    
  } else if (the_codon == "TTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTT"])
    
  } else if (the_codon == "TTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTG","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCC","TCA","TCG","AGT","AGC")])
  } else if (the_codon == "TCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCC","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCC","AGT","AGC")])
    
  } else if (the_codon == "TAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAC")])
    
  } else if (the_codon == "TAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAT")])
    
  } else if (the_codon == "TGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGC")])
    
  } else if (the_codon == "TGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGT")])
    
  } else if (the_codon == "TGG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "CTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTC","CTA","CTG")])
    
  } else if (the_codon == "CTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTA","CTG")])
    
  } else if (the_codon == "CTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTG")])
    
  } else if (the_codon == "CTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTA")])
    
  } else if (the_codon == "CCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCC","CCA","CCG")])
    
  } else if (the_codon == "CCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCA","CCG")])
    
    
  } else if (the_codon == "CCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCG")])
    
  } else if (the_codon == "CCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCA")])
    
  } else if (the_codon == "CAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAC")])
    
  } else if (the_codon == "CAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAT")])
    
  } else if (the_codon == "CAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAG")])
    
  } else if (the_codon == "CAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAA")])
    
  } else if (the_codon == "CGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGC","CGA","CGG")])
    
  } else if (the_codon == "CGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGG")])
    
  } else if (the_codon == "CGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGC","CGG")])
    
  } else if (the_codon == "CGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGC")])
    
  } else if (the_codon == "ATT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATA")])
    
  } else if (the_codon == "ATC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATT","ATA")])
    
  } else if (the_codon == "ATA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATT")])
    
  } else if (the_codon == "ATG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "ACT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACC","ACA","ACG")])
    
    
  } else if (the_codon == "ACC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACA","ACG")])
    
  } else if (the_codon == "ACA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACG")])
    
    
  } else if (the_codon == "ACG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACA")])
    
  } else if (the_codon == "AAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAC")])
    
  } else if (the_codon == "AAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAT")])
    
  } else if (the_codon == "AAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAG")])
    
  } else if (the_codon == "AAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAA")])
    
  } else if (the_codon == "AGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGC","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGT","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGG")])
    
  } else if (the_codon == "AGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGA")])
    
  } else if (the_codon == "GTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTG")])
    
  } else if (the_codon == "GTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTT","GTA","GTG")])
    
  } else if (the_codon == "GTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTT","GTG")])
    
  } else if (the_codon == "GTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTT")])
    
  } else if (the_codon == "GCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCG")])
    
  } else if (the_codon == "GCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCT","GCA","GCG")])
    
  } else if (the_codon == "GCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCT","GCG")])
    
  } else if (the_codon == "GCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCT")])
    
  } else if (the_codon == "GAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAC")])
    
  } else if (the_codon == "GAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAT")])
    
  } else if (the_codon == "GAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAG")])
    
  } else if (the_codon == "GAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAA")])
    
  } else if (the_codon == "GGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGG")])
    
  } else if (the_codon == "GGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGT","GGA","GGG")])
    
  } else if (the_codon == "GGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGT","GGG")])
    
    
  } else if (the_codon == "GGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGT")])
  }
  return((nb_unique_syn_mut_codons/3))
}

calculate_epitope_related_sites_nb_ss <- function(start_pos,end_pos){
  Nb_syn_sites_peptide <- 0
  for (pos_in_gene in seq(from =start_pos,to = end_pos,by = 3)){
    current_codon_gene <- substr(x = genome_refseq,start = pos_in_gene,stop=pos_in_gene+2)
    Nb_syn_sites_peptide <- Nb_syn_sites_peptide + calculate_third_of_possible_ns_codon(current_codon_gene)
  }
  return(Nb_syn_sites_peptide)
}

v_nb_ss_epitope_related_coding_regions <- NULL 
for (i in 1:nrow(df_epitopes)){
  if (!paste0(df_epitopes$Genomic_start[i],"-",df_epitopes$Genomic_End[i])%in%names(v_nb_ss_epitope_related_coding_regions)){
    v_nb_ss_epitope_related_coding_regions <- c(v_nb_ss_epitope_related_coding_regions,calculate_epitope_related_sites_nb_ss(start_pos = df_epitopes$Genomic_start[i],end_pos = df_epitopes$Genomic_End[i]))
    names(v_nb_ss_epitope_related_coding_regions)[length(v_nb_ss_epitope_related_coding_regions)] <- paste0(df_epitopes$Genomic_start[i],"-",df_epitopes$Genomic_End[i])
  }
}

get_position_mutation <- function(the_mut_name){
  return(as.integer(substr(x = the_mut_name,start = 2,stop = as.vector(regexpr(pattern = ";",text = the_mut_name,fixed = T))-2)))
}

v_group_patient <- unique(df_epitopes[,c("patient_ID","Group")])$Group
names(v_group_patient) <- unique(df_epitopes[,c("patient_ID","Group")])$patient_ID

v_seq_peptide <- unique(df_epitopes[,c("peptide_id","Peptide")])$Peptide
names(v_seq_peptide) <- unique(df_epitopes[,c("peptide_id","Peptide")])$peptide_id

#known epitope sites
v_unique_epitope_positions <- sort(unique(c(v_position_peptide_group1,v_position_peptide_group2,v_position_peptide_group3)))
#Non-epitope sites in analyzed ORFs
v_non_epitope_sites <- sort(unique(setdiff(c(v_start_orfs["orf1a"]:v_end_orfs["orf1a"],v_start_orfs["orf1b"]:v_end_orfs["orf1b"],v_start_orfs["S"]:v_end_orfs["S"],v_start_orfs["E"]:v_end_orfs["E"],v_start_orfs["M"]:v_end_orfs["M"],v_start_orfs["N"]:v_end_orfs["N"]),v_unique_epitope_positions)))#sort(unique(setdiff(setdiff(1:nchar(genome_refseq),c(v_start_orfs["ORF3a"]:v_end_orfs["ORF3a"],v_start_orfs["ORF3b"]:v_end_orfs["ORF3b"],v_start_orfs["ORF3c"]:v_end_orfs["ORF3c"],v_start_orfs["ORF6"]:v_end_orfs["ORF6"],v_start_orfs["ORF7a"]:v_end_orfs["ORF7a"],v_start_orfs["ORF7b"]:v_end_orfs["ORF7b"],v_start_orfs["ORF8"]:v_end_orfs["ORF8"],v_start_orfs["ORF10"]:v_end_orfs["ORF10"])),v_unique_epitope_positions)))

v_length_epitope_sites_vs_others <- c(length(v_unique_epitope_positions),length(v_non_epitope_sites))
names(v_length_epitope_sites_vs_others) <- c("TRUE","FALSE")

v_length_top100_epitope_hotspots_vs_others <- c(100,length(v_unique_epitope_positions)-100)
names(v_length_top100_epitope_hotspots_vs_others) <- c("TRUE","FALSE")

ggplot(data = df_epitopes[,c("Genomic_start","Mapping_region","RFU")]) + geom_line(mapping = aes(x=Genomic_start,y=RFU,col=Mapping_region)) + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Position")+ guides(col=guide_legend(title="Genomic Region")) 
ggsave(filename = "Epitopes_RFU_across_SARS_CoV_2_genome.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)

#NCBI dataset
df_variants_NCBI_SRA_amplicon$is_epitope_related <- df_variants_NCBI_SRA_amplicon$Position%in% sort(v_unique_epitope_positions)
df_variants_NCBI_SRA_amplicon$is_in_top100_epitope <- df_variants_NCBI_SRA_amplicon$Position %in% df_epitopes_top100_hotspots$Position
df_variants_NCBI_SRA_amplicon$is_exclusive_to_group1 <- (df_variants_NCBI_SRA_amplicon$Position%in%v_position_peptide_group1)&(!df_variants_NCBI_SRA_amplicon$Position%in%c(v_position_peptide_group2,v_position_peptide_group3))
df_variants_NCBI_SRA_amplicon$is_exclusive_to_group2 <- (df_variants_NCBI_SRA_amplicon$Position%in%v_position_peptide_group1)&(!df_variants_NCBI_SRA_amplicon$Position%in%c(v_position_peptide_group1,v_position_peptide_group3))
df_variants_NCBI_SRA_amplicon$is_exclusive_to_group3 <- (df_variants_NCBI_SRA_amplicon$Position%in%v_position_peptide_group1)&(!df_variants_NCBI_SRA_amplicon$Position%in%c(v_position_peptide_group1,v_position_peptide_group2))
df_variants_NCBI_SRA_amplicon$is_shared_epitope_position <- df_variants_NCBI_SRA_amplicon$Position%in%v_shared_epitope_positions
lst_samples_NCBI_SRA_amplicon <- sort(unique(df_variants_NCBI_SRA_amplicon$Sample))
#determine what's the minimum coverage required for a site to have at least 80% power for detecting at least 10 copies of SNVs at >=5%
min_cov <- 1
p <- 0
min_nb_reads_supporting_snv <- 5

while (p<0.8){
  p <- 1 - pbinom(q = min_nb_reads_supporting_snv, size = min_cov, prob = 0.05)
  min_cov <- min_cov + 1
  if (min_cov%% 10){
    print(paste0("current min cov :", min_cov))
  }
}

df_variants_site_enough_covered_NCBI_SRA_amplicon <- subset(df_variants_NCBI_SRA_amplicon, total_depth>=min_cov)
v_unique_epitope_positions_with_enough_coverage_NCBI_SRA_amplicon <- intersect(v_unique_epitope_positions,df_variants_site_enough_covered_NCBI_SRA_amplicon$Position)

v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon <- unique(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,mutation_type=="Non-Synonymous")$mutation_name)
v_lst_synonymous_mutations_NCBI_SRA_amplicon <- unique(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,mutation_type=="Synonymous")$mutation_name)

#compare mutation and substitution rate in epitopes vs outside
df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon <- data.frame(Sample=rep(lst_samples_NCBI_SRA_amplicon,length(c(TRUE,FALSE))),is_epitope_related=rep(c(TRUE,FALSE),each=length(lst_samples_NCBI_SRA_amplicon)),nb_mutations=0,nb_fixed_mutations=0,mut_rate=0,subst_rate=0,stringsAsFactors = F)
nb_cores <- nb_cpus
lst_splits <- split(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon), ceiling(seq_along(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon))/(nrow(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon)/nb_cores)))
the_f_parallel_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  df_metrics_current_subset <- df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon[the_vec,]
  count_iter <- 0
  for (the_i in 1:nrow(df_metrics_current_subset)){
    df_depth_NCBI_SRA_amplicon_current_sample <- read.csv2(file = paste0(depth_data_wp,"depth_report_NCBI_SRA_amplicon/df_depth_NCBI_SRA_amplicon_",df_metrics_current_subset$Sample[the_i],".csv"),sep = ",",header = F,stringsAsFactors = FALSE)
    colnames(df_depth_NCBI_SRA_amplicon_current_sample) <- c("sample","position","depth")
    df_depth_NCBI_SRA_amplicon_current_sample$ORF <- vapply(X = df_depth_NCBI_SRA_amplicon_current_sample$position,FUN = find_ORF_of_mutation,FUN.VALUE = c(""))
    df_depth_NCBI_SRA_amplicon_current_sample <- unique(df_depth_NCBI_SRA_amplicon_current_sample)
    v_currentsample_positions_enough_covered <- subset(df_depth_NCBI_SRA_amplicon_current_sample,(depth>=min_cov))$position
    if (df_metrics_current_subset$is_epitope_related[the_i]){
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_unique_epitope_positions,v_currentsample_positions_enough_covered)
    }else{
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_non_epitope_sites,v_currentsample_positions_enough_covered)
    }
    df_metrics_current_subset$nb_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(is_fixed=="No")&(!Position%in%v_positions_inter)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_fixed_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(is_fixed=="Yes")&(is_prevalence_above_transmission_threshold)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] <- length(v_the_sites_with_enough_cov_for_current_category)
    df_metrics_current_subset$mut_rate[the_i] <- df_metrics_current_subset$nb_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$subst_rate[the_i] <- df_metrics_current_subset$nb_fixed_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$Nb_ss[the_i] <- sum(unname(vapply(X = v_the_sites_with_enough_cov_for_current_category,FUN = calculate_nb_ss_position_in_genome,FUN.VALUE = c(0))),na.rm=T)
    df_metrics_current_subset$Nb_nss[the_i] <- df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] - df_metrics_current_subset$Nb_ss[the_i]
    df_metrics_current_subset$within_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$within_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    
    if (the_i%%100==0){
      print(paste0("[Epitope sites vs others (NCBI)] Core ",i_cl,": Step ",the_i," done out of ",nrow(df_metrics_current_subset),"!"))
    }
  }
  df_metrics_current_subset$pN <- df_metrics_current_subset$within_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$pS <- df_metrics_current_subset$within_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$pN_pS <- ifelse(test=df_metrics_current_subset$pS==0,yes=NA,no=df_metrics_current_subset$pN/df_metrics_current_subset$pS)
  df_metrics_current_subset$dN <- df_metrics_current_subset$between_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$dS <- df_metrics_current_subset$between_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$dN_dS <- ifelse(test=df_metrics_current_subset$dS==0,yes=NA,no=df_metrics_current_subset$dN/df_metrics_current_subset$dS)
  df_metrics_current_subset$alpha_MK_Test <- ifelse(test=df_metrics_current_subset$dN_dS==0,yes=NA,no=(1-((df_metrics_current_subset$pN_pS)/(df_metrics_current_subset$dN_dS))))
  
  return(df_metrics_current_subset)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_Evo_rates_analyses.txt"))
registerDoParallel(cl)
df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon(i_cl)
stopCluster(cl)
# #saveRDS(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,paste0(output_workspace,"df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.rds"))

#ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_epitope_related),levels=c("TRUE","FALSE")),y = mut_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("Epitope sites?") + ylab("Within-host mutation rate (Count / Length)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$mut_rate)+1e-2,1e-2),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$mut_rate)+1e-2)) + stat_compare_means(method = "wilcox")
#ggsave(filename = "Mutation_rate_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
#ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_epitope_related),levels=c("TRUE","FALSE")),y = subst_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("Epitope sites?") + ylab("Substitution rate (Count / Length)") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$subst_rate)+1e-4,1e-4),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$subst_rate)+1e-4)) + stat_compare_means(method = "wilcox")
#ggsave(filename = "Substitution_rate_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)

#ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_epitope_related),levels=c("TRUE","FALSE")),y = pN_pS,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("Epitope sites?") + ylab("pN/pS") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$pN_pS)+1e-2,1e-2),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$pN_pS)+1e-2))
#ggsave(filename = "pN_pS_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
#ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_epitope_related),levels=c("TRUE","FALSE")),y = dN_dS,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("Epitope sites?") + ylab("dN/dS") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$dN_dS)+1e-4,1e-4),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$dN_dS)+1e-4))
#ggsave(filename = "dN_dS_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
#ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_epitope_related),levels=c("TRUE","FALSE")),y = alpha_MK_Test,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("Epitope sites?") + ylab("McDonald-Kreitman test \U003B1") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$alpha_MK_Test)+1e-2,1e-2),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon$alpha_MK_Test)+1e-2))
#ggsave(filename = "MK_test_alpha_in_epitope_vs_out_of_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)

#compare mutation and substitution rate in epitopes vs outside (split by genomic region)
v_orfs_of_interest <- c("orf1a","orf1b","S","E","M","N")
df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon <- data.frame(Sample=rep(lst_samples_NCBI_SRA_amplicon,each=length(c(TRUE,FALSE))*length(v_orfs_of_interest)),is_epitope_related=rep(c(TRUE,FALSE),each=length(lst_samples_NCBI_SRA_amplicon)*length(v_orfs_of_interest)),ORF=rep(v_orfs_of_interest,length(lst_samples_NCBI_SRA_amplicon)*length(c(TRUE,FALSE))),nb_mutations=0,nb_fixed_mutations=0,mut_rate=0,subst_rate=0,stringsAsFactors = F)
nb_cores <- nb_cpus
lst_splits <- split(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon), ceiling(seq_along(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon))/(nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon)/nb_cores)))
the_f_parallel_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  df_metrics_current_subset <- df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon[the_vec,]
  count_iter <- 0
  for (the_i in 1:nrow(df_metrics_current_subset)){
    df_depth_NCBI_SRA_amplicon_current_sample <- read.csv2(file = paste0(depth_data_wp,"depth_report_NCBI_SRA_amplicon/df_depth_NCBI_SRA_amplicon_",df_metrics_current_subset$Sample[the_i],".csv"),sep = ",",header = F,stringsAsFactors = FALSE)
    colnames(df_depth_NCBI_SRA_amplicon_current_sample) <- c("sample","position","depth")
    df_depth_NCBI_SRA_amplicon_current_sample$ORF <- vapply(X = df_depth_NCBI_SRA_amplicon_current_sample$position,FUN = find_ORF_of_mutation,FUN.VALUE = c(""))
    df_depth_NCBI_SRA_amplicon_current_sample <- unique(df_depth_NCBI_SRA_amplicon_current_sample)
    v_currentsample_positions_enough_covered <- subset(df_depth_NCBI_SRA_amplicon_current_sample,(depth>=min_cov)&(ORF==df_metrics_current_subset$ORF[the_i]))$position
    if (df_metrics_current_subset$is_epitope_related[the_i]){
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_unique_epitope_positions,v_currentsample_positions_enough_covered)
    }else{
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_non_epitope_sites,v_currentsample_positions_enough_covered)
    }
    df_metrics_current_subset$nb_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(is_fixed=="No")&(!Position%in%v_positions_inter)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_fixed_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(is_fixed=="Yes")&(is_prevalence_above_transmission_threshold)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] <- length(v_the_sites_with_enough_cov_for_current_category)
    df_metrics_current_subset$mut_rate[the_i] <- df_metrics_current_subset$nb_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$subst_rate[the_i] <- df_metrics_current_subset$nb_fixed_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$Nb_ss[the_i] <- sum(unname(vapply(X = v_the_sites_with_enough_cov_for_current_category,FUN = calculate_nb_ss_position_in_genome,FUN.VALUE = c(0))),na.rm=T)
    df_metrics_current_subset$Nb_nss[the_i] <- df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] - df_metrics_current_subset$Nb_ss[the_i]
    df_metrics_current_subset$within_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$within_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(ORF==df_metrics_current_subset$ORF[the_i])&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    
    if (the_i%%100==0){
      print(paste0("[Epitope sites vs others by ORF (NCBI)] Core ",i_cl,": Step ",the_i," done out of ",nrow(df_metrics_current_subset),"!"))
    }
  }
  df_metrics_current_subset$pN <- df_metrics_current_subset$within_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$pS <- df_metrics_current_subset$within_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$pN_pS <- ifelse(test=df_metrics_current_subset$pS==0,yes=NA,no=df_metrics_current_subset$pN/df_metrics_current_subset$pS)
  df_metrics_current_subset$dN <- df_metrics_current_subset$between_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$dS <- df_metrics_current_subset$between_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$dN_dS <- ifelse(test=df_metrics_current_subset$dS==0,yes=NA,no=df_metrics_current_subset$dN/df_metrics_current_subset$dS)
  df_metrics_current_subset$alpha_MK_Test <- ifelse(test=df_metrics_current_subset$dN_dS==0,yes=NA,no=(1-((df_metrics_current_subset$pN_pS)/(df_metrics_current_subset$dN_dS))))
  
  return(df_metrics_current_subset)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_Evo_rates_analyses.txt"))
registerDoParallel(cl)
df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon(i_cl)
stopCluster(cl)
# #saveRDS(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,paste0(output_workspace,"df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.rds"))

df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$label <- ifelse(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$is_epitope_related,yes="Epitope sites",no="Other sites")
#within-host mutation rate by ORF (NCBI_SRA_amplicon)
ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,aes(x=factor(ORF,levels=v_orfs_of_interest),y = mut_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("ORF") + ylab("Within-host mutation rate (Count / Length)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$mut_rate)+1e-2,1e-2),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$mut_rate)+1e-2)) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1)
ggsave(filename = "Mutation_rate_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$mut_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$mut_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$mut_rate)$p.value
#substitution rate by ORF (NCBI_SRA_amplicon)
ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,aes(x=factor(ORF,levels=v_orfs_of_interest),y = subst_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("ORF") + ylab("Substitution rate (Count / Length)") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$subst_rate)+1e-4,1e-4),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$subst_rate)+1e-4)) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1)
ggsave(filename = "Substitution_rate_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$subst_rate);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$subst_rate,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$subst_rate)$p.value

#pN/pS by ORF (NCBI_SRA_amplicon)
ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,aes(x=factor(ORF,levels=v_orfs_of_interest),y = pN_pS,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("ORF") + ylab("pN/pS") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$pN_pS)+1e-2,1e-2),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$pN_pS)+1e-2))
ggsave(filename = "pN_pS_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$pN_pS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$pN_pS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$pN_pS)$p.value
#dN/dS by ORF (NCBI_SRA_amplicon)
ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,aes(x=factor(ORF,levels=v_orfs_of_interest),y = dN_dS,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("ORF") + ylab("dN/dS") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$dN_dS)+1e-4,1e-4),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$dN_dS)+1e-4))
ggsave(filename = "dN_dS_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$dN_dS);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$dN_dS,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$dN_dS)$p.value

#M-K test alpha by ORF (NCBI_SRA_amplicon)
ggplot(data = df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,aes(x=factor(ORF,levels=v_orfs_of_interest),y = alpha_MK_Test,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("ORF") + ylab("McDonald-Kreitman test \U003B1") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) #+ scale_y_continuous(breaks=seq(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$alpha_MK_Test)+1e-4,1e-4),limits = c(0,max(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon$alpha_MK_Test)+1e-4))
ggsave(filename = "MK_test_alpha_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1a")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="orf1b")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="S")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="E")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="M")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(is_epitope_related))$alpha_MK_Test,subset(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,(ORF=="N")&(!is_epitope_related))$alpha_MK_Test)$p.value

#compare mutation and substitution rate in epitopes vs outside (S protein domains)
df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon <- data.frame(Sample=rep(lst_samples_NCBI_SRA_amplicon,each=length(c(TRUE,FALSE))*length(v_S_protein_domains)),is_epitope_related=rep(c(TRUE,FALSE),each=length(lst_samples_NCBI_SRA_amplicon)*length(v_S_protein_domains)),S_protein_domain=rep(v_S_protein_domains,length(lst_samples_NCBI_SRA_amplicon)*length(c(TRUE,FALSE))),nb_mutations=0,nb_fixed_mutations=0,mut_rate=0,subst_rate=0,stringsAsFactors = F)
nb_cores <- nb_cpus
lst_splits <- split(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon), ceiling(seq_along(1:nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon))/(nrow(df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon)/nb_cores)))
the_f_parallel_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  df_metrics_current_subset <- df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon[the_vec,]
  count_iter <- 0
  for (the_i in 1:nrow(df_metrics_current_subset)){
    df_depth_NCBI_SRA_amplicon_current_sample <- read.csv2(file = paste0(depth_data_wp,"depth_report_NCBI_SRA_amplicon/df_depth_NCBI_SRA_amplicon_",df_metrics_current_subset$Sample[the_i],".csv"),sep = ",",header = F,stringsAsFactors = FALSE)
    colnames(df_depth_NCBI_SRA_amplicon_current_sample) <- c("sample","position","depth")
    df_depth_NCBI_SRA_amplicon_current_sample$ORF <- vapply(X = df_depth_NCBI_SRA_amplicon_current_sample$position,FUN = find_ORF_of_mutation,FUN.VALUE = c(""))
    df_depth_NCBI_SRA_amplicon_current_sample <- subset(df_depth_NCBI_SRA_amplicon_current_sample,ORF=="S")
    df_depth_NCBI_SRA_amplicon_current_sample$pos_in_protein <- ceiling((df_depth_NCBI_SRA_amplicon_current_sample$position - v_start_genes["S"] + 1)/3)
    df_depth_NCBI_SRA_amplicon_current_sample$S_protein_domain <- vapply(X = df_depth_NCBI_SRA_amplicon_current_sample$pos_in_protein,FUN = find_S_protein_domain_of_mutation,FUN.VALUE = c(""))
    df_depth_NCBI_SRA_amplicon_current_sample <- unique(df_depth_NCBI_SRA_amplicon_current_sample)
    v_currentsample_positions_enough_covered <- subset(df_depth_NCBI_SRA_amplicon_current_sample,(depth>=min_cov)&(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i]))$position
    if (df_metrics_current_subset$is_epitope_related[the_i]){
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_unique_epitope_positions,v_currentsample_positions_enough_covered)
    }else{
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_non_epitope_sites,v_currentsample_positions_enough_covered)
    }
    df_metrics_current_subset$nb_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(is_fixed=="No")&(!Position%in%v_positions_inter)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_fixed_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(is_fixed=="Yes")&(is_prevalence_above_transmission_threshold)&(is_epitope_related==df_metrics_current_subset$is_epitope_related[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] <- length(v_the_sites_with_enough_cov_for_current_category)
    df_metrics_current_subset$mut_rate[the_i] <- df_metrics_current_subset$nb_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$subst_rate[the_i] <- df_metrics_current_subset$nb_fixed_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$Nb_ss[the_i] <- sum(unname(vapply(X = v_the_sites_with_enough_cov_for_current_category,FUN = calculate_nb_ss_position_in_genome,FUN.VALUE = c(0))),na.rm=T)
    df_metrics_current_subset$Nb_nss[the_i] <- df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] - df_metrics_current_subset$Nb_ss[the_i]
    df_metrics_current_subset$within_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$within_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(S_protein_domain==df_metrics_current_subset$S_protein_domain[the_i])&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    
    if (the_i%%100==0){
      print(paste0("[Epitope sites vs others by S_protein_domain (NCBI)] Core ",i_cl,": Step ",the_i," done out of ",nrow(df_metrics_current_subset),"!"))
    }
  }
  df_metrics_current_subset$pN <- df_metrics_current_subset$within_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$pS <- df_metrics_current_subset$within_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$pN_pS <- ifelse(test=df_metrics_current_subset$pS==0,yes=NA,no=df_metrics_current_subset$pN/df_metrics_current_subset$pS)
  df_metrics_current_subset$dN <- df_metrics_current_subset$between_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$dS <- df_metrics_current_subset$between_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$dN_dS <- ifelse(test=df_metrics_current_subset$dS==0,yes=NA,no=df_metrics_current_subset$dN/df_metrics_current_subset$dS)
  df_metrics_current_subset$alpha_MK_Test <- ifelse(test=df_metrics_current_subset$dN_dS==0,yes=NA,no=(1-((df_metrics_current_subset$pN_pS)/(df_metrics_current_subset$dN_dS))))
  
  return(df_metrics_current_subset)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_Evo_rates_analyses.txt"))
registerDoParallel(cl)
df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon(i_cl)
stopCluster(cl)
saveRDS(df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon,paste0(output_workspace,"df_metrics_in_epitope_vs_out_of_epitopes_by_ORF_NCBI_SRA_amplicon.rds"))

df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$label <- ifelse(df_metrics_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$is_epitope_related,yes="Epitope sites",no="Other sites")
#within-host mutation rate by S_protein_domain (NCBI_SRA_amplicon)
ggplot(data = df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,aes(x=factor(S_protein_domain,levels=v_S_protein_domains),y = mut_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_point() + xlab("S protein domain") + ylab("Within-host mutation rate (Count / Length)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + scale_y_continuous(breaks=seq(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$mut_rate)+1e-2,1e-2),limits = c(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$mut_rate)+1e-2)) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1)
ggsave(filename = "Mutation_rate_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$mut_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$mut_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$mut_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$mut_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$mut_rate)$p.value
#substitution rate by S_protein_domain (NCBI_SRA_amplicon)
ggplot(data = df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,aes(x=factor(S_protein_domain,levels=v_S_protein_domains),y = subst_rate,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("S protein domain") + ylab("Substitution rate (Count / Length)") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_y_continuous(breaks=seq(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$subst_rate)+1e-4,1e-4),limits = c(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$subst_rate)+1e-4)) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1)
ggsave(filename = "Substitution_rate_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$subst_rate)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$subst_rate);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$subst_rate); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$subst_rate,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$subst_rate)$p.value

# pN/pS by S_protein_domain (NCBI_SRA_amplicon)
ggplot(data = df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,aes(x=factor(S_protein_domain,levels=v_S_protein_domains),y = pN_pS,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("S protein domain") + ylab("pN/pS") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) # + scale_y_continuous(breaks=seq(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$pN_pS)+1e-2,1e-2),limits = c(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$pN_pS)+1e-2))
ggsave(filename = "pN_pS_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$pN_pS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$pN_pS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$pN_pS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$pN_pS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$pN_pS)$p.value
#dN/dS by S_protein_domain (NCBI_SRA_amplicon)
ggplot(data = df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,aes(x=factor(S_protein_domain,levels=v_S_protein_domains),y = dN_dS,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("S protein domain") + ylab("dN/dS") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) # + scale_y_continuous(breaks=seq(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$dN_dS)+1e-4,1e-4),limits = c(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$dN_dS)+1e-4))
ggsave(filename = "dN_dS_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$dN_dS)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$dN_dS);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$dN_dS); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$dN_dS,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$dN_dS)$p.value

#McDonald-Kreitman test alpha by S_protein_domain (NCBI_SRA_amplicon)
ggplot(data = df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,aes(x=factor(S_protein_domain,levels=v_S_protein_domains),y = alpha_MK_Test,fill=as.character(is_epitope_related))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("S protein domain") + ylab("McDonald-Kreitman test alpha") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "kruskal") + facet_wrap(~label,ncol=1) # + scale_y_continuous(breaks=seq(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$alpha_MK_Test)+1e-4,1e-4),limits = c(0,max(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon$alpha_MK_Test)+1e-4))
ggsave(filename = "alpha_MK_Test_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="NTD")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="RBD")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD1")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD2")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CR")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR1")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="CH-BH")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="SD3")&(!is_epitope_related))$alpha_MK_Test)$p.value
mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$alpha_MK_Test);mean(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$alpha_MK_Test); wilcox.test(subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(is_epitope_related))$alpha_MK_Test,subset(df_rates_in_epitope_vs_out_of_epitopes_by_S_protein_domain_NCBI_SRA_amplicon,(S_protein_domain=="HR2-TM-CT")&(!is_epitope_related))$alpha_MK_Test)$p.value

##############
#Identify peptides that are significantly enriched in one group
#df_epitope_frequency_per_group <- aggregate(df_epitopes$patient_ID,by=list(peptide_id=df_epitopes$peptide_id,Group=df_epitopes$Group),FUN=function(x) return(length(x)))
# mtx_epitope_frequency_per_group <- (as.matrix((reshape2::acast(as.data.frame(table(df_epitopes$peptide_id,df_epitopes$Group)), Var1~Var2, value.var="Freq"))))
# # range(p.adjust(vapply(X = rownames(mtx_epitope_frequency_per_group),FUN = function(the_pep) chisq.test(mtx_epitope_frequency_per_group[the_pep,])$p.value,FUN.VALUE = c(0.0)),method = "fdr"))
# df_data_indval_pa <- as.data.frame(t(ifelse(as.matrix((reshape2::acast(as.data.frame(table(df_epitopes$peptide_id,df_epitopes$patient_ID)), Var1~Var2, value.var="Freq")))>0, yes = 1,no=0)))
# res_indvalpa_analysis <- multipatt(x = df_data_indval_pa,cluster = v_group_patient[rownames(df_data_indval_pa)],func = "IndVal.g",duleg=TRUE, max.order = 1,control = how(nperm = 9999),print.perm = TRUE)
# df_res_indvalpa_analysis <- as.data.frame(res_indvalpa_analysis$sign)
# df_signif_res_indval_pa <- subset(df_res_indvalpa_analysis,p.value<=0.05)
# df_signif_res_indval_pa$Genomic_start <- unname(vapply(X = rownames(df_signif_res_indval_pa),FUN = function(x) subset(df_epitopes,peptide_id==x)$Genomic_start[1],FUN.VALUE = c(0)))
# df_signif_res_indval_pa$Genomic_End <- unname(vapply(X = rownames(df_signif_res_indval_pa),FUN = function(x) subset(df_epitopes,peptide_id==x)$Genomic_End[1],FUN.VALUE = c(0)))
# df_signif_res_indval_pa$mutation_rate_NCBI_SRA_amplicon <- unname(vapply(X = 1:nrow(df_signif_res_indval_pa),FUN = function(i) length(unique(subset(df_variants_NCBI_SRA_amplicon,(Position>=(df_signif_res_indval_pa$Genomic_start[i]))&(Position<=(df_signif_res_indval_pa$Genomic_End[i])))$mutation_name))/(df_signif_res_indval_pa$Genomic_End[i]-df_signif_res_indval_pa$Genomic_start[i]+1),FUN.VALUE = c(0)))
# df_signif_res_indval_pa$substitution_rate_NCBI_SRA_amplicon <- unname(vapply(X = 1:nrow(df_signif_res_indval_pa),FUN = function(i) length(unique(subset(df_variants_NCBI_SRA_amplicon,(Position>=(df_signif_res_indval_pa$Genomic_start[i]))&(Position<=(df_signif_res_indval_pa$Genomic_End[i]))&(is_fixed=="Yes"))$mutation_name))/(df_signif_res_indval_pa$Genomic_End[i]-df_signif_res_indval_pa$Genomic_start[i]+1),FUN.VALUE = c(0)))
# names(df_signif_res_indval_pa)[4:6] <- c("Group_with_strongest_association","stat","p_value_association")
# df_signif_res_indval_pa$sequence <- v_seq_peptide[rownames(df_signif_res_indval_pa)]
# df_signif_res_indval_pa_to_save <- df_signif_res_indval_pa
# df_signif_res_indval_pa_to_save$peptide_id <- rownames(df_signif_res_indval_pa_to_save)
# rownames(df_signif_res_indval_pa_to_save) <- NULL
# df_signif_res_indval_pa_to_save <- df_signif_res_indval_pa_to_save[,c("peptide_id","sequence","Group_with_strongest_association","p_value_association","Genomic_start","Genomic_End","mutation_rate_NCBI_SRA_amplicon","substitution_rate_NCBI_SRA_amplicon")]
# #write.table(x=df_signif_res_indval_pa_to_save,file = paste0(output_workspace,"Epitopes_with_significant_association_to_a_specific_group.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

#import the 4 ancient hCoVs epitopes that mapped to SARS-CoV-2 genome
df_epitopes_HKU1 <- read.csv2(file = paste0(output_workspace,"Epitopes_HKU1_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes_HKU1$Genomic_start <- as.integer(df_epitopes_HKU1$Genomic_start)
df_epitopes_HKU1$Genomic_End <- as.integer(df_epitopes_HKU1$Genomic_End)
df_epitopes_HKU1$Mapping_region <- vapply(X = df_epitopes_HKU1$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes_HKU1$peptide_id <- NULL
df_epitopes_HKU1$hCoVs_overlap <- NULL
df_epitopes_HKU1 <- subset(df_epitopes_HKU1,(!is.na(Genomic_start))&(!is.na(Genomic_start)))
df_epitopes_HKU1$Mapping_region <- factor(as.character(df_epitopes_HKU1$Mapping_region),intersect(v_orfs,df_epitopes_HKU1$Mapping_region))
df_epitopes_HKU1 <- subset(df_epitopes_HKU1,vapply(X = 1:nrow(df_epitopes_HKU1),FUN = function(i) (return(grepl(pattern = df_epitopes_HKU1$Mapping_region[i],x = df_epitopes_HKU1$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes_HKU1$Group <- as.character(df_epitopes_HKU1$Group)
df_epitopes_HKU1$RFU <- as.numeric(df_epitopes_HKU1$RFU)
df_epitopes_HKU1 <- subset(df_epitopes_HKU1, RFU>=1000)

df_epitopes_OC43 <- read.csv2(file = paste0(output_workspace,"Epitopes_OC43_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes_OC43$Genomic_start <- as.integer(df_epitopes_OC43$Genomic_start)
df_epitopes_OC43$Genomic_End <- as.integer(df_epitopes_OC43$Genomic_End)
df_epitopes_OC43$Mapping_region <- vapply(X = df_epitopes_OC43$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes_OC43$peptide_id <- NULL
df_epitopes_OC43$hCoVs_overlap <- NULL
df_epitopes_OC43 <- subset(df_epitopes_OC43,(!is.na(Genomic_start))&(!is.na(Genomic_start)))
df_epitopes_OC43$Mapping_region <- factor(as.character(df_epitopes_OC43$Mapping_region),intersect(v_orfs,df_epitopes_OC43$Mapping_region))
df_epitopes_OC43 <- subset(df_epitopes_OC43,vapply(X = 1:nrow(df_epitopes_OC43),FUN = function(i) (return(grepl(pattern = df_epitopes_OC43$Mapping_region[i],x = df_epitopes_OC43$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes_OC43$Group <- as.character(df_epitopes_OC43$Group)
df_epitopes_OC43$RFU <- as.numeric(df_epitopes_OC43$RFU)
df_epitopes_OC43 <- subset(df_epitopes_OC43, RFU>=1000)

df_epitopes_NL63 <- read.csv2(file = paste0(output_workspace,"Epitopes_NL63_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes_NL63$Genomic_start <- as.integer(df_epitopes_NL63$Genomic_start)
df_epitopes_NL63$Genomic_End <- as.integer(df_epitopes_NL63$Genomic_End)
df_epitopes_NL63$Mapping_region <- vapply(X = df_epitopes_NL63$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes_NL63$peptide_id <- NULL
df_epitopes_NL63$hCoVs_overlap <- NULL
df_epitopes_NL63 <- subset(df_epitopes_NL63,(!is.na(Genomic_start))&(!is.na(Genomic_start)))
df_epitopes_NL63$Mapping_region <- factor(as.character(df_epitopes_NL63$Mapping_region),intersect(v_orfs,df_epitopes_NL63$Mapping_region))
df_epitopes_NL63 <- subset(df_epitopes_NL63,vapply(X = 1:nrow(df_epitopes_NL63),FUN = function(i) (return(grepl(pattern = df_epitopes_NL63$Mapping_region[i],x = df_epitopes_NL63$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes_NL63$Group <- as.character(df_epitopes_NL63$Group)
df_epitopes_NL63$RFU <- as.numeric(df_epitopes_NL63$RFU)
df_epitopes_NL63 <- subset(df_epitopes_NL63, RFU>=1000)

df_epitopes_229E <- read.csv2(file = paste0(output_workspace,"Epitopes_229E_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes_229E$Genomic_start <- as.integer(df_epitopes_229E$Genomic_start)
df_epitopes_229E$Genomic_End <- as.integer(df_epitopes_229E$Genomic_End)
df_epitopes_229E$Mapping_region <- vapply(X = df_epitopes_229E$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes_229E$peptide_id <- NULL
df_epitopes_229E$hCoVs_overlap <- NULL
df_epitopes_229E <- subset(df_epitopes_229E,(!is.na(Genomic_start))&(!is.na(Genomic_start)))
df_epitopes_229E$Mapping_region <- factor(as.character(df_epitopes_229E$Mapping_region),intersect(v_orfs,df_epitopes_229E$Mapping_region))
df_epitopes_229E <- subset(df_epitopes_229E,vapply(X = 1:nrow(df_epitopes_229E),FUN = function(i) (return(grepl(pattern = df_epitopes_229E$Mapping_region[i],x = df_epitopes_229E$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes_229E$Group <- as.character(df_epitopes_229E$Group)
df_epitopes_229E$RFU <- as.numeric(df_epitopes_229E$RFU)
df_epitopes_229E <- subset(df_epitopes_229E, RFU>=1000)

#hCoVs overlap
v_more_ancient_hcovs <- c("HKU1","OC43","NL63","229E")
#lst_df_epitopes_more_ancient_hcovs <- list(HKU1=df_epitopes_HKU1,OC43=df_epitopes_OC43,NL63=df_epitopes_NL63,"229E"=df_epitopes_229E)
df_epitope_overlap_ancient_hcovs <- rbind(df_epitopes_HKU1[,c("Peptide","Genomic_start","Mapping_region","patient_ID")],df_epitopes_OC43[,c("Peptide","Genomic_start","Mapping_region","patient_ID")],df_epitopes_NL63[,c("Peptide","Genomic_start","Mapping_region","patient_ID")],df_epitopes_229E[,c("Peptide","Genomic_start","Mapping_region","patient_ID")])
df_epitope_overlap_ancient_hcovs <- unique(df_epitope_overlap_ancient_hcovs)
df_epitope_overlap_ancient_hcovs <- aggregate(df_epitope_overlap_ancient_hcovs$Peptide,by=list(df_epitope_overlap_ancient_hcovs$Peptide,df_epitope_overlap_ancient_hcovs$Genomic_start,df_epitope_overlap_ancient_hcovs$Mapping_region),FUN=function(x) length(x))
names(df_epitope_overlap_ancient_hcovs) <- c("Peptide","Genomic_start","Mapping_region","Count")
df_epitope_overlap_ancient_hcovs$hcovs_overlap <- paste0("SARS-CoV-2/",vapply(X = df_epitope_overlap_ancient_hcovs$Peptide,FUN = function(the_pep) return(paste0(v_more_ancient_hcovs[c(the_pep%in%df_epitopes_HKU1$Peptide,the_pep%in%df_epitopes_OC43$Peptide,the_pep%in%df_epitopes_NL63$Peptide,the_pep%in%df_epitopes_229E$Peptide)],collapse="/")),FUN.VALUE = c("")))
df_epitope_overlap_ancient_hcovs$hcovs_overlap <- ifelse(test = df_epitope_overlap_ancient_hcovs$hcovs_overlap==paste0(c("SARS-CoV-2",v_more_ancient_hcovs),collapse = "/"),yes="All",no = df_epitope_overlap_ancient_hcovs$hcovs_overlap)

#cairo_ps(filename = paste0(output_workspace,"hCoVs_overlap_for_ALL_Epitopes.eps"), height = 7.8, width=5.9,fallback_resolution = 1200)
#grid.arrange(ggplot(data = df_epitope_overlap_ancient_hcovs) + geom_linerange(mapping = aes(x = Genomic_start,ymin=0, ymax = Count)) + geom_point(mapping = aes(x=Genomic_start,y=Count,col=hcovs_overlap),size=2)+ theme_bw()  + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),legend.position = "none",axis.text.y = element_text(size=12)) + ylab("Prevalence in the dataset (n=15)") + xlab("") + guides(col=guide_legend(title="hCoVs overlap")) + scale_y_continuous(limits = c(-1.5,15),breaks=seq(0,15,1))+ scale_x_continuous(limits = c(0,nchar(genome_refseq)),breaks=seq(0,nchar(genome_refseq),5000)),#ggplot(data = df_epitope_overlap_ancient_hcovs) + geom_smooth(mapping = aes(x=Genomic_start,y=Count,col=hcovs_overlap),method = "loess",size=0.1)+ theme_bw()  + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position = "none") + ylab("LOESS") + xlab("Position")+ guides(col=guide_legend(title="hCoVs overlap")) + scale_y_continuous(limits = c(0,15),breaks=seq(0,15,1))+ scale_x_continuous(limits = c(0,nchar(genome_refseq)),breaks=seq(0,nchar(genome_refseq),5000)), ncol=1)
#dev.off()
#cairo_ps(filename = paste0(output_workspace,"hCoVs_overlap_for_Epitopes_that_are_present_in_most_of_the_samples.eps"), height = 7.8, width=5.9,fallback_resolution = 1200)
#grid.arrange(ggplot(data = subset(df_epitope_overlap_ancient_hcovs,Count>(max(df_epitope_overlap_ancient_hcovs$Count)/2))) + geom_linerange(mapping = aes(x = Genomic_start,ymin=0, ymax = Count)) + geom_point(mapping = aes(x=Genomic_start,y=Count,col=hcovs_overlap),size=3)+ theme_bw()  + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position = "none") + ylab("Prevalence in the dataset (n=15)") + xlab("") + guides(col=guide_legend(title="hCoVs overlap")) + scale_y_continuous(limits = c(-1.5,15),breaks=seq(0,15,1))+ scale_x_continuous(limits = c(0,nchar(genome_refseq)),breaks=seq(0,nchar(genome_refseq),5000)),#ggplot(data =subset(df_epitope_overlap_ancient_hcovs,Count>(max(df_epitope_overlap_ancient_hcovs$Count)/2))) + geom_smooth(mapping = aes(x=Genomic_start,y=Count,col=hcovs_overlap),method = "loess",size=0.1)+ theme_bw()  + theme(title =  element_text(size=12),legend.text = element_text(size = 8),legend.position = "none",axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("LOESS") + xlab("Position")+ guides(col=guide_legend(title="hCoVs overlap")) + scale_y_continuous(limits = c(0,15),breaks=seq(0,15,1))+ scale_x_continuous(limits = c(0,nchar(genome_refseq)),breaks=seq(0,nchar(genome_refseq),5000)), ncol=1)
#dev.off()

df_epitopes_ancient_hcovs <- rbind(df_epitopes_HKU1,df_epitopes_OC43,df_epitopes_NL63,df_epitopes_229E)
df_all_epitopes_ancient_hcovs <- unique(data.frame(Peptide=c(read.csv2(file = paste0(output_workspace,"Epitopes_HKU1_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Peptide,read.csv2(file = paste0(output_workspace,"Epitopes_OC43_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Peptide,read.csv2(file = paste0(output_workspace,"Epitopes_NL63_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Peptide, read.csv2(file = paste0(output_workspace,"Epitopes_229E_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Peptide),Genomic_start=as.integer(c(read.csv2(file = paste0(output_workspace,"Epitopes_HKU1_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Genomic_start,read.csv2(file = paste0(output_workspace,"Epitopes_OC43_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Genomic_start,read.csv2(file = paste0(output_workspace,"Epitopes_NL63_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Genomic_start, read.csv2(file = paste0(output_workspace,"Epitopes_229E_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)$Genomic_start),stringsAsFactors = FALSE)))
df_all_epitopes_ancient_hcovs$Mapping_region <- ifelse(test = is.na(df_all_epitopes_ancient_hcovs$Genomic_start),yes = "NA",no = vapply(X = df_all_epitopes_ancient_hcovs$Genomic_start,FUN = find_ORF_of_mutation,FUN.VALUE = c("")))
v_lst_id_epitope_seq_ancient_hcovs_in_order <- 1:length(sort(unique(df_all_epitopes_ancient_hcovs$Peptide),decreasing = FALSE))
names(v_lst_id_epitope_seq_ancient_hcovs_in_order) <- sort(unique(df_all_epitopes_ancient_hcovs$Peptide),decreasing = FALSE)
df_epitopes_ancient_hcovs$peptide_id <- paste0("EhCov_",df_epitopes_ancient_hcovs$Mapping_region,"_",unname(v_lst_id_epitope_seq_ancient_hcovs_in_order[df_epitopes_ancient_hcovs$Peptide]))
df_epitope_overlap_ancient_hcovs$peptide_id <- paste0("EhCov_",df_epitope_overlap_ancient_hcovs$Mapping_region,"_",unname(v_lst_id_epitope_seq_ancient_hcovs_in_order[df_epitope_overlap_ancient_hcovs$Peptide]))
df_all_epitopes_ancient_hcovs$peptide_id <- paste0("EhCov_",as.character(df_all_epitopes_ancient_hcovs$Mapping_region),"_",unname(v_lst_id_epitope_seq_ancient_hcovs_in_order[df_all_epitopes_ancient_hcovs$Peptide]))
rownames(df_all_epitopes_ancient_hcovs) <- df_all_epitopes_ancient_hcovs$peptide_id

df_sars_cov_2_epitopes <- readRDS(file = paste0(output_workspace,"df_sars_cov_2_epitopes.rds"))
df_sars_cov_2_epitopes <- subset(df_sars_cov_2_epitopes, peptide_id %in% df_epitopes$peptide_id) #filter for RFU >= 1000

#RFU across ORFs and groups
#p <- ggboxplot(data = df_epitopes_ancient_hcovs, x = "Mapping_region", y="RFU",color = "Group",add = "jitter")+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Genomic region") +scale_y_continuous(limits = c(0,max(df_epitopes_ancient_hcovs$RFU)+10000),breaks=seq(0,max(df_epitopes_ancient_hcovs$RFU)+10000,10000))
#facet(p +  stat_compare_means(), facet.by = "Group", ncol = 1)
#ggsave(filename = "Ancient_hCoVs_RFU_by_ORF_and_Groups.png", path=output_workspace, width = 20, height = 20, units = "cm",dpi = 1200)

#RFU across groups
#ggboxplot(data = df_epitopes_ancient_hcovs, x = "Group", y="RFU",color = "Group",add = "jitter") +  stat_compare_means() + theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Group")
#ggsave(filename = "Ancient_hCoVs_RFU_by_Groups.png", path=output_workspace, width = 20, height = 15, units = "cm",dpi = 1200)

#RFU across Antibody and groups
#p <- ggboxplot(data = df_epitopes_ancient_hcovs, x = "Antibody", y="RFU",color = "Group",add = "jitter")+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("RFU") + xlab("Antibody") +scale_y_continuous(limits = c(0,max(df_epitopes_ancient_hcovs$RFU)+10000),breaks=seq(0,max(df_epitopes_ancient_hcovs$RFU)+10000,10000))
#facet(p +  stat_compare_means(), facet.by = "Group", ncol = 1)
#ggsave(filename = "Ancient_hCoVs_RFU_by_Antibody_and_Groups.png", path=output_workspace, width = 20, height = 20, units = "cm",dpi = 1200)

#create Jalview annotation file
df_positions_epitope_proteome_genome <- NULL
df_positions_epitope_proteome_genome$Genomic_Position <- 1:nchar(genome_refseq)
df_positions_epitope_proteome_genome <- as.data.frame(df_positions_epitope_proteome_genome)
df_positions_epitope_proteome_genome$ORF <- unname(vapply(X = df_positions_epitope_proteome_genome$Genomic_Position,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c("")))
df_positions_epitope_proteome_genome <- subset(df_positions_epitope_proteome_genome,!is.na(ORF))
df_positions_epitope_proteome_genome$pos_in_protein <- unname(ceiling((df_positions_epitope_proteome_genome$Genomic_Position - v_start_orfs[df_positions_epitope_proteome_genome$ORF] + 1)/3))
df_positions_epitope_proteome_genome$is_epitope_related <- df_positions_epitope_proteome_genome$Genomic_Position%in% v_unique_epitope_positions
for (current_orf in c("orf1a","orf1b","S","E","M","N")){
  v_pos_in_prot_to_highlight <- sort(unique(subset(df_positions_epitope_proteome_genome,(ORF==current_orf)&(is_epitope_related))$pos_in_protein))
  df_align_colour_current_orf <- data.frame(charact_align = unlist(strsplit(x = seqinr::getSequence(object = toupper(read.fasta(paste0(output_workspace,"hCoVs_alignments/",paste0(current_orf,"_aligned.fasta")),seqtype = "AA",as.string = TRUE,forceDNAtolower = FALSE)$`SARS-CoV-2`),as.string = TRUE)[[1]],split="")),stringsAsFactors = F)
  df_align_colour_current_orf$pos_in_align <- 1:(nrow(df_align_colour_current_orf))
  df_align_colour_current_orf <- subset(df_align_colour_current_orf,charact_align!="-")
  df_align_colour_current_orf$pos_in_prot <- 1:nrow(df_align_colour_current_orf)
  current_v_pos_in_alignment_to_highlight <- subset(df_align_colour_current_orf,pos_in_prot%in%v_pos_in_prot_to_highlight)$pos_in_align
  sink(paste0(output_workspace,"hCoVs_alignments/",current_orf,"_annotations.txt"))
  cat("JALVIEW_ANNOTATION\n")
  i <- 1
  for (current_pos in current_v_pos_in_alignment_to_highlight){
    cat(paste0("SEQUENCE_GROUP\tGroup_",i,"\t",current_pos,"\t",current_pos,"\t*\n"),append = T)
    i <- i + 1
  }
  i <- 1
  for (current_pos in current_v_pos_in_alignment_to_highlight){
    cat(paste0("PROPERTIES\tGroup_",i,"\toutlineColour=red\n"),append = T)
    i <- i + 1
  }
  sink()
}

# #pairwise alignments between SARS-CoV-2 epitopes and 4 endemic hCoVs epitopes
# nb_cores <- nb_cpus
# lst_splits <- split(1:nrow(df_sars_cov_2_epitopes), ceiling(seq_along(1:nrow(df_sars_cov_2_epitopes))/(nrow(df_sars_cov_2_epitopes)/nb_cores)))
# the_f_parallel <- function(i_cl){
#   the_vec<- lst_splits[[i_cl]]
#   df_pairwise_align_score <- NULL
#   data("PAM120")
#   count_iter <- 0
#   for (i_epitope_sars_cov_2 in the_vec){
#     df_pairwise_align_score <- rbind(df_pairwise_align_score,data.frame(id_sars_cov_2_epitope=rownames(df_sars_cov_2_epitopes)[i_epitope_sars_cov_2],id_EhCovs_epitope=rownames(df_all_epitopes_ancient_hcovs),score=vapply(X = rownames(df_all_epitopes_ancient_hcovs),FUN = function(current_ep_Ehcovs) return(pairwiseAlignment(pattern = AAString(df_sars_cov_2_epitopes$Peptide[i_epitope_sars_cov_2]),subject = AAString(df_all_epitopes_ancient_hcovs[current_ep_Ehcovs,"Peptide"]),  substitutionMatrix = PAM120,type="global",gapOpening = 6, gapExtension = 4,scoreOnly=T)[1]),FUN.VALUE = c(0)),stringsAsFactors = FALSE))
#     count_iter <- count_iter + 1
#     if (i_epitope_sars_cov_2%%10==0){
#       print(paste0("Core ",i_cl,": ",count_iter," SARS-CoV-2 epitopes fully compared out of ",length(the_vec)))
#     }
#   }
#   return(df_pairwise_align_score)
# }
# cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_pairwise_align_sc2_vs_EhCoVs.txt"))
# registerDoParallel(cl)
# df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel(i_cl)
# stopCluster(cl)
# #saveRDS(object = df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes,file = paste0(output_workspace,"df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes.rds"))

#compare proportion of shared epitopes
df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes <- readRDS(file = paste0(output_workspace,"df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes.rds"))
df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes <- subset(df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes,(id_sars_cov_2_epitope%in%df_epitopes$peptide_id)&(id_EhCovs_epitope%in%df_epitopes_ancient_hcovs$peptide_id))
##ggplot(data=df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes) + geom_density(mapping=aes(x=score)) + xlab("PAM120 (gapOpening = 6, gapExtension = 4)")

# png(filename = paste0(output_workspace,"Distribution_SC2_vs_EhCoVs_epitopes_PAM_score.png"),width = 8000,height = 8000,units="px",res=600)
# plot(density(df_scores_pairwise_alignments_sars_cov_2_vs_Ehcovs_epitopes$score),xlab="Epitopes pairwise alignment score SARS-CoV-2 vs EhCoVs\n(PAM120, gapOpening = 6 and gapExtension = 4)",ylab="Density",main="")
# dev.off()

#import conservation scores into one dataframe
df_conservation_scores <- NULL
for (current_orf in c("orf1a","orf1b","S","E","M","N")){
  df_to_add <- read.csv2(file = paste0(output_workspace,"hCoVs_alignments/",current_orf,"_conservation_scores.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
  df_to_add$ORF <- current_orf
  df_to_add$pos_in_alignment <- 1:nrow(df_to_add)
  df_to_add$char_sc2_in_align <- unlist(strsplit(x = seqinr::getSequence(object = toupper(read.fasta(paste0(output_workspace,"hCoVs_alignments/",paste0(current_orf,"_aligned.fasta")),seqtype = "AA",as.string = TRUE,forceDNAtolower = FALSE)$`SARS-CoV-2`),as.string = TRUE)[[1]],split=""))
  df_to_add <- subset(df_to_add, char_sc2_in_align!="-")
  df_to_add$pos_in_protein <- 1:nrow(df_to_add)
  df_to_add <- df_to_add[,c("ORF","pos_in_alignment","pos_in_protein","char_sc2_in_align","pcp_conservation_score")]
  df_conservation_scores <- rbind(df_conservation_scores,df_to_add)
}
df_conservation_scores$pcp_conservation_score <- as.numeric(df_conservation_scores$pcp_conservation_score)
#calculate average conservation scores and number of perfectly conserved sites
df_epitopes$avg_pcp_conservation_score <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) mean(subset(df_conservation_scores,(pos_in_protein>=ceiling((df_epitopes$Genomic_start[i]-v_start_orfs[as.character(df_epitopes$Mapping_region[i])]+1)/3))&(pos_in_protein<=ceiling((df_epitopes$Genomic_End[i]-unname(v_start_orfs[as.character(df_epitopes$Mapping_region[i])])+1)/3))&(ORF==as.character(df_epitopes$Mapping_region[i])))$pcp_conservation_score,na.rm=T),FUN.VALUE = c(0.0)))
df_epitopes$nb_perfectly_conserved_residues <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) sum(subset(df_conservation_scores,(pos_in_protein>=ceiling((df_epitopes$Genomic_start[i]-v_start_orfs[as.character(df_epitopes$Mapping_region[i])]+1)/3))&(pos_in_protein<=ceiling((df_epitopes$Genomic_End[i]-unname(v_start_orfs[as.character(df_epitopes$Mapping_region[i])])+1)/3))&(ORF==as.character(df_epitopes$Mapping_region[i])))$pcp_conservation_score==11.0),FUN.VALUE = c(0)))
#df_epitopes$Mapping_region <- factor(df_epitopes$Mapping_region,levels=names(palette_orfs_epitopes))

#Figures avg_pcp_conservation_score and nb_perfectly_conserved_residues_per_epitope across groups and ORFS
list_comp_groups <- list(c("1", "2"), c("1", "3"), c("2", "3"))

#ggplot(data = df_epitopes,mapping = aes(x = as.factor(Group), y=avg_pcp_conservation_score)) + geom_violin(aes(fill=as.factor(Group)))+ geom_boxplot(width=0.075,fill="white") + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("Average epitope residues P.C.P. conservation score") + xlab("Patient group") +stat_compare_means(comparisons = list_comp_groups,method="wilcox")+ stat_compare_means(method="kruskal") + scale_y_continuous(breaks = seq(0,11,1),limits = c(0,15)) + scale_fill_manual(values = palette_patient_groups)
#ggsave(filename = "Avg_epitope_residues_pcp_conservation_scores_across_groups.png", path=output_workspace, width = 40, height = 20, units = "cm",dpi = 1200)

#ggplot(data = df_epitopes,mapping = aes(x = as.factor(Group), y=nb_perfectly_conserved_residues)) + geom_violin(aes(fill=as.factor(Group))) + geom_boxplot(width=0.075)+ theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("Number of perfectly conserved residues per epitope") + xlab("Patient group") +stat_compare_means(comparisons = list_comp_groups,method="wilcox")+ stat_compare_means(method="kruskal") + scale_y_continuous(breaks = seq(0,15,1),limits = c(0,20))+ scale_fill_manual(values = palette_patient_groups)
#ggsave(filename = "Nb_perfectly_conserved_residues_per_epitope_across_groups.png", path=output_workspace, width = 40, height = 20, units = "cm",dpi = 1200)

#ggplot(data = df_epitopes,mapping = aes(x = Mapping_region, y=avg_pcp_conservation_score))+ geom_violin(aes(fill=Mapping_region))+ geom_boxplot(width=0.075)+ theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("Average epitope residues P.C.P. conservation score") + xlab("ORF") + stat_compare_means() + scale_y_continuous(breaks = seq(0,11,1),limits = c(0,15))+ scale_fill_manual(values = palette_orfs_epitopes)
#ggviolin(data = df_epitopes,x="Mapping_region",y="avg_pcp_conservation_score",fill = "Mapping_region",add = "boxplot",add.params = list(fill = "white",width=0.075))+ theme_bw()+theme(legend.position="none",title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12)) + ylab("Average epitope residues P.C.P. conservation score") + xlab("Patient group") + stat_compare_means() + scale_y_continuous(breaks = seq(0,11,1),limits = c(0,15))+ scale_fill_manual(values = palette_orfs_epitopes)
#ggsave(filename = "Avg_epitope_residues_pcp_conservation_scores_across_ORFs.png", path=output_workspace, width = 40, height = 20, units = "cm",dpi = 1200)

#ggplot(data = df_epitopes,mapping = aes(x = factor(Mapping_region,levels=names(palette_orfs_epitopes)), y=nb_perfectly_conserved_residues))+ geom_violin(aes(fill=factor(Mapping_region,levels=names(palette_orfs_epitopes))))+ geom_boxplot(width=0.075)+ theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("Number of perfectly conserved residues per epitope") + xlab("ORF") + stat_compare_means() + scale_y_continuous(breaks = seq(0,15,1),limits = c(0,20))+ scale_fill_manual(values = palette_orfs_epitopes)
#ggsave(filename = "Nb_perfectly_conserved_residues_per_epitope_across_ORFs.png", path=output_workspace, width = 40, height = 20, units = "cm",dpi = 1200)

#vapply(X = 1:3, FUN = function(x) nrow(subset(df_epitopes,(RFU>100)&(nb_perfectly_conserved_residues>7)&(Group==x))),FUN.VALUE = c(0))

#patients average RFU sars-cov-2 vs average RFU ehcovs
df_avg_immune_responses_sc2_vs_ehcovs <- data.frame(patient=unique(df_epitopes$patient_ID),stringsAsFactors = F)
df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU <- unname(vapply(X=df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) mean(subset(df_epitopes,patient_ID==x)$RFU,na.rm=T),FUN.VALUE = c(0.0)))
df_avg_immune_responses_sc2_vs_ehcovs$avg_ehcovs_RFU <- unname(vapply(X=df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) mean(subset(df_epitopes_ancient_hcovs,patient_ID==x)$RFU,na.rm=T),FUN.VALUE = c(0.0)))
df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes <- unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(avg_pcp_conservation_score>=6))$peptide_id)),FUN.VALUE = c(0)))#unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes_ancient_hcovs,patient_ID==x)$peptide_id)),FUN.VALUE = c(0)))
df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes <- unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x))$peptide_id)),FUN.VALUE = c(0)))
df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes_IgA <- unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(Antibody=="IgA"))$peptide_id)),FUN.VALUE = c(0)))
df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes_IgG <- unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(Antibody=="IgG"))$peptide_id)),FUN.VALUE = c(0)))
df_avg_immune_responses_sc2_vs_ehcovs$Group <- v_group_patient[df_avg_immune_responses_sc2_vs_ehcovs$patient]
df_avg_immune_responses_sc2_vs_ehcovs$pos_neg <- ifelse(test = df_avg_immune_responses_sc2_vs_ehcovs$Group%in%c(1,2),yes="positive",no="negative")
#ggplotregression(fit = lm(formula = y~x,data = data.frame(x=df_avg_immune_responses_sc2_vs_ehcovs$avg_ehcovs_RFU,y=df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU)),ggsave_path = output_workspace,the_filename = paste0("avg_RFU_sc2_vs_avg_RFU_ehCoVs_across_patients.png"),xlabl = paste0("average immune response (RFU) to\nthe EhCoVs epitopes that mapped to SARS-CoV-2 genome"),ylabl = "average immune response (RFU) \nto SARS-CoV-2 epitopes")
#ggplotRegression_export_eps(fit = lm(formula = y~x,data = data.frame(y=df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU,x=df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes)),ggsave_path = output_workspace,the_filename = paste0("avg_RFU_sc2_vs_nb_cross_reactive_epitopes_across_patients.eps"),xlabl = paste0("Number of cross-reactive epitopes"),ylabl = "average immune response to SARS-CoV-2 (RFU)")
#Nb cross-reactive epitopes vs patients' Group
#ggplot(data = df_avg_immune_responses_sc2_vs_ehcovs,mapping = aes(x = as.factor(Group), y=nb_cross_reactive_epitopes))+ geom_violin(aes(fill=as.factor(Group)))+ geom_boxplot(width=0.05)+ theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("Number of cross-reactive epitopes") + xlab("Patient group") + stat_compare_means(method = "wilcox",comparisons = list_comp_groups) + scale_fill_manual(values = palette_patient_groups)
#ggsave(filename = "Nb_cross-reactive_epitopes_vs_patient_Group.png", path=output_workspace, width = 17, height = 15, units = "cm",dpi = 1200)
#Avg RFU vs Nb epitopes
#ggplotRegression_export_eps(fit = lm(formula = y~x,data = data.frame(y=df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU,x=df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes)),ggsave_path = output_workspace,the_filename = paste0("avg_RFU_sc2_vs_nb_epitopes_across_patients.eps"),xlabl = paste0("Number of epitopes"),ylabl = "average immune response to SARS-CoV-2 (RFU)")
#Avg RFU vs Nb epitopes (IgA)
#ggplotRegression_export_eps(fit = lm(formula = y~x,data = data.frame(y=df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU,x=df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes_IgA)),ggsave_path = output_workspace,the_filename = paste0("avg_RFU_sc2_vs_nb_epitopes_IgA_across_patients.eps"),xlabl = paste0("Number of epitopes"),ylabl = "average immune response to SARS-CoV-2 (RFU)")
#Avg RFU vs Nb epitopes (IgG)
#ggplotRegression_export_eps(fit = lm(formula = y~x,data = data.frame(y=df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU,x=df_avg_immune_responses_sc2_vs_ehcovs$nb_epitopes_IgG)),ggsave_path = output_workspace,the_filename = paste0("avg_RFU_sc2_vs_nb_epitopes_IgG_across_patients.eps"),xlabl = paste0("Number of epitopes"),ylabl = "average immune response to SARS-CoV-2 (RFU)")
#Avg immune response vs Nb epitopes
#ggplot(df_avg_immune_responses_sc2_vs_ehcovs) +
#  geom_smooth(aes(y = log10(avg_sc2_RFU),x = nb_epitopes,col="Overall"), method=lm, se=FALSE) + geom_smooth(aes(y = log10(avg_sc2_RFU),x = nb_epitopes_IgA,col="IgA"), method=lm, se=FALSE) + geom_point(aes(y = log10(avg_sc2_RFU),x = nb_epitopes_IgA,col="IgA")) + geom_smooth(aes(y = log10(avg_sc2_RFU),x = nb_epitopes_IgG,col="IgG"), method=lm, se=FALSE) + geom_point(aes(y = log10(avg_sc2_RFU),x = nb_epitopes_IgG,col="IgG")) + scale_color_manual(values = c("Overall"="black","IgA"="green3","IgG"="tan2")) +
#  facet_wrap(~pos_neg) +
# labs(x = "Number of epitopes", y = "log10(Average RFU)",col="Antibody") + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position = "right")
#ggsave(filename = "Avg_RFU_vs_Nb_epitopes_by_Antibody_and_Covid19_Status.eps", path=output_workspace, width = 20/2.54, height = 10/2.54,dpi = 1200,device = cairo_ps)
#ggsave(filename = "Avg_RFU_vs_Nb_epitopes_by_Antibody_and_Covid19_Status.png", path=output_workspace, width = 20, height = 10, units = "cm",dpi = 1200)


#Proportion of sc2 perfectly conserved sites that are epitopes in one the 4 hCoVs
get_pos_in_prot_ehcovs_epitopes <- function(the_orf){
  df_current_epitopes <- subset(df_epitopes_ancient_hcovs,Mapping_region==the_orf)
  v_out <- NULL
  if(nrow(df_current_epitopes)==0){
    return(NA)
  }
  for (j in 1:nrow(df_current_epitopes)){
    v_out <- c(v_out,(ceiling((df_current_epitopes$Genomic_start[j]-v_start_orfs[as.character(df_current_epitopes$Mapping_region)[j]]+1)/3)):(ceiling((df_current_epitopes$Genomic_End[j]-v_start_orfs[as.character(df_current_epitopes$Mapping_region)[j]]+1)/3)))
  }
  return(sort(unique(v_out)))
}
lst_orfs_pos_in_prot_ehcovs_epitopes <- NULL
for (current_orf in names(palette_orfs_epitopes)){
  lst_orfs_pos_in_prot_ehcovs_epitopes <- c(lst_orfs_pos_in_prot_ehcovs_epitopes,list(c(get_pos_in_prot_ehcovs_epitopes(the_orf = current_orf))))
}

names(lst_orfs_pos_in_prot_ehcovs_epitopes) <- names(palette_orfs_epitopes)
df_sc2_perfectly_conserved_sites <- subset(df_conservation_scores,pcp_conservation_score==11)
proportion_sc2_perfectly_conserved_sites_associated_to_cross_reactivity <- sum(unname(vapply(X = 1:nrow(df_sc2_perfectly_conserved_sites),FUN = function(i) ifelse(test = df_sc2_perfectly_conserved_sites$ORF[i] %in% df_epitopes_ancient_hcovs$Mapping_region ,yes=df_sc2_perfectly_conserved_sites$pos_in_protein[i]%in%(lst_orfs_pos_in_prot_ehcovs_epitopes[[df_sc2_perfectly_conserved_sites$ORF[i]]]),no = F),FUN.VALUE = c(T))))/nrow(df_sc2_perfectly_conserved_sites)

#Epitope cross-reactivity profile (Heatmap)
#explore the use of df_epitope_overlap_ancient_hcovs
mtx_pres_abs_cross_reactivity_immunity_related_epitopes <- matrix(0,ncol=length(c("HKU1","NL63","OC43","229E")),nrow=nrow(df_epitope_overlap_ancient_hcovs))
rownames(mtx_pres_abs_cross_reactivity_immunity_related_epitopes) <- df_epitope_overlap_ancient_hcovs$peptide_id
colnames(mtx_pres_abs_cross_reactivity_immunity_related_epitopes) <- c("HKU1","NL63","OC43","229E")
for (i in 1:nrow(mtx_pres_abs_cross_reactivity_immunity_related_epitopes)){
  mtx_pres_abs_cross_reactivity_immunity_related_epitopes[i,] <- as.integer(colnames(mtx_pres_abs_cross_reactivity_immunity_related_epitopes) %in% strsplit(x = df_epitope_overlap_ancient_hcovs$hcovs_overlap[i],"/")[[1]])
}
#png(filename = paste0(output_workspace,"Heatmap_mtx_pres_abs_sc2_cross_reactive_epitopes_across_EhCoVs.png"),width = 6000,height = 8000,units="px",res=600)
#heatmap.2(x = mtx_pres_abs_cross_reactivity_immunity_related_epitopes, distfun = function(x) dist(x, method = "binary"), hclustfun = function(d) hclust(d,method = "ward.D2"),main = "Presence (black) / absence (white)",trace="none",scale="none", labRow = NA,cexRow = 1,cexCol = 1,xlab = "Endemic human coronavirus",ylab="SARS-CoV-2 epitope",col = colorRampPalette(c("white","black"), space = "rgb")(100),key = FALSE)
#dev.off()

#Explore different cross-reactivity cut-offs
df_cross_reactivity_cutoffs <- data.frame(cutoff_value=c(seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25),seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,1)),cutoff_type=c(rep("Cut-off on the average\np.c.p conservation score",length(seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25))),rep("Cut-off on the\nbest match PAM score",length(seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,1)))))
df_cross_reactivity_cutoffs$p_value <- unname(vapply(X = 1:nrow(df_cross_reactivity_cutoffs),FUN = function(i) ifelse(test=df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score",yes=wilcox.test(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score<df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU)$p.value,no = wilcox.test(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope<df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU)$p.value),FUN.VALUE = c(0.0)))
df_cross_reactivity_cutoffs$normalized_bin_size_difference <- unname(vapply(X = 1:nrow(df_cross_reactivity_cutoffs),FUN = function(i) ifelse(test=df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score",yes=abs(nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score<df_cross_reactivity_cutoffs$cutoff_value[i]))- nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i]))),no = abs(nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope<df_cross_reactivity_cutoffs$cutoff_value[i]))-nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i]))) ),FUN.VALUE = c(0)))
df_cross_reactivity_cutoffs$normalized_bin_size_difference <- 5*scale(df_cross_reactivity_cutoffs$normalized_bin_size_difference,center = F,scale=T)
df_cross_reactivity_cutoffs$normalized_bin_size_difference <- unname(vapply(X = 1:nrow(df_cross_reactivity_cutoffs),FUN = function(i) ifelse(test=df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score",yes=abs(nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score<df_cross_reactivity_cutoffs$cutoff_value[i]))- nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i]))),no = abs(nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope<df_cross_reactivity_cutoffs$cutoff_value[i]))-nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i]))) ),FUN.VALUE = c(0)))
df_cross_reactivity_cutoffs$sign_bin_size_difference <- unname(vapply(X = 1:nrow(df_cross_reactivity_cutoffs),FUN = function(i) ifelse(test=df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score",yes=sign(nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i]))- nrow(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score<df_cross_reactivity_cutoffs$cutoff_value[i]))),no = sign(nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i]))-nrow(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope<df_cross_reactivity_cutoffs$cutoff_value[i]))) ),FUN.VALUE = c(0)))
df_cross_reactivity_cutoffs$sign_mean_difference <- unname(vapply(X = 1:nrow(df_cross_reactivity_cutoffs),FUN = function(i) ifelse(test=df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score",yes=sign(mean(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,na.rm=T)- mean(subset(df_sars_cov_2_epitopes,avg_pcp_conservation_score<df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,na.rm=T)),no = sign(mean(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,na.rm=T)-mean(subset(df_sars_cov_2_epitopes,PAM_score_best_match_with_a_Ehcovs_epitope<df_cross_reactivity_cutoffs$cutoff_value[i])$avg_RFU,na.rm=T)) ),FUN.VALUE = c(0)))
#assess the effect of the cut-offs on avg_RFU~nb_cross_reactive_epitopes
df_cross_reactivity_cutoffs$p_value_avg_RFU_vs_nb_cross_reactive_epitopes <- NA
df_cross_reactivity_cutoffs$sign_slope_avg_RFU_vs_nb_cross_reactive_epitopes <- NA
for (i in 1:nrow(df_cross_reactivity_cutoffs)){
  if (df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score"){
    cp_df_avg_immune_responses_sc2_vs_ehcovs <- df_avg_immune_responses_sc2_vs_ehcovs
    cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes <- unname(vapply(X = cp_df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i]))$peptide_id)),FUN.VALUE = c(0)))
    coeficients_summary_fit <- coefficients(summary(lm(cp_df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU~cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes)))
    df_cross_reactivity_cutoffs$p_value_avg_RFU_vs_nb_cross_reactive_epitopes[i] <- coeficients_summary_fit[2,4]
    df_cross_reactivity_cutoffs$sign_slope_avg_RFU_vs_nb_cross_reactive_epitopes[i] <- sign(coeficients_summary_fit[2,1])
  }else{
    cp_df_avg_immune_responses_sc2_vs_ehcovs <- df_avg_immune_responses_sc2_vs_ehcovs
    cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes <- unname(vapply(X = cp_df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) sum((subset(df_sars_cov_2_epitopes,peptide_id%in%(subset(df_epitopes,(patient_ID==x))$peptide_id))$PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i])),FUN.VALUE = c(0)))
    coeficients_summary_fit <- coefficients(summary(lm(cp_df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU~cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes)))
    df_cross_reactivity_cutoffs$p_value_avg_RFU_vs_nb_cross_reactive_epitopes[i] <- coeficients_summary_fit[2,4]
    df_cross_reactivity_cutoffs$sign_slope_avg_RFU_vs_nb_cross_reactive_epitopes[i] <- sign(coeficients_summary_fit[2,1])
  }
}
remove(cp_df_avg_immune_responses_sc2_vs_ehcovs)
#assess the effect of the cut-offs on nb_cross_reactive_epitopes_per_patient~Group
df_cross_reactivity_cutoffs$p_value_nb_cross_reactive_epitopes_per_patient_vs_Group <- NA
for (i in 1:nrow(df_cross_reactivity_cutoffs)){
  if (df_cross_reactivity_cutoffs$cutoff_type[i]=="Cut-off on the average\np.c.p conservation score"){
    cp_df_avg_immune_responses_sc2_vs_ehcovs <- df_avg_immune_responses_sc2_vs_ehcovs
    cp_df_avg_immune_responses_sc2_vs_ehcovs$Group <- as.factor(cp_df_avg_immune_responses_sc2_vs_ehcovs$Group)
    cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes <- unname(vapply(X = cp_df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(avg_pcp_conservation_score>=df_cross_reactivity_cutoffs$cutoff_value[i]))$peptide_id)),FUN.VALUE = c(0)))
    current_fit <-lm(cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes~cp_df_avg_immune_responses_sc2_vs_ehcovs$Group)
    coeficients_summary_fit <- coefficients(summary(current_fit))
    df_cross_reactivity_cutoffs$p_value_nb_cross_reactive_epitopes_per_patient_vs_Group[i] <- broom::glance(current_fit)$p.value
  }else{
    cp_df_avg_immune_responses_sc2_vs_ehcovs <- df_avg_immune_responses_sc2_vs_ehcovs
    cp_df_avg_immune_responses_sc2_vs_ehcovs$Group <- as.factor(cp_df_avg_immune_responses_sc2_vs_ehcovs$Group)
    cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes <- unname(vapply(X = cp_df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) sum((subset(df_sars_cov_2_epitopes,peptide_id%in%(subset(df_epitopes,(patient_ID==x))$peptide_id))$PAM_score_best_match_with_a_Ehcovs_epitope>=df_cross_reactivity_cutoffs$cutoff_value[i])),FUN.VALUE = c(0)))
    current_fit <-lm(cp_df_avg_immune_responses_sc2_vs_ehcovs$nb_cross_reactive_epitopes~cp_df_avg_immune_responses_sc2_vs_ehcovs$Group)
    coeficients_summary_fit <- coefficients(summary(current_fit))
    df_cross_reactivity_cutoffs$p_value_nb_cross_reactive_epitopes_per_patient_vs_Group[i] <- broom::glance(current_fit)$p.value
  }
}
remove(cp_df_avg_immune_responses_sc2_vs_ehcovs)
#Plots illustrating the cross-reactivity Cut-offs' effects
#RFU cross-reactive vs RFU SC2-specific
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the average\np.c.p conservation score"),mapping = aes(x = cutoff_value, y=-log10(p_value)))+ geom_point(mapping=aes(size=normalized_bin_size_difference,col=as.character(sign_mean_difference))) + scale_color_manual(values = c("-1"="red","1"="blue")) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\nRFU cross-reactive vs SARS-CoV-2-specific epitopes") + xlab("Cut-off value for defining cross-reactive epitopes\n(Average p.c.p. conservation score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25),limits = range(seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25))) +geom_vline(xintercept = 6)
#ggsave(filename = "pvalue_avg_RFU_cross-reactive_epitopes_vs_SC2_specific_eptopes_at_different_avg_pcp_conservation_score_cut-offs.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the\nbest match PAM score"),mapping = aes(x = cutoff_value, y=-log10(p_value))) + geom_point(mapping=aes(size=normalized_bin_size_difference,col=as.character(sign_mean_difference))) + scale_color_manual(values = c("-1"="red","1"="blue")) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\nRFU cross-reactive vs SARS-CoV-2-specific epitopes") + xlab("Cut-off value for defining cross-reactive epitopes\n(Best match PAM score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2),limits = range(seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2)))
#ggsave(filename = "pvalue_avg_RFU_cross-reactive_epitopes_vs_SC2_specific_eptopes_at_different_cut-offs_on_best_match_PAM_score.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)
#avg_RFU vs Nb_cross_reactive_epitopes
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the average\np.c.p conservation score"),mapping = aes(x = cutoff_value, y=-log10(p_value_avg_RFU_vs_nb_cross_reactive_epitopes)))+ geom_point(mapping=aes(size=normalized_bin_size_difference,col=as.character(sign_slope_avg_RFU_vs_nb_cross_reactive_epitopes))) + scale_color_manual(values = c("-1"="red","1"="blue")) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\navreage RFU ~ Number of cross-reactive epitopes") + xlab("Cut-off value for defining cross-reactive epitopes\n(Average p.c.p. conservation score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25),limits = range(seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25))) +geom_vline(xintercept = 6)
#ggsave(filename = "pvalue_avg_RFU_vs_nb_cross_reactive_epitopes_at_different_avg_pcp_conservation_score_cut-offs.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the\nbest match PAM score"),mapping = aes(x = cutoff_value, y=-log10(p_value_avg_RFU_vs_nb_cross_reactive_epitopes))) + geom_point(mapping=aes(size=normalized_bin_size_difference,col=as.character(sign_slope_avg_RFU_vs_nb_cross_reactive_epitopes))) + scale_color_manual(values = c("-1"="red","1"="blue")) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\navreage RFU ~ Number of cross-reactive epitopes") + xlab("Cut-off value for defining cross-reactive epitopes\n(Best match PAM score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2),limits = range(seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2)))
#ggsave(filename = "pvalue_avg_RFU_vs_nb_cross_reactive_epitopes_at_different_cut-offs_on_best_match_PAM_score.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)
#nb_cross_reactive_epitopes_per_patient~Group
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the average\np.c.p conservation score"),mapping = aes(x = cutoff_value, y=-log10(p_value_nb_cross_reactive_epitopes_per_patient_vs_Group)))+ geom_point(mapping=aes(size=normalized_bin_size_difference)) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\nNumber of cross-reactive epitopes per patient vs Group") + xlab("Cut-off value for defining cross-reactive epitopes\n(Average p.c.p. conservation score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25),limits = range(seq(min(floor(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)+1,max(ceiling(df_sars_cov_2_epitopes$avg_pcp_conservation_score),na.rm=T)-1,0.25))) +geom_vline(xintercept = 6)
#ggsave(filename = "pvalue_nb_cross_reactive_epitopes_per_patient_vs_Group_at_different_avg_pcp_conservation_score_cut-offs.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)
#ggplot(data = subset(df_cross_reactivity_cutoffs,cutoff_type=="Cut-off on the\nbest match PAM score"),mapping = aes(x = cutoff_value, y=-log10(p_value_nb_cross_reactive_epitopes_per_patient_vs_Group))) + geom_point(mapping=aes(size=normalized_bin_size_difference)) + theme_bw()+ theme(title =  element_text(size=12),legend.text = element_text(size = 8),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),legend.position="none") + ylab("-log10(p-value)\nNumber of cross-reactive epitopes per patient vs Group") + xlab("Cut-off value for defining cross-reactive epitopes\n(Best match PAM score)") + geom_hline(yintercept = -log10(0.05),lty=2,col="red") +scale_x_continuous(breaks=seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2),limits = range(seq(min(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)+1,max(df_sars_cov_2_epitopes$PAM_score_best_match_with_a_Ehcovs_epitope,na.rm=T)-1,2)))
#ggsave(filename = "pvalue_nb_cross_reactive_epitopes_per_patient_vs_Group_at_different_cut-offs_on_best_match_PAM_score.eps", path=output_workspace, width = 17, height = 15, units = "cm",device = cairo_ps,dpi = 1200)

#Evolution rates Cross-reactive vs SC2-specific epitopes
df_high_confidence_epitope_metrics <- readRDS(file = paste0(output_workspace,"df_high_confidence_epitope_metrics.rds"))
df_epitopes$is_cross_reactive <- NA
df_epitopes$is_cross_reactive <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_epitopes$Mapping_region[i])&(pos_in_protein>=df_epitopes$protein_start[i])&(pos_in_protein<=df_epitopes$protein_end[i])&(is_cross_reactive)))>=5,FUN.VALUE = c(F)))
df_epitopes$is_cross_reactive <- ifelse(test=df_epitopes$is_cross_reactive,yes="Cross-reactive",no="SC2-specific") #ifelse(test=df_epitopes$avg_pcp_conservation_score>=6,yes="Cross-reactive",no="SC2-specific")
v_position_Cross_reactive_epitopes <- NULL
for (i in 1:nrow(subset(df_epitopes,(!is.na(is_cross_reactive))&(is_cross_reactive == "Cross-reactive")))){
  v_position_Cross_reactive_epitopes <- c(v_position_Cross_reactive_epitopes,((subset(df_epitopes,(!is.na(is_cross_reactive))&(is_cross_reactive == "Cross-reactive"))$Genomic_start[i]):(subset(df_epitopes,(!is.na(is_cross_reactive))&(is_cross_reactive == "Cross-reactive"))$Genomic_End[i])))
}
v_position_Cross_reactive_epitopes<- sort(unique(v_position_Cross_reactive_epitopes))
v_position_SC2_specific_epitopes <- setdiff(v_unique_epitope_positions,v_position_Cross_reactive_epitopes)
v_length_cross_reactive_vs_sc2_specific_epitopes_positions <- c("Cross-reactive"=length(v_position_Cross_reactive_epitopes),"SC2-specific"=length(v_position_SC2_specific_epitopes))
v_coverage_cross_reactive_vs_sc2_specific_epitopes_positions <- c("Cross-reactive"=mean(subset(df_variants_NCBI_SRA_amplicon,Position%in%v_position_Cross_reactive_epitopes)$total_depth,na.rm=T),"SC2-specific"=mean(subset(df_variants_NCBI_SRA_amplicon,Position%in%v_position_SC2_specific_epitopes)$total_depth,na.rm=T))

df_variants_NCBI_SRA_amplicon$is_cross_reactive_position <- df_variants_NCBI_SRA_amplicon$Position%in%v_position_Cross_reactive_epitopes
df_variants_site_enough_covered_NCBI_SRA_amplicon <- subset(df_variants_NCBI_SRA_amplicon, total_depth>=min_cov)

#compare mutation and substitution rate cross_reactive_vs_sc2_specific epitopes
df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon <- data.frame(Sample=rep(lst_samples_NCBI_SRA_amplicon,length(c(TRUE,FALSE))),is_cross_reactive_position=rep(c(TRUE,FALSE),each=length(lst_samples_NCBI_SRA_amplicon)),nb_mutations=0,nb_fixed_mutations=0,mut_rate=0,subst_rate=0,stringsAsFactors = F)
nb_cores <- nb_cpus
lst_splits <- split(1:nrow(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon), ceiling(seq_along(1:nrow(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon))/(nrow(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon)/nb_cores)))
the_f_parallel_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  df_metrics_current_subset <- df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon[the_vec,]
  count_iter <- 0
  for (the_i in 1:nrow(df_metrics_current_subset)){
    df_depth_NCBI_SRA_amplicon_current_sample <- read.csv2(file = paste0(depth_data_wp,"depth_report_NCBI_SRA_amplicon/df_depth_NCBI_SRA_amplicon_",df_metrics_current_subset$Sample[the_i],".csv"),sep = ",",header = F,stringsAsFactors = FALSE)
    colnames(df_depth_NCBI_SRA_amplicon_current_sample) <- c("sample","position","depth")
    df_depth_NCBI_SRA_amplicon_current_sample$ORF <- vapply(X = df_depth_NCBI_SRA_amplicon_current_sample$position,FUN = find_ORF_of_mutation,FUN.VALUE = c(""))
    df_depth_NCBI_SRA_amplicon_current_sample <- unique(df_depth_NCBI_SRA_amplicon_current_sample)
    v_currentsample_positions_enough_covered <- subset(df_depth_NCBI_SRA_amplicon_current_sample,(depth>=min_cov))$position
    if (df_metrics_current_subset$is_cross_reactive_position[the_i]){
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_position_Cross_reactive_epitopes,v_currentsample_positions_enough_covered)
    }else{
      v_the_sites_with_enough_cov_for_current_category <- intersect(v_position_SC2_specific_epitopes,v_currentsample_positions_enough_covered)
    }
    df_metrics_current_subset$nb_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(is_fixed=="No")&(!Position%in%v_positions_inter)&(is_cross_reactive_position==df_metrics_current_subset$is_cross_reactive_position[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_fixed_mutations[the_i] <- nrow(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(is_fixed=="Yes")&(is_prevalence_above_transmission_threshold)&(is_cross_reactive_position==df_metrics_current_subset$is_cross_reactive_position[the_i])&(Sample==df_metrics_current_subset$Sample[the_i])))
    df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] <- length(v_the_sites_with_enough_cov_for_current_category)
    df_metrics_current_subset$mut_rate[the_i] <- df_metrics_current_subset$nb_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$subst_rate[the_i] <- df_metrics_current_subset$nb_fixed_mutations[the_i]/df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i]
    df_metrics_current_subset$Nb_ss[the_i] <- sum(unname(vapply(X = v_the_sites_with_enough_cov_for_current_category,FUN = calculate_nb_ss_position_in_genome,FUN.VALUE = c(0))),na.rm=T)
    df_metrics_current_subset$Nb_nss[the_i] <- df_metrics_current_subset$nb_sites_with_enough_cov_for_this_category[the_i] - df_metrics_current_subset$Nb_ss[the_i]
    df_metrics_current_subset$within_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$within_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(!Position%in%v_positions_inter)&(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq<0.75)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_nsm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(!(is_synonymous))&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_nonsynonymous_mutations_NCBI_SRA_amplicon))
    df_metrics_current_subset$between_host_Nb_sm[the_i] <- length(intersect(subset(df_variants_site_enough_covered_NCBI_SRA_amplicon,(Position %in% v_the_sites_with_enough_cov_for_current_category)&(!is.na(is_synonymous))&(is_synonymous)&(VarFreq>=0.75)&(is_prevalence_above_transmission_threshold)&(Sample==df_metrics_current_subset$Sample[the_i]))$mutation_name,v_lst_synonymous_mutations_NCBI_SRA_amplicon))
    
    if (the_i%%100==0){
      print(paste0("[Cross-reactive vs SC2-secific epitopes (NCBI)] Core ",i_cl,": Step ",the_i," done out of ",nrow(df_metrics_current_subset),"!"))
    }
  }
  df_metrics_current_subset$pN <- df_metrics_current_subset$within_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$pS <- df_metrics_current_subset$within_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$pN_pS <- ifelse(test=df_metrics_current_subset$pS==0,yes=NA,no=df_metrics_current_subset$pN/df_metrics_current_subset$pS)
  df_metrics_current_subset$dN <- df_metrics_current_subset$between_host_Nb_nsm/df_metrics_current_subset$Nb_nss
  df_metrics_current_subset$dS <- df_metrics_current_subset$between_host_Nb_sm/df_metrics_current_subset$Nb_ss
  df_metrics_current_subset$dN_dS <- ifelse(test=df_metrics_current_subset$dS==0,yes=NA,no=df_metrics_current_subset$dN/df_metrics_current_subset$dS)
  df_metrics_current_subset$alpha_MK_Test <- ifelse(test=df_metrics_current_subset$dN_dS==0,yes=NA,no=(1-((df_metrics_current_subset$pN_pS)/(df_metrics_current_subset$dN_dS))))
  
  return(df_metrics_current_subset)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_Evo_rates_analyses.txt"))
registerDoParallel(cl)
df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel","Biostrings"))  %dopar% the_f_parallel_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon(i_cl)
stopCluster(cl)
saveRDS(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,paste0(output_workspace,"df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.rds"))

ggplot(data = df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_cross_reactive_position),levels=c("TRUE","FALSE")),y = mut_rate,fill=as.character(is_cross_reactive_position))) + geom_violin() + geom_point() + xlab("Cross-reactive epitope sites?") + ylab("Within-host mutation rate (Count / Length)") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + scale_y_continuous(breaks=seq(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$mut_rate)+1e-2,1e-2),limits = c(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$mut_rate)+1e-2)) + stat_compare_means(method = "wilcox")
ggsave(filename = "Mutation_rate_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
ggplot(data = df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_cross_reactive_position),levels=c("TRUE","FALSE")),y = subst_rate,fill=as.character(is_cross_reactive_position))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("Cross-reactive epitope sites?") + ylab("Substitution rate (Count / Length)") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_y_continuous(breaks=seq(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$subst_rate)+1e-4,1e-4),limits = c(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$subst_rate)+1e-4)) + stat_compare_means(method = "wilcox")
ggsave(filename = "Substitution_rate_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)

ggplot(data = df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_cross_reactive_position),levels=c("TRUE","FALSE")),y = pN_pS,fill=as.character(is_cross_reactive_position))) + geom_violin() + geom_point() + xlab("Cross-reactive epitope sites?") + ylab("pN/pS") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$pN_pS)+1e-2,1e-2),limits = c(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$pN_pS)+1e-2))
ggsave(filename = "pN_pS_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
ggplot(data = df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_cross_reactive_position),levels=c("TRUE","FALSE")),y = dN_dS,fill=as.character(is_cross_reactive_position))) + geom_violin() + geom_boxplot(width=0.075,fill="white") + geom_point() + xlab("Cross-reactive epitope sites?") + ylab("dN/dS") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$dN_dS)+1e-4,1e-4),limits = c(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$dN_dS)+1e-4))
ggsave(filename = "dN_dS_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)
ggplot(data = df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon,aes(x=factor(as.character(is_cross_reactive_position),levels=c("TRUE","FALSE")),y = alpha_MK_Test,fill=as.character(is_cross_reactive_position))) + geom_violin() + geom_point() + xlab("Cross-reactive epitope sites?") + ylab("McDonald-Kreitman test \U003B1") + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.position = "none") + scale_fill_manual(values = c("TRUE"="tan2","FALSE"="grey60")) + stat_compare_means(method = "wilcox") #+ scale_y_continuous(breaks=seq(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$alpha_MK_Test)+1e-2,1e-2),limits = c(0,max(df_metrics_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon$alpha_MK_Test)+1e-2))
ggsave(filename = "MK_test_alpha_cross_reactive_vs_sc2_specific_epitopes_NCBI_SRA_amplicon.png", path=output_workspace, width = 15, height = 15, units = "cm",dpi = 1200)

#Proportion of epitope sites with mutation
proportion_epitope_sites_with_mutation_NCBI_SRA_amplicon <- length(unique(subset(df_variants_NCBI_SRA_amplicon,is_epitope_related)$Position))/length(v_unique_epitope_positions)

#add SNVs position in ORF protein
df_variants_NCBI_SRA_amplicon$pos_in_ORF_protein_seq <- vapply(X = 1:nrow(df_variants_NCBI_SRA_amplicon),FUN = function(i) ceiling((df_variants_NCBI_SRA_amplicon$Position[i] - v_start_orfs[df_variants_NCBI_SRA_amplicon$ORF[i]] + 1)/3),FUN.VALUE = c(0))

#EhCoVs prevalence
df_epitopes_HKU1$protein_start <- unname(ceiling((df_epitopes_HKU1$Genomic_start - v_start_orfs[as.character(df_epitopes_HKU1$Mapping_region)] + 1)/3))
df_epitopes_HKU1$protein_end <- unname(ceiling((df_epitopes_HKU1$Genomic_End - v_start_orfs[as.character(df_epitopes_HKU1$Mapping_region)] + 1)/3))
df_epitopes_NL63$protein_start <- unname(ceiling((df_epitopes_NL63$Genomic_start - v_start_orfs[as.character(df_epitopes_NL63$Mapping_region)] + 1)/3))
df_epitopes_NL63$protein_end <- unname(ceiling((df_epitopes_NL63$Genomic_End - v_start_orfs[as.character(df_epitopes_NL63$Mapping_region)] + 1)/3))
df_epitopes_OC43$protein_start <- unname(ceiling((df_epitopes_OC43$Genomic_start - v_start_orfs[as.character(df_epitopes_OC43$Mapping_region)] + 1)/3))
df_epitopes_OC43$protein_end <- unname(ceiling((df_epitopes_OC43$Genomic_End - v_start_orfs[as.character(df_epitopes_OC43$Mapping_region)] + 1)/3))
df_epitopes_229E$protein_start <- unname(ceiling((df_epitopes_229E$Genomic_start - v_start_orfs[as.character(df_epitopes_229E$Mapping_region)] + 1)/3))
df_epitopes_229E$protein_end <- unname(ceiling((df_epitopes_229E$Genomic_End - v_start_orfs[as.character(df_epitopes_229E$Mapping_region)] + 1)/3))

#Positive vs negative patients RFU boxplots
df_epitopes$pos_neg <- v_patients_to_binary_group[df_epitopes$patient_ID]
#cross-reactive epitopes definition
df_epitopes$is_cross_reactive <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_epitopes$Mapping_region[i])&(pos_in_protein>=df_epitopes$protein_start[i])&(pos_in_protein<=df_epitopes$protein_end[i])&(is_cross_reactive)))>=5,FUN.VALUE = c(F)))
df_proteins_nb_cr_and_SC2_specific_epitopes_current_study <- data.frame(ORF=c("orf1a","orf1b","S","E","M","N"), Number_of_cross_reactive_epitopes = NA, Number_of_SC2_specific_epitopes = NA)
df_proteins_nb_cr_and_SC2_specific_epitopes_current_study$Number_of_cross_reactive_epitopes <- unname(vapply(X = c("orf1a","orf1b","S","E","M","N"),FUN=function(x) nrow(subset(df_epitopes,(Mapping_region==x)&(is_cross_reactive))),FUN.VALUE = c(0)))
df_proteins_nb_cr_and_SC2_specific_epitopes_current_study$Number_of_SC2_specific_epitopes <- unname(vapply(X = c("orf1a","orf1b","S","E","M","N"),FUN=function(x) nrow(subset(df_epitopes,(Mapping_region==x)&(!is_cross_reactive))),FUN.VALUE = c(0)))
df_epitopes$is_cross_reactive <- ifelse(test=df_epitopes$is_cross_reactive,yes="Cross-reactive",no="SC2-specific") #ifelse(test=df_epitopes$avg_pcp_conservation_score>=6,yes="Cross-reactive",no="SC2-specific")

#effect of threshold nb cross-reactive epitope sites on the number of cross-reactive epitopes (sensitivity analysis)
df_nb_cr_epitopes <- NULL
df_slope_correlation_avg_immune_resp_vs_nb_cr_epitopes_current_study <- NULL
for (n in 1:15){
  df_epitopes$is_cross_reactive <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_epitopes$Mapping_region[i])&(pos_in_protein>=df_epitopes$protein_start[i])&(pos_in_protein<=df_epitopes$protein_end[i])&(is_cross_reactive)))>=n,FUN.VALUE = c(F)))
  current_fit_avg_immune_resp_vs_nb_epitopes <- lm(formula = scale(df_avg_immune_responses_sc2_vs_ehcovs$avg_sc2_RFU,T,T)~scale(unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(is_cross_reactive))$peptide_id)),FUN.VALUE = c(0))),T,T))
  df_nb_cr_epitopes <- rbind(df_nb_cr_epitopes,data.frame(dataset="Current study",threshold_nb_cr_sites=n,Sample=df_avg_immune_responses_sc2_vs_ehcovs$patient,nb_cr_epitopes=unname(vapply(X = df_avg_immune_responses_sc2_vs_ehcovs$patient,FUN = function(x) length(unique(subset(df_epitopes,(patient_ID==x)&(is_cross_reactive))$peptide_id)),FUN.VALUE = c(0))),stringsAsFactors = F))
  df_slope_correlation_avg_immune_resp_vs_nb_cr_epitopes_current_study <- rbind(df_slope_correlation_avg_immune_resp_vs_nb_cr_epitopes_current_study,data.frame(dataset="Current study",threshold_nb_cr_sites=n,slope=ifelse(test=broom::glance(current_fit_avg_immune_resp_vs_nb_epitopes)$p.value<0.05,yes = unname(coefficients(current_fit_avg_immune_resp_vs_nb_epitopes)[2]),no=NA),stringsAsFactors = F))
}

df_epitopes$is_cross_reactive <- unname(vapply(X = 1:nrow(df_epitopes),FUN = function(i) nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_epitopes$Mapping_region[i])&(pos_in_protein>=df_epitopes$protein_start[i])&(pos_in_protein<=df_epitopes$protein_end[i])&(is_cross_reactive)))>=5,FUN.VALUE = c(F)))
df_epitopes$is_cross_reactive <- ifelse(test=df_epitopes$is_cross_reactive,yes="Cross-reactive",no="SC2-specific") #ifelse(test=df_epitopes$avg_pcp_conservation_score>=6,yes="Cross-reactive",no="SC2-specific")


# ##########################Mutations of interest####################
#mutation emergence (which wave)
df_variants_NCBI_SRA_amplicon$wave <- ifelse(is.na(df_variants_NCBI_SRA_amplicon$collection_date),yes=NA,no=ifelse(test=df_variants_NCBI_SRA_amplicon$collection_date < "2020-07",yes=1,no=ifelse(test=df_variants_NCBI_SRA_amplicon$collection_date < "2021-03",yes=2,no=3)))
df_variants_NCBI_SRA_amplicon$short_label_mut <- paste0(df_variants_NCBI_SRA_amplicon$gene,":",df_variants_NCBI_SRA_amplicon$old_aa,df_variants_NCBI_SRA_amplicon$pos_in_protein,df_variants_NCBI_SRA_amplicon$new_aa)
v_lst_short_label_mut_NCBI_SRA_amplicon <- sort(unique(df_variants_NCBI_SRA_amplicon$short_label_mut))
v_wave_of_emergence_of_mutations_NCBI_SRA_amplicon <- vapply(X = v_lst_short_label_mut_NCBI_SRA_amplicon, FUN = function(the_mut) ifelse(all(is.na(subset(df_variants_NCBI_SRA_amplicon,short_label_mut==the_mut)$collection_date)),yes=NA,no=min(subset(df_variants_NCBI_SRA_amplicon,short_label_mut==the_mut)$wave,na.rm=T)),FUN.VALUE = c(0))
v_prevalence_mutations <- vapply(X = v_lst_short_label_mut_NCBI_SRA_amplicon, FUN = function(the_mut) length(unique(subset(df_variants_NCBI_SRA_amplicon,short_label_mut==the_mut)$Sample)),FUN.VALUE = c(0))

#S protein mutations of concern or under investigation
v_S_region_mutations_of_interest_in_epitope_sites <- intersect(paste0("S:",v_S_region_mutations_of_interest),subset(df_variants_NCBI_SRA_amplicon,is_epitope_related)$short_label_mut)
df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon <- data.frame(mutation=unique(v_S_region_mutations_of_interest_in_epitope_sites),prevalence_mutation=v_prevalence_mutations[v_S_region_mutations_of_interest_in_epitope_sites],Wave_first_time_observed=v_wave_of_emergence_of_mutations_NCBI_SRA_amplicon[v_S_region_mutations_of_interest_in_epitope_sites],avg_RFU=NA,avg_RFU_epitope_site_in_positive_patients_current_study=NA,avg_RFU_epitope_site_in_negative_patients_current_study=NA,is_detected_as_cross_reactive_in_current_study=NA,prevalence_epitope_site_in_positive_patients_current_study=NA,prevalence_epitope_site_in_negative_patients_current_study=NA,month_first_time_detected=NA,lst_candidate_country_first_apparition=NA,stringsAsFactors = F)
for (i in 1:nrow(df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon)){
  pos_in_S_current_mut <- as.integer(substr(df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$mutation[i],gregexpr(pattern = ":",text = df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$mutation[i],fixed=T)[[1]][1]+2,nchar(df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$mutation[i])-1))
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$avg_RFU[i] <- mean(subset(df_epitopes,(!is.na(Mapping_region))&(Mapping_region=="S")&(pos_in_S_current_mut>=protein_start)&(pos_in_S_current_mut<=protein_end))$RFU, na.rm=T)
  
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$prevalence_epitope_site_in_positive_patients_current_study[i] <- ifelse(nrow(subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut)))>0,subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut))$prevalence_in_positive_patients,0)
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$prevalence_epitope_site_in_negative_patients_current_study[i] <- ifelse(nrow(subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut)))>0,subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut))$prevalence_in_negative_patients,0)
  if (df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$prevalence_epitope_site_in_positive_patients_current_study[i]>0){
    df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$avg_RFU_epitope_site_in_positive_patients_current_study[i] <- subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut))$avg_RFU_in_positive_patients
  }
  if (df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$prevalence_epitope_site_in_negative_patients_current_study[i]>0){
    df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$avg_RFU_epitope_site_in_negative_patients_current_study[i] <- subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut))$avg_RFU_in_negative_patients
  }
  
  v <- subset(df_high_confidence_epitope_metrics,(ORF=="S")&(pos_in_protein==pos_in_S_current_mut))$is_cross_reactive
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$is_detected_as_cross_reactive_in_current_study[i] <- ifelse(length(v)>0,v,FALSE)
  
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$month_first_time_detected[i] <- min(subset(df_variants_NCBI_SRA_amplicon,short_label_mut==df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$mutation[i])$collection_date,na.rm=T)
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$month_first_time_detected[i] <- ifelse(df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$month_first_time_detected[i]==Inf,yes=NA,no=df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$month_first_time_detected[i])
  df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$lst_candidate_country_first_apparition[i] <- paste0(names(table(subset(df_variants_NCBI_SRA_amplicon,(short_label_mut==df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$mutation[i])&(collection_date==df_metrics_S_region_mutations_of_interest_in_epitopes_NCBI_SRA_amplicon$month_first_time_detected[i]))$seq_country)),collapse="/")
  
}

#list of non-synonymous mutations
v_lst_ns_mut <- sort(unique(subset(df_variants_NCBI_SRA_amplicon,(!is.na(mutation_type))&(mutation_type=="Non-Synonymous"))$mutation_name))

#metrics_ALL_MISSENSE_mutations_sites (only select N-S mutations in epitope sites)
df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon <- data.frame(mutation=unique(subset(df_variants_NCBI_SRA_amplicon,(mutation_type=="Non-Synonymous")&(is_epitope_related))$mutation_name),complete_mut_name=NA,protein=NA,pos_in_prot=NA,avg_RFU=NA,prevalence_mutation=NA,Wave_first_time_observed=NA,avg_RFU_epitope_site_in_positive_patients_current_study=NA,avg_RFU_epitope_site_in_negative_patients_current_study=NA,is_detected_as_cross_reactive_in_current_study=NA,prevalence_epitope_site_in_positive_patients_current_study=NA,prevalence_epitope_site_in_negative_patients_current_study=NA,month_first_time_detected=NA,lst_candidate_country_first_apparition=NA,stringsAsFactors = F)
for (i in 1:nrow(df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon)){
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$complete_mut_name[i] <- df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i]
  current_mut_short_label <- subset(df_variants_NCBI_SRA_amplicon,mutation_name==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i])$short_label_mut[1]
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i] <- subset(df_variants_NCBI_SRA_amplicon,mutation_name==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i])$ORF[1]
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i] <- subset(df_variants_NCBI_SRA_amplicon,mutation_name==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i])$pos_in_ORF_protein_seq[1]
  #pos_in_prot in actually pos_in_ORF from last line
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$avg_RFU[i] <- mean(subset(df_epitopes,(!is.na(Mapping_region))&(Mapping_region==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]>=protein_start)&(df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]<=protein_end))$RFU, na.rm=T)
  #mutation name becomes mut_short_label from here on
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i] <- current_mut_short_label
  
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$prevalence_mutation[i] <- v_prevalence_mutations[current_mut_short_label]
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$Wave_first_time_observed[i] <- v_wave_of_emergence_of_mutations_NCBI_SRA_amplicon[current_mut_short_label]
  
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$prevalence_epitope_site_in_positive_patients_current_study[i] <- ifelse(nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i])))>0,subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]))$prevalence_in_positive_patients,0)
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$prevalence_epitope_site_in_negative_patients_current_study[i] <- ifelse(nrow(subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i])))>0,subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]))$prevalence_in_negative_patients,0)
  
  if (df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$prevalence_epitope_site_in_positive_patients_current_study[i]>0){
    df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$avg_RFU_epitope_site_in_positive_patients_current_study[i] <- subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]))$avg_RFU_in_positive_patients
  }
  if (df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$prevalence_epitope_site_in_negative_patients_current_study[i]>0){
    df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$avg_RFU_epitope_site_in_negative_patients_current_study[i] <- subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]))$avg_RFU_in_negative_patients
  }
  
  v <- subset(df_high_confidence_epitope_metrics,(ORF==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$protein[i])&(pos_in_protein==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$pos_in_prot[i]))$is_cross_reactive
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$is_detected_as_cross_reactive_in_current_study[i] <- ifelse(length(v)>0,v,FALSE)
  
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$month_first_time_detected[i] <- min(subset(df_variants_NCBI_SRA_amplicon,short_label_mut==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i])$collection_date,na.rm=T)
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$month_first_time_detected[i] <- ifelse(df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$month_first_time_detected[i]==Inf,yes=NA,no=df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$month_first_time_detected[i])
  df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$lst_candidate_country_first_apparition[i] <- paste0(names(table(subset(df_variants_NCBI_SRA_amplicon,(short_label_mut==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$mutation[i])&(collection_date==df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$month_first_time_detected[i]))$seq_country)),collapse="/")
  
}
#measure difference in avg_RFU between Covid19- and Covid19+
df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$Difference_avg_RFU_epitope_sites_neg_vs_pos_patients_current_study <- df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$avg_RFU_epitope_site_in_negative_patients_current_study - df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon$avg_RFU_epitope_site_in_positive_patients_current_study
write.table(x=df_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon,file = paste0(output_workspace,"Table_immunological_metrics_ALL_MISSENSE_mutations_sites_NCBI_SRA_amlicon.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

#function that converts genomic mutation into ORF protein seq mutation
convert_genomic_mut_to_ORF_prot_mut <- function(the_genomic_mut){
  ref_nucl <- substr(the_genomic_mut,1,1)
  the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
  new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(ref_codon)!=3){
    stop("codon length should be 3!")
  }
  pos_in_ORF_prot_seq <- unname(ceiling((the_position - v_start_orfs[the_orf] + 1)/3))
  ref_aa <- translate_seq(the_codon = ref_codon)
  new_aa <- translate_seq(the_codon = mut_codon)
  return(as.character(paste0(the_orf,":",ref_aa,pos_in_ORF_prot_seq,new_aa)))
}#function that converts genomic mutation into ORF protein seq mutation
convert_genomic_mut_to_ORF_prot_mut <- function(the_genomic_mut){
  ref_nucl <- substr(the_genomic_mut,1,1)
  the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
  new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
  the_orf <- find_ORF_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(ref_codon)!=3){
    stop("codon length should be 3!")
  }
  pos_in_ORF_prot_seq <- unname(ceiling((the_position - v_start_orfs[the_orf] + 1)/3))
  ref_aa <- translate_seq(the_codon = ref_codon)
  new_aa <- translate_seq(the_codon = mut_codon)
  return(as.character(paste0(the_orf,":",ref_aa,pos_in_ORF_prot_seq,new_aa)))
}
#function that converts genomic mutation into protein mutation
convert_genomic_mut_to_prot_mut <- function(the_genomic_mut){
  ref_nucl <- substr(the_genomic_mut,1,1)
  the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
  new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
  the_orf <- find_ORF_of_mutation(the_position)
  the_gene <- find_gene_of_mutation(the_position)
  if (is.na(the_orf)||(grepl(pattern = "UTR",x = the_orf,fixed = TRUE))){
    return(NA)
  }else{
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
  }
  if (nchar(ref_codon)!=3){
    stop("codon length should be 3!")
  }
  pos_in_prot <- unname(ceiling((the_position - v_start_genes[the_gene] + 1)/3))
  ref_aa <- translate_seq(the_codon = ref_codon)
  new_aa <- translate_seq(the_codon = mut_codon)
  return(as.character(paste0(the_gene,":",ref_aa,pos_in_prot,new_aa)))
}

#function that converts genomic mutation into ORF protein seq mutation
find_ORF9C_prot_mut_from_N_prot_mut <- function(the_N_prot_mut){
  if (!"short_genomic_mut"%in%names(df_variants_NCBI_SRA_amplicon)){
    df_variants_NCBI_SRA_amplicon$short_genomic_mut <- paste0(df_variants_NCBI_SRA_amplicon$Ref,df_variants_NCBI_SRA_amplicon$Position,df_variants_NCBI_SRA_amplicon$VarAllele)
  }
  delim_pos <- gregexpr(pattern = ":",text = the_N_prot_mut,fixed = T)[[1]][1]
  the_N_prot_mut <- substr(the_N_prot_mut,delim_pos+1,nchar(the_N_prot_mut))
  old_aa <- substr(the_N_prot_mut,1,1)
  pos_in_N_prot <- as.integer(substr(the_N_prot_mut,2,nchar(the_N_prot_mut)-1))
  new_aa <-substr(the_N_prot_mut,nchar(the_N_prot_mut),nchar(the_N_prot_mut))
  genomic_position <- v_start_orfs["N"] + (pos_in_N_prot*3) - 1
  pos_in_codon <- ((genomic_position - v_start_orfs["N"] + 1)%%3)+(3*as.integer(((genomic_position - v_start_orfs["N"] + 1)%%3)==0))
  if (pos_in_codon==1){
    v_genomic_positions <- genomic_position:(genomic_position+2)
  }else if (pos_in_codon==2){
    v_genomic_positions <- (genomic_position-1):(genomic_position+1)
  }else if (pos_in_codon==3){
    v_genomic_positions <- (genomic_position-2):(genomic_position)
  }else{
    stop("Codon position must be between 1 and 3!!!")
  }
  possible_nucls <- c("A","T","C","G")
  the_orf <- "N"
  candidate_genomic_muts_N_prot <- NULL
  for (the_position in v_genomic_positions){
    pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
    if (pos_in_codon==1){
      ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
      ref_nucl <- substr(ref_codon,1,1)
    }else if (pos_in_codon==2){
      ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
      ref_nucl <- substr(ref_codon,2,2)
    }else if (pos_in_codon==3){
      ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
      ref_nucl <- substr(ref_codon,3,3)
    }else{
      stop("Codon position must be between 1 and 3!!!")
    }
    for (current_allele in setdiff(possible_nucls,ref_nucl)){
      if (pos_in_codon==1){
        current_mut_codon <- paste0(current_allele,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
        mut_codon_aa <- translate_seq(the_codon = current_mut_codon)
      }else if (pos_in_codon==2){
        current_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),current_allele,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
        mut_codon_aa <- translate_seq(the_codon = current_mut_codon)
      }else if (pos_in_codon==3){
        current_mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),current_allele)
        mut_codon_aa <- translate_seq(the_codon = current_mut_codon)
      }else{
        stop("Codon position must be between 1 and 3!!!")
      }
      if (mut_codon_aa == new_aa){
        candidate_genomic_muts_N_prot <- c(candidate_genomic_muts_N_prot,(paste0(ref_nucl,the_position,current_allele)))
      }
    }
  }
  
  if (length(candidate_genomic_muts_N_prot)==0){
    return(NA)
  }else{
    #if only one of the candidate mutation has been observed, return it eslse return all candidates
    if (sum(vapply(X = candidate_genomic_muts_N_prot,FUN = function(x) (x%in%df_variants_NCBI_SRA_amplicon$short_genomic_mut),FUN.VALUE = F))==1){
      candidate_genomic_muts_N_prot <- candidate_genomic_muts_N_prot[vapply(X = candidate_genomic_muts_N_prot,FUN = function(x) (x%in%df_variants_NCBI_SRA_amplicon$short_genomic_mut),FUN.VALUE = F)]
    }
    #print(candidate_genomic_muts_N_prot)
    v_orf9c_prot_muts <- NULL
    for (the_genomic_mut in candidate_genomic_muts_N_prot){
      ref_nucl <- substr(the_genomic_mut,1,1)
      the_position <- as.integer(substr(the_genomic_mut,2,nchar(the_genomic_mut)-1))
      new_nucl <- substr(the_genomic_mut,nchar(the_genomic_mut),nchar(the_genomic_mut))
      the_orf <- "ORF9c"
      the_gene <- "ORF9c"
      pos_in_codon <- ((the_position - v_start_orfs[the_orf] + 1)%%3)+(3*as.integer(((the_position - v_start_orfs[the_orf] + 1)%%3)==0))
      if (pos_in_codon==1){
        ref_codon <- substr(x = genome_refseq,start = the_position,stop = the_position+2)
        mut_codon <- paste0(new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+2))
      }else if (pos_in_codon==2){
        ref_codon <- substr(x = genome_refseq,start = the_position-1,stop = the_position+1)
        mut_codon <- paste0(substr(x = genome_refseq,start = the_position-1,stop = the_position-1),new_nucl,substr(x = genome_refseq,start = the_position+1,stop = the_position+1))
      }else if (pos_in_codon==3){
        ref_codon <- substr(x = genome_refseq,start = the_position-2,stop = the_position)
        mut_codon <- paste0(substr(x = genome_refseq,start = the_position-2,stop = the_position-1),new_nucl)
      }else{
        stop("Codon position must be between 1 and 3!!!")
      }
      if (nchar(ref_codon)!=3){
        stop("codon length should be 3!")
      }
      pos_in_ORF_prot_seq <- unname(ceiling((the_position - v_start_orfs[the_orf] + 1)/3))
      ref_aa <- translate_seq(the_codon = ref_codon)
      new_aa <- translate_seq(the_codon = mut_codon)
      v_orf9c_prot_muts <- c(v_orf9c_prot_muts,as.character(paste0(the_orf,":",ref_aa,pos_in_ORF_prot_seq,new_aa)))
      #if only one of the candidate mutation has been observed, return it eslse return all candidates
      if (sum(vapply(X = v_orf9c_prot_muts,FUN = function(x) (x%in%df_variants_NCBI_SRA_amplicon$label_mut_ORF_prot_effect),FUN.VALUE = F))==1){
        return(v_orf9c_prot_muts[vapply(X = v_orf9c_prot_muts,FUN = function(x) (x%in%df_variants_NCBI_SRA_amplicon$label_mut_ORF_prot_effect),FUN.VALUE = F)])
      }else{
        return(v_orf9c_prot_muts)
      }
    }
  }
}

#df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$label_mut_prot_effect <- unname(vapply(X = df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages$mutation_name,FUN = convert_genomic_mut_to_prot_mut,FUN.VALUE = c("")))
#Define lineages of interest
v_lineages_of_interest <- c("Alpha (B.1.1.7)","Beta (B.1.351)","Gamma (P.1)","Epsilon (B.1.427+B.1.429)","Kappa+Delta (B.1.617.X)","Eta (B.1.525)","Iota (B.1.526)","Lambda (C.37)","Omicron (B.1.1.529)","Theta (P.3)","Zeta (P.2)","A.2.5","B.1.1.318","B.1.1.519","B.1.466.2 ","B.1.621","B.1.214.2","AV.1","AT.1","C.36.3","R.1","R.2")
#Alpha (B.1.1.7)
v_Alpha_signature_muts <- c("N:M1Stop","N:D3L","N:R203K","N:G204R","N:S235F","ORF1a:T1001I","ORF1a:A1708D","ORF1a:I2230T","ORF1b:P314L","ORF8:R52I","ORF8:Y73C","S:N501Y","S:A570D","S:D614G","S:P681H","S:T716I","S:S982A","S:D1118H")
#Beta (B.1.351)
v_Beta_signature_muts <- c("E:P71L","N:T205I","ORF1a:T265I","ORF1a:K1655N","ORF1a:K3353R","ORF3a:Q57H","ORF3a:S171L","S:D80A","S:D215G","S:L242H","S:K417N","S:D614G","S:A701V")
#Gamma (P.1)
v_Gamma_signature_muts <- c("N:P80R","N:R203K","N:G204R","ORF1a:S1188L","ORF1a:K1795Q","ORF1b:P314L","ORF1b:E1264D","ORF3a:S253P","ORF8:E92K","S:L18F","S:T20N","S:P26S","S:D138Y","S:R190S","S:K417T","S:E484K","S:N501Y","S:D614G","S:H655Y","S:T1027I","S:V1176F")
#Epsilon (B.1.427/B.1.429)
v_Epsilon_signature_muts <- c("N:T205I","ORF1a:T265I","ORF1a:I4205V","ORF1a:S3158T","ORF1b:P314L","ORF1b:P976L","ORF1b:D1183Y","ORF3a:Q57H","S:S13I","S:W152C","S:L452R","S:D614G")
#Kappa/Delta (B.1.617.X)
v_Kappa_Delta_signature_muts <- c("M:I82T","N:D63G","N:R203M","N:D377Y","ORF1b:P314L","ORF3a:S26L","ORF7a:V82A","S:T19R","S:L452R","S:T478K","S:D614G","S:P681R","S:D950N")
#Eta (B.1.525)
v_Eta_signature_muts <- c("M:I82T","N:S2Y","N:A12G","N:T205I","ORF1a:T2007I","ORF1b:P314F","S:A67V","S:E484K","S:D614G","S:Q677H","S:F888L")
#Iota (B.1.526)
v_Iota_signature_muts <- c("ORF1a:T265I","ORF1a:L3201P","ORF1b:P314L","ORF3a:P42L","ORF3a:Q57H","ORF8:T11I","S:D614G")
#Lambda (C.37)
v_lambda_signature_muts <- c("N:P13L","N:R203K","N:G204R","N:G214C","ORF1a:T1246I","ORF1a:P2287S","ORF1a:F2387V","ORF1a:L3201P","ORF1a:T3255I","ORF1a:G3278S","ORF1b:P314L","S:G75V","S:T76I","S:L452Q","S:F490S","S:D614G","S:T859N")
#Theta (P.3)
v_P3_signature_muts <- c("N:R203K","N:G204R","ORF1a:D1554G","ORF1a:S2625F","ORF1a:D2980N","ORF1a:L3201P","ORF1a:D3681E","ORF1a:L3930F","ORF1b:P314L","ORF1b:A1291V","ORF8:K2Q","S:E484K","S:N501Y","S:D614G","S:P681H","S:E1092K","S:H1101Y","S:V1176F")
#Zeta (P.2)
v_Zeta_signature_muts <- c("N:A119S","N:R203K","N:G204R","N:M234I","ORF1a:L3468V","ORF1a:L3930F","ORF1b:P314L","S:E484K","S:D614G","S:V1176F")
#A.2.5
v_A_2_5_signature_muts <- c("N:S197L","N:M234I","N:P365S","N:P383L","ORF1a:L4F","ORF1a:K1657E","ORF1a:F3071Y","ORF1a:T3255I","ORF1a:H3580Q","ORF1b:P1000L","ORF3a:S74F","ORF3a:G196V","ORF8:L84S","S:L452R","S:D614G")
#B.1.1.318
v_B_1_1_318_signature_muts <- c("M:I82T","N:R203K","N:G204R","N:A208G","ORF1a:E1196V","ORF1a:K2511N","ORF1a:T2936I","ORF1a:A3209V","ORF1a:T3284I","ORF1b:P314L","ORF1b:V2371M","ORF8:F3Stop","S:T95I","S:E484K","S:D614G","S:P681H","S:D796H")
#B.1.1.519
v_B_1_1_519_signature_muts <- c("N:R203K","N:G204R","ORF1a:P959S","ORF1a:T3255I","ORF1a:I3618V","ORF1a:T4175I","ORF1b:P314L","S:T478K","S:D614G","S:P681H","S:T732A")
#B.1.466.2
v_B_1_466_2_signature_muts <- c("N:T205I","ORF1a:T1168I","ORF1a:P1640L","ORF1b:P314L","ORF3a:Q57H","S:N439K","S:D614G")
#B.1.621
v_B_1_621_signature_muts <- c("N:T205I","ORF1a:T1055A","ORF1a:T1538I","ORF1a:T3255I","ORF1a:Q3729R","ORF1b:P314L","ORF1b:P1342S","ORF3a:Q57H","ORF3a:N257Stop","ORF8:T11K","ORF8:P38S","S:T95I","S:R346K","S:E484K","S:N501Y","S:D614G","S:P681H","S:D950N")
#B.1.214.2
v_B_1_214_2_signature_muts <- c("N:T205I","ORF1a:I1398V","ORF1a:T1881I","ORF1a:A4016V","ORF1b:P314L","S:Q414K","S:D614G","S:T716I")
#AV.1
v_AV_1_signature_muts <- c("M:A63T","M:H125Y","N:M1Stop","N:I157V","N:R203K","N:G204R","ORF1a:D75G","ORF1a:G519S","ORF1a:A591V","ORF1a:H1160Y","ORF1a:P1640L","ORF1a:K1745I","ORF1a:N2405S","ORF1a:A3209V","ORF1b:P314L","ORF1b:S425A","ORF1b:T1404M","ORF1b:A1643V","ORF3a:N257Stop","S:D80G","S:T95I","S:G142D","S:E484K","S:D614G","S:P681H","S:I1130V","S:D1139H")
#AT.1
v_AT_1_signature_muts <- c("M:L16I","N:R203K","N:G204R","ORF1a:S376L","ORF1a:V1006F","ORF1a:T2247N","ORF1a:T3255I","ORF1a:Q3729R","ORF1a:S4119T","ORF1b:P314L","ORF1b:T1173N","ORF1b:V1905L","ORF1b:A2431V","ORF3a:L95M","S:P9L","S:C136Y","S:D215G","S:H245P","S:E484K","S:D614G","S:E780K")
#C.36.3
v_C_36_3_signature_muts <- c("M:I82T","N:M1Stop","N:R203K","N:G204R","N:G212V","ORF1a:E102K","ORF1a:A859V","ORF1a:T1246I","ORF1a:D1639N","ORF1a:P2287S","ORF1a:D2980N","ORF1a:D3222N","ORF1a:G3278S","ORF1a:S3687L","ORF1a:L3691S","ORF1a:T4090I","ORF1b:P314L","ORF1b:D1028Y","ORF7b:A43S","S:S12F","S:W152R","S:R346S","S:L452R","S:D614G","S:Q677H","S:A899S")
#R.1
v_R_1_signature_muts <- c("M:F28L","N:S187L","N:R203K","N:G204R","N:Q418H","ORF1b:P314L","ORF1b:G1362R","ORF1b:P1936H","S:W152L","S:E484K","S:D614G","S:G769V")
#R.2
v_R_2_signature_muts <- c("N:R203K","N:G204R","ORF1a:A1049V","ORF1a:K1202N","ORF1a:Q1592R","ORF1a:D4085G","ORF1b:P314L","S:D614G","S:Q677H","S:T732S")
#Omicron (B.1.1.529)
v_Omicron_B_1_1_529_signature_muts <- c("S:Y505H","N:G204R","N:R203K","ORF1a:A2710T","ORF1a:I3758V","ORF1a:K856R","ORF1a:T3255I","ORF1b:I1566V","ORF1b:P314L","S:D614G","S:E484A","S:G339D","S:G446S","S:G496S","S:H655Y","S:K417N","S:N440K","S:N501Y","S:N679K","S:P681H","S:Q493R","S:Q498R","S:S371L","S:S373P","S:S375F","S:S477N","S:T478K","S:T547K","ORF1a:L2084I","S:A67V","S:Y145D","E:T9I","M:A63T","M:D3G","M:Q19E","S:L212I","S:T95I","ORF1a:P3395H")

#Lineages of interest (high confidence epitope analysis)
df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- readRDS(file = paste0(output_workspace,"Table_df_all_mutations_prevalence_in_lineages.rds"))
df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages <- subset(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages,!lineage%in%c("B.1.1.7","B.1.351","P.1","B.1.427","B.1.429","B.1.617.X","B.1.525","B.1.526","C.37","P.3","P.2","A.2.5","B.1.160","B.1.177","B.1.1.318","B.1.1.519","B.1.466.2 ","B.1.621","B.1.214.2","AV.1","AT.1","C.36.3","R.1","R.2","B.1.1.529"))
df_signature_mutations <- subset(df_prevalence_mut_of_interest_in_NCBI_WORLDWIDE_lineages,(!is.na(lineage))&(prevalence>=0.9))
df_signature_mutations$label_mut_prot_effect <- NA
df_signature_mutations$is_not_sense <- NA
for (i in 1:nrow(df_signature_mutations)){
  df_signature_mutations$label_mut_ORF_prot_effect[i] <- convert_genomic_mut_to_ORF_prot_mut(the_genomic_mut = df_signature_mutations$mutation_name[i])
  df_signature_mutations$label_mut_prot_effect[i] <- convert_genomic_mut_to_prot_mut(the_genomic_mut = df_signature_mutations$mutation_name[i])
  # print(i)
  pos_delim_char <- gregexpr(pattern = ":",text = df_signature_mutations$label_mut_prot_effect[i],fixed = T)[[1]][1]
  current_old_aa <- substr(df_signature_mutations$label_mut_prot_effect[i],pos_delim_char+1,pos_delim_char+1)
  if (grepl(pattern = "Stop",x = df_signature_mutations$label_mut_prot_effect[i],fixed = T)){
    current_new_aa <- "Stop"
  }else{
    current_new_aa <- substr(df_signature_mutations$label_mut_prot_effect[i],nchar(df_signature_mutations$label_mut_prot_effect[i]),nchar(df_signature_mutations$label_mut_prot_effect[i]))
  }
  #determine if mutation is not synonymous (so either non-synonymous (Missense) or nonsense)
  df_signature_mutations$is_not_sense[i] <- (!is.na(current_old_aa))&(current_old_aa!=current_new_aa)
} 
df_signature_mutations <- subset(df_signature_mutations,(!is.na(label_mut_ORF_prot_effect))&(!is.na(is_not_sense))&(is_not_sense))[,c("label_mut_ORF_prot_effect","lineage")]
df_signature_mutations <- df_signature_mutations[,c("label_mut_ORF_prot_effect","lineage")]
df_signature_mutations_to_add <- data.frame(label_mut_ORF_prot_effect=c(v_Alpha_signature_muts,v_Beta_signature_muts,v_Gamma_signature_muts,v_Epsilon_signature_muts,v_Kappa_Delta_signature_muts,v_Eta_signature_muts,v_Iota_signature_muts,v_lambda_signature_muts,v_Omicron_B_1_1_529_signature_muts,v_P3_signature_muts,v_Zeta_signature_muts,v_A_2_5_signature_muts,v_B_1_1_318_signature_muts,v_B_1_1_519_signature_muts,v_B_1_466_2_signature_muts,v_B_1_621_signature_muts,v_B_1_214_2_signature_muts,v_AV_1_signature_muts,v_AT_1_signature_muts,v_C_36_3_signature_muts,v_R_1_signature_muts,v_R_2_signature_muts),lineage=rep(x = v_lineages_of_interest,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         times=c(length(v_Alpha_signature_muts),length(v_Beta_signature_muts),length(v_Gamma_signature_muts),length(v_Epsilon_signature_muts),length(v_Kappa_Delta_signature_muts),length(v_Eta_signature_muts),length(v_Iota_signature_muts),length(v_lambda_signature_muts),length(v_Omicron_B_1_1_529_signature_muts),length(v_P3_signature_muts),length(v_Zeta_signature_muts),length(v_A_2_5_signature_muts),length(v_B_1_1_318_signature_muts),length(v_B_1_1_519_signature_muts),length(v_B_1_466_2_signature_muts),length(v_B_1_621_signature_muts),length(v_B_1_214_2_signature_muts),length(v_AV_1_signature_muts),length(v_AT_1_signature_muts),length(v_C_36_3_signature_muts),length(v_R_1_signature_muts),length(v_R_2_signature_muts))),stringsAsFactors = F)
df_signature_mutations <- unique(rbind(df_signature_mutations,df_signature_mutations_to_add))
df_signature_mutations <- subset(df_signature_mutations,vapply(X = df_signature_mutations$label_mut_ORF_prot_effect,FUN = function (x) !grepl(pattern = "Stop",x = x,fixed = T),FUN.VALUE = F))
v_start_orfs_cap <- v_start_orfs
names(v_start_orfs_cap) <- toupper(names(v_start_orfs))
for (i in 1:nrow(df_signature_mutations)){
  current_pos_delim_char <- gregexpr(pattern = ":",text = df_signature_mutations$label_mut_ORF_prot_effect[i],fixed = T)[[1]][1]
  df_signature_mutations$ORF[i] <- toupper(substr(df_signature_mutations$label_mut_ORF_prot_effect[i],1,current_pos_delim_char-1))
  current_pos_in_ORF <- as.integer(substr(df_signature_mutations$label_mut_ORF_prot_effect[i],current_pos_delim_char+2,nchar(df_signature_mutations$label_mut_ORF_prot_effect[i])-1))
  genomic_position <- unname(v_start_orfs_cap[df_signature_mutations$ORF[i]] + (current_pos_in_ORF*3) - 1)
  pos_in_codon <- unname((genomic_position - v_start_orfs_cap[df_signature_mutations$ORF[i]] + 1)%%3)+(3*as.integer(((genomic_position - v_start_orfs_cap[df_signature_mutations$ORF[i]] + 1)%%3)==0))
  if (pos_in_codon==1){
    v_genomic_positions <- genomic_position:(genomic_position+2)
  }else if (pos_in_codon==2){
    v_genomic_positions <- (genomic_position-1):(genomic_position+1)
  }else if (pos_in_codon==3){
    v_genomic_positions <- (genomic_position-2):(genomic_position)
  }else{
    stop("Codon position must be between 1 and 3!!!")
  }
  
  df_signature_mutations$is_in_epitope_site[i] <- any(vapply(X = v_genomic_positions,FUN = function(x) x %in%v_unique_epitope_positions,FUN.VALUE = F)) 
  df_signature_mutations$is_in_cr_epitope_site[i] <- any(vapply(X = v_genomic_positions,FUN = function(x) x %in%v_position_Cross_reactive_epitopes,FUN.VALUE = F))
}
df_variants_NCBI_SRA_amplicon$label_mut_ORF_prot_effect <- paste0(toupper(df_variants_NCBI_SRA_amplicon$ORF),":",df_variants_NCBI_SRA_amplicon$old_aa,df_variants_NCBI_SRA_amplicon$pos_in_ORF_protein_seq,df_variants_NCBI_SRA_amplicon$new_aa)

for (i in 1:nrow(df_signature_mutations)){
  if (df_signature_mutations$ORF[i]=="N"){
    current_pos_delim_char <- gregexpr(pattern = ":",text = df_signature_mutations$label_mut_ORF_prot_effect[i],fixed = T)[[1]][1]
    #print((substr(df_signature_mutations$label_mut_ORF_prot_effect[i],current_pos_delim_char+2,nchar(df_signature_mutations$label_mut_ORF_prot_effect[i])-1)))
    current_pos_in_ORF <- as.integer(substr(df_signature_mutations$label_mut_ORF_prot_effect[i],current_pos_delim_char+2,nchar(df_signature_mutations$label_mut_ORF_prot_effect[i])-1))
    if (current_pos_in_ORF%in%(unname(ceiling((v_start_orfs["ORF9c"] - v_start_orfs["N"] + 1)/3)):unname(ceiling((v_end_orfs["ORF9c"] - v_start_orfs["N"] + 1)/3)))){
      the_line_to_add <- df_signature_mutations[i,]
      the_line_to_add$ORF <- "ORF9c"
      the_line_to_add$label_mut_ORF_prot_effect <- find_ORF9C_prot_mut_from_N_prot_mut(the_N_prot_mut = df_signature_mutations$label_mut_ORF_prot_effect[i])
      df_signature_mutations <- rbind(df_signature_mutations,the_line_to_add)
    }
  }
}
df_variants_NCBI_SRA_amplicon$label_mut_ORF_prot_effect <- paste0(toupper(df_variants_NCBI_SRA_amplicon$ORF),":",df_variants_NCBI_SRA_amplicon$old_aa,df_variants_NCBI_SRA_amplicon$pos_in_ORF_protein_seq,df_variants_NCBI_SRA_amplicon$new_aa)
df_signature_mutations$is_nonsyn <- NA
df_signature_mutations$ORF_prot_site <- NA
for (i in 1:nrow(df_signature_mutations)){
  the_prot_mut <- df_signature_mutations$label_mut_ORF_prot_effect[i]
  delim_pos <- gregexpr(pattern = ":",text = the_prot_mut,fixed = T)[[1]][1]
  the_prot_mut <- substr(the_prot_mut,delim_pos+1,nchar(the_prot_mut))
  old_aa <- substr(the_prot_mut,1,1)
  new_aa <-substr(the_prot_mut,nchar(the_prot_mut),nchar(the_prot_mut))
  df_signature_mutations$is_nonsyn[i] <- (!is.na(old_aa))&(old_aa!=new_aa)
  df_signature_mutations$ORF_prot_site[i] <- substr(the_prot_mut,2,nchar(the_prot_mut)-1)
}
df_signature_mutations <- subset(df_signature_mutations,is_nonsyn)
df_signature_mutations$is_nonsyn <- NULL 
v_mut_name_to_label_mut_ORF_prot_effect <- unique(df_variants_NCBI_SRA_amplicon[,c("mutation_name","label_mut_ORF_prot_effect")])$label_mut_ORF_prot_effect
names(v_mut_name_to_label_mut_ORF_prot_effect) <- unique(df_variants_NCBI_SRA_amplicon[,c("mutation_name","label_mut_ORF_prot_effect")])$mutation_name
saveRDS(object = df_signature_mutations,file = paste0(output_workspace,"Table_Missense_and_Nonsense_signature_mutations_prevalence_in_SC2_lineages_consensus_sequences_as_of_2021_01_16_plus_VOCs.rds"))

#VOC nb mutations in epitope sites and c-r epitope sites https://outbreak.info/situation-reports and https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/ as of 2021-07-06
df_signature_mutations_for_VOCs <- subset(df_signature_mutations,lineage%in%v_lineages_of_interest)
df_signature_mutations_for_VOCs$ORF <- ifelse(df_signature_mutations_for_VOCs$ORF=="orf1a","ORF1A",df_signature_mutations_for_VOCs$ORF)
df_signature_mutations_for_VOCs$ORF <- ifelse(df_signature_mutations_for_VOCs$ORF=="orf1b","ORF1B",df_signature_mutations_for_VOCs$ORF)
df_signature_mutations_for_VOCs$ORF <- ifelse(df_signature_mutations_for_VOCs$ORF=="ORF1a","ORF1A",df_signature_mutations_for_VOCs$ORF)
df_signature_mutations_for_VOCs$ORF <- ifelse(df_signature_mutations_for_VOCs$ORF=="ORF1b","ORF1B",df_signature_mutations_for_VOCs$ORF)
df_VOC_number_of_sig_muts_in_epitope_sites<- aggregate(df_signature_mutations_for_VOCs$is_in_epitope_site,by=list(VOC=df_signature_mutations_for_VOCs$lineage),FUN = function(x) sum(x,na.rm=T))
df_VOC_number_of_sig_muts_in_epitope_sites <- subset(df_VOC_number_of_sig_muts_in_epitope_sites,x>0)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF <- aggregate(df_signature_mutations_for_VOCs$is_in_epitope_site,by=list(VOC=df_signature_mutations_for_VOCs$lineage,ORF=df_signature_mutations_for_VOCs$ORF),FUN = function(x) sum(x,na.rm=T))
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF <- subset(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF,x>0)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="orf1a","ORF1A",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="orf1b","ORF1B",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="ORF1a","ORF1A",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="ORF1b","ORF1B",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
palette_orfs_epitopes2 <- palette_orfs_epitopes
names(palette_orfs_epitopes2) <- ifelse(names(palette_orfs_epitopes2)=="orf1a","ORF1A",names(palette_orfs_epitopes2))
names(palette_orfs_epitopes2) <- ifelse(names(palette_orfs_epitopes2)=="orf1b","ORF1B",names(palette_orfs_epitopes2))
names(palette_orfs_epitopes2) <- ifelse(names(palette_orfs_epitopes2)=="ORF1a","ORF1A",names(palette_orfs_epitopes2))
names(palette_orfs_epitopes2) <- ifelse(names(palette_orfs_epitopes2)=="ORF1b","ORF1B",names(palette_orfs_epitopes2))
palette_orfs_epitopes2 <- c(palette_orfs_epitopes2,c("N sites over-\nlapping ORF9c"="black"))
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="ORF9c","N sites over-\nlapping ORF9c",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF=="N sites overlapping ORF9c","N sites over-\nlapping ORF9c",df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF)
ggplot(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF) + geom_bar(aes(factor(VOC,levels=df_VOC_number_of_sig_muts_in_epitope_sites$VOC[order(df_VOC_number_of_sig_muts_in_epitope_sites$x,decreasing = T)]), x,fill=factor(ORF,levels = names(sort(vapply(X = names(palette_orfs_epitopes2),FUN = function(the_orf) max(subset(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF,ORF==the_orf)$x), FUN.VALUE = c(0)))))), stat="identity",position = "stack") + xlab("SARS-CoV-2 variant of concern") + ylab("Number of nonsynonymous signature\nmutations located at epitope sites") + scale_y_continuous(limits = c(0,max(df_VOC_number_of_sig_muts_in_epitope_sites$x,na.rm=T)+2),breaks=seq(0,max(df_VOC_number_of_sig_muts_in_epitope_sites$x,na.rm=T)+1,4)) + theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust = 1,size=11/2),axis.text.y = element_text(size=14/2),axis.title.y = element_text(size=16/2),axis.title.x = element_text(size=16/2),legend.title = element_text(size=13/2),legend.text = element_text(size=10/2),legend.key.size = unit(0.4, 'cm')) + scale_fill_manual(values = palette_orfs_epitopes2) + guides(fill=guide_legend(title="ORF"))
ggsave(filename = "VOCs_nb_NS_sig_mutations_in_epitope_sites_by_ORF.eps", path=output_workspace, width = 9.13, height = 9, units = "cm",device="eps")
#save table signature mutations with data about affected differential epitope sites
df_signature_mutations_for_VOCs$is_epitope_site_differential <- NA
for (i in 1:nrow(df_signature_mutations_for_VOCs)){
  the_subset_df <- subset(df_high_confidence_epitope_metrics,(ORF==df_signature_mutations_for_VOCs$ORF[i])&(pos_in_protein==df_signature_mutations_for_VOCs$ORF_prot_site[i]))
  if (nrow(the_subset_df)==0){
    df_signature_mutations_for_VOCs$is_epitope_site_differential[i] <- FALSE
  }else{
    if ((((the_subset_df$avg_RFU_in_positive_patients/the_subset_df$avg_RFU_in_negative_patients)>=2))&(the_subset_df$avg_RFU_in_negative_patients>0)){
      df_signature_mutations_for_VOCs$is_epitope_site_differential[i] <- TRUE
    }else{
      df_signature_mutations_for_VOCs$is_epitope_site_differential[i] <- FALSE
    }
  }
}
write.table(x=df_signature_mutations_for_VOCs,file = paste0(output_workspace,"Table_signature_mutations_VOCs_VUIs.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)
# v_orfs_length2 <- v_orfs_length
# names(v_orfs_length2) <- ifelse(names(v_orfs_length2)=="orf1a","ORF1A",names(v_orfs_length2))
# names(v_orfs_length2) <- ifelse(names(v_orfs_length2)=="orf1b","ORF1B",names(v_orfs_length2))
# names(v_orfs_length2) <- ifelse(names(v_orfs_length2)=="ORF1a","ORF1A",names(v_orfs_length2))
# names(v_orfs_length2) <- ifelse(names(v_orfs_length2)=="ORF1b","ORF1B",names(v_orfs_length2))
# names(v_orfs_length2) <- ifelse(names(v_orfs_length2)=="ORF9c","N sites overlapping ORF9c",names(v_orfs_length2))
# df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$density <- df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$x/(v_orfs_length2[df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF])
v_nb_epitopes_per_ORF <- vapply(X=names(v_orfs_length),FUN=function(x) length(intersect((v_start_orfs[x]:v_end_orfs[x]),v_unique_epitope_positions)),FUN.VALUE=c(0))
names(v_nb_epitopes_per_ORF) <- ifelse(names(v_nb_epitopes_per_ORF)=="orf1a","ORF1A",names(v_nb_epitopes_per_ORF))
names(v_nb_epitopes_per_ORF) <- ifelse(names(v_nb_epitopes_per_ORF)=="orf1b","ORF1B",names(v_nb_epitopes_per_ORF))
names(v_nb_epitopes_per_ORF) <- ifelse(names(v_nb_epitopes_per_ORF)=="ORF1a","ORF1A",names(v_nb_epitopes_per_ORF))
names(v_nb_epitopes_per_ORF) <- ifelse(names(v_nb_epitopes_per_ORF)=="ORF1b","ORF1B",names(v_nb_epitopes_per_ORF))
names(v_nb_epitopes_per_ORF) <- ifelse(names(v_nb_epitopes_per_ORF)=="ORF9c","N sites overlapping ORF9c",names(v_nb_epitopes_per_ORF))
df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$density <- df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$x/(v_nb_epitopes_per_ORF[df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF$ORF])
ggplot(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF) + geom_bar(aes(factor(VOC,levels=df_VOC_number_of_sig_muts_in_epitope_sites$VOC[order(df_VOC_number_of_sig_muts_in_epitope_sites$x,decreasing = T)]), density,fill=factor(ORF,levels = names(sort(vapply(X = names(palette_orfs_epitopes2),FUN = function(the_orf) max(subset(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF,ORF==the_orf)$x), FUN.VALUE = c(0)))))), stat="identity",position = "stack") + xlab("SARS-CoV-2 variant of concern") + ylab("Density of nonsynonymous signature\nmutations located at epitope sites") + scale_y_continuous(limits = c(0,0.035),breaks=seq(0,0.035,0.005)) + theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=10),axis.title.y = element_text(size=13),axis.title.x = element_text(size=13)) + scale_fill_manual(values = palette_orfs_epitopes2) + guides(fill=guide_legend(title="ORF"))
ggsave(filename = "VOCs_density_NS_sig_mutations_in_epitope_sites_by_ORF.eps", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200,device="eps")

df_VOC_number_of_sig_muts_in_cr_epitope_sites <- aggregate(df_signature_mutations_for_VOCs$is_in_cr_epitope_site,by=list(VOC=df_signature_mutations_for_VOCs$lineage),FUN = function(x) sum(x,na.rm=T))
df_VOC_number_of_sig_muts_in_cr_epitope_sites <- subset(df_VOC_number_of_sig_muts_in_cr_epitope_sites,x>0)
df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF <- aggregate(df_signature_mutations_for_VOCs$is_in_cr_epitope_site,by=list(VOC=df_signature_mutations_for_VOCs$lineage,ORF=df_signature_mutations_for_VOCs$ORF),FUN = function(x) sum(x,na.rm=T))
df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF <- subset(df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF,x>0)
df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF$ORF <- ifelse(df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF$ORF=="ORF9c","N sites overlapping ORF9c",df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF$ORF)
ggplot(df_VOC_number_of_sig_muts_in_cr_epitope_sites_by_ORF) + geom_bar(aes(factor(VOC,levels=df_VOC_number_of_sig_muts_in_cr_epitope_sites$VOC[order(df_VOC_number_of_sig_muts_in_cr_epitope_sites$x,decreasing = T)]), x,fill=factor(ORF,levels = names(sort(vapply(X = names(palette_orfs_epitopes2),FUN = function(the_orf) max(subset(df_VOC_number_of_sig_muts_in_epitope_sites_by_ORF,ORF==the_orf)$x), FUN.VALUE = c(0)))))), stat="identity",position = "stack") + xlab("SARS-CoV-2 variant of concern") + ylab("Number of nonsynonymous signature\nmutations located at cross-reactive epitope sites") + scale_y_continuous(limits = c(0,max(df_VOC_number_of_sig_muts_in_cr_epitope_sites$x,na.rm=T)+1),breaks=seq(0,max(df_VOC_number_of_sig_muts_in_cr_epitope_sites$x,na.rm=T)+1,1)) + theme_bw() + theme(axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=10),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12)) + scale_fill_manual(values = palette_orfs_epitopes2) + guides(fill=guide_legend(title="ORF"))
ggsave(filename = "VOCs_nb_NS_sig_mutations_in_Cross-reactive_epitope_sites_by_ORF.eps", path=output_workspace, width = 18.3, height = 15, units = "cm",dpi = 1200,device="eps")

df_lineage_number_of_sig_muts_in_epitope_sites <- aggregate(df_signature_mutations$is_in_epitope_site,by=list(lineage=df_signature_mutations$lineage),FUN = function(x) sum(x,na.rm=T))
df_lineage_number_of_sig_muts_in_epitope_sites$lineage_type <- ifelse(test = df_lineage_number_of_sig_muts_in_epitope_sites$lineage %in% v_lineages_of_interest,yes="VOCs and VUIs",no="Others")
ggplot(df_lineage_number_of_sig_muts_in_epitope_sites,aes(x=as.factor(lineage_type),y=x,fill=as.factor(lineage_type))) + geom_violin() + geom_jitter() + geom_boxplot(width=0.075,fill="white") + xlab("SARS-CoV-2 lineage type") + ylab("Number of nonsynoymous signature\nmutations located at epitope sites") + scale_y_continuous(limits = c(0,max(df_lineage_number_of_sig_muts_in_epitope_sites$x,na.rm=T)+1),breaks=seq(0,max(df_lineage_number_of_sig_muts_in_epitope_sites$x,na.rm=T)+1,4)) + theme_bw() + theme(axis.text.x = element_text(size=16/2),axis.text.y = element_text(size=14/2),axis.title.y = element_text(size=16/2),axis.title.x = element_text(size=16/2),legend.position="None",legend.title = element_text(size=13/2),legend.text = element_text(size=13/2))+ scale_fill_manual(values = c("VOCs and VUIs"="tan2","Others"="grey60")) + stat_compare_means(method = "wilcox")
ggsave(filename = "VOCs_VS_OTHER_SC2_lineages_nb_NS_sig_mutations_in_epitope_sites.eps", path=output_workspace, width = 9.26, height = 9.42, units = "cm",dpi = 1200,device="eps")

#RFU >=100
df_epitopes_RFU_100 <- read.csv2(file = paste0(output_workspace,"Epitopes_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_epitopes_RFU_100$Genomic_start <- as.integer(df_epitopes_RFU_100$Genomic_start)
df_epitopes_RFU_100$Genomic_End <- as.integer(df_epitopes_RFU_100$Genomic_End)

df_epitopes_RFU_100$Mapping_region <- vapply(X = df_epitopes_RFU_100$Genomic_End,FUN = function(x) return(find_ORF_of_mutation(the_site_position = x)),FUN.VALUE = c(""))
df_epitopes_RFU_100$Mapping_region <- factor(as.character(df_epitopes_RFU_100$Mapping_region),intersect(v_orfs,df_epitopes_RFU_100$Mapping_region))
df_epitopes_RFU_100$peptide_id <- paste0(df_epitopes_RFU_100$Mapping_region,"_",unname(v_lst_id_peptide_seq_in_order[df_epitopes_RFU_100$Peptide]))
df_epitopes_RFU_100 <- df_epitopes_RFU_100[,c("peptide_id",names(read.csv2(file = paste0(output_workspace,"Epitopes_mapped.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)))]
#exclude possible annotation mistakes
df_epitopes_RFU_100 <- subset(df_epitopes_RFU_100,vapply(X = 1:nrow(df_epitopes_RFU_100),FUN = function(i) (return(grepl(pattern = df_epitopes_RFU_100$Mapping_region[i],x = df_epitopes_RFU_100$Annotated_region[i],fixed = TRUE) )),FUN.VALUE = c(FALSE) ))
df_epitopes_RFU_100$Group <- as.character(df_epitopes_RFU_100$Group)
df_epitopes_RFU_100$RFU <- as.numeric(df_epitopes_RFU_100$RFU)
df_epitopes_RFU_100 <- subset(df_epitopes_RFU_100, RFU>=100)
df_epitopes_RFU_100$protein_start <- unname(ceiling((df_epitopes_RFU_100$Genomic_start - v_start_orfs[as.character(df_epitopes_RFU_100$Mapping_region)] + 1)/3))
df_epitopes_RFU_100$protein_end <- unname(ceiling((df_epitopes_RFU_100$Genomic_End - v_start_orfs[as.character(df_epitopes_RFU_100$Mapping_region)] + 1)/3))

#Get p-values differential RFU
df_high_confidence_epitope_metrics_RFU_100 <- df_high_confidence_epitope_metrics
df_high_confidence_epitope_metrics_RFU_100$pvalue <- NA
df_high_confidence_epitope_metrics_RFU_100$adjusted_pvalue <- NA
for (i in 1:nrow(df_high_confidence_epitope_metrics_RFU_100)){
  v_avg_RFU_pos <- subset(df_epitopes_RFU_100,(patient_ID%in%v_covid_pos_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_100$ORF[i]))$RFU
  v_avg_RFU_neg <- subset(df_epitopes_RFU_100,(patient_ID%in%v_covid_neg_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_100$ORF[i]))$RFU
  df_high_confidence_epitope_metrics_RFU_100$avg_RFU_in_positive_patients[i] <- mean(v_avg_RFU_pos,na.rm=T)
  df_high_confidence_epitope_metrics_RFU_100$avg_RFU_in_negative_patients[i] <- mean(v_avg_RFU_neg,na.rm=T)
  df_high_confidence_epitope_metrics_RFU_100$prevalence_in_positive_patients[i] <- length(unique(subset(df_epitopes_RFU_100,(patient_ID%in%v_covid_pos_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_100$ORF[i]))$patient_ID))
  df_high_confidence_epitope_metrics_RFU_100$prevalence_in_negative_patients[i] <- length(unique(subset(df_epitopes_RFU_100,(patient_ID%in%v_covid_neg_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_100$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_100$ORF[i]))$patient_ID))
  df_high_confidence_epitope_metrics_RFU_100$ratio[i] <- ifelse(test = df_high_confidence_epitope_metrics_RFU_100$avg_RFU_in_negative_patients[i]==0,yes = NA, no =df_high_confidence_epitope_metrics_RFU_100$avg_RFU_in_positive_patients[i]/df_high_confidence_epitope_metrics_RFU_100$avg_RFU_in_negative_patients[i])
  if (df_high_confidence_epitope_metrics_RFU_100$prevalence_in_negative_patients[i]==0 | df_high_confidence_epitope_metrics_RFU_100$prevalence_in_positive_patients[i]==0){
    next()
  }
  df_high_confidence_epitope_metrics_RFU_100$pvalue[i] <- ks.test(x=v_avg_RFU_pos,y=v_avg_RFU_neg)$p.value
}
df_high_confidence_epitope_metrics_RFU_100$adjusted_pvalue <- p.adjust(p=df_high_confidence_epitope_metrics_RFU_100$pvalue,method="fdr")
write.table(x=df_high_confidence_epitope_metrics_RFU_100[,c("ORF","pos_in_protein","avg_RFU_in_positive_patients","avg_RFU_in_negative_patients","prevalence_in_positive_patients","prevalence_in_negative_patients","ratio","pvalue","adjusted_pvalue")],file = paste0(output_workspace,"df_high_confidence_epitope_metrics_RFU_100.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

#Get p-values differential RFU (>=1000)
df_high_confidence_epitope_metrics_RFU_1000 <- df_high_confidence_epitope_metrics
df_high_confidence_epitope_metrics_RFU_1000$pvalue <- NA
df_high_confidence_epitope_metrics_RFU_1000$adjusted_pvalue <- NA
for (i in 1:nrow(df_high_confidence_epitope_metrics_RFU_1000)){
  v_avg_RFU_pos <- subset(df_epitopes,(patient_ID%in%v_covid_pos_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_1000$ORF[i]))$RFU
  v_avg_RFU_neg <- subset(df_epitopes,(patient_ID%in%v_covid_neg_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_1000$ORF[i]))$RFU
  df_high_confidence_epitope_metrics_RFU_1000$avg_RFU_in_positive_patients[i] <- mean(v_avg_RFU_pos,na.rm=T)
  df_high_confidence_epitope_metrics_RFU_1000$avg_RFU_in_negative_patients[i] <- mean(v_avg_RFU_neg,na.rm=T)
  df_high_confidence_epitope_metrics_RFU_1000$ratio[i] <- ifelse(test = df_high_confidence_epitope_metrics_RFU_1000$avg_RFU_in_negative_patients[i]==0,yes = NA, no =df_high_confidence_epitope_metrics_RFU_1000$avg_RFU_in_positive_patients[i]/df_high_confidence_epitope_metrics_RFU_1000$avg_RFU_in_negative_patients[i])
  df_high_confidence_epitope_metrics_RFU_1000$prevalence_in_positive_patients[i] <- length(unique(subset(df_epitopes,(patient_ID%in%v_covid_pos_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_1000$ORF[i]))$patient_ID))
  df_high_confidence_epitope_metrics_RFU_1000$prevalence_in_negative_patients[i] <- length(unique(subset(df_epitopes,(patient_ID%in%v_covid_neg_patients)&(protein_start<=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(protein_end>=df_high_confidence_epitope_metrics_RFU_1000$pos_in_protein[i])&(Mapping_region==df_high_confidence_epitope_metrics_RFU_1000$ORF[i]))$patient_ID))
  if (df_high_confidence_epitope_metrics_RFU_1000$prevalence_in_negative_patients[i]==0 | df_high_confidence_epitope_metrics_RFU_1000$prevalence_in_positive_patients[i]==0){
    next()
  }
  df_high_confidence_epitope_metrics_RFU_1000$pvalue[i] <- ks.test(x=v_avg_RFU_pos,y=v_avg_RFU_neg)$p.value
}
df_high_confidence_epitope_metrics_RFU_1000$adjusted_pvalue <- p.adjust(p=df_high_confidence_epitope_metrics_RFU_1000$pvalue,method="fdr")
write.table(x=df_high_confidence_epitope_metrics_RFU_1000[,c("ORF","pos_in_protein","avg_RFU_in_positive_patients","avg_RFU_in_negative_patients","prevalence_in_positive_patients","prevalence_in_negative_patients","ratio","pvalue","adjusted_pvalue")],file = paste0(output_workspace,"df_high_confidence_epitope_metrics_RFU_1000.csv"),sep = ",",na = "NA",row.names = FALSE,col.names = TRUE)

#save Rsession
library("session")
save.session(file = paste0(output_workspace,"high_confidence_epitope_analysis_RSession.Rda"))
