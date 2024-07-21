library(tidyr)
library(RNOmni)

# Major data input
trait_pid <- read.csv("trait_rsid_phecode.csv") # phenotype and target common loci info
trait_qtl_gene <- read.csv("trait_1_ukb22418_qtl_gene_signals.csv") # associated gene and QTL info
trait_value <- read.csv("trait_phenotype.csv") # trait value 
cov22 <- read.csv("cov.csv") # covariates
common_gt_all <- read.csv("trait_common_gt.csv",sep = "\t") # cRV genotype
ukb_wes_sid <- read.csv("ukb23157_sid.csv",header = F) # UKB WES eid 

# Merge tables above
pheindex <- separate(trait_pid,"phecode",into = c("phecode","irnt"),sep = "_")
trait_qtl_gene <- unique(trait_qtl_gene[,c(3,6)])
pheindex <- merge(pheindex,trait_qtl_gene,by="ID")
all_phecode <- unique(trait_pid$phecode) 
all_phecode_index <- paste0("p",all_phecode,"_i0")
trait_value_1 <- trait_value[,c("eid",all_phecode_index)]
phenotype <- merge(trait_value_1,cov22,by='eid')

# Function for exchanging common variant genotype (0 & 2) 
# Exchange when non-carrier count < 100,000
exchange_02 <- function(col_name){
  if (('0' %in% col_name)&(table(col_name)['0']<100000)){
    col_name[col_name==2] <- 3
    col_name[col_name==0] <- 2
    col_name[col_name==3] <- 0
  } else if ('0' %in% col_name == FALSE) {
    col_name[col_name==2] <- 0
  }
  return(col_name)
}

# Variant interaction test for one trait
# loop for all cRVs associated with the trait
for (i in seq(ncol(common_gt_all)-1)) {
  common_single_gt <- as.data.frame(cbind(common_gt_all[,1],exchange_02(common_gt_all[,i+1])))
  colnames(common_single_gt) <- c(colnames(common_gt_all)[1],colnames(common_gt_all)[i+1])
  rsid <- colnames(common_gt_all)[i+1]
  cv_index <- pheindex[which(pheindex$ID==rsid),]
  # loop for all variant pairs
  for (n in seq(nrow(cv_index))) {
    # normalize trait value
    phenodes <- cv_index[n,2]
    fieldid <- paste0("p",as.character(phenodes),"_i0")
    phe_data <- phenotype[,c(1,which(names(phenotype) == fieldid),58:79)]
    phe_model <- phe_data[complete.cases(phe_data),]
    phe_model[,2] <- RankNorm(phe_model[,2])
    # input PV genotype
    rare_gene <- cv_index[n,3]
    lof_url <- paste0("./gt/",rare_gene,"_rare_lof_gt.tsv")
    mis_url <- paste0("./gt/",rare_gene,"_rare_mis_gt.tsv")
    if (file.exists(lof_url)){
      rare_gt <- read.delim(lof_url)
      rare_gt['sid'] <- ukb_wes_sid
      gt_lof <- merge(common_single_gt,rare_gt,by='sid')
      phe_lof <- merge(phe_model,gt_lof,by='sid')
      # filter less than 5 carriers with double variants
      remain_index <- c()
      for (m in seq(ncol(gt_lof)-2)) {
        dat1 <- phe_lof[which(phe_lof[,25]>0),]
        num <- nrow(dat1[which(dat1[,m+25]>0),])
        num2 <- nrow(phe_lof[which((phe_lof[,25]==0)&(phe_lof[,m+25]>0)),])
        if (num>4 & num2>4 & num<nrow(dat1)) {
          remain_index <- append(remain_index,m+25)
        }
      }
      model_lof <- phe_lof[,c(c(1:25),remain_index)]
      model_lof[,25] <- as.numeric(model_lof[,25])
      # modeling for LoF PV
      if (ncol(model_lof)>25) {
        sample_size <- c()
        model_lof_res <- data.frame()
        # loop for all LoF
        for (j in seq(ncol(model_lof)-25)) {
          n_common = nrow(model_lof[which(model_lof[,25]>0),])
          n_rare = nrow(model_lof[which(model_lof[,25+j]>0),])
          n_inter = nrow(model_lof[which((model_lof[,25+j]>0)&(model_lof[,25]>0)),])
          sam <- c(n_common,n_rare,n_inter)
          sample_size <- append(sample_size,sam)
          gmodel <- lm(model_lof[,2]~model_lof[,25]*model_lof[,25+j]+model_lof[,25]+model_lof[,25+j]+p31+p21022+p31*p21022*p21022+p31*p31*p21022+
                          p22009_a1+p22009_a2+p22009_a3+p22009_a4+p22009_a5+p22009_a6+p22009_a7+p22009_a8+p22009_a9+p22009_a10+
                          PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data = model_lof)
          coef_table <- coef(summary(gmodel))
          coef_table <- coef_table[c(2,3,26),]
          row.names(coef_table)<-c(paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"common",sep = '_'),
                                   paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"rare",sep = '_'),
                                   paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"interaction",sep = '_'))
          model_lof_res <- rbind(model_lof_res,coef_table)
        }
        model_lof_res['sample_size'] <- sample_size
        model_lof_res['variant_type'] <- "LoF"
        model_lof_res['phenotype'] <- phenodes
        url_lof_res <- paste0("./res/",rsid,"_",rare_gene,"_",phenodes,"_lof.csv")
        write.csv(model_lof_res,url_lof_res)
      } else {print("no suitable variant")}
    }
    if (file.exists(mis_url)){
      rare_gt <- read.delim(mis_url)
      rare_gt['sid'] <- ukb_wes_sid
      gt_lof <- merge(common_single_gt,rare_gt,by='sid')
      phe_lof <- merge(phe_model,gt_lof,by='sid')
      # filter less than 5 carriers with double variants
      remain_index <- c()
      for (m in seq(ncol(gt_lof)-2)) {
        dat1 <- phe_lof[which(phe_lof[,25]>0),]
        num <- nrow(dat1[which(dat1[,m+25]>0),])
        num2 <- nrow(phe_lof[which((phe_lof[,25]==0)&(phe_lof[,m+25]>0)),])
        if (num>4 & num2>4 & num<nrow(dat1)) {
          remain_index <- append(remain_index,m+25)
        }
      }
      model_lof <- phe_lof[,c(c(1:25),remain_index)]
      model_lof[,25] <- as.numeric(model_lof[,25])
      # modeling for missense PV
      if (ncol(model_lof)>25) {
        sample_size <- c()
        model_lof_res <- data.frame()
        # loop for all Mis
        for (j in seq(ncol(model_lof)-25)) {
          n_common = nrow(model_lof[which(model_lof[,25]>0),])
          n_rare = nrow(model_lof[which(model_lof[,25+j]>0),])
          n_inter = nrow(model_lof[which((model_lof[,25+j]>0)&(model_lof[,25]>0)),])
          sam <- c(n_common,n_rare,n_inter)
          sample_size <- append(sample_size,sam)
          gmodel <- lm(model_lof[,2]~model_lof[,25]*model_lof[,25+j]+model_lof[,25]+model_lof[,25+j]+p31+p21022+p31*p21022*p21022+p31*p31*p21022+
                          p22009_a1+p22009_a2+p22009_a3+p22009_a4+p22009_a5+p22009_a6+p22009_a7+p22009_a8+p22009_a9+p22009_a10+
                          PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data = model_lof)
          coef_table <- coef(summary(gmodel))
          coef_table <- coef_table[c(2,3,26),]
          row.names(coef_table)<-c(paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"common",sep = '_'),
                                   paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"rare",sep = '_'),
                                   paste(colnames(model_lof)[25],colnames(model_lof)[25+j],"interaction",sep = '_'))
          model_lof_res <- rbind(model_lof_res,coef_table)
        }
        model_lof_res['sample_size'] <- sample_size
        model_lof_res['variant_type'] <- "Mis"
        model_lof_res['phenotype'] <- phenodes
        url_lof_res <- paste0("./res/",rsid,"_",rare_gene,"_",phenodes,"_mis.csv")
        write.csv(model_lof_res,url_lof_res)
      } else {print("no suitable variant")}
    }
  }
}

