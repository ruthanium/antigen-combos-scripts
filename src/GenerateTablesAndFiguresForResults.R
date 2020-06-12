# *************************************
# Author: Ruth Dannenfelser
# Date: March 15, 2019
#
# Print out result tables and 
# generate top plots for easy viewing
#
#**************************************
library(data.table)
library(ggplot2)
library(cluster)
library(gridExtra)
library(grid)
setwd("/Genomics/ogtr04/rd6/BioGPS/")
setwd("/home/rd6/BioGPS/") # for lumos

# array data
arraymat = readRDS("scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")
cans = unique(arraymat[arraymat$type == "cancer",]$tissue.cancer)

# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-array-vs-normals-res-allTM.RData"))
  res = res[order(db),]
  res[,dbrank := seq(1, length(res$db))]
  res = res[order(-ch),]
  res[,chrank := seq(1, length(res$db))]
  res[,avgrank := apply(res[,c("dbrank", "chrank"), with=F], 1, mean)]
  res[,minrank := apply(res[,c("dbrank", "chrank"), with=F], 1, min)]
  res = res[order(minrank),]
  
  fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-array.txt", sep="\t")) 
  fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "array-top2000.txt", sep="\t")) 
  
  # # original - make plots for the top figures
  # # problem - often dominated by one single antigen
  # res = res[1:50,]
  # for (p in res$pair)
  # {
  #   ps = unlist(strsplit(p, ":"))
  #   genes = c(ps[1], ps[2])
  #   rm(ps)
  #   
  #   tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  #   tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
  #   tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
  #   names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
  #   genes[1] = gsub("-", "", genes[1])
  #   genes[2] = gsub("-", "", genes[2])
  #   
  #   gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
  #     xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + theme_minimal() + theme(legend.position="bottom")
  #   
  #   ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-array.png"), device="png")
  # }
  
  # new try and print out top 50 with unique potentials
  res = res[1:2000,]
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
      tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
        xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-array.png"), device="png")
      
    }
    
  }
}



# RNASeq data
rnaseqmat = readRDS("scripts/clustering-method/rnaseq-tm-clin.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "cancer",]$tissue.cancer)

# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-res-allTM.RData"))
  res = res[order(db),]
  res[,dbrank := seq(1, length(res$db))]
  res = res[order(-ch),]
  res[,chrank := seq(1, length(res$db))]
  res[,avgrank := apply(res[,c("dbrank", "chrank"), with=F], 1, mean)]
  res[,minrank := apply(res[,c("dbrank", "chrank"), with=F], 1, min)]
  res = res[order(minrank),]
  
  fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq.txt", sep="\t")) 
  fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-top2000.txt", sep="\t")) 
  
  # print out top 50 with unique potentials in the top 2000
  res = res[1:2000,]
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
      tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
        xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-rnaseq.png"), device="png")
      
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-array.png"), device="png")
        
      }
      
    }
    
  }
}
  
  
  
  
  
  
  
  
  
  
  
# ----------------------------------
# array data combined wtih RNAseq
# ----------------------------------
arraymat = readRDS("scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")
rnaseqmat = readRDS("scripts/clustering-method/rnaseq-tm-clin.RData")
rnaseqmat = rnaseqmat[rnaseqmat$tissue.cancer != "V7460",]
rnaseqmat = rnaseqmat[rnaseqmat$tissue.cancer != "V7057",]
cans = unique(arraymat[arraymat$type == "cancer",]$tissue.cancer)
cans2 = unique(rnaseqmat[rnaseqmat$type == "cancer",]$tissue.cancer)

rnaseq2array = fread("scripts/clustering-method/other-supporting-data/rnaseq-array-cancer-names.txt", sep="\t", header=T)

# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "array-top2000.txt", sep="\t"))
  if (can %in% rnaseq2array$array)
  {
    rres = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", rnaseq2array[rnaseq2array$array == can,]$rnaseq), "-rnaseq-top2000.txt", sep="\t"))
    
    # print out top 50 with unique potentials
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        
        if (genes[1] %in% names(rnaseqmat) & genes[2] %in% names(rnaseqmat))
        {
          tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
          tmp$target = ifelse(tmp$tissue.cancer == rnaseq2array[rnaseq2array$array == can,]$rnaseq, "target", "other")
          tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
          names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
          genes[1] = gsub("-", "", genes[1])
          genes[2] = gsub("-", "", genes[2])
          
          gg2 = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
            xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + 
            theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
          
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.pdf"))
          grid.arrange(gg, gg2, ncol=2)
          dev.off()
        }
        else
        {
          ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.png"), device="png")
        }
        
      }
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.png"), device="png")
        
      }
    }
  }
}

# would be nice to generate box plots for the distribution of samples as well
for (can in cans)
{
  print (can);
  res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "array-top2000.txt", sep="\t"))
  if (can %in% rnaseq2array$array)
  {
    rres = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", rnaseq2array[rnaseq2array$array == can,]$rnaseq), "-rnaseq-top2000.txt", sep="\t"))
    
    # print out top 50 with unique potentials
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
          geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("array UPC") + labs(title=genes[1]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
          geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("array UPC") + labs(title=genes[2]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        
        if (genes[1] %in% names(rnaseqmat) & genes[2] %in% names(rnaseqmat))
        {
          tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
          tmp$target = ifelse(tmp$tissue.cancer == rnaseq2array[rnaseq2array$array == can,]$rnaseq, "target", "other")
          tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
          names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
          genes[1] = gsub("-", "", genes[1])
          genes[2] = gsub("-", "", genes[2])
          
          gg2a = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
            geom_boxplot() + theme_minimal()  + xlab("") + coord_flip() + ylab("RNAseq UPC") + labs(title=genes[1]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
          gg2b = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
            geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("RNAseq UPC") + labs(title=genes[2]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
          
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
          grid.arrange(gga, ggb, ncol=2)
          dev.off()
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box2.pdf"))
          grid.arrange(gg2a, gg2b, ncol=2)
          dev.off()
        }
        else
        {
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
          grid.arrange(gga, ggb, ncol=2)
          dev.off()
        }
        
      }
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg2a = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
          geom_boxplot() + theme_minimal() + coord_flip() + xlab("") + ylab("array UPC") + labs(title=genes[1]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        gg2b = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
          geom_boxplot() + theme_minimal() + coord_flip() + xlab("") + ylab("array UPC") + labs(title=genes[2]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
        grid.arrange(gga, ggb, ncol=2)
        dev.off()
        
      }
    }
  }
}


















# ----------------------------------
# tissueGPS
# ----------------------------------
rnaseqmat = readRDS("scripts/clustering-method/rnaseq-tm-clin.RData")
rnaseqmat = rnaseqmat[rnaseqmat$tissue.cancer != "V7460",]
rnaseqmat = rnaseqmat[rnaseqmat$tissue.cancer != "V7057",]
rnaseqmat = rnaseqmat[rnaseqmat$type == "tissue",]
tissues = unique(rnaseqmat$tissue.cancer)

# for each cancer print the output to a table
for (tis in tissues)
{
  print (tis);
  
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", tis), "-tissueGPS-res-allTM.RData"))
  res = res[order(db),]
  res[,dbrank := seq(1, length(res$db))]
  res = res[order(-ch),]
  res[,chrank := seq(1, length(res$db))]
  res[,avgrank := apply(res[,c("dbrank", "chrank"), with=F], 1, mean)]
  res[,minrank := apply(res[,c("dbrank", "chrank"), with=F], 1, min)]
  res = res[order(minrank),]
  
  fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", tis), "-gtex.txt", sep="\t")) 
  fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", tis), "-gtex-top2000.txt", sep="\t")) 
  
  res = res[1:2000]
  
  # print out top 50 with unique potentials
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {          
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == tis, "target", "other")
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gmain = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) +
        xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + #ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("GTEx UPC") + ylim(0,1) + labs(title=genes[1]) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("GTEx UPC") + ylim(0,1) + labs(title=genes[2]) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      lay <- rbind(c(1,1,2,3),c(1,1,2,3))
      oGG2 <- arrangeGrob(gmain, gga, ggb, layout_matrix=lay, respect=T, top=textGrob(paste0(tis, "-", genes[1], ":", genes[2]), gp=gpar(fontsize=16)));
      plot(oGG2)
      ggsave(file=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", tis), "-", rank, "-", genes[1], "-", genes[2], "-tissueGPS-box.pdf"), oGG2, width = 10, height=5, units = "in")
      
    }
  }
}

# ----------------------------------
# combat TPM rnaseq 
# ----------------------------------
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)

ending = "-rnaseq-combat"

# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), ending, ".RData"))
  res = res[order(db),]
  res[,dbrank := seq(1, length(res$db))]
  res = res[order(-ch),]
  res[,chrank := seq(1, length(res$db))]
  res[,avgrank := apply(res[,c("dbrank", "chrank"), with=F], 1, mean)]
  res[,minrank := apply(res[,c("dbrank", "chrank"), with=F], 1, min)]
  res = res[order(minrank),]

  fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), ending, ".txt"), sep="\t")
  fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), ending, "-top2000.txt"), sep="\t")

  res = res[1:2000]
  
  # shortcut when tables are already made
  #res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-combat-top2000.txt"))
  
  # print out top 50 with unique potentials
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {          
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == can, "cancer", "normal")
      tmp = tmp[tmp$type=="tissue" | (tmp$type=="tumor" & tmp$tissue.cancer== can),]
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gmain = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) +
        scale_color_manual(values=c("#e74c3c", "#34495e", "#bdc3c7")) + #ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)")  + labs(title=genes[1]) +
        scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)") + labs(title=genes[2]) +
        scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      lay <- rbind(c(1,1,2,3),c(1,1,2,3))
      oGG2 <- arrangeGrob(gmain, gga, ggb, layout_matrix=lay, respect=T, top=textGrob(paste0(can, "-", genes[1], ":", genes[2]), gp=gpar(fontsize=16)));
      plot(oGG2)
      ggsave(file=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], ending, "-box.pdf"), oGG2, width = 10, height=5, units = "in")
      
    }
  }
}



# ----------------------------------
# array data combined wtih combat 
# corrected RNAseq runs
# ----------------------------------
arraymat = readRDS("scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin.RData")
cans = unique(arraymat[arraymat$type == "cancer",]$tissue.cancer)
cans2 = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)

rnaseq2array = fread("scripts/clustering-method/other-supporting-data/rnaseq-array-cancer-names.txt", sep="\t", header=T)



# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "array-top2000.txt", sep="\t"))
  if (can %in% rnaseq2array$array)
  {
    rres = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", rnaseq2array[rnaseq2array$array == can,]$rnaseq), "-rnaseq-top2000.txt", sep="\t"))
    
    # print out top 50 with unique potentials
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        
        if (genes[1] %in% names(rnaseqmat) & genes[2] %in% names(rnaseqmat))
        {
          tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
          tmp$target = ifelse(tmp$tissue.cancer == rnaseq2array[rnaseq2array$array == can,]$rnaseq, "target", "other")
          tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
          names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
          genes[1] = gsub("-", "", genes[1])
          genes[2] = gsub("-", "", genes[2])
          
          gg2 = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
            xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + 
            theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
          
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.pdf"))
          grid.arrange(gg, gg2, ncol=2)
          dev.off()
        }
        else
        {
          ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.png"), device="png")
        }
        
      }
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array.png"), device="png")
        
      }
    }
  }
}

# would be nice to generate box plots for the distribution of samples as well
for (can in cans)
{
  print (can);
  res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "array-top2000.txt", sep="\t"))
  if (can %in% rnaseq2array$array)
  {
    rres = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", rnaseq2array[rnaseq2array$array == can,]$rnaseq), "-rnaseq-top2000.txt", sep="\t"))
    
    # print out top 50 with unique potentials
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
          geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("array UPC") + labs(title=genes[1]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
          geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("array UPC") + labs(title=genes[2]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        
        if (genes[1] %in% names(rnaseqmat) & genes[2] %in% names(rnaseqmat))
        {
          tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
          tmp$target = ifelse(tmp$tissue.cancer == rnaseq2array[rnaseq2array$array == can,]$rnaseq, "target", "other")
          tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
          names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
          genes[1] = gsub("-", "", genes[1])
          genes[2] = gsub("-", "", genes[2])
          
          gg2a = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
            geom_boxplot() + theme_minimal()  + xlab("") + coord_flip() + ylab("RNAseq UPC") + labs(title=genes[1]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
          gg2b = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
            geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("RNAseq UPC") + labs(title=genes[2]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
          
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
          grid.arrange(gga, ggb, ncol=2)
          dev.off()
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box2.pdf"))
          grid.arrange(gg2a, gg2b, ncol=2)
          dev.off()
        }
        else
        {
          pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
          grid.arrange(gga, ggb, ncol=2)
          dev.off()
        }
        
      }
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg2a = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
          geom_boxplot() + theme_minimal() + coord_flip() + xlab("") + ylab("array UPC") + labs(title=genes[1]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        gg2b = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
          geom_boxplot() + theme_minimal() + coord_flip() + xlab("") + ylab("array UPC") + labs(title=genes[2]) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        pdf(paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-array-box.pdf"))
        grid.arrange(gga, ggb, ncol=2)
        dev.off()
        
      }
    }
  }
}




# ----------------------------------
# combat TPM rnaseq - test with
# silhouette plots
# ----------------------------------
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)
sil = fread("scripts/clustering-method/uveal-melanoma-silhouette.txt", sep="\t")
can = c("Uveal Melanoma")
names(sil) = c("pair", "sil")
# for each cancer print the output to a table
res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat.RData"))
res = merge(res, sil, by="pair")

res = res[order(db),]
res[,dbrank := seq(1, length(res$db))]
res[,':='(logdb=log(db), logch=log(ch))]
minlogdb = min(res$logdb)
maxlogdb = max(res$logdb)
res = res[,dbscaled := 1-((logdb - minlogdb) / (maxlogdb-minlogdb)) ]
minlogch = min(res$logch)
maxlogch = max(res$logch)
res = res[,chscaled := ((logch - minlogch) / (maxlogch-minlogch)) ]
minsil = min(res$sil)
maxsil = max(res$sil)
res = res[,silscaled := ((sil - minsil) / (maxsil-minsil)) ]

res = res[order(-ch),]
res[,chrank := seq(1, length(res$db))]
res = res[order(-sil),]
res[,silrank := seq(1, length(res$sil))]
res[,avgrank := apply(res[,c("dbrank", "chrank", "sil"), with=F], 1, mean)]

res[,minrank := apply(res[,c("dbrank", "chrank"), with=F], 1, min)]
res = res[order(avgrank),]

fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-combat.txt"), sep="\t")
fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-combat-top2000.txt"), sep="\t")

res = res[1:2000]

# shortcut when tables are already made
res = fread(paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-combat-top2000.txt"))

# print out top 50 with unique potentials
upperlim = 50;
curr = 0;
rank = 0;
one = character()
two = character()
for (p in res$pair)
{
  rank = rank + 1
  if (curr == upperlim)
  {
    break;
  }
  
  ps = unlist(strsplit(p, ":"))
  genes = c(ps[1], ps[2])
  rm(ps)
  if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
  {
    curr = curr + 1;
    if (genes[1] %in% one)
    {
      two = append(two, genes[1])
      one = one[!one %in% c(genes[1])]
    }
    else
    {          
      one = append(one, genes[1])
    }
    if (genes[2] %in% one)
    {
      two = append(two, genes[2])
      one = one[!one %in% c(genes[2])]
    }
    else
    {
      one = append(one, genes[2])
    }
    
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == can, "cancer", "normal")
    tmp = tmp[tmp$type=="tissue" | (tmp$type=="tumor" & tmp$tissue.cancer== can),]
    names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
    genes[1] = gsub("-", "", genes[1])
    genes[2] = gsub("-", "", genes[2])
    
    gmain = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) +
      scale_color_manual(values=c("#e74c3c", "#34495e", "#bdc3c7")) + #ggtitle(paste0("rank ", rank)) + 
      theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    
    gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
      geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)")  + labs(title=genes[1]) +
      scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
      geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)") + labs(title=genes[2]) +
      scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    lay <- rbind(c(1,1,2,3),c(1,1,2,3))
    oGG2 <- arrangeGrob(gmain, gga, ggb, layout_matrix=lay, respect=T, top=textGrob(paste0(can, "-", genes[1], ":", genes[2]), gp=gpar(fontsize=16)));
    plot(oGG2)
    ggsave(file=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], "-rnaseq-combat-box.pdf"), oGG2, width = 10, height=5, units = "in")
    
  }
}






# ----------------------------------
# combat TPM rnaseq - distance and DB
# clinicals only
# ----------------------------------
clin.dists = readRDS("scripts/clustering-method/rdata/singles-distances.RData")
# test on one cancer
can = "Colon Adenocarcinoma"
res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-clin.RData"))
res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
res = res[order(-score),]

# plot multiple of the top perfomers
p <- list()
for(i in 1:9){
  
  genes = unlist(strsplit(as.character(res$pair[i]), ":"))
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
  tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
  names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
  genes[1] = gsub("-", "", genes[1])
  genes[2] = gsub("-", "", genes[2])
  
  p[[i]] <- ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + #facet_wrap(~target) +
    scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
}
do.call(grid.arrange,p)

# to do
# combine with boxplot and export for browsing

# best clinical, best clinical-clinical, best clinical-novel across cancers

# want to show the distributions of the top 100 across cancers
first = 1
bigres = data.table()
for (can in cans)
{
  print(can)
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-clin.RData"))
  res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
  res = res[order(-score),]
  res[, type := can]
  if (first == 1)
  {
    bigres = res
    first = 2 
  }
  else
  {
    bigres = rbind(bigres, res)
  }
}
saveRDS(bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin.RData")
#bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin.RData")

# do the same combination for all of the pairs - combine into one giant pile
first = 1
bigres = data.table()
for (can in cans)
{
  print(can)
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2.RData"))
  res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
  res = res[order(-score),]
  res[, type := can]
  if (first == 1)
  {
    bigres = res
    first = 2 
  }
  else
  {
    bigres = rbind(bigres, res)
  }
}
saveRDS(bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2.RData")

# do the same for triples now - one giant pile
first = 1
restriples = data.table()
for (can in cans)
{
  print(can)
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-clin-triples.RData"))
  res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
  res = res[order(-score),]
  res[, type := can]
  if (first == 1)
  {
    restriples = res
    first = 2 
  }
  else
  {
    restriples = rbind(restriples, res)
  }
}
saveRDS(restriples, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin-triples.RData")




# do the same for triples now - one giant pile (latest set of triples created in Dec 2019/Jan2020)
first = 1
restriples = data.table()
for (can in cans)
{
  print(can)
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-triples-S-reduction.RData"))
  #res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  #res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  #res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
  #res = res[order(-score),]
  res[, type := can]
  if (first == 1)
  {
    restriples = res
    first = 2 
  }
  else
  {
    restriples = rbind(restriples, res)
  }
}
#saveRDS(restriples, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-triples-S-reduction.RData")
# scale across cancers and not per tumor type
restriples[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
restriples[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
restriples[,score := apply(restriples[,c("dbscaled", "distscaled"), with=F], 1, min)]
restriples = restriples[order(-score),]

saveRDS(restriples, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-triples-S-reduction.RData", version = 2) # make backwards compatible






# add annotations to the bigres of what type of pairs
clinSmall = fread("data/clin-greg.txt")
clinSmall = clinSmall$V2

tisSpec = fread("data/tissue-specific-surface-markers.txt", header=F)
tisSpec = tisSpec$V1

fib = fread("data/fibroblast-surface-markers.txt", header=F)
fib = fib$V1

bigres[, c("gene1", "gene2") := tstrsplit(pair, ":", fixed=T)]
bigres[, c("g1lab", "g2lab") := "N"]

# order of preference - Tissue specific, then clinical, then fibroblast - ordered b/c
# some might override other categories
bigres[, "g1lab" := ifelse(gene1 %in% tisSpec, "T", g1lab)]
bigres[, "g2lab" := ifelse(gene2 %in% tisSpec, "T", g2lab)]

bigres[, "g1lab" := ifelse(gene1 %in% clinSmall, "C", g1lab)]
bigres[, "g2lab" := ifelse(gene2 %in% clinSmall, "C", g2lab)]

bigres[, "g1lab" := ifelse(gene1 %in% fib, "F", g1lab)]
bigres[, "g2lab" := ifelse(gene2 %in% fib, "F", g2lab)]

bigres[, "ttype" := paste0(g1lab, ":", g2lab) ]
bigres[, "ttype" := ifelse(ttype == "N:T", "T:N", ttype)]
bigres[, "ttype" := ifelse(ttype == "N:F", "F:N", ttype)]
bigres[, "ttype" := ifelse(ttype == "N:C", "C:N", ttype)]

bigres[, "ttype" := ifelse(ttype == "F:C", "C:F", ttype)]
bigres[, "ttype" := ifelse(ttype == "T:C", "C:T", ttype)]
bigres[, "ttype" := ifelse(ttype == "F:T", "T:F", ttype)]

#saveRDS(bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2.RData")

# plot showing the total distance distribution of pairs across cancer types
ggplot(bigres, aes(x=type, y=dist.man)) + geom_violin() + theme_minimal() + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()

# same plot as above but showing for the 26 clinicals across cancer types
ggplot(clin.dists, aes(x=can, y=dist)) + geom_violin() + theme_minimal() + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + ylim(0, 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()



# ----------------------------------
# combat TPM rnaseq - DB + dist
# ----------------------------------
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "cancer",]$tissue.cancer)

ending = "-rnaseq-combat-sketch-0.2"

# for each cancer print the output to a table
for (can in cans)
{
  print (can);
  
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), ending, ".RData"))
  res[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  res[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, min)]
  res = res[order(-score),]
  
  fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/sketch02/", gsub(" ", "-", can), ending, ".txt"), sep="\t")
  fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/sketch02/", gsub(" ", "-", can), ending, "-top2000.txt"), sep="\t")
  
  res = res[1:2000]
  
  # shortcut when tables are already made
  #res = fread(paste0("scripts/clustering-method/res-tables/allTM/sketch02/", gsub(" ", "-", can), "-rnaseq-combat-top2000.txt"))
  
  # print out top 50 with unique potentials
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {          
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == can, "cancer", "normal")
      tmp = tmp[tmp$type=="tissue" | (tmp$type=="tumor" & tmp$tissue.cancer== can),]
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gmain = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) +
        scale_color_manual(values=c("#e74c3c", "#34495e", "#bdc3c7")) + #ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      gga = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[1])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)")  + labs(title=genes[1]) +
        scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      ggb = ggplot(tmp, aes_string(x=paste0("reorder(tissue.cancer, ",genes[1],", FUN=mean)"), y= genes[2])) + 
        geom_boxplot() + theme_minimal() + xlab("") + coord_flip() + ylab("log(TPM)") + labs(title=genes[2]) +
        scale_x_discrete(label=function(x) abbreviate(x, minlength=7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      lay <- rbind(c(1,1,2,3),c(1,1,2,3))
      oGG2 <- arrangeGrob(gmain, gga, ggb, layout_matrix=lay, respect=T, top=textGrob(paste0(can, "-", genes[1], ":", genes[2]), gp=gpar(fontsize=16)));
      plot(oGG2)
      ggsave(file=paste0("scripts/clustering-method/res-figures/allTM/sketch02/", gsub(" ", "-", can), "-", rank, "-", genes[1], "-", genes[2], ending, "-box.pdf"), oGG2, width = 10, height=5, units = "in")
      
    }
  }
}







# TO DO - kmeans default library, feed clusters
# iterate once inside no outside to get cluster assignments
# for top 5 pairings with each clinical and bottom 5 pairings with each clinical


# calculate the difference between the singles and the doubles
# get the top pair for each clinical marker
res[, c("gene1", "gene2") := tstrsplit(pair, ":", fixed=TRUE)]

pairs = character()
dbs = numeric()
dists = numeric()
for (gene in unique(clin.dists$gene))
{
  tmpp = res[gene1 == gene | gene2 == gene, ] 
  pairs = append(pairs, as.character(tmpp[1,]$pair))
  dbs = append(dbs, tmpp[1,]$db)
  dists = append(dists, tmpp[1,]$dist.man)
}
topPerClin = data.table(pair=pairs, db=dbs, dist.man=dists)
topPerClin[, c("gene1", "gene2") := tstrsplit(pair, ":", fixed=TRUE)]
topPerClin[, clin := ifelse(topPerClin$gene1 %in% unique(clin.dists$gene), gene1, gene2)]
tcan = can
cd = clin.dists[clin.dists$can == tcan,]
rm(tcan)
topPerClin = merge(topPerClin, cd, by.x="clin", by.y="gene")
edist = numeric()
for (p in topPerClin$pair)
{
  tmp = rnaseqmat[rnaseqmat$type == "tissue" | rnaseqmat$tissue.cancer == can,]
  genes = unlist(strsplit(p, ":"))
  tmp = tmp[,c(genes[1], genes[2], "tissue.cancer"), with=F]
  tmp[, "target" := ifelse(tmp$tissue.cancer == can, "target", "other")]
  
  cm = colMeans(tmp[target == "other",1:2, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="euclidian")[1]
  edist = append(edist, cent)
}
topPerClin$eDist = edist 

#arraymat = readRDS("scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin.RData")
#cans = unique(arraymat[arraymat$type == "cancer",]$tissue.cancer)
cans = unique(rnaseqmat[rnaseqmat$type == "cancer",]$tissue.cancer)

rnaseq2array = fread("scripts/clustering-method/other-supporting-data/rnaseq-array-cancer-names.txt", sep="\t", header=T)

for (can in cans)
{
  print (can);
  res = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-clin.RData"))
  res[,dbscaled := (log(db) - min(log(db)))/(max(log(db)) - min(log(db)))] # make score between 0 and 1
  res[,distscaled := (dist.man - min(dist.man))/(max(dist.man) - min(dist.man))]
  res[,score := apply(res[,c("dbscaled", "distscaled"), with=F], 1, mean)]
  res = res[order(db),]
  
  #fwrite(res, paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq.txt", sep="\t")) 
  #fwrite(res[1:2000], paste0("scripts/clustering-method/res-tables/allTM/", gsub(" ", "-", can), "-rnaseq-top2000.txt", sep="\t")) 
  
  # print out top 50 with unique potentials in the top 2000
  res = res[1:2000,]
  upperlim = 50;
  curr = 0;
  rank = 0;
  one = character()
  two = character()
  for (p in res$pair)
  {
    rank = rank + 1
    if (curr == upperlim)
    {
      break;
    }
    
    ps = unlist(strsplit(p, ":"))
    genes = c(ps[1], ps[2])
    rm(ps)
    if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
    {
      curr = curr + 1;
      if (genes[1] %in% one)
      {
        two = append(two, genes[1])
        one = one[!one %in% c(genes[1])]
      }
      else
      {
        one = append(one, genes[1])
      }
      if (genes[2] %in% one)
      {
        two = append(two, genes[2])
        one = one[!one %in% c(genes[2])]
      }
      else
      {
        one = append(one, genes[2])
      }
      
      tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
      tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
      tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
      names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
      genes[1] = gsub("-", "", genes[1])
      genes[2] = gsub("-", "", genes[2])
      
      gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
        xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
        theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
      
      ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-rnaseq.png"), device="png")
      
    }
  }
  else
  {
    upperlim = 50;
    curr = 0;
    rank = 0;
    one = character()
    two = character()
    for (p in res$pair)
    {
      rank = rank + 1
      if (curr == upperlim)
      {
        break;
      }
      
      ps = unlist(strsplit(p, ":"))
      genes = c(ps[1], ps[2])
      rm(ps)
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        curr = curr + 1;
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        tmp = as.data.frame(arraymat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
        tmp$target = ifelse(tmp$tissue.cancer == can, "target", "other")
        tmp = tmp[tmp$target == "target" | (tmp$type == "tissue"),]
        names(tmp) = gsub("-", "", names(tmp)) # rare case where gene name contains a hypen (e.g. ERVW-1)
        genes[1] = gsub("-", "", genes[1])
        genes[2] = gsub("-", "", genes[2])
        
        gg = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + facet_wrap(~target) +
          xlim(0,1) + ylim(0,1) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + ggtitle(paste0("rank ", rank)) + 
          theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
        
        ggsave(gg,filename=paste0("scripts/clustering-method/res-figures/allTM/", gsub(" ", "-", can), "-", genes[1], "-", genes[2], "-array.png"), device="png")
        
      }
      
    }
    
  }
}