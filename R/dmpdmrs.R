### TODO:
# Make vulcano plot for dmrs/dmps parameters (better visualize where to set the cut points) Remeber to use shiny app vulcano plot.



library(ggplot2)
transparent_theme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))


### DMPS:

summary_dmps <- function(DMPextr_object,dir="./results/dmps/",name="raw",write=F){
  require(data.table)
  dt<-data.table::as.data.table(DMPextr_object)
  dt_summary<-dt[,.(
    "Hyper" = sum(Type=="Hyper"),Hypo=sum(Type=="Hypo")
    # "min/max" = paste(range(round(diff_meanMeth,2)),collapse=" / "),
    # "p<0.05" = sum(adj.P.Val<0.05),
    # "min5" = paste(head(round(sort(diff_meanMeth),2),5),collapse=";"),
    # "max5" = paste(tail(round(sort(diff_meanMeth),2),5),collapse=";"),
    
  ),by=c("Contrast")
  ]
  if(write){
    data.table::fwrite(dt,paste0(dir,name,"_dmp_raw.csv.gz"))
    data.table::fwrite(dt_summary,paste0(dir,name,"_dmp_summary.txt"))
  }
  
  return(dt_summary)
}


filter_dmps <- function(dmps, p.value = 0.01, mDiff = 0.3, s=F){
  require(data.table)
  dmps<-data.table::as.data.table(dmps)
  dmps_f <- dmps[adj.P.Val <= p.value & 
                abs(diff_meanMeth) >= mDiff, ]
  dmps_s<-summary_dmps(dmps_f)
  
  if(s==T){out<-dmps_s[,.SD,.SDcols=c("Hypo","Hyper","Contrast")]}else{out<-dmps_f}
  return(out)
}

apply_filter_dmps<-function(dmps,dev = "png",p.value=seq(0.00,0.1,.01),
                            mDiff=seq(0.15,0.5,.05),path="analysis/intermediate/dmps/"){
  require(ggplot2)
  dir.create(path)
  params<-expand.grid(p.value,mDiff,T)
  # names(params)<-c("p.val","mDiff","s")
  # params|>purrr::map(function(x)filter_dmps(dmps,x[1],x[2],x[3]))
  # p2<-split(params, seq(nrow(params)))
  # res2<-p2 |> purrr::map(\(x)filter_dmps(dmps,unlist(x[1]),unlist(x[2]),unlist(x[3])))
  
  res1<-with(params, Map(function(a,b,c) {
    dt<-filter_dmps(dmps, a,b,c)
    dt$p.val<-a
    dt$mDiff<-b
    dt
  }, Var1, Var2, Var3))
  pdata<-rbindlist(res1)
  pdata[,All:=Hypo+Hyper]
  pd<-melt(pdata,measure.vars=c("Hyper","Hypo","All"))
  
  pd$variable<-factor(pd$variable)
  if(NROW(pd)>0 & length(levels(pd$variable))>1){
    plt_list<-list()
    plt_list[["dmp_params"]] <-ggplot2::ggplot(data=pd,aes(x=mDiff,y=value,group=p.val,color=p.val))+#,color=Contrast,group=varaible))+
      geom_line(aes(linetype=factor(p.val)))+
      facet_grid(variable~Contrast)#,margin="variable")
    lapply(1:length(plt_list),function(x)
      ggsave(plot = plt_list[[x]],
             filename = paste0(names(plt_list[x]),".",dev),
             path = path,
             device = dev))
  }
}

#DMRs
summary_dmrs <- function(dmrs,path="/results/dmrs/",write=T){
  dmrs[,Type:=ifelse(meandiff>0,"Hyper","Hypo")]
  dmrs.l<-dmrs[,list(Hyper.DMRS=sum(Type=="Hyper"),Hypo.DMRS=sum(Type=="Hypo")),by=c("Contrast")]
  genes.l<-dmrs[,list(Hyper.Genes=length(unique(unlist(strsplit(overlapping.genes[Type=="Hyper"],",")))),Hypo.Genes=length(unique(unlist(strsplit(overlapping.genes[Type=="Hypo"],","))))),by=c("Contrast")]
  summary<-merge(dmrs.l,genes.l)
  if(write)data.table::fwrite(summary,path)
  return(summary)
}

filter_dmrs <- function(dmrs, p.value = 0.01, mDiff = 0.2, min.cpg=5,s=F){
  require(data.table)             
  dmrs<-data.table::as.data.table(dmrs)
  out <- dmrs[HMFDR <= p.value & abs(meandiff) >= mDiff & no.cpgs >= min.cpg,]
  if(s)out<-summary_dmrs(out,write=F)
  return(out)
}

apply_filter_dmrs<-function(dmrs,plots=T,dev="png",p.value=seq(0.001,0.11,.025),
                            mDiff=seq(0.15,0.5,.05),min.cpg=seq(3,6,1),path="analysis/intermediate/dmrs/"){
  require(ggplot2)
  dir.create(path)
  params<-expand.grid(p.value,mDiff,min.cpg,s=T)
 
  
  res1<-with(params, Map(function(a,b,c,d) {
    dt<-filter_dmrs(dmrs, a,b,c,d)
    dt$p.val<-a
    dt$mDiff<-b
    dt$min.cpg<-c
    dt
  }, Var1, Var2, Var3,s))
  pdata<-rbindlist(res1)
  colA<-names(pdata)[endsWith(names(pdata),"DMRS")]
  colB<-names(pdata)[endsWith(names(pdata),"Genes")]
  pd<-melt(pdata,measure.vars=list(colA,colB),value.name = c("DMRS","Genes"),
           variable.name = "type",variable.factor = T)
  levels(pd$type)<-c("Hyper","Hypo")
  if(plots){
    plt_list<-list()
  
    pd2<-make_ribbon_dt(dt = pd,var = "min.cpg")
  
    merge(pd,pd2)->pd3
    (plt_list[["p_Genes"]]<-ggplot2::ggplot(data=pd3,aes(x=mDiff,y=DMRS,group=factor(type),color=factor(type)))+#,color=Contrast,group=varaible))+
      # geom_line(aes(linetype=factor(p.val)))+
      geom_line(aes(y=Genes))+
        # theme(axis.title.y.right = element_blank(),
        #       axis.text.y.right = element_blank())+
        theme_minimal()+
      # geom_line(aes(y=DMRS))+
      # geom_ribbon(aes(ymin=min.min.cpg.Genes,ymax=max.min.cpg.Genes))+
      facet_grid(Contrast~min.cpg)#,margins="min.cpg")
    )
    (plt_list[["p_r_Genes"]]<-ggplot2::ggplot(data=pd3,aes(x=mDiff,y=DMRS,group=factor(type),color=factor(type)))+#,color=Contrast,group=varaible))+
        # geom_line(aes(linetype=factor(p.val)))+
        geom_line(aes(y=Genes))+
        # geom_line(aes(y=DMRS))+
        geom_ribbon(aes(ymin=min.min.cpg.Genes,ymax=max.min.cpg.Genes,))+
        # theme_minimal()+
        facet_grid(Contrast~.,margins=F)
        
    )
    
    (plt_list[["p_DMRS"]]<-ggplot2::ggplot(data=pd3,aes(x=mDiff,y=DMRS,group=factor(type),color=factor(type)))+#,color=Contrast,group=varaible))+
        # geom_line(aes(linetype=factor(p.val)))+
        geom_line()+
        # geom_line(aes(y=DMRS))+
        # geom_ribbon(aes(ymin=min.min.cpg.Genes,ymax=max.min.cpg.Genes))+
        facet_grid(Contrast~min.cpg,margins="min.cpg")
    )
    (plt_list[["p_r_DMRS"]]<-ggplot2::ggplot(data=pd3,aes(x=mDiff,y=DMRS,group=factor(type),color=factor(type)))+#,color=Contrast,group=varaible))+
        # geom_line(aes(linetype=factor(p.val)))+
        # geom_line(aes(y=DMRS,color=factor(min.cpg)))+
        # geom_line(aes(y=DMRS))+
        geom_ribbon(aes(ymin=min.min.cpg.Genes,ymax=max.min.cpg.Genes))+
        # theme_minimal()+
        facet_grid(Contrast~.,margins=F)
      
    )
    
    lapply(1:length(plt_list),function(x)
      ggsave(plot = plt_list[[x]],
             filename = paste0(names(plt_list[x]),".",dev),
             path = path,
             device = dev))
  
  }
  return(pd)

}


make_ribbon_names<-function(var){
  
  varnames<-apply(expand.grid(c("min","max"),var,c("DMRS","Genes")),1,function(x)paste(x,collapse="."))
  return(varnames)
}

make_ribbon_dt<-function(var,dt){
  n<-names(dt)
  sdcols<-setdiff(n,c("Genes","DMRS",var))
  # make_ribbon(min.cpg,dt)
  
  dt2<-dt[,
          make_ribbon_names(var):=.(
            min(DMRS),
            max(DMRS),
            min(Genes),
            max(Genes)
            
          ),by=sdcols]
  return(dt2)
}

betasdmps <- function(betas,dmps,rgSet,anno_cols=c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")){
  require(data.table)
  dmps<-data.table::setDT(as.data.frame(dmps))
  ss<-rgSet@colData
  pvals <- minfi::detectionP(rgSet)
  p<-pvals[dmps$Name,]
  colnames(p)<-paste0(ss[colnames(pvals),"Sample_Name"],".Detection_Pval")
  b<-betas[dmps$Name,]
  colnames(b)<-paste0(ss[colnames(betas),"Sample_Name"],".AVG_betas")
  d<-dmps[,.SD,.SDcols=anno_cols]
  df<-cbind(d,b,p)
  df<-data.table::setDT(as.data.frame(df))
  data.table::setorder(df,chr,pos)
  # df[order(chr,pos)]
  return(df)
}



#' For each contrast extract a set of DMPs and add gene annotation and methylation values
#'
#'
#' @title Extract DMPs, annotation and methylation difference for each contrast
#'
#' @return data.table
#' @author Izar de Villasante
#' @export
#' @import minfi
#' @import data.table
#' @import limma
#' @param beta_normalized normalized betavalues, as produce by minfi::getBeta(grSet_noob)),
#'  where colnames(beta_normalized) == metadata$sample_Name
#' @param ContrastsDM list of contrasts as returned by limma::makeContrasts()
#' which will pass to limma topTable as input
#' @param mDiff absolute mean methylation difference between groups to filter by
#' @param ann annotation dataset from manifest with metadata such as gene info,
#' CGI, RefGene, etc. see topTable genelist arg.
#' @param writeOut save result as .csv default = TRUE.
#' @param writedir
#'
#' @inheritParams limma::topTable
#' @examples
#'
#' betas<-readRDS("data/beta_noob.rds")
#' fit<-readRDS("data/fit2.rds")
#' ann<-readRDS("data/ann.rds")
#' DMPann <- DMPextr(fit = fit,                       # linear contrast model
#'                   ContrastsDM = ContrastsDM,          # contrasts
#'                   p.value = 0.01,                      # filter significantly different probes
#'                   beta_normalized = beta_noob,        # extract mean group betas
#'                   mDiff = 0.5,                        # select mean methylation differences
#'                   ann = ann,                          # annotate positions (CGI, RefGene, etc)
#'                   writeOut = FALSE                    # write output to file


DMPextr <- function(
    fit, ContrastsDM=colnames(fit$contrasts), p.value, beta_normalized, mDiff, ann,
    writedir = "analysis/DMP_", writeOut = TRUE,ncores=NULL){
  require(foreach)
  
  ann<-data.table::as.data.table(ann,keep.rownames = "rn")
  data.table::setkey(ann,"rn")
  # ann      = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  annSub   = ann[rownames(beta_normalized)]

  if(is.null(ncores))ncores=length(ContrastsDM)
  
  cl<- parallel::makePSOCKcluster(ncores,outfile="")
  # parallel::clusterEvalQ(cl,{
  #   requireNamespace(c("limma","data.table"))
  #   data.table::setDTthreads(0)
  # })
  doParallel::registerDoParallel(cl)
  
  message("Processing ", length(ContrastsDM)," contrasts. Using ",ncores," cores.")
  res<-foreach::foreach(i=itertools::isplitIndices(length(ContrastsDM), chunks=ncores),#isplitIndices(1400,chunks=ncores),
                        .combine='rbind',
                        .multicombine = F,.packages = c("limma","data.table"),
                        # .export = c("guessArrayTypes",".default.epic.annotation"),
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{
    DMP_1 <- limma::topTable(fit,
                             num = Inf,
                             coef = i ,
                             genelist = annSub,
                             p.value = p.value  # = p.adj
    )
    
    if(nrow(DMP_1)<2){
      warning(paste("No DMP found for contrast:", ContrastsDM[i], sep=" "))
      next
    }
    else{
      DMP_1$Type <- "Hyper"
      DMP_1$Type[which(DMP_1$logFC < 0)] <- "Hypo"  # invert direction (see above0)
      #
      # DMP_1 <- DMP_1[ , c("chr", "pos", "strand", "Name",
      #                     "Type", "P.Value" , "adj.P.Val" ,
      #                     "Islands_Name", "Relation_to_Island",
      #                     "UCSC_RefGene_Name", "UCSC_RefGene_Accession","UCSC_RefGene_Group",
      #                     "Phantom4_Enhancers","Phantom5_Enhancers","X450k_Enhancer",
      #                     "Regulatory_Feature_Name", "Regulatory_Feature_Group",
      #                     "GencodeBasicV12_NAME","GencodeBasicV12_Accession", "GencodeBasicV12_Group",
      #                     "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group"
      # )]
      DMP_1$Contrast  <- ContrastsDM[i]
      
      # Extract the methylation values for the respective DMP and samples
      # Get contrast i:
      c1 <-fit$contrasts[,ContrastsDM[i]]
      c1 <-c1[c1!=0]
      # Variable names in contrast i:
      vars1 <-names(c1)
      # design matrix for contrasts:
      design <-fit$design[,vars1]
      design <-t(design) * c1
      
      # betas
      cg         <- which(rownames(beta_normalized) %in% rownames(DMP_1))
      betas <-beta_normalized[cg,]
      # diff  _mean:
      DMP_1$diff_meanMeth<-rowSums(apply(design,1,function(x) apply(betas%*%diag(x),1,function(y)mean(y[y!=0]))))
      
      # filter for absolute methylation difference
      DMP_1 <- DMP_1[which(abs(DMP_1$diff_meanMeth)>= mDiff), ]
      
      # save DMPs
      
      
      # write output file
      if(writeOut == TRUE){
        
        cat(paste("writing analysis/DMP_", ContrastsDM[i], ".csv\n",sep =""))
        data.table::fwrite(DMP_1, file = paste(writedir, ContrastsDM[i], ".csv", sep =""))
      }
  
      return(DMP_1)
    }
  }
  
  parallel::stopCluster(cl)
  return(res)
}




plotDMP <- function(DMPann,names,path=NULL){
  if(NROW(DMPann) >0){
    library(ggplot2)
    DMPresults <- data.frame(table(DMPann[ ,c("Contrast","Type")]))
    # plot DMPs (hypo/hyper)
    g1<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
      geom_bar(position="dodge", stat= "identity")+
      theme_bw()+
      scale_fill_manual(values=c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))+
      labs(x = "", y = "count", fill='Methylation')+
      ggtitle('Differently methylated probes')
    
    # ggplot2::ggsave(g1,paste0("analysis/DMPplots/",names,".png"))
    
    # plot with facets
    g2<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
      geom_bar(position="dodge", stat= "identity")+
      facet_wrap(.~Type, scales = "free_x") +
      theme_bw()+
      scale_fill_manual(values=c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))+
      labs(x = "", y = "count", fill='Methylation')+
      ggtitle('Differently methylated probes')
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_p0.01_m0.3.facet.png")
    
    # plot proportion of DMPs in CGI
    DMP_annCGI <- data.frame(DMPann[ ,c("Contrast","Type", "Relation_to_Island")])
    g3<-ggplot2::ggplot(DMP_annCGI, aes(Contrast, fill = Relation_to_Island)) +
      facet_wrap(.~Type, scales = "free_x") +
      geom_bar(position ="fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ylab("DMPs") +
      xlab("")
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_annCGI.png")
    
    #table(DMPann[ ,c("Relation_to_Island","Type", "Contrast")])
    
    # plot proportion of DMPs in genomic elements
    DMPann$UCSC_RefGene_Group[which(DMPann$UCSC_RefGene_Group == "")] <- "."
    DMP_annGenomic<-DMPann[,.(
      UCSC_RefGene_Group_short = unlist(lapply(strsplit(UCSC_RefGene_Group, ";"),'['))),
      by = c("Contrast","Type")]
    g4<-ggplot2::ggplot(DMP_annGenomic, aes(Contrast, fill = UCSC_RefGene_Group_short)) +
      facet_wrap(.~Type, scales = "free_x") +
      geom_bar(position = "fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(fill = "RefGene") +
      ylab("DMPs") +
      xlab("")
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_annGenomic.png")
    
    # table(DMPann[ ,c("UCSC_RefGene_Group_short","Type", "Contrast")])
    plt_list<-list(g1,g2,g3,g4)
    n <- c("DMP_annCGI.png", "DMP_annGenomic.png", "DMP_count_facet.png", "DMP_count.png")
    
    if(!is.null(path)){
      sapply(1:length(plt_list),function(x){
        ggplot2::ggsave(
          filename = n[x],
          plot = plt_list[[x]],
          device = NULL,
          path = path,
          scale = 1,
          width = 4,
          height = 5,
          units = c("in"),
          dpi = 600,
          limitsize = TRUE,
          bg = NULL
        )
      })
    } 
    return(plt_list)
  }else{
    warning("empty data.frame please, try again modifying filter params.")
    
  }
  
}



find_dmrs<-function(object, model, fdr = 0.05, p.value = "fdr", betacutoff = 0.3, min.cpg=5){
  # FDR threshold used to define DMRS is indexed at the rate of that of DMPs. This rate is defined at fdr. Should only use fdr not p.value
  require(DMRcate)
  require(S4Vectors)
  contrasts <- colnames(model$contrasts)
  future::plan("multisession")
  results <- future.apply::future_lapply(contrasts,future.seed = T, function(x){
    
    myAnnotation <- DMRcate::cpg.annotate(
      object = object,                    
      datatype = "array", 
      what = "Beta", 
      fdr = fdr,
      analysis.type = "differential", 
      design = model$design, 
      contrasts = TRUE,
      cont.matrix = model$contrasts, 
      coef = x, 
      arraytype = "EPIC"
    )
    
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one 
        # R expression in the "try" part then you'll have to 
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression 
        # in case the "try" part was completed successfully
        
        # message("This is the 'try' part")
        DMRs <- DMRcate::dmrcate(myAnnotation, 
                                 pcutoff = p.value,
                                 betacutoff= betacutoff,
                                 min.cpgs =min.cpg,
                                 C=2
        )
        
      },
      error=function(cond) {
        
        message(cond)
        # Choose a return value in case of error
        return(NA)
      },
      warning=function(cond) {
        
        message(cond)
        # Choose a return value in case of warning
        return(NULL)
      },
      finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
        message(paste("Processed contrast:", x))
        message("Next.")
      }
    )
    
    
    if(!identical(out@no.cpgs, integer(0))){
      results.ranges <- DMRcate::extractRanges(out)
      results.ranges$Contrast = x
      # data.frame(Contrast = x, results.ranges)
    }else{
      results.ranges<-NULL
    }
    
    return(results.ranges)
  }
  )
  
  # return(data.table::rbindlist(results) )
  
  while (is.list(results)) {
    results<-suppressWarnings(do.call("c",results))
  }
  return(results)
}
summarize<-function(dmps,dmrs,path="results/"){
  sdmps<-summary_dmps(dmps,write = F)
  sdmrs<-summary_dmrs(dmrs,write=F)
  s<-merge(sdmps,sdmrs)
  data.table::fwrite(s,path)
  return(s)
}
