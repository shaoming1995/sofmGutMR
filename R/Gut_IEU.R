#' @title 肠道菌群MR工具
#' @param savefile 设置一个保存结果的文件夹名字
#' @param PATH 切分的肠道菌群暴露文件位置
#' @param GWASID 结局的GWAS id号
#' @param kb 聚类的距离
#' @param r2 相关系数
#' @param outname 结局的英文名称或者缩写
#'
#' @export
Gut_IEU<-function(savefile,PATH,GWASID,outname,kb,r2){
  library(TwoSampleMR)
  dir.create(savefile)
  filename<-data.frame(dir(PATH))
  A_temp<-c()
  B_temp<-c()
  C_temp<-c()
for (i in filename[,1]){
  ipath<-paste0(PATH,"/",i)
  exp_temp<-read.csv(ipath,header = T)
  test2 <- (try(exp_temp1 <- clump_data(exp_temp, clump_kb = kb,
                                       clump_r2 = r2)))
  if (class(test2) != "try-error") {
   OUT<-extract_outcome_data(snps=exp_temp1$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
   if(dim(OUT)[[1]]!=0){
   OUT$id.outcome<-outname
   OUT$outcome<-outname
   OUT<-OUT[!duplicated(OUT$SNP),]
   exp_temp_out<-merge(exp_temp1,OUT,by="SNP",all=F)
   exp_temp_out$eaf.exposure<-NA
   #分别取出暴露与结局的数据
   exp<- exp_temp_out[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure",
                         "pval.exposure","id.exposure","exposure", "eaf.exposure")]
   out<- exp_temp_out[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome",
                         "pval.outcome","id.outcome","outcome","eaf.outcome")]
   #调控方向
   dat <-harmonise_data(exposure_dat=exp,outcome_dat=out,action=2)
   res <- mr(dat)
   res$cluster <- 1    #0表达富集成功
   mr_OR<-generate_odds_ratios(res)
   mr_OR$or<-round(mr_OR$or,3)
   mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
   mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
   mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
   het <- mr_heterogeneity(dat)
   ple <- mr_pleiotropy_test(dat)
   A_temp<-rbind(res,A_temp)
   B_temp<-rbind(het,B_temp)
   C_temp<-rbind(ple,C_temp)
   print(paste0("当前运行到",i,"文件"))
   Aname<-paste0(savefile,"/","肠道菌群与",outname,"的MR结果.csv")
   Bname<-paste0(savefile,"/","肠道菌群与",outname,"的异质性结果.csv")
   Cname<-paste0(savefile,"/","肠道菌群与",outname,"的多效性结果.csv")
   write.csv(A_temp,Aname,row.names = F)
   write.csv(B_temp,Bname,row.names = F)
   write.csv(C_temp,Cname,row.names = F)}
   else{
     cat(i,"与",outname,"未找到工具变量不进行计算")}}
  else{
    cat(i,"由于网络502问题未完成clump")
  OUT<-extract_outcome_data(snps=exp_temp$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
  if(dim(OUT)[[1]]!=0){
  OUT$id.outcome<-outname
  OUT$outcome<-outname
  OUT<-OUT[!duplicated(OUT$SNP),]
  exp_temp_out<-merge(exp_temp,OUT,by="SNP",all=F)
  exp_temp_out$eaf.exposure<-exp_temp_out$eaf.outcome
  #分别取出暴露与结局的数据
  exp<- exp_temp_out[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure",
                 "pval.exposure","id.exposure","exposure", "eaf.exposure")]
  out<- exp_temp_out[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome",
                 "pval.outcome","id.outcome","outcome","eaf.outcome")]
    #调控方向
  dat <-harmonise_data(exposure_dat=exp,outcome_dat=out,action=2)
      res <- mr(dat)
      res$cluster <- 0    #0表达没有富集成功
      mr_OR<-generate_odds_ratios(res)
      mr_OR$or<-round(mr_OR$or,3)
      mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
      mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
      mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
      het <- mr_heterogeneity(dat)
      ple <- mr_pleiotropy_test(dat)
      A_temp<-rbind(res,A_temp)
      B_temp<-rbind(het,B_temp)
      C_temp<-rbind(ple,C_temp)
      print(paste0("当前运行到",i,"文件"))
      Aname<-paste0(savefile,"/","肠道菌群与",outname,"的MR结果.csv")
      Bname<-paste0(savefile,"/","肠道菌群与",outname,"的异质性结果.csv")
      Cname<-paste0(savefile,"/","肠道菌群与",outname,"的多效性结果.csv")
      write.csv(A_temp,Aname,row.names = F)
      write.csv(B_temp,Bname,row.names = F)
      write.csv(C_temp,Cname,row.names = F)}else{
        cat(i,"与",outname,"未找到工具变量不进行计算")}
   }
  }
    cat("当前分析已完成请前往肠道菌群MR结果文件夹下查看结果")
}

