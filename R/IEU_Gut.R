#' @title 肠道菌群MR工具
#' @param name 学号
#' @param key 密码
#' @param savefile 设置一个保存结果的文件夹名字
#' @param PATH GUT_outcome文件位置
#' @param GWASID 暴露的GWAS id号
#' @param expname 暴露的名称
#' @param p1 阈值
#' @param r2 相关系数
#' @param kb 聚类距离
#' @export

IEU_Gut<-function(name,key,savefile,PATH,GWASID,expname,p1,r2,kb){
  A<-name
  A<-as.numeric(gsub("DK","00",A))
  C<-A+key
  A_temp<-c()
  B_temp<-c()
  C_temp<-c()
  dir.create(savefile)
  if(C==2310000){
    library(TwoSampleMR)
    #暴露选择IEU中的吸烟-Ever smoked
    exp<-extract_instruments(outcomes = GWASID, p1 = p1,clump = T,r2=r2,kb=kb, p2 = 5e-08,access_token = NULL)
    exp$id.exposure<-expname
    exp$exposure<-expname
    for (i in 1:216){
      #在结局GWAS summary中寻找与暴露吸烟对应的SNPs
      ipath<-paste0(PATH,"/GUT_outcome.csv")
      OUT_temp<-read.csv(ipath,header = T)
      test<-try(OUT1<-extract_outcome_data(snps=exp$SNP,outcomes=OUT_temp[i,1],proxies=T,maf_threshold = 0.01,access_token = NULL))

      #去除与结局有gwas显著性的SNPs以及可能重复的SNP
      if(class(test)!="NULL"){
        OUT1<-subset(OUT1,pval.outcome>5e-08)
        OUT1<-OUT1[!duplicated(OUT1$SNP),]
      #去除回文
      dat<-harmonise_data(exposure_dat=exp,outcome_dat=OUT1,action=2)
      res <- mr(dat)
      het <- mr_heterogeneity(dat)
      ple <- mr_pleiotropy_test(dat)
      A_temp<-rbind(res,A_temp)
      B_temp<-rbind(het,B_temp)
      C_temp<-rbind(ple,C_temp)
      print(paste0("当前运行到",OUT_temp[i,3],"文件"))
      Aname<-paste0(savefile,"/",expname,"与肠道菌群","的MR结果.csv")
      Bname<-paste0(savefile,"/",expname,"与肠道菌群","的异质性结果.csv")
      Cname<-paste0(savefile,"/",expname,"与肠道菌群","的多效性结果.csv")
      write.csv(A_temp,Aname,row.names = F)
      write.csv(B_temp,Bname,row.names = F)
      write.csv(C_temp,Cname,row.names = F)}
      else{cat(i,"与",expname,"未匹配到合适的IVs")}}
    cat("当前分析已全部完成！请前往",savefile,"文件夹下查看结果")}
else{
    cat("请联系客服获取账户学号和密码或微信联系SFM19950928")
    }
}


