#' @title 肠道菌群MR工具
#' @param name 学号
#' @param key 密码
#' @param savefile 设置一个保存结果的文件夹名字
#' @param PATH 切分的肠道菌群暴露文件位置
#' @param GWASsummay 本地结局的预处理好的GWAS summay
#' @param kb 聚类距离
#' @param r2 相关系数
#' @param outname 本地结局的名称
#' @export
Gut_local<-function(name,key,savefile,PATH,GWASsummay,outname,kb,r2){
  A<-name
  A<-as.numeric(gsub("DK","00",A))
  C<-A+key
  if(C==2310000){
  library(TwoSampleMR)
  dir.create(savefile)
  filename<-data.frame(dir(PATH))
  A_temp<-c()
  B_temp<-c()
  C_temp<-c()
  for (i in filename[,1]){
      ipath<-paste0(PATH,"\\",i)
      exp_temp<-read.csv(ipath,header = T)
      exp_temp<-clump_data(exp_temp,clump_kb = kb,clump_r2 = r2)
      total<-merge(GWASsummay,exp_temp,by="SNP")
      if(dim(total)[[1]]!=0){
      total$eaf.exposure<-total$eaf.outcome
      #分别取出暴露与结局的数据
      exp3<-total[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure",
                     "pval.exposure","id.exposure","exposure", "eaf.exposure")]
      out3<-total[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome",
                     "pval.outcome","id.outcome","outcome","eaf.outcome")]
      #去除回文
      dat<-harmonise_data(exposure_dat=exp3,outcome_dat=out3,action=2)
      res <- mr(dat)
      het <- mr_heterogeneity(dat)
      ple <- mr_pleiotropy_test(dat)
      A_temp<-rbind(res,A_temp)
      B_temp<-rbind(het,B_temp)
      C_temp<-rbind(ple,C_temp)
      print(paste0("当前运行到",i,"文件"))
      Aname<-paste0(savefile,"\\","肠道菌群与",outname,"的MR结果.csv")
      Bname<-paste0(savefile,"\\","肠道菌群与",outname,"的异质性结果.csv")
      Cname<-paste0(savefile,"\\","肠道菌群与",outname,"的多效性结果.csv")
      write.csv(A_temp,Aname,row.names = F)
      write.csv(B_temp,Bname,row.names = F)
      write.csv(C_temp,Cname,row.names = F)}else{
        cat(i,"与",outname,"未找到工具变量不进行计算")
      }}
  cat("当前分析已全部完成！请前往",savefile,"文件夹下查看结果")}
  else{
    cat("请联系客服获取账户学号和密码或微信联系SFM19950928")
  }
}


