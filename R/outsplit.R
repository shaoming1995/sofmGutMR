#' @title 肠道菌群MR工具
#' @param name 学号
#' @param key 密码
#' @param outfilename 肠道菌群作为结局下载的SNP数据文件位置
#' @export
outsplit<-function(name,key,outfilename){
  A<-name
  A<-as.numeric(gsub("DK","00",A))
  C<-A+key
  if(C==2310000){out<-read.csv(outfilename)
  head(out)
  colnames(out)<-c("id.outcome","chr","pos","SNP","other_allele.outcome","effect_allele.outcome","beta.outcome","se.outcome","meta","pval.outcome","N","nn")
  out$outcome<-out$id.outcome
  data_class<-unique(out$id.outcome)
  dir.create("切分好的肠道菌群结局文件")
  for (i in data_class){
    A<-subset(out,id.outcome==i)
    cfilename = paste0("切分好的肠道菌群结局文件/",i,".csv")
    write.csv(A,cfilename) }
  cat("肠道菌群结局已经完成切分")}else{
    cat("请联系客服获取账户学号和密码或微信联系SFM19950928")}
}



