#' @title 肠道菌群MR工具
#' @param name 学号
#' @param key 密码
#' @param expfilename 肠道菌群作为暴露下载的SNP数据文件位置
#' @export
expsplit<-function(name,key,expfilename){
A<-name
A<-as.numeric(gsub("DK","00",A))
C<-A+key
if(C==2310000){
exp<-read.csv(expfilename)
head(exp)
colnames(exp)<-c("id.exposure","chr","pos","SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","meta","pval.exposure","N","nn")
exp$exposure<-exp$id.exposure
data_class<-unique(exp$id.exposure)
dir.create("切分好的肠道菌群暴露文件")
for (i in data_class){
  A<-subset(exp,id.exposure==i)
  cfilename = paste0("切分好的肠道菌群暴露文件\\",i,".csv")
  write.csv(A,cfilename) }
cat("肠道菌群暴露已经完成切分")}else{
  cat("请联系客服获取账户学号和密码或微信联系SFM19950928")}
}

