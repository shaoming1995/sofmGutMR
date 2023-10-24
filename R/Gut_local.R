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
    A <- name
    A <- as.numeric(gsub("DK", "00", A))
    C <- A + key
    if (C == 2310000) {
      library(TwoSampleMR)
      dir.create(savefile)
      filename <- data.frame(dir(PATH))
      A_temp <- c()
      B_temp <- c()
      C_temp <- c()
      D_temp <- c()
      for (i in filename[, 1]) {
        ipath <- paste0(PATH, "\\", i)
        exp_temp <- read.csv(ipath, header = T)
        test2 <- (try(exp_temp <- clump_data(exp_temp, clump_kb = kb,
                                             clump_r2 = r2)))
        if (class(test2) == "try-error") {
          total <- merge(GWASsummay, exp_temp, by = "SNP")
          total$eaf.exposure <- NA   #修改这里，让从本地和IEU数据库获得的菌群数据一致
          exp3 <- total[, c("SNP", "effect_allele.exposure",
                            "other_allele.exposure", "beta.exposure", "se.exposure",
                            "pval.exposure", "id.exposure", "exposure",
                            "eaf.exposure")]
          out3 <- total[, c("SNP", "effect_allele.outcome",
                            "other_allele.outcome", "eaf.outcome", "beta.outcome",
                            "se.outcome", "pval.outcome", "id.outcome",
                            "outcome", "eaf.outcome")]
          dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                                action = 2)
          data_h<-dat%>%subset(dat$mr_keep==TRUE)
          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure)*(data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure","SNP","effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure","Fvalue","pval.exposure",
                                       "beta.outcome","se.outcome", "pval.outcome")] #导出SNP表格用
          data_h_TableS1$cluster <- 0  #0表达没有富集成功
          res <- mr(data_h)
          res$cluster <- 0    #0表达没有富集成功
          mr_OR<-generate_odds_ratios(res)
          mr_OR$or<-round(mr_OR$or,3)
          mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
          mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
          mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")") #生成成一个类似1.475(0.806-2.701)的格式，方便复制用
        }
        else {
          total <- merge(GWASsummay, exp_temp, by = "SNP")
          total$eaf.exposure <- NA
          exp3 <- total[, c("SNP", "effect_allele.exposure",
                            "other_allele.exposure", "beta.exposure", "se.exposure",
                            "pval.exposure", "id.exposure", "exposure",
                            "eaf.exposure")]
          out3 <- total[, c("SNP", "effect_allele.outcome",
                            "other_allele.outcome", "eaf.outcome", "beta.outcome",
                            "se.outcome", "pval.outcome", "id.outcome",
                            "outcome", "eaf.outcome")]
          dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                                action = 2)
          data_h<-dat%>%subset(dat$mr_keep==TRUE)
          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure)*(data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure","SNP","effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure","Fvalue","pval.exposure",
                                       "beta.outcome","se.outcome", "pval.outcome")]
          data_h_TableS1$cluster <- 1
          res <- mr(data_h)
          res$cluster <- 1
          mr_OR<-generate_odds_ratios(res)
          mr_OR$or<-round(mr_OR$or,3)
          mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
          mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
          mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
        }
        if (dim(res)[[1]] != 0) {
          het <- mr_heterogeneity(dat)
          ple <- mr_pleiotropy_test(dat)
          A_temp <- rbind(mr_OR, A_temp)
          B_temp <- rbind(het, B_temp)
          C_temp <- rbind(ple, C_temp)
          D_temp <- rbind(data_h_TableS1, D_temp)
          print(paste0("当前运行到", i, "文件"))
          Aname <- paste0(savefile, "\\", "肠道菌群与",
                          outname, "的MR结果.csv")
          Bname <- paste0(savefile, "\\", "肠道菌群与",
                          outname, "的异质性结果.csv")
          Cname <- paste0(savefile, "\\", "肠道菌群与",
                          outname, "的多效性结果.csv")
          Dname <- paste0(savefile, "\\", "肠道菌群与",
                          outname, "的SNPs情况.csv")
          write.csv(A_temp, Aname, row.names = F)
          write.csv(B_temp, Bname, row.names = F)
          write.csv(C_temp, Cname, row.names = F)
          write.csv(D_temp, Dname, row.names = F)
        }
        else {
          cat("请前往切分好的肠道菌群暴露文件下删除",
              i, "文件再次重新运行")
        }
      }
      cat("当前分析已全部完成！请前往", savefile,
          "文件夹下查看结果")
    }
    else {
      cat("请联系客服获取账户学号和密码或微信联系SFM19950928")
    }
 }


