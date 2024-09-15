clinical_file="Data/clinical_data_7_features(325samples)_DFAGCNS_ying.csv" 
meta <- read.table(clinical_file, sep=',', row.names=1, header=T, quote="", comment="", check.names=F)   #meta：临床数据
meta$subtype=as.numeric(factor(meta$subtype))
meta$race=as.numeric(factor(meta$race))
colnames(meta)
meta$OS.time=as.numeric(meta[,6])/30   #将生存时间天数变成月数


#年龄和年龄分组
meta$age_group=ifelse(meta$age>median(meta$age),'older','younger')
# colnames(meta)[7] <- "survival_risk"
colnames(meta)[6] <- "time"  #改变某一列的名字

library(survival)
library(survminer)
#年龄
# sfit1 <- survfit(Surv(time, OS)~age_group, data=meta)   #这里用年龄分组先画一个为例，~gender可以换成其他的列，只要是有分组的就可以
# ggsurvplot(sfit1, conf.int=F, pval=TRUE)
# 
# #画的好看些
# ggsurvplot(sfit1,palette = c("#FF4500", "#2E9FDF"),
#            risk.table =TRUE,pval =TRUE,
#            conf.int =TRUE,xlab ="Time in months", 
#            ggtheme =theme_light(), 
#            ncensor.plot = TRUE)

#按风险画
sfit2 <-survfit(Surv(time, OS)~Risk, data=meta)
ggsurvplot(sfit2, conf.int=F, pval=TRUE)

#画的好看些
# ggsurvplot(sfit2,palette = c("#FF4500", "#2E9FDF"),
#            risk.table =TRUE,pval =TRUE,
#            conf.int =TRUE,xlab ="Time in months", 
#            ggtheme =theme_light(), #网格线
#            surv.median.line = "hv", # 添加中位生存时间线
#            legend.labs=c("1", "0"), 
#            legend.title="Risk_group",  #改图例名称 
#            pval.coord = c(152, 0.8),#p值位置坐标
#            pval.size =6,#p值字体大小
#            break.x.by=30 ,#x轴刻度的间距.例如30月一个刻度
#            ncensor.plot = FALSE)
dev.off()
ggsurvplot(sfit2,palette = c("#FF4500", "#2E9FDF"),  #ED1010红色  #185192蓝色  原来的代码：sfit2,palette = c("#FF4500", "#2E9FDF"),"#FF4500", "#2E9FDF" 
           risk.table =TRUE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           surv.median.line = "hv", # 添加中位生存时间线
           legend.labs=c("High", "Low"), 
           legend.title="Risk_group",  #改图例名称 
           pval.coord = c(152, 0.95),#p值位置坐标
           pval.size =6,#p值字体大小
           break.x.by=30 ,#x轴刻度的间距.例如30月一个刻度
           # ggtheme =theme_light(), #网格线
           ncensor.plot = FALSE)
# ggsurvplot(sfit2, # 创建的拟合对象
#            # data = meta,  # 指定变量数据来源
#            conf.int = TRUE, # 显示置信区间
#            pval = TRUE, # 添加P值
#            risk.table = FALSE, # 绘制累计风险曲线
#            xlab ="Time in months", 
#            # surv.median.line = "hv", # 添加中位生存时间线
#            # add.all = TRUE, # 添加总患者生存曲线
#            palette = "hue")  # 自定义调色板


# ggsave(("Pic/name.tiff")plot =plot=gg$plot, width = 18, height = 8.5,units = "in"))
#------------------------------------------------------------------------------
#单个基因的生存曲线


library(survival)
library(survminer)

exprSet_file <- "Data/four-omics_WGCNAafter_(sample-label)genes.csv"
clinical_file <- "Data/clinical_data_7_features(325samples).csv"
exprSet <- read.table(exprSet_file, sep=',', row.names=1, header=T, quote="", comment="", check.names=F)
meta <- read.table(clinical_file, sep=',', row.names=1, header=T, quote="", comment="", check.names=F)   #meta：临床数据

meta$subtype=as.numeric(factor(meta$subtype))   #将subtype转变成数字
meta$race=as.numeric(factor(meta$race))      #将race转变成数字
colnames(meta)    #展示meta的列
meta$OS.time=as.numeric(meta[,6])/30   #将生存时间天数变成月数


#年龄和年龄分组
meta$age_group=ifelse(meta$age>median(meta$age),'older','younger')    #在meta的最后一列增加age_group
colnames(meta)[7] <- "survival_risk"
colnames(meta)[6] <- "time"  #改变某一列的名字


#单个基因
g=rownames(exprSet)[1]
meta$gene=ifelse(as.integer(exprSet[g,])> median(as.integer(exprSet[g,])),'high','low')
#上面一行代码实现的是根据(exprSet[g,])是否大于中位数median来判断hign和low，这样就把病人分成了2组
sfit3=survfit(Surv(time,OS)~gene,data=meta)
ggsurvplot(sfit3,pval =TRUE, data = meta, risk.table = TRUE)

#多个基因:多个基因输出结果
gs=rownames(exprSet)[1:4]   #取4个基因为例
splots <- lapply(gs, function(g){
  meta$gene=ifelse(as.integer(exprSet[g,]) > median(as.integer(exprSet[g,])),'high','low')
  sfit4=survfit(Surv(time, OS)~gene, data=meta)
  ggsurvplot(sfit4,pval =TRUE, data = meta, risk.table = TRUE)
}) 
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 2, risk.table.height = 0.4)



### 2.logrank批量生存分析
#这里是根据logrankfile把每个基因p值都挑出来，然后根据需求挑出那些p<0.05或者p<0.01的基因
logrankfile = "OV_log_rank_p.Rdata"
if(!file.exists(logrankfile)){
  mySurv=with(meta,Surv(time, OS))
  log_rank_p <- apply(exprSet , 1 , function(gene){
    # gene=exprSet[1,]
    meta$group=ifelse(gene>median(gene),'high','low')  
    data.survdiff=survdiff(mySurv~group,data=meta)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    return(p.val)
  })
  log_rank_p=sort(log_rank_p)
  save(log_rank_p,file = logrankfile)
}
load(logrankfile)
table(log_rank_p<0.01) 
table(log_rank_p<0.05) 
#可以根据names把这些p<0.05/p<0.01的基因挑出来
lr = log_rank_p[log_rank_p<0.05]
names(lr)   #就可以看到p<0.05有哪些基因

