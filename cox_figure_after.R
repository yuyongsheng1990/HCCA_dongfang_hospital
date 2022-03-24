# Title     : 肝门胆管癌病理分期c_index和BS折线图比较
# Objective : cox regression based R
# Created by: Yongsheng Yu
# Created on: 2022/3/16

library(xlsx)
library(randomForestSRC)
library(prodlim)
library(pec)
library(survival)
library(Matrix)
library(foreach, warn.conflicts = FALSE)
library(glmnet)
library(Hmisc, warn.conflicts = FALSE)
library(lattice)
library(Formula)
library(ggplot2)
library(SparseM, warn.conflicts = FALSE)
library(rms)
library(MASS)

#-----------------------------------------ML+backward stepwise----------------------------------------------------------
## 读入机器学习变量
df_after_ML_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital//data/result/cox/df_1.7.2.5_after_ML_tran_cox.xlsx", 1, encoding = "UTF-8")
df_after_ML_tran <- df_after_ML_tran[,-1]
names(df_after_ML_tran)
df_after_ML_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital//data/result/cox/df_1.7.2.5_after_ML_test_cox.xlsx", 1, encoding = "UTF-8")
df_after_ML_test <- df_after_ML_test[,-1]
names(df_after_ML_test)

# 删除训练集空值
dim(df_after_ML_tran) # 数据维度
df_after_ML_tran <-df_after_ML_tran[!(is.na(df_after_ML_tran$target_survival_month)),]
dim(df_after_ML_tran)
# 删除测试集空值
dim(df_after_ML_test) # 数据维度
df_after_ML_test <-df_after_ML_test[!(is.na(df_after_ML_test$target_survival_month)),]
dim(df_after_ML_test)
## 机器学习提取的cox回归
names(df_after_ML_tran)
b.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~surgery_result+AJCC_8+tumor_CA19.9+MSKCC+HBcAb+surgery_plasm+HBsAg+surgery_bleeding+smoke+emaciation+drinking+DB+tumor_AFP+tumor_CEA+PTCD_ERCP+biliary_disease+gender,data=df_after_ML_tran, x=TRUE,y=T)
summary(b.cox)
## backward stepwise筛选机器ML变量
AICBWCOX.obj <- stepAIC(b.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
# cox回归：ml+backward stepwise
b.d.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ surgery_result+AJCC_8+tumor_CA19.9+tumor_CEA, data=df_after_ML_tran, x=TRUE) 
summary(b.d.cox)
# 95%CI：c_index +/-1.96*se
CI.min <- 0.786-1.96*0.021
CI.min
CI.max <- 0.786+1.96*0.021
CI.max

# -------------------------------------模型筛选能力对比筛选-----------------------------------------------
## 读入全部变量
df_after_all_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.6_after_all_tran_cox.xlsx", 1, encoding = "UTF-8")
df_after_all_tran <- df_after_all_tran[,-1]
names(df_after_all_tran)
df_after_all_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.6_after_all_test_cox.xlsx", 1, encoding = "UTF-8")
df_after_all_test <- df_after_all_test[,-1]
names(df_after_all_test)

# 删除训练集空值
dim(df_after_all_tran) # 数据维度
df_after_all_tran <-df_after_all_tran[!(is.na(df_after_all_tran$target_survival_month)),]
dim(df_after_all_tran)
# 删除测试集空值
dim(df_after_all_test) # 数据维度
df_after_all_test <-df_after_all_test[!(is.na(df_after_all_test$target_survival_month)),]
dim(df_after_all_test)

# 全部变量cox回归
all.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ ., data=df_after_all_tran, x=TRUE,y=T) 
## backward stepwise筛选全部变量
AICBWCOX.obj <- stepAIC(all.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
# cox回归：全部变量+backward stepwise
all.d.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ gender+DB+tumor_CA19.9+AJCC_8+MSKCC+surgery_result, data=df_after_all_tran, x=TRUE) 
summary(all.d.cox)

# 比较模型c_index
surv.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ .)
pec.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ 1)
## compare C-index
set.seed(4) # 设置随机数种子，不同种子有不同效果
bcvCindex  <- pec::cindex(list("Machine Learning"=b.cox,"Backward Stepwise"=all.d.cox,"Machine Learning and Backward Stepwise"=b.d.cox),
                          formula=pec.f, newdata=df_after_all_test,
                          splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(0,88,1))
print(bcvCindex)
op<-par(mfrow=c(1,1))
plot(bcvCindex,Legend.cex=1, ylim = c(0.6,1),xlim = c(0,80), xlab="Time(month)")


## compare Brier score
set.seed(4)
bcvBS <- pec(list("Machine Learning"=b.cox,"Backward Stepwise"=all.d.cox,"Machine Learning and Backward Stepwise"=b.d.cox), 
             newdata = df_after_all_test, formula = pec.f, cens.model = "marginal", 
             splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F)

print(bcvBS)
plot(bcvBS,Legend.cex=1,xlab = "Time(month)",ylab="Barier Score")

# -----------------------------------------术后风险评分直方图-------------------------------------------
## 画出风险得分的直方图，带正态分布拟合曲线
## 术后风险评分_tran
df_after_staging_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.8.3_术后分期_训练集.xlsx", 1, encoding = "UTF-8")
df_after_staging_tran <- df_after_staging_tran[,-1]
names(df_after_staging_tran)

w=c(df_after_staging_tran$risk_score_3Y)
hist(w,breaks=50,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_after_staging_tran$risk_score_3Y),1)),xlab = "The Risk Scores of Training Set",xlim=range(0,8))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*160,col="blue",lwd=2)#添加一条正态曲线

## 术后风险评分_tran
df_after_staging_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.8.3_术后分期_测试集.xlsx", 1, encoding = "UTF-8")
df_after_staging_test <- df_after_staging_test[,-1]
names(df_after_staging_test)

w=c(df_after_staging_test$risk_score_3Y)
hist(w,breaks=20,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_after_staging_test$risk_score_3Y),1)),xlab = "The Risk Scores of Test Set",xlim=range(0,8))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*40,col="blue",lwd=2)#添加一条正态曲线

#-------------------------------------术后分期生存分析曲线图-------------------------------------------------

# 删除训练集空值
dim(df_after_staging_tran) # 数据维度
df_after_staging_tran <-df_after_staging_tran[!(is.na(df_after_staging_tran$target_survival_month)),]
dim(df_after_staging_tran)

# 画出术后训练集生存曲线
fit_us <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_after_3Y, data = df_after_staging_tran)
class(fit_us)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_us, data = df_after_staging_tran,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("red","black","green","dodgerblue","deepskyblue4"),
           xlab="The Survival Time of Training Set",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4","5"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# 删除测试集空值
dim(df_after_staging_test) # 数据维度
df_after_staging_test <-df_after_staging_test[!(is.na(df_after_staging_test$target_survival_month)),]
dim(df_after_staging_test)

# 画出术后测试集生存曲线
fit_us <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_after_3Y, data = df_after_staging_test)
class(fit_us)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_us, data = df_after_staging_test,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("red","black","green","dodgerblue","deepskyblue4"),
           xlab="The Survival Time of Test Set",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4","5"))

# ---------------------------------术后病理分期比较----------------------------------------------------

## compare 术前分期C-index
set.seed(142)

# 术后分期生存曲线
fit_score <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_after_3Y, data = df_after_staging_tran)
class(fit_score)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_score, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","green","dodgerblue","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4","5"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# 上一版分期系统
fit_existing <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_after_existing, data = df_after_staging_tran)
class(fit_existing)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_existing, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Former Risk Score System",
           legend.labs = c("1","2","3","4","5"),
           legend.adj = 1,
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# AJCC_8分期
fit_ajcc <- survfit(Surv(target_survival_month, target_3Y_death)~ AJCC_8, data = df_after_staging_tran)
class(fit_ajcc)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_ajcc, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "AJCC_8",
           legend.labs = c("0","2","3","4"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# 删除训练集MSKCC分期、Gazzaniga_T分期空值
df_after_staging_tran <-df_after_staging_tran[!(is.na(df_after_staging_tran$MSKCC)),]
dim(df_after_staging_tran)
# 删除测试集MSKCC分期空值
df_after_staging_test <-df_after_staging_test[!(is.na(df_after_staging_test$MSKCC)),]
dim(df_after_staging_test)

# Gazzaniga_T分期
fit_gazz <- survfit(Surv(target_survival_month, target_3Y_death)~ Gazzaniga_T, data = df_after_staging_tran)
class(fit_gazz)
ggsurvplot(fit_gazz, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Gazzaniga_T",
           legend.labs = c("1","2","3","4"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# MSKCC分期
fit_mskcc <- survfit(Surv(target_survival_month, target_3Y_death)~ MSKCC, data = df_after_staging_tran)
class(fit_mskcc)
ggsurvplot(fit_mskcc, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "MSKCC",
           legend.labs = c("1","2","3"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# Blumgart_T分期
fit_blumgart <- survfit(Surv(target_survival_month, target_3Y_death)~ Blumgart_T, data = df_after_staging_tran)
class(fit_blumgart)
ggsurvplot(fit_blumgart, data = df_after_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Blumgart_T",
           legend.labs = c("1","2","3","4"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

##画出多个评分系统的c_index和BS比较图
pec.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ 1)

set.seed(22) # 设置随机数种子，不同的随机数种子，模型效果不同
compare_Cindex  <- pec::cindex(list("Risk Score System"=fit_score,"Former Risk Score System"=fit_existing,"AJCC_8"=fit_ajcc, "Gazzaniga_T"=fit_gazz,"MSKCC"=fit_mskcc,"Blumgart_T"=fit_blumgart),
                               formula=pec.f, newdata=df_after_staging_test,
                               splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(1,80,1)) # BootCv集成方式会导致结果随机，可能报错，也可能不报错,要不就用"noPlan"
print(compare_Cindex)
op<-par(mfrow=c(1,1))
plot(compare_Cindex,legend.cex=0.8,ylim = c(0.5,1),xlim = c(0,80),xlab = "Time(month)")
x # 通过数据查看最好的随机数种子是多少

## compare 术前分期BS
set.seed(22) # 
compare_bs <- pec(list("Risk Score System"=fit_score,"Former Risk Score System"=fit_existing,"AJCC_8"=fit_ajcc, "Gazzaniga_T"=fit_gazz,"MSKCC"=fit_mskcc,"Blumgart_T"=fit_blumgart),
                  newdata = df_after_staging_test, formula = pec.f, cens.model = "marginal",
                  splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F) # BootCv集成方式会导致结果随机，可能报错，也可能不报错
print(compare_bs)
plot(compare_bs,legend.cex=0.8,legend.adj=0,xlim = c(0,80),xlab="Time(month)",ylab="Barier Score") # legend.cex调整标签大小缩放比例;legend.col调整颜色;legend.adj标签文字相对文职
x # 通过数据查看最好的随机数种子是多少

