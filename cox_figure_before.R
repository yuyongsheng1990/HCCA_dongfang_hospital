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
df_before_ML_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.5_before_ML_tran_cox.xlsx", 1, encoding = "UTF-8")
df_before_ML_tran <- df_before_ML_tran[,-1]
names(df_before_ML_tran)
dim(df_before_ML_tran)
df_before_ML_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.5_before_ML_test_cox.xlsx", 1, encoding = "UTF-8")
df_before_ML_test <- df_before_ML_test[,-1]
names(df_before_ML_test)

# 删除训练集空值
dim(df_before_ML_tran) # 数据维度
df_before_ML_tran <-df_before_ML_tran[!(is.na(df_before_ML_tran$target_survival_month)),]
dim(df_before_ML_tran)
# 删除测试集空值
dim(df_before_ML_test) # 数据维度
df_before_ML_test <-df_before_ML_test[!(is.na(df_before_ML_test$target_survival_month)),]
dim(df_before_ML_test)
## 机器学习提取的cox回归
b.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~tumor_CA19.9+Gazzaniga_T+Blumgart_T+MSKCC+jaundice+drinking+biliary_disease+emaciation+tumor_CEA+HBsAg+DB+tumor_AFP+Bismuth_C+cardio_disease,
               data=df_before_ML_tran, x=TRUE,y=T) 
summary(b.cox)
## backward stepwise筛选机器ML变量
AICBWCOX.obj <- stepAIC(b.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
# cox回归：ml+backward stepwise
b.d.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ tumor_CA19.9+MSKCC+biliary_disease+tumor_CEA+DB, data=df_before_ML_tran, x=TRUE) 
summary(b.d.cox)
# 95%CI：c_index +/-1.96*se
CI.min <- 0.746-1.96*0.022
CI.min
CI.max <- 0.746+1.96*0.022
CI.max

# -------------------------------------模型筛选能力对比筛选-----------------------------------------------
## 读入全部变量
df_before_all_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.6_before_all_tran_cox.xlsx", 1, encoding = "UTF-8")
df_before_all_tran <- df_before_all_tran[,-1]
names(df_before_all_tran)
df_before_all_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.6_before_all_test_cox.xlsx", 1, encoding = "UTF-8")
df_before_all_test <- df_before_all_test[,-1]
names(df_before_all_test)

# 删除训练集空值
dim(df_before_all_tran) # 数据维度
df_before_all_tran <-df_before_all_tran[!(is.na(df_before_all_tran$target_survival_month)),]
dim(df_before_all_tran)
# 删除测试集空值
dim(df_before_all_test) # 数据维度
df_before_all_test <-df_before_all_test[!(is.na(df_before_all_test$target_survival_month)),]
dim(df_before_all_test)

# 全部变量cox回归
all.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ ., data=df_before_all_tran, x=TRUE,y=T) 
## backward stepwise筛选全部变量
AICBWCOX.obj <- stepAIC(all.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
# cox回归：全部变量+backward stepwise
all.d.cox <- coxph(Surv(target_survival_month, target_3Y_death) ~ DB+tumor_CEA+tumor_CA19.9+Blumgart_T, data=df_before_all_tran, x=TRUE) 
summary(all.d.cox)

# 比较模型c_index
surv.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ .)
pec.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ 1)
## compare C-index
set.seed(3)
bcvCindex  <- pec::cindex(list("Machine Learning"=b.cox,"All Variables and Backward Stepwise"=all.d.cox,"Machine Learning and Backward Stepwise"=b.d.cox),
                          formula=pec.f, newdata=df_before_all_test,
                          splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(1,88,1))
print(bcvCindex)
op<-par(mfrow=c(1,1))
plot(bcvCindex,Legend.cex=1, ylim = c(0.6,1),xlim=range(0,80),xlab = "Time(month)", legend.title=element_blank())  # element_blank()移除所有表题
legend("topright", inset = 0.05, c("Machine Learning","All Variables and Backward Stepwise","Machine Learning and Backward Stepwise"),
        lty = c(1), col = c('black','red','green'), box.lty=0) # 添加线条标签; legend去除外边框线box.lty=0
# mtext("Chick1-diet1", side = 2,line = -6,las = 1, col = "red") # 添加线条标注
x # 获得随机数种子22

## compare Brier score
set.seed(3)
bcvBS <- pec(list("Machine Learning"=b.cox,"Backward Stepwise"=all.d.cox,"Machine Learning and Backward Stepwise"=b.d.cox), 
                          newdata = df_before_all_test, formula = pec.f, cens.model = "marginal", 
                          splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F)
print(bcvBS)
plot(bcvBS,Legend.cex=1, ylim = c(0,0.3),xlim=range(0,80),xlab = "Time(month)",ylab = "Barier Score")

# -----------------------------------------术前风险评分直方图-------------------------------------------
## 画出风险得分的直方图，带正态分布拟合曲线
## 术前风险评分_tran
df_before_staging_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.8.3_术前分期_训练集.xlsx", 1, encoding = "UTF-8")
df_before_staging_tran <- df_before_staging_tran[,-1]
names(df_before_staging_tran)

w=c(df_before_staging_tran$risk_score_3Y)
hist(w,breaks=50,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_before_staging_tran$risk_score_3Y),1)),xlab = "The Risk Scores of Training Set",xlim=range(0,3))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*80,col="blue",lwd=2)#添加一条正态曲线

## 术前风险评分_tran
df_before_staging_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.8.3_术前分期_测试集.xlsx", 1, encoding = "UTF-8")
df_before_staging_test <- df_before_staging_test[,-1]
names(df_before_staging_test)

w=c(df_before_staging_test$risk_score_3Y)
hist(w,breaks=20,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_before_staging_test$risk_score_3Y),1)),xlab = "The Risk Scores of Test Set",xlim=range(0,3))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*20,col="blue",lwd=2)#添加一条正态曲线

#-------------------------------------术前分期生存分析曲线图(with 置信区间)-------------------------------------------------
# 删除训练集生存时长空值
dim(df_before_staging_tran) # 数据维度
df_before_staging_tran <-df_before_staging_tran[!(is.na(df_before_staging_tran$target_survival_month)),]
dim(df_before_staging_tran)

# 删除测试集生存时长空值
dim(df_before_staging_test) # 数据维度
df_before_staging_test <-df_before_staging_test[!(is.na(df_before_staging_test$target_survival_month)),]
dim(df_before_staging_test)

# 画出术前训练集生存曲线
fit_us <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_before_3Y, data = df_before_staging_tran)
class(fit_us)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_us, data = df_before_staging_tran,
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

# 画出术前测试集生存曲线
fit_us <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_before_3Y, data = df_before_staging_test)
class(fit_us)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_us, data = df_before_staging_test,
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
           legend.labs = c("1","2","3","4","5"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间

# ---------------------------------术前分期生存曲线----------------------------------------------------

## compare 术前分期C-index
set.seed(3)

# 术前分期生存曲线
fit_score <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_before_3Y, data = df_before_staging_tran)
class(fit_score)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_score, data = df_before_staging_tran,
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
'''
# 现有方法分期
fit_existing <- survfit(Surv(target_survival_month, target_3Y_death)~ staging_before_existing, data = df_before_staging_tran)
class(fit_existing)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit_existing, data = df_before_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","blue2"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Former Risk Score System",
           legend.labs = c("1","2","3","4"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间
'''

# Gazzaniga_T分期
fit_gazz <- survfit(Surv(target_survival_month, target_3Y_death)~ Gazzaniga_T, data = df_before_staging_tran)
class(fit_gazz)
ggsurvplot(fit_gazz, data = df_before_staging_tran,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold"),
           xlab="Survival Time(month)",
           legend = c(0.90,0.85),
           legend.title = "Gazzaniga_T",
           legend.labs = c("1","2","3","4"),
           break.x.by=20, # 设置x轴刻度间距
           xlim=c(0,80)) # 设置x轴刻度区间
# MSKCC分期
fit_mskcc <- survfit(Surv(target_survival_month, target_3Y_death)~ MSKCC, data = df_before_staging_tran)
class(fit_mskcc)
ggsurvplot(fit_mskcc, data = df_before_staging_tran,
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
fit_blumgart <- survfit(Surv(target_survival_month, target_3Y_death)~ Blumgart_T, data = df_before_staging_tran)
class(fit_blumgart)
ggsurvplot(fit_blumgart, data = df_before_staging_tran,
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

#------------------------------比较术前分期c_index和BS----------------------------------------
# 删除训练集空值
dim(df_before_staging_tran) # 数据维度
df_before_staging_tran <-df_before_staging_tran[!(is.na(df_before_staging_tran$Gazzaniga_T)),]
dim(df_before_staging_tran)

# 删除测试集空值
dim(df_before_staging_test) # 数据维度
df_before_staging_test <-df_before_staging_test[!(is.na(df_before_staging_test$Gazzaniga_T)),]
dim(df_before_staging_test)

##画出多个评分系统的c_index和BS比较图
pec.f <- as.formula(Surv(target_survival_month, target_3Y_death) ~ 1)
set.seed(3)
compare_Cindex  <- pec::cindex(list("Risk Score System"=fit_score,"Gazzaniga_T"=fit_gazz,"MSKCC"=fit_mskcc,"Blumgart_T"=fit_blumgart),
                          formula=pec.f, newdata=df_before_staging_test,
                          splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(1,88,1)) # BootCv集成方式会导致结果随机，可能报错，也可能不报错,要不就用"noPlan"
print(compare_Cindex)
op<-par(mfrow=c(1,1))
plot(compare_Cindex,legend.cex=1,ylim = c(0.5,0.9),xlim = c(0,80),xlab = "Time(month)")
# 查看随机数种子 x

## compare 术前分期BS
set.seed(3)
compare_bs <- pec(list("Risk Score System"=fit_score,"Gazzaniga_T"=fit_gazz,"MSKCC"=fit_mskcc,"Blumgart_T"=fit_blumgart),
                  newdata = df_before_staging_test, formula = pec.f, cens.model = "marginal", 
                  splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F) # BootCv集成方式会导致结果随机，可能报错，也可能不报错,要不就用"noPlan"
print(compare_bs)
plot(compare_bs,legend.cex=0.8,xlim = c(0,80), xlab = "Time(month)", ylab = "Barier Score")



