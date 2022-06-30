library(xlsx)
library(randomForestSRC)
library(prodlim)
library(pec)
library(survival)
library(Matrix)
library(foreach)
library(glmnet)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(SparseM)
library(rms)
library(MASS)

df_train <- read.xlsx("D:/aaa诺道医学/东方肝胆/df_importance_three_year_tran.xlsx", 1, encoding = "UTF-8")
df_train <- df_train[,-1]
names(df_train)
df_test <- read.xlsx("D:/aaa诺道医学/东方肝胆/df_importance_three_year_test.xlsx", 1, encoding = "UTF-8")
df_test <- df_test[,-1]
names(df_test)


## 机器学习提取的重要变量
b.cox <- coxph(Surv(Survival.time.M., Survival_state_three_year) ~ALP.U.l.+γ.GT.U.l.+N+T+albumin.g.l.+Diameter.cm.+AST.U.l.+DBIL.umol.l.+TBIL.umol.l.+PA.mg.l.+ALT.U.l.+AFP.μg.l.+CEA.μg.l.+CA19.9.U.ml.+Age,
               data=df_train, x=TRUE,y=T) 
summary(b.cox)
## backward stepwise筛选变量
AICBWCOX.obj <- stepAIC(b.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
## 选取上述模型中显著的变量
C.cox <- coxph(Surv(Survival.time.M., Survival_state_three_year) ~ T+N+CEA.μg.l.+CA19.9.U.ml.+AFP.μg.l.+PA.mg.l., data=df_train, x=TRUE) 
summary(C.cox)

## 构建风险评分系统，计算每个病人的风险评分，找出最小值和最大值
df_train$score=0.186*df_train$T+0.656*df_train$N+0.147*df_train$CEA.μg.l.+0.120*df_train$CA19.9.U.ml.+0.055*df_train$AFP.μg.l.-0.187*df_train$PA.mg.l.
range(df_train$score)
df_train$score=(1.2+df_train$score)*10
range(df_train$score)
describe(df_train$score)

df_test$score=0.186*df_test$T+0.656*df_test$N+0.147*df_test$CEA.μg.l.+0.120*df_test$CA19.9.U.ml.+0.055*df_test$AFP.μg.l.-0.187*df_test$PA.mg.l.
range(df_test$score)
df_test$score=(1.2+df_test$score)*10
range(df_test$score)
describe(df_test$score)


## 全部变量用于cox回归
df_all_tran <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.2_df_before_all_tran.xlsx", 1, encoding = "UTF-8")
df_all_tran <- df_all_tran[,-1]
df_all_tran <- df_all_tran[,-1]
df_all_tran <- df_all_tran[,-1]
names(df_all_tran)
df_all_tran <- df_all_tran[,-2]
##删除T.LCSGJ.
df_all_tran <- df_all_tran[,-11]
names(df_all_tran)
df_all_test <- read.xlsx("D:/PycharmProjects/HCCA_dongfang_hospital/data/result/cox/df_1.7.2.2_df_before_all_test.xlsx", 1, encoding = "UTF-8")
df_all_test <- df_all_test[,-1]
df_all_test <- df_all_test[,-1]
df_all_test <- df_all_test[,-1]
names(df_all_test)
df_all_test <- df_all_test[,-2]
##删除T.LCSGJ.
df_all_test <- df_all_test[,-11]
names(df_all_test)
all.cox <- coxph(Surv(target_3Y_survival, target_3Y_death) ~ ., data=df_all_tran, x=TRUE) 


## 所有变量经过backward stepwise cox回归
## backward stepwise筛选变量
AICBWCOX.obj <- stepAIC(all.cox,direction= c("backward"), steps = 1000)
summary(AICBWCOX.obj)
## 选取上述模型中显著的变量
d.cox <- coxph(Surv(Survival.time.M., Survival_state_three_year) ~ Sex+Age+History.of.stone+Smoking+HBV+T+N+M+CA19.9.U.ml.+PA.mg.l.+CEA.μg.l.+DBIL.umol.l.+TBIL.umol.l.+Excision+BloodType_A, data=df_all_tran, x=TRUE) 
summary(d.cox)


surv.f <- as.formula(Surv(Survival.time.M., Survival_state_three_year) ~ .)
pec.f <- as.formula(Surv(Survival.time.M., Survival_state_three_year) ~ 1)
## compare C-index
set.seed(142)
bcvCindex  <- pec::cindex(list("Machine Learning"=b.cox,"Backward Stepwise"=d.cox,"Machine Learning and Backward Stepwise"=C.cox),
                          formula=pec.f, data=df_all_test,
                          splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(1,80,1))
print(bcvCindex)
op<-par(mfrow=c(1,1))
plot(bcvCindex,Legend.cex=1.5)


## 训练集c-index
f <- cph(Surv(Survival.time.M., Survival_state_three_year) ~ T+N+CEA.μg.l.+CA19.9.U.ml.+AFP.μg.l.+PA.mg.l., x=T, y=T, surv=T, data=df_train)
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(Survival.time.M., Survival_state_three_year) ~ predict(f), data = df_train)

## compare Brier score
set.seed(1234)
prederror.df_top20 <- pec(list("Machine Learning"=b.cox,"Backward Stepwise"=d.cox,"Machine Learning and Backward Stepwise"=C.cox), 
                          data = df_all_test, formula = pec.f, cens.model = "marginal", 
                          splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F)

print(prederror.df_top20)
plot(prederror.df_top20,Legend.cex=1.5)

#op<-par(mfrow=c(2,2))
## 画出风险得分的直方图，带正态分布拟合曲线
w=c(df_train$score)
hist(w,breaks=20,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_train$score),1)),xlab = "The Risk Scores of Training Set",xlim=range(0,40))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*1112,col="blue",lwd=2)#添加一条正态曲线

w=c(df_test$score)
hist(w,breaks=20,freq=T,col = "gray",main=paste("The Median of Risk Socre：",round(median(df_test$score),1)),xlab = "The Risk Scores of Test Set",xlim=range(0,40))
x=seq(min(w),max(w),by=0.001)#做一组序列，用于绘制normal curve的x坐标
y=dnorm(x,mean(w),sd(w))#求x的正态分布函数值
lines(x,y*278,col="blue",lwd=2)#添加一条正态曲线

## 画出不同风险等级分组的生存曲线及logrank检验p值
myfun <- function(x){
  if(x<=10){x=1}
  else if((x>10)&(x<=20)){x=2}
  else if((x>20)&(x<=30)){x=3}
  else{x=4}
}

aaa <- lapply(c(df_train$score),myfun)
risk_group <-c()
for(i in 1:(length(aaa))){
  risk_group <- c(risk_group, as.data.frame(aaa[i])[1,1])
}
df_train$risk_group=risk_group
df_train$risk_group
names(df_train)

bbb <- lapply(c(df_test$score),myfun)
risk_group <-c()
for(i in 1:(length(bbb))){
  risk_group <- c(risk_group, as.data.frame(bbb[i])[1,1])
}
df_test$risk_group=risk_group
df_test$risk_group
names(df_test)

#df_train_1<-df_train[which((df_train$risk_group==3)|(df_train$risk_group==4)),]
fit <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data = df_train)
class(fit)
library(survminer)
library(ggpubr)
library(magrittr)
ggsurvplot(fit, data = df_train,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("skyblue3","cornflowerblue","dodgerblue","deepskyblue4"),
           xlab="The Survival Time of Training Set",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4"))

#df_test_1<-df_test[which((df_test$risk_group==3)|(df_test$risk_group==4)),]
fit_test <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data = df_test)
class(fit_test)
ggsurvplot(fit_test, data = df_test,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("skyblue3","cornflowerblue","dodgerblue","deepskyblue4"),
           xlab="The Survival Time of Test Set",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4"))


## 计算风险评分系统的c_index（基于测试集）
## 生成包括目标变量和系统指标变量的数据集
ccc <- read.xlsx("D:/aaa诺道医学/东方肝胆/df_all_test.xlsx", 1, encoding = "UTF-8")
ccc <- ccc[,-1]
names(ccc)
names(df_test)
df_test_compare <- cbind(df_test[c("Survival_state_three_year","Survival.time.M.","risk_group")],ccc[c("TNM","TNM.LCSGJ.")])
names(df_test_compare)

fit_score <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data=df_test_compare)
c_index <- summary(fit_score)$concordance
c_index
## 可以计算p值
#fit_score <- cph(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, x=T, y=T, surv=T, data=df_test_compare)
#validate(fit_score, method="boot", B=1000, dxy=T)
#rcorrcens(Surv(Survival.time.M., Survival_state_three_year) ~ predict(fit_score), data = df_test_compare)

fit_AJCC <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ TNM, data=df_test_compare)
c_index <- summary(fit_AJCC)$concordance
c_index
#fit_AJCC <- cph(Surv(Survival.time.M., Survival_state_three_year)~ TNM, x=T, y=T, surv=T, data=df_test_compare)
#validate(fit_AJCC, method="boot", B=1000, dxy=T)
#rcorrcens(Surv(Survival.time.M., Survival_state_three_year) ~ predict(fit_AJCC), data = df_test_compare)

fit_LCSGJ <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ TNM.LCSGJ., data=df_test_compare)
c_index <- summary(fit_LCSGJ)$concordance
c_index
#fit_LCSGJ <- cph(Surv(Survival.time.M., Survival_state_three_year)~ TNM.LCSGJ., x=T, y=T, surv=T, data=df_test_compare)
#validate(fit_LCSGJ, method="boot", B=1000, dxy=T)
#rcorrcens(Surv(Survival.time.M., Survival_state_three_year) ~ predict(fit_LCSGJ), data = df_test_compare)


## 画生存曲线图，不带置信区间和检验
fit_score <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data = df_test_compare)
class(fit_score)
ggsurvplot(fit_score, data = df_test_compare,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4"))

fit_AJCC <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ TNM, data = df_test_compare)
class(fit_AJCC)
ggsurvplot(fit_AJCC, data = df_test_compare,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "AJCC 8th Stage",
           legend.labs = c("1A","2","3A","3B","4"))
fit_LCSGJ <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ TNM.LCSGJ., data = df_test_compare)
class(fit_LCSGJ)
ggsurvplot(fit_LCSGJ, data = df_test_compare,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "LCSGJ Stage",
           legend.labs = c("1","2","3","4A","4B"))


##画出多个评分系统的c_index和BS比较图
pec.f <- as.formula(Surv(Survival.time.M., Survival_state_three_year) ~ 1)
## compare C-index
set.seed(142)
compare_Cindex  <- pec::cindex(list("Risk Score System"=fit_score,"AJCC 8th Stage"=fit_AJCC,"LCSGJ Stage"=fit_LCSGJ),
                          formula=pec.f, data=df_test_compare,
                          splitMethod="BootCv",B=3,confInt = TRUE,eval.times=seq(1,80,1))
print(compare_Cindex)
op<-par(mfrow=c(1,1))
plot(compare_Cindex)

## compare Brier score
set.seed(1234)
compare_bs <- pec(list("Risk Score System"=fit_score,"AJCC 8th Stage"=fit_AJCC,"LCSGJ Stage"=fit_LCSGJ), 
                          data = df_test_compare, formula = pec.f, cens.model = "marginal", 
                          splitMethod = "BootCv",B=3, confInt = TRUE, confLevel = 0.95, reference = F)

print(compare_bs)
plot(compare_bs)


## 验证集情况
df_val <- read.xlsx("D:/aaa诺道医学/东方肝胆/df_val.xlsx", 1, encoding = "UTF-8")
names(df_val)

df_val$score=0.186*df_val$T+0.656*df_val$N+0.147*df_val$CEA.μg.l.+0.120*df_val$CA19.9.U.ml.+0.055*df_val$AFP.μg.l.-0.187*df_val$PA.mg.l.
range(df_val$score)
df_val$score=(1.2+df_val$score)*10
range(df_val$score)
describe(df_val$score)

vvv <- lapply(c(df_val$score),myfun)
risk_group <-c()
for(i in 1:(length(vvv))){
  risk_group <- c(risk_group, as.data.frame(vvv[i])[1,1])
}
df_val$risk_group=risk_group
df_val$risk_group
names(df_val)

fit_score <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data=df_val)
c_index <- summary(fit_score)$concordance
c_index

fit_AJCC <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ TNM, data=df_val)
c_index <- summary(fit_AJCC)$concordance
c_index

fit_LCSGJ <-coxph(Surv(Survival.time.M., Survival_state_three_year)~ TNM.LCSGJ., data=df_val)
c_index <- summary(fit_LCSGJ)$concordance
c_index

## 画生存曲线图，不带置信区间和检验
fit_score <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ risk_group, data = df_val)
class(fit_score)
ggsurvplot(fit_score, data = df_val,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "Risk Score System",
           legend.labs = c("1","2","3","4"))

fit_AJCC <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ TNM, data = df_val)
class(fit_AJCC)
ggsurvplot(fit_AJCC, data = df_val,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "AJCC 8th Stage",
           legend.labs = c("1A","2","3A","3B","4"))
fit_LCSGJ <- survfit(Surv(Survival.time.M., Survival_state_three_year)~ TNM.LCSGJ., data = df_val)
class(fit_LCSGJ)
ggsurvplot(fit_LCSGJ, data = df_val,
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = F,
           size = 1,
           linetype = "solid",
           palette = c("springgreen2","red","dodgerblue","gold","blue2"),
           xlab="Survival Time",
           legend = c(0.90,0.85),
           legend.title = "LCSGJ Stage",
           legend.labs = c("1","2","3","4A","4B"))


##画出多个评分系统的c_index和BS比较图
pec.f <- as.formula(Surv(Survival.time.M., Survival_state_three_year) ~ 1)
## compare C-index
set.seed(142)
compare_Cindex  <- pec::cindex(list("Risk Score System"=fit_score,"AJCC 8th Stage"=fit_AJCC,"LCSGJ Stage"=fit_LCSGJ),
                               formula=pec.f, data=df_val,
                               splitMethod="BootCv",B=1,confInt = TRUE,eval.times=seq(1,80,1))
print(compare_Cindex)
op<-par(mfrow=c(1,1))
plot(compare_Cindex)

## compare Brier score
set.seed(1234)
compare_bs <- pec(list("Risk Score System"=fit_score,"AJCC 8th Stage"=fit_AJCC,"LCSGJ Stage"=fit_LCSGJ), 
                  data = df_val, formula = pec.f, cens.model = "marginal", 
                  splitMethod = "BootCv",B=1, confInt = TRUE, confLevel = 0.95, reference = F)

print(compare_bs)
plot(compare_bs)

