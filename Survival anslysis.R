library(Rmisc)
library(survival)
library(survminer)
library(lubridate)
library(foreign)
library(tableone)
library(broom)
library(dplyr)
library(knitr)
library(rms)
library(ComparisonSurv)
library(labelled)

survival <- read.dta("./survival.dta")

#create pvalue function
pvalue <- function(x, ...) UseMethod("pvalue")
pvalue.survdiff <- function (x, ...) 
{
  if (length(x$n) == 1) {
    df <- 1
    pval <- pchisq(x$chisq, 1, lower.tail = FALSE)
  } else {
    if (is.matrix(x$obs)) {
      otmp <- rowSums(x$obs)
      etmp <- rowSums(x$exp)
    } else {
      otmp <- x$obs
      etmp <- x$exp
    }
    df <- sum(etmp > 0) - 1
    pval <- pchisq(x$chisq, df, lower.tail = FALSE)
  }
  list(chisq = x$chisq, p.value = pval, df = df)
}


#event = dead, prog
#time = t_prog / t_dead
#adding survival object
survival$survObj_dead <- with(survival, Surv(t_dead, dead==1))
survival$survObj_prog <- with(survival, Surv(t_prog, prog==1))

## Descriptive study
#t_dead description
mean(survival$t_dead)
range(survival$t_dead)
# Use table one package
myvars <- c("stage", "sex", "age", "prog", "dead", "t_prog","t_dead" )
catvars <- c("stage", "sex", "prog", "dead")
tab1 <- CreateTableOne(vars = myvars, data=survival, factorVars = catvars, strata="tx")
tab1Mat <- print(tab1,  quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1Mat, file="Table1.csv")

###Survival plot
prog_fit <- survfit(Surv(t_prog, prog)~tx, data=survival)
ggsurvplot(
    fit = prog_fit, 
    legend.labs = c("No Treatment", "Treatment"), 
    legend.title = "Treatment group",
    xlab = "years",
    ylab="Overall progression-free probability", 
    title= "KM curve for Progression Free Rate", 
    pval=T, 
    risk.table = T
)

summary(prog_fit, times=c(1,2,3,5))
Fixpoint.test(survival$t_prog, survival$prog, survival$tx, 5)

#log-rank
survdiff(Surv(t_prog, prog)~tx, data=survival) %>% pvalue()
#cox ph
coxph(Surv(t_prog, prog)~tx, data=survival)

#using dead as event
dead_fit <-survfit(Surv(t_dead, dead)~tx, data=survival)
ggsurvplot(
    fit = dead_fit, 
    legend.labs = c("No Treatment", "Treatment"), 
    legend.title = "Treatment group",
    xlab = "years",
    ylab="Overall survival probability",
    title= "KM curve for Overall Survival", 
    pval=T, risk.table = T
)
summary(dead_fit, times = c(1,2,3,5))
Fixpoint.test(survival$t_dead, survival$dead, survival$tx, 5)

#log-rank
survdiff(Surv(t_dead, dead)~tx, data=survival) %>% pvalue()

### Coxph for progression free survival ###
#label data first
survival$tx <- as.factor(survival$tx)
survival$sex <- as.factor(survival$sex)
survival$stage<- as.factor(survival$stage)
var_label(survival)[1:4]<- c("Treatment", "Stage", "Sex", "Age")

mdrcox <- coxph(Surv(t_prog, prog)~tx +sex+stage+age, data=survival )
gtsummary::tbl_regression(
  mdrcox, exponentiate=T)

##checking coxph assumption
test <- cox.zph(mdrcox, global=TRUE)
test
ggcoxzph(test)
survplot(fit=npsurv(survObj_prog~factor(tx), survival), 
         loglog=TRUE, logt=TRUE,
         xlim = c(0,3), ylim=c(-1.5,0.5),
         conf="none",
         label.curves=list(keys="lines"))


##Coxph for overall  durvival
mdrcox <- coxph(survObj_dead~tx+sex+stage+age, data=survival)
gtsummary::tbl_regression(mdrcox, exponentiate=T)
test <- cox.zph(mdrcox, global=TRUE)
test
ggcoxzph(test)
#significance regarding stage was noted. (p=0.008)

#loglog method
survplot(fit=npsurv(survObj_dead~factor(stage), survival), 
         loglog=TRUE, logt=TRUE,
         xlim = c(0,3),
         label.curves=list(keys="lines"))

#stratified with sex
sex_plot <- list()
#sex =1
stra_survival <- survival %>% filter(sex==1)
stra_survModel <-survfit(survObj_dead~tx, data=stra_survival)
sex_plot[[1]] <- ggsurvplot(fit=stra_survModel, pval=T)
mdrcox <- coxph(survObj_dead~tx+stage+age, data=stra_survival)
cox.zph(mdrcox)
gtsummary::tbl_regression(mdrcox)

#sex= 0
stra_survival <- survival %>% filter(sex==0)
stra_survModel <-survfit(survObj_dead~tx, data=stra_survival)
sex_plot[[2]] <- ggsurvplot(fit=stra_survModel, pval=T)
mdrcox <- coxph(survObj_dead~tx+stage+age, data=stra_survival)
cox.zph(mdrcox) # not pass, stage p< 0.05
gtsummary::tbl_regression(mdrcox)

#stratified with stage
stage_plot <- list()
#Stage 1
stra_survival <- survival %>% filter(stage==1)
stra_survModel <-survfit(survObj_dead~tx, data=stra_survival)
stage_plot[[1]] <- ggsurvplot(
  fit = stra_survModel, pval=T, title="Stage 1")
survdiff(survObj_dead~tx, data=stra_survival)
mdrcox <- coxph(survObj_dead~tx+sex+age, data=stra_survival)
cox.zph(mdrcox)
gtsummary::tbl_regression(mdrcox, exponentiate=T)

#Stage 2
stra_survival <- survival %>% filter(stage==2)
stra_survModel <-survfit(survObj_dead~tx, data=stra_survival)
stage_plot[[2]] <- ggsurvplot(fit = stra_survModel, pval=T, 
                              title="Stage 2")
survdiff(survObj_dead~tx, data=stra_survival)
mdrcox <- coxph(survObj_dead~tx+sex+age, data=stra_survival)
cox.zph(mdrcox)
gtsummary::tbl_regression(mdrcox)

#Stage 3
stra_survival <- survival %>% filter(stage==3)
stra_survModel <-survfit(survObj_dead~tx, data=stra_survival)
stage_plot[[3]] <- ggsurvplot(fit = stra_survModel, pval=T,
                              title="Stage 3")
survdiff(survObj_dead~tx, data=stra_survival)
mdrcox <- coxph(survObj_dead~tx+sex+age, data=stra_survival)
cox.zph(mdrcox)
gtsummary::tbl_regression(mdrcox)

#combining the plot
arrange_ggsurvplots(stage_plot, ncol=3, nrow=1)
