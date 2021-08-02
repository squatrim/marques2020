# survival plot
survPlot <- function (data) {
  surv_data <- data %>%
    filter(IDH_codel.subtype == "IDHwt") %>% # select only IDH wildtype  
    dplyr::select(survival_months,vital_status,FOSL1) %>% 
    # mutate(cutoff_group = ifelse(FOSL1 >= quantile(FOSL1, probs=c(0.30, 0.5, 0.7), na.rm = TRUE)[3],"high","low"))
    mutate(cutoff_group = case_when(FOSL1 >= quantile(FOSL1, probs=c(0.30, 0.5, 0.7), na.rm = TRUE)[3] ~ "high",
                                    FOSL1 <= quantile(FOSL1, probs=c(0.30, 0.5, 0.7), na.rm = TRUE)[1] ~ "low"))  
   surv_data <-  na.omit(surv_data)
  
  my.Surv <- with(surv_data,Surv(time = survival_months, event = vital_status == 1))
  expr.surv <- survminer::surv_fit(my.Surv ~ cutoff_group, data = surv_data)
  smax <- max(surv_data[ ,"survival_months"], na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  log.rank <- survdiff(my.Surv ~ cutoff_group, rho = 0, data = surv_data)
  mantle.cox <- survdiff(my.Surv ~ cutoff_group, rho = 1, data = surv_data)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ cutoff_group, data=surv_data))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 6)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 6)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  legend.labs = c(sprintf("%s High \n(n=%s, median=%s)", "FOSL1", surv$records[1], surv$median[1]),
                  sprintf("%s Low \n(n=%s, median=%s)", "FOSL1", surv$records[2], surv$median[2]))
  
  surv_plot <- survminer::ggsurvplot(fit = expr.surv, censor = F, conf.int = F, legend = c(0.55,0.9), 
                                     surv.scale = "percent", ylab = "Surviving", xlab = "Survival time (Months)",
                                     legend.labs = legend.labs, legend.title = "", xlim = c(0,smax),
                                     font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(), 
                                     risk.table.y.text = F, risk.table.y.text.col = T, risk.table.height = 0.4, data = surv_data)
  
  surv_plot$plot <-  surv_plot$plot + annotate("text", x = xmax, y = c(0.65,0.575,0.505), size = 8/3,
                                               label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                         sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                                         sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p)))
  surv_plot <- ggpar(surv_plot,
                     font.main = c(10,"black","plain"),
                     font.x = c(10,"black","plain"),
                     font.y = c(10,"black","plain"),
                     font.tickslab = c(8,"black","plain"))
  surv_plot
}
