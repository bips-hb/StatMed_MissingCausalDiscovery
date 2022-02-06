library(ggplot2)

load("res_CurvePlot.RData")

rate_A_50_00 <- apply(res_A_50_00, 1, function(x) {mean(x<0.05)})
rate_A_50_10 <- apply(res_A_50_10, 1, function(x) {mean(x<0.05)})
rate_A_50_30 <- apply(res_A_50_30, 1, function(x) {mean(x<0.05)})
rate_A_50_50 <- apply(res_A_50_50, 1, function(x) {mean(x<0.05)})
rate_A_50_70 <- apply(res_A_50_70, 1, function(x) {mean(x<0.05)})

rate_B_50_00 <- apply(res_B_50_00, 1, function(x) {mean(x<0.05)})
rate_B_50_10 <- apply(res_B_50_10, 1, function(x) {mean(x<0.05)})
rate_B_50_30 <- apply(res_B_50_30, 1, function(x) {mean(x<0.05)})
rate_B_50_50 <- apply(res_B_50_50, 1, function(x) {mean(x<0.05)})
rate_B_50_70 <- apply(res_B_50_70, 1, function(x) {mean(x<0.05)})

rate_C_50_00 <- apply(res_C_50_00, 1, function(x) {mean(x<0.05)})
rate_C_50_10 <- apply(res_C_50_10, 1, function(x) {mean(x<0.05)})
rate_C_50_30 <- apply(res_C_50_30, 1, function(x) {mean(x<0.05)})
rate_C_50_50 <- apply(res_C_50_50, 1, function(x) {mean(x<0.05)})
rate_C_50_70 <- apply(res_C_50_70, 1, function(x) {mean(x<0.05)})

rate_D_50_00 <- apply(res_D_50_00, 1, function(x) {mean(x<0.05)})
rate_D_50_10 <- apply(res_D_50_10, 1, function(x) {mean(x<0.05)})
rate_D_50_30 <- apply(res_D_50_30, 1, function(x) {mean(x<0.05)})
rate_D_50_50 <- apply(res_D_50_50, 1, function(x) {mean(x<0.05)})
rate_D_50_70 <- apply(res_D_50_70, 1, function(x) {mean(x<0.05)})

rate_A_500_00 <- apply(res_A_500_00, 1, function(x) {mean(x<0.05)})
rate_A_500_10 <- apply(res_A_500_10, 1, function(x) {mean(x<0.05)})
rate_A_500_30 <- apply(res_A_500_30, 1, function(x) {mean(x<0.05)})
rate_A_500_50 <- apply(res_A_500_50, 1, function(x) {mean(x<0.05)})
rate_A_500_70 <- apply(res_A_500_70, 1, function(x) {mean(x<0.05)})

rate_B_500_00 <- apply(res_B_500_00, 1, function(x) {mean(x<0.05)})
rate_B_500_10 <- apply(res_B_500_10, 1, function(x) {mean(x<0.05)})
rate_B_500_30 <- apply(res_B_500_30, 1, function(x) {mean(x<0.05)})
rate_B_500_50 <- apply(res_B_500_50, 1, function(x) {mean(x<0.05)})
rate_B_500_70 <- apply(res_B_500_70, 1, function(x) {mean(x<0.05)})

rate_C_500_00 <- apply(res_C_500_00, 1, function(x) {mean(x<0.05)})
rate_C_500_10 <- apply(res_C_500_10, 1, function(x) {mean(x<0.05)})
rate_C_500_30 <- apply(res_C_500_30, 1, function(x) {mean(x<0.05)})
rate_C_500_50 <- apply(res_C_500_50, 1, function(x) {mean(x<0.05)})
rate_C_500_70 <- apply(res_C_500_70, 1, function(x) {mean(x<0.05)})

rate_D_500_00 <- apply(res_D_500_00, 1, function(x) {mean(x<0.05)})
rate_D_500_10 <- apply(res_D_500_10, 1, function(x) {mean(x<0.05)})
rate_D_500_30 <- apply(res_D_500_30, 1, function(x) {mean(x<0.05)})
rate_D_500_50 <- apply(res_D_500_50, 1, function(x) {mean(x<0.05)})
rate_D_500_70 <- apply(res_D_500_70, 1, function(x) {mean(x<0.05)})

y <- 100*c(rate_A_50_00, rate_A_50_10, rate_A_50_30, rate_A_50_50, rate_A_50_70,
       rate_B_50_00, rate_B_50_10, rate_B_50_30, rate_B_50_50, rate_B_50_70,
       rate_C_50_00, rate_C_50_10, rate_C_50_30, rate_C_50_50, rate_C_50_70,
       rate_D_50_00, rate_D_50_10, rate_D_50_30, rate_D_50_50, rate_D_50_70,
      rate_A_500_00, rate_A_500_10, rate_A_500_30, rate_A_500_50, rate_A_500_70,
      rate_B_500_00, rate_B_500_10, rate_B_500_30, rate_B_500_50, rate_B_500_70,
      rate_C_500_00, rate_C_500_10, rate_C_500_30, rate_C_500_50, rate_C_500_70,
      rate_D_500_00, rate_D_500_10, rate_D_500_30, rate_D_500_50, rate_D_500_70)



Method <- rep(c("Test-wise deletion","Multiple imputation"), 40)
q <- rep(c(0, 10, 30, 50, 70), 8, each=2)
n <- rep(c("n=50", "n=500"), each=40)
scenario <- rep(c("Scenario A","Scenario B","Scenario C","Scenario D"), 2,
                each=10)

plot_dat <- data.frame(y, Method, q, n, scenario)

cairo_ps("CurvePlot.eps",
    height=4, width=6)
  ggplot(data=plot_dat, aes(x=q, y=y, color=Method, lty=Method)) +
    geom_point() + geom_line() + facet_grid(cols=vars(scenario), rows=vars(n)) +
    xlab("Percent missing") + ylab("Power (%)") +  expand_limits(x = 0.8) +
    scale_x_continuous(breaks=c(0,10,30,50,70), minor_breaks = NULL) +
    scale_y_continuous(limits=c(0, 100)) +
    theme_bw() +
    theme(legend.position = c(0.75, 0.83),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
dev.off()


### Supplementary Material

load("res_CurveNull.RData")

rate_A_50_00 <- apply(res_A_50_00, 1, function(x) {mean(x<0.05)})
rate_A_50_10 <- apply(res_A_50_10, 1, function(x) {mean(x<0.05)})
rate_A_50_30 <- apply(res_A_50_30, 1, function(x) {mean(x<0.05)})
rate_A_50_50 <- apply(res_A_50_50, 1, function(x) {mean(x<0.05)})
rate_A_50_70 <- apply(res_A_50_70, 1, function(x) {mean(x<0.05)})

rate_B_50_00 <- apply(res_B_50_00, 1, function(x) {mean(x<0.05)})
rate_B_50_10 <- apply(res_B_50_10, 1, function(x) {mean(x<0.05)})
rate_B_50_30 <- apply(res_B_50_30, 1, function(x) {mean(x<0.05)})
rate_B_50_50 <- apply(res_B_50_50, 1, function(x) {mean(x<0.05)})
rate_B_50_70 <- apply(res_B_50_70, 1, function(x) {mean(x<0.05)})

rate_C_50_00 <- apply(res_C_50_00, 1, function(x) {mean(x<0.05)})
rate_C_50_10 <- apply(res_C_50_10, 1, function(x) {mean(x<0.05)})
rate_C_50_30 <- apply(res_C_50_30, 1, function(x) {mean(x<0.05)})
rate_C_50_50 <- apply(res_C_50_50, 1, function(x) {mean(x<0.05)})
rate_C_50_70 <- apply(res_C_50_70, 1, function(x) {mean(x<0.05)})

rate_D_50_00 <- apply(res_D_50_00, 1, function(x) {mean(x<0.05)})
rate_D_50_10 <- apply(res_D_50_10, 1, function(x) {mean(x<0.05)})
rate_D_50_30 <- apply(res_D_50_30, 1, function(x) {mean(x<0.05)})
rate_D_50_50 <- apply(res_D_50_50, 1, function(x) {mean(x<0.05)})
rate_D_50_70 <- apply(res_D_50_70, 1, function(x) {mean(x<0.05)})

rate_A_500_00 <- apply(res_A_500_00, 1, function(x) {mean(x<0.05)})
rate_A_500_10 <- apply(res_A_500_10, 1, function(x) {mean(x<0.05)})
rate_A_500_30 <- apply(res_A_500_30, 1, function(x) {mean(x<0.05)})
rate_A_500_50 <- apply(res_A_500_50, 1, function(x) {mean(x<0.05)})
rate_A_500_70 <- apply(res_A_500_70, 1, function(x) {mean(x<0.05)})

rate_B_500_00 <- apply(res_B_500_00, 1, function(x) {mean(x<0.05)})
rate_B_500_10 <- apply(res_B_500_10, 1, function(x) {mean(x<0.05)})
rate_B_500_30 <- apply(res_B_500_30, 1, function(x) {mean(x<0.05)})
rate_B_500_50 <- apply(res_B_500_50, 1, function(x) {mean(x<0.05)})
rate_B_500_70 <- apply(res_B_500_70, 1, function(x) {mean(x<0.05)})

rate_C_500_00 <- apply(res_C_500_00, 1, function(x) {mean(x<0.05)})
rate_C_500_10 <- apply(res_C_500_10, 1, function(x) {mean(x<0.05)})
rate_C_500_30 <- apply(res_C_500_30, 1, function(x) {mean(x<0.05)})
rate_C_500_50 <- apply(res_C_500_50, 1, function(x) {mean(x<0.05)})
rate_C_500_70 <- apply(res_C_500_70, 1, function(x) {mean(x<0.05)})

rate_D_500_00 <- apply(res_D_500_00, 1, function(x) {mean(x<0.05)})
rate_D_500_10 <- apply(res_D_500_10, 1, function(x) {mean(x<0.05)})
rate_D_500_30 <- apply(res_D_500_30, 1, function(x) {mean(x<0.05)})
rate_D_500_50 <- apply(res_D_500_50, 1, function(x) {mean(x<0.05)})
rate_D_500_70 <- apply(res_D_500_70, 1, function(x) {mean(x<0.05)})

y <- 100*c(rate_A_50_00, rate_A_50_10, rate_A_50_30, rate_A_50_50, rate_A_50_70,
       rate_B_50_00, rate_B_50_10, rate_B_50_30, rate_B_50_50, rate_B_50_70,
       rate_C_50_00, rate_C_50_10, rate_C_50_30, rate_C_50_50, rate_C_50_70,
       rate_D_50_00, rate_D_50_10, rate_D_50_30, rate_D_50_50, rate_D_50_70,
      rate_A_500_00, rate_A_500_10, rate_A_500_30, rate_A_500_50, rate_A_500_70,
      rate_B_500_00, rate_B_500_10, rate_B_500_30, rate_B_500_50, rate_B_500_70,
      rate_C_500_00, rate_C_500_10, rate_C_500_30, rate_C_500_50, rate_C_500_70,
      rate_D_500_00, rate_D_500_10, rate_D_500_30, rate_D_500_50, rate_D_500_70)



Method <- rep(c("Test-wise deletion","Multiple imputation"), 40)
q <- rep(c(0, 10, 30, 50, 70), 8, each=2)
n <- rep(c("n=50", "n=500"), each=40)
scenario <- rep(c("Scenario A","Scenario B","Scenario C","Scenario D"), 2,
                each=10)

plot_dat <- data.frame(y, Method, q, n, scenario)

cairo_ps("CurveNull.eps",
    height=4, width=6)
  ggplot(data=plot_dat, aes(x=q, y=y, color=Method, lty=Method)) +
    geom_point() + geom_line() + facet_grid(cols=vars(scenario), rows=vars(n)) +
    xlab("Percent missing") + ylab("Type I error rate (%)") +  expand_limits(x = 0.8) +
    scale_x_continuous(breaks=c(0,10,30,50,70), minor_breaks = NULL) +
    scale_y_continuous(limits=c(0, 10)) +
    theme_bw() +
    theme(legend.position = c(0.75, 0.50),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
dev.off()


