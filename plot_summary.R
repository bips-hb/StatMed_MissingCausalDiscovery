library(ggplot2)

### load summary data
load("summary_ECOLI.RData")
load("summary_MAGIC.RData")
load("summary_ASIA.RData")
load("summary_SACHS.RData")
load("summary_HEALTHCARE.RData")
load("summary_MEHRA.RData")
load("summary_ECOLI_large.RData")

rel <- function(GRAPH, comp_string) {
  su <- eval(parse(text=paste("sum", GRAPH, 100, "MCAR", sep="_")))
  s_100_MCAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 100, "MAR", sep="_")))
  s_100_MAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 100, "MNAR", sep="_")))
  s_100_MNAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  
  su <- eval(parse(text=paste("sum", GRAPH, 1000, "MCAR", sep="_")))
  s_1000_MCAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 1000, "MAR", sep="_")))
  s_1000_MAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 1000, "MNAR", sep="_")))
  s_1000_MNAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  
  su <- eval(parse(text=paste("sum", GRAPH, 5000, "MCAR", sep="_")))
  s_5000_MCAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 5000, "MAR", sep="_")))
  s_5000_MAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  su <- eval(parse(text=paste("sum", GRAPH, 5000, "MNAR", sep="_")))
  s_5000_MNAR <- ( su["-- without correction","h"] -
                    su[comp_string,"h"] ) / su["complete","h"]
  
  return(c(s_100_MCAR, s_100_MAR, s_100_MNAR,
           s_1000_MCAR, s_1000_MAR, s_1000_MNAR,
           s_5000_MCAR, s_5000_MAR, s_5000_MNAR))
}


ECOLI <- rel("ECOLI", "-- linear models")
MAGIC <- rel("MAGIC", "-- linear models")
ASIA <- rel("ASIA", "-- saturated logistic")
SACHS <- rel("SACHS", "-- 2-way logistic")
HEALTHCARE <- rel("HEALTHCARE", "-- CALIBER")
MEHRA <- rel("MEHRA", "-- CALIBER")
ECOLI_large <- rel("ECOLI_large", "-- linear models")

y <- c(ECOLI, MAGIC, ASIA, SACHS, HEALTHCARE, MEHRA, ECOLI_large)
#logy <- (-1)^(y<0) * log(abs(y)+1)

n <- factor(rep(c("100", "1 000", "5 000"), each=3, times=7), levels=rev(c("100", "1 000", "5 000")))
m <- factor(rep(c("MCAR","MAR","MNAR"), 21), levels=rev(c("MCAR","MAR","MNAR")))
d <- factor(rep(c("ECOLI","MAGIC","ASIA","SACHS","HEALTHCARE","MEHRA","ECOLI_large"), each=9))
s <- factor(rep(c("Gaussian","discrete","mixed","Gaussian"), c(18,18,18,9)), levels=c("Gaussian","discrete","mixed"))
l <- factor(rep(c("7-12 nodes","46 nodes"), c(54, 9)))

ggdat <- data.frame(y, n, m, d, s, l)


cairo_ps(file="plot_summary.eps", width=8, height=6)
  ggplot(ggdat, aes(x=m, y=y)) +
    geom_hline(yintercept=0, color="gray50", lwd=1) +
    scale_shape_manual(name="Variable scale",
                       values=c("Gaussian"=21, "discrete"=22, "mixed"=24)) +
    scale_fill_manual(name="Sample size",
                values=c("100"="gray90", "1 000"="gray60", "5 000"="black")) +
    geom_point(size=5, alpha=0.85,
               position=position_nudge(x=-0.25*(n=="100")+0.25*(n=="5 000")+
                              0.03*(s=="Gaussian")-0.03*(s=="mixed")),
               aes(shape=s, fill=n)) +
    scale_colour_manual(name="Graph size",
                        values=c("7-12 nodes"="transparent", "46 nodes"="black")) +
    geom_point(aes(colour=l), size=5, shape=4,
            position=position_nudge(x=-0.25*(n=="100")+0.25*(n=="5 000")+
                              0.03*(s=="Gaussian")-0.03*(s=="mixed"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks = element_blank(),
          axis.text = element_text(size=14),
          legend.position = c(0.6, 0.2),
          legend.box = "horizontal") +
    scale_y_discrete(limits=c(-1.5, 0, 2.5),
            labels=c(expression("" %<-% "        test-wise\n    deletion better"), "0",
                     expression("multiple imputation better    " %->% ""))) +
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    coord_flip()
dev.off()


  
