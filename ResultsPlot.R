library(ggplot2)

### load all the summary tables
load("summary_ECOLI.RData")
load("summary_MAGIC.RData")
load("summary_ASIA.RData")
load("summary_SACHS.RData")
load("summary_HEALTHCARE.RData")
load("summary_MEHRA.RData")
load("summary_ECOLI_large.RData")


### ECOLI bar plots
pre_ECOLI_100 <- sum_ECOLI_100_MCAR[c(1,2,3,6,7,8,9,10),"pre.tdr"]
rec_ECOLI_100 <- sum_ECOLI_100_MCAR[c(1,2,3,6,7,8,9,10),"rec.tpr"]
edg_ECOLI_100 <- sum_ECOLI_100_MCAR[c(1,2,3,6,7,8,9,10),"e"]
y_ECOLI_100 <- c(rep(c(-1,1), each=8), pre_ECOLI_100, -rec_ECOLI_100)
method_ECOLI <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle",
                         "MI", "MI random forests", "MI CALIBER",
                         "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI random forests", "MI CALIBER",
                       "Mean imputation")))

ggdat_ECOLI_100 <- data.frame(y=y_ECOLI_100, method=method_ECOLI,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_ECOLI_100 <- data.frame(e=edg_ECOLI_100, method=method_ECOLI)

cairo_ps(file="plot_ECOLI_100.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_ECOLI_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3), c(24,8)), size=3) +
    geom_text(data=edges_ECOLI_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ECOLI")*" (17 edges),"*italic(" n = 100           "),
                            " #E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_ECOLI_5000 <- sum_ECOLI_5000_MCAR[c(1,2,3,6,7,8,9,10),"pre.tdr"]
rec_ECOLI_5000 <- sum_ECOLI_5000_MCAR[c(1,2,3,6,7,8,9,10),"rec.tpr"]
edg_ECOLI_5000 <- sum_ECOLI_5000_MCAR[c(1,2,3,6,7,8,9,10),"e"]
y_ECOLI_5000 <- c(rep(c(-1,1), each=8), pre_ECOLI_5000, -rec_ECOLI_5000)
method_ECOLI <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle",
                         "MI", "MI random forests", "MI CALIBER",
                         "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI random forests", "MI CALIBER",
                       "Mean imputation")))

ggdat_ECOLI_5000 <- data.frame(y=y_ECOLI_5000, method=method_ECOLI,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_ECOLI_5000 <- data.frame(e=edg_ECOLI_5000, method=method_ECOLI)

cairo_ps(file="plot_ECOLI_5000.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_ECOLI_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3), c(24,8)), size=3) +
    geom_text(data=edges_ECOLI_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ECOLI")*" (17 edges),"*italic(" n = 5 000             "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()




### MAGIC bar plots
pre_MAGIC_100 <- sum_MAGIC_100_MCAR[c(1,2,3,6,7,8,9,10),"pre.tdr"]
rec_MAGIC_100 <- sum_MAGIC_100_MCAR[c(1,2,3,6,7,8,9,10),"rec.tpr"]
edg_MAGIC_100 <- sum_MAGIC_100_MCAR[c(1,2,3,6,7,8,9,10),"e"]
y_MAGIC_100 <- c(rep(c(-1,1), each=8), pre_MAGIC_100, -rec_MAGIC_100)
method_MAGIC <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle", "MI",
                         "MI random forests", "MI CALIBER", "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI random forests", "MI CALIBER",
                       "Mean imputation")))

ggdat_MAGIC_100 <- data.frame(y=y_MAGIC_100, method=method_MAGIC,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_MAGIC_100 <- data.frame(e=edg_MAGIC_100, method=method_MAGIC)

cairo_ps(file="plot_MAGIC_100.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_MAGIC_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3,0,-0.3), c(24,1,1,6)), size=3) +
    geom_text(data=edges_MAGIC_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("MAGIC")*" (7 edges),"*italic(" n = 100            "),
                            " #E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_MAGIC_5000 <- sum_MAGIC_5000_MCAR[c(1,2,3,6,7,8,9,10),"pre.tdr"]
rec_MAGIC_5000 <- sum_MAGIC_5000_MCAR[c(1,2,3,6,7,8,9,10),"rec.tpr"]
edg_MAGIC_5000 <- sum_MAGIC_5000_MCAR[c(1,2,3,6,7,8,9,10),"e"]
y_MAGIC_5000 <- c(rep(c(-1,1), each=8), pre_MAGIC_5000, -rec_MAGIC_5000)
method_MAGIC <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle", "MI",
                         "MI random forests", "MI CALIBER", "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI random forests", "MI CALIBER",
                       "Mean imputation")))

ggdat_MAGIC_5000 <- data.frame(y=y_MAGIC_5000, method=method_MAGIC,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_MAGIC_5000 <- data.frame(e=edg_MAGIC_5000, method=method_MAGIC)

cairo_ps(file="plot_MAGIC_5000.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_MAGIC_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3), c(24,8)), size=3) +
    geom_text(data=edges_MAGIC_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("MAGIC")*" (7 edges),"*italic(" n = 5 000              "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()









### ASIA bar plots
pre_ASIA_100 <- sum_ASIA_100_MCAR[c(1,2,3,6,7,8,9,10,11),"pre.tdr"]
rec_ASIA_100 <- sum_ASIA_100_MCAR[c(1,2,3,6,7,8,9,10,11),"rec.tpr"]
edg_ASIA_100 <- sum_ASIA_100_MCAR[c(1,2,3,6,7,8,9,10,11),"e"]
y_ASIA_100 <- c(rep(c(-1,1), each=9), pre_ASIA_100, -rec_ASIA_100)
method_ASIA <- factor(c("Full data", "List-wise deletion", "Test-wise deletion",
                        "MI oracle", "MI", "MI main effects",
                        "MI random forests", "MI CALIBER", "Mode imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI main effects",
                       "MI random forests", "MI CALIBER", "Mode imputation")))

ggdat_ASIA_100 <- data.frame(y=y_ASIA_100, method=method_ASIA,
                              w=rep(c("plus", "minus", "prec", "rec"), each=9))
edges_ASIA_100 <- data.frame(e=edg_ASIA_100, method=method_ASIA)

cairo_ps(file="plot_ASIA_100.eps",
         width=3.8, height=1.9)
  ggplot(ggdat_ASIA_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=18),
              hjust=rep(c(1.3,-0.3,0,-0.3), c(27,1,1,7)), size=3) +
    geom_text(data=edges_ASIA_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ASIA")*" (8 edges),"*italic(" n = 100               "),
                            " #E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_ASIA_5000 <- sum_ASIA_5000_MCAR[c(1,2,3,6,7,8,9,10,11),"pre.tdr"]
rec_ASIA_5000 <- sum_ASIA_5000_MCAR[c(1,2,3,6,7,8,9,10,11),"rec.tpr"]
edg_ASIA_5000 <- sum_ASIA_5000_MCAR[c(1,2,3,6,7,8,9,10,11),"e"]
y_ASIA_5000 <- c(rep(c(-1,1), each=9), pre_ASIA_5000, -rec_ASIA_5000)
method_ASIA <- factor(c("Full data", "List-wise deletion", "Test-wise deletion",
                        "MI oracle", "MI", "MI main effects",
                        "MI random forests", "MI CALIBER", "Mode imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle", "MI", "MI main effects",
                       "MI random forests", "MI CALIBER", "Mode imputation")))

ggdat_ASIA_5000 <- data.frame(y=y_ASIA_5000, method=method_ASIA,
                              w=rep(c("plus", "minus", "prec", "rec"), each=9))
edges_ASIA_5000 <- data.frame(e=edg_ASIA_5000, method=method_ASIA)

cairo_ps(file="plot_ASIA_5000.eps",
         width=3.8, height=1.9)
  ggplot(ggdat_ASIA_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=18),
              hjust=rep(c(1.3,-0.3), c(27,9)), size=3) +
    geom_text(data=edges_ASIA_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ASIA")*" (8 edges),"*italic(" n = 5 000                 "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()







### SACHS bar plots
pre_SACHS_100 <- sum_SACHS_100_MCAR[1:9,"pre.tdr"]
rec_SACHS_100 <- sum_SACHS_100_MCAR[1:9,"rec.tpr"]
edg_SACHS_100 <- sum_SACHS_100_MCAR[1:9,"e"]
y_SACHS_100 <- c(rep(c(-1,1), each=9), pre_SACHS_100, -rec_SACHS_100)
method_SACHS <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle 3-way", "MI 2-way",
                         "MI main effects", "MI random forests", "MI CALIBER",
                         "Mode imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                         "MI oracle 3-way", "MI 2-way", "MI main effects",
                         "MI random forests", "MI CALIBER", "Mode imputation")))

ggdat_SACHS_100 <- data.frame(y=y_SACHS_100, method=method_SACHS,
                              w=rep(c("plus", "minus", "prec", "rec"), each=9))
edges_SACHS_100 <- data.frame(e=edg_SACHS_100, method=method_SACHS)

cairo_ps(file="plot_SACHS_100.eps",
         width=3.8, height=1.9)
  ggplot(ggdat_SACHS_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=18),
              hjust=1.3, size=3) +
    geom_text(data=edges_SACHS_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("SACHS")*" (17 edges)"*italic(", n = 100         "),
                            " #E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_SACHS_5000 <- sum_SACHS_5000_MCAR[1:9,"pre.tdr"]
rec_SACHS_5000 <- sum_SACHS_5000_MCAR[1:9,"rec.tpr"]
edg_SACHS_5000 <- sum_SACHS_5000_MCAR[1:9,"e"]
y_SACHS_5000 <- c(rep(c(-1,1), each=9), pre_SACHS_5000, -rec_SACHS_5000)
method_SACHS <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI oracle 3-way", "MI 2-way",
                         "MI main effects", "MI random forests", "MI CALIBER",
                         "Mode imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                         "MI oracle 3-way", "MI 2-way", "MI main effects",
                         "MI random forests", "MI CALIBER", "Mode imputation")))

ggdat_SACHS_5000 <- data.frame(y=y_SACHS_5000, method=method_SACHS,
                              w=rep(c("plus", "minus", "prec", "rec"), each=9))
edges_SACHS_5000 <- data.frame(e=edg_SACHS_5000, method=method_SACHS)

cairo_ps(file="plot_SACHS_5000.eps",
         width=3.8, height=1.9)
  ggplot(ggdat_SACHS_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=18),
              hjust=rep(c(1.3,-0.3), c(27,9)), size=3) +
    geom_text(data=edges_SACHS_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("SACHS")*" (17 edges)"*italic(", n = 5 000          "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()







### HEALTHCARE bar plots
pre_HEALTHCARE_100 <- sum_HEALTHCARE_100_MCAR[c(1,2,3,4,6,7,8,9),"pre.tdr"]
rec_HEALTHCARE_100 <- sum_HEALTHCARE_100_MCAR[c(1,2,3,4,6,7,8,9),"rec.tpr"]
edg_HEALTHCARE_100 <- sum_HEALTHCARE_100_MCAR[c(1,2,3,4,6,7,8,9),"e"]
y_HEALTHCARE_100 <- c(rep(c(-1,1), each=8), pre_HEALTHCARE_100, -rec_HEALTHCARE_100)
method_HEALTHCARE <- factor(c("Full data", "List-wise deletion",
                              "Test-wise deletion", "MI oracle 2-way",
                              "MI main effects", "MI random forests",
                              "MI CALIBER", "Mean/mode imp."),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle 2-way", "MI main effects",
                       "MI random forests", "MI CALIBER", "Mean/mode imp.")))

ggdat_HEALTHCARE_100 <- data.frame(y=y_HEALTHCARE_100, method=method_HEALTHCARE,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_HEALTHCARE_100 <- data.frame(e=edg_HEALTHCARE_100, method=method_HEALTHCARE)

cairo_ps(file="plot_HEALTHCARE_100.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_HEALTHCARE_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3,0,-0.3), c(24,5,1,2)), size=3) +
    geom_text(data=edges_HEALTHCARE_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic(" HEALTHCARE")*" (9 edges)"*italic(", n = 100 "),
                            "#E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_HEALTHCARE_5000 <- sum_HEALTHCARE_5000_MCAR[c(1,2,3,4,6,7,8,9),"pre.tdr"]
rec_HEALTHCARE_5000 <- sum_HEALTHCARE_5000_MCAR[c(1,2,3,4,6,7,8,9),"rec.tpr"]
edg_HEALTHCARE_5000 <- sum_HEALTHCARE_5000_MCAR[c(1,2,3,4,6,7,8,9),"e"]
y_HEALTHCARE_5000 <- c(rep(c(-1,1), each=8), pre_HEALTHCARE_5000, -rec_HEALTHCARE_5000)
method_HEALTHCARE <- factor(c("Full data", "List-wise deletion",
                              "Test-wise deletion", "MI oracle 2-way",
                              "MI main effects", "MI random forests",
                              "MI CALIBER", "Mean/mode imp."),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI oracle 2-way", "MI main effects",
                       "MI random forests", "MI CALIBER", "Mean/mode imp.")))

ggdat_HEALTHCARE_5000 <- data.frame(y=y_HEALTHCARE_5000, method=method_HEALTHCARE,
                              w=rep(c("plus", "minus", "prec", "rec"), each=8))
edges_HEALTHCARE_5000 <- data.frame(e=edg_HEALTHCARE_5000, method=method_HEALTHCARE)

cairo_ps(file="plot_HEALTHCARE_5000.eps",
         width=3.8, height=1.8)
  ggplot(ggdat_HEALTHCARE_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=16),
              hjust=rep(c(1.3,-0.3), c(24,8)), size=3) +
    geom_text(data=edges_HEALTHCARE_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("HEALTHCARE")*" (9 edges)"*italic(", n = 5 000  "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()







### MEHRA bar plots
pre_MEHRA_100 <- sum_MEHRA_100_MCAR[c(1,2,3,6,7,8,9),"pre.tdr"]
rec_MEHRA_100 <- sum_MEHRA_100_MCAR[c(1,2,3,6,7,8,9),"rec.tpr"]
edg_MEHRA_100 <- sum_MEHRA_100_MCAR[c(1,2,3,6,7,8,9),"e"]
y_MEHRA_100 <- c(rep(c(-1,1), each=7), pre_MEHRA_100, -rec_MEHRA_100)
method_MEHRA <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI main effects",
                         "MI random forests", "MI CALIBER", "Mean/mode imp."),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI main effects", "MI random forests",
                   "MI CALIBER", "Mean/mode imp.")))

ggdat_MEHRA_100 <- data.frame(y=y_MEHRA_100, method=method_MEHRA,
                              w=rep(c("plus", "minus", "prec", "rec"), each=7))
edges_MEHRA_100 <- data.frame(e=edg_MEHRA_100, method=method_MEHRA)

cairo_ps(file="plot_MEHRA_100.eps",
         width=3.8, height=1.65)
  ggplot(ggdat_MEHRA_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=14),
              hjust=1.3, size=3) +
    geom_text(data=edges_MEHRA_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("MEHRA")*" (15 edges)"*italic(", n = 100         "),
                            " #E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_MEHRA_5000 <- sum_MEHRA_5000_MCAR[c(1,2,3,6,7,8,9),"pre.tdr"]
rec_MEHRA_5000 <- sum_MEHRA_5000_MCAR[c(1,2,3,6,7,8,9),"rec.tpr"]
edg_MEHRA_5000 <- sum_MEHRA_5000_MCAR[c(1,2,3,6,7,8,9),"e"]
y_MEHRA_5000 <- c(rep(c(-1,1), each=7), pre_MEHRA_5000, -rec_MEHRA_5000)
method_MEHRA <- factor(c("Full data", "List-wise deletion",
                         "Test-wise deletion", "MI main effects",
                         "MI random forests", "MI CALIBER", "Mean/mode imp."),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                       "MI main effects", "MI random forests",
                   "MI CALIBER", "Mean/mode imp.")))

ggdat_MEHRA_5000 <- data.frame(y=y_MEHRA_5000, method=method_MEHRA,
                              w=rep(c("plus", "minus", "prec", "rec"), each=7))
edges_MEHRA_5000 <- data.frame(e=edg_MEHRA_5000, method=method_MEHRA)

cairo_ps(file="plot_MEHRA_5000.eps",
         width=3.8, height=1.65)
  ggplot(ggdat_MEHRA_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=14),
              hjust=1.3, size=3) +
    geom_text(data=edges_MEHRA_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("MEHRA")*" (15 edges)"*italic(", n = 5 000           "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()







### ECOLI_large bar plots
pre_ECOLI_large_100 <- sum_ECOLI_large_100_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"pre.tdr"]
rec_ECOLI_large_100 <- sum_ECOLI_large_100_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"rec.tpr"]
edg_ECOLI_large_100 <- sum_ECOLI_large_100_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"e"]
y_ECOLI_large_100 <- c(rep(c(-1,1), each=11),
                       pre_ECOLI_large_100, -rec_ECOLI_large_100)
method_ECOLI_large <- factor(c("Full data", "List-wise deletion",
                               "Test-wise deletion", "MI oracle", "MI",
                               "MI random forests", "MI CALIBER", "Two-step A",
                               "Two-step B", "Two-step C", "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                          "MI oracle", "MI", "MI random forests",
                          "MI CALIBER", "Two-step A", "Two-step B",
                         "Two-step C", "Mean imputation")))

ggdat_ECOLI_large_100 <- data.frame(y=y_ECOLI_large_100,
                                    method=method_ECOLI_large,
                              w=rep(c("plus", "minus", "prec", "rec"), each=11))
edges_ECOLI_large_100 <- data.frame(e=edg_ECOLI_large_100,
                                    method=method_ECOLI_large)

cairo_ps(file="plot_ECOLI_large_100.eps",
         width=3.8, height=2.2)
  ggplot(ggdat_ECOLI_large_100, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=22),
              hjust=rep(c(1.3,-0.3), c(32,12)), size=3) +
    geom_text(data=edges_ECOLI_large_100,
              aes(x=method, label=round(e, digit=1),
                  y=-1.1, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray")) +
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1.2,1) +
    theme(axis.text.y = element_text(size=9),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ECOLI_large")*" (70 edges),"*italic(" n = 100  "),
                            "#E    % Recall              % Precision"))) +
   coord_flip(clip="off")
dev.off()



pre_ECOLI_large_5000 <- sum_ECOLI_large_5000_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"pre.tdr"]
rec_ECOLI_large_5000 <- sum_ECOLI_large_5000_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"rec.tpr"]
edg_ECOLI_large_5000 <- sum_ECOLI_large_5000_MCAR[c(1,2,3,6,7,8,9,10,11,12,13),"e"]
y_ECOLI_large_5000 <- c(rep(c(-1,1), each=11),
                        pre_ECOLI_large_5000, -rec_ECOLI_large_5000)
method_ECOLI_large <- factor(c("Full data", "List-wise deletion",
                               "Test-wise deletion", "MI oracle", "MI",
                               "MI random forests", "MI CALIBER", "Two-step A",
                               "Two-step B", "Two-step C", "Mean imputation"),
          levels=rev(c("Full data", "List-wise deletion", "Test-wise deletion",
                          "MI oracle", "MI", "MI random forests",
                          "MI CALIBER", "Two-step A", "Two-step B",
                         "Two-step C", "Mean imputation")))

ggdat_ECOLI_large_5000 <- data.frame(y=y_ECOLI_large_5000,
                                     method=method_ECOLI_large,
                              w=rep(c("plus", "minus", "prec", "rec"), each=11))
edges_ECOLI_large_5000 <- data.frame(e=edg_ECOLI_large_5000,
                                     method=method_ECOLI_large)

cairo_ps(file="plot_ECOLI_large_5000.eps",
         width=3.8, height=2.2)
  ggplot(ggdat_ECOLI_large_5000, aes(x=method, y=y, fill=w)) +
    geom_bar(stat="identity", position="identity", colour="black", width=1) +
    geom_text(aes(label=round(abs(y*100))),
              colour=rep(c("transparent","black"), each=22),
              hjust=rep(c(1.3,-0.3), c(32,12)), size=3) +
    geom_text(data=edges_ECOLI_large_5000,
              aes(x=method, label=round(e, digit=1),
                  y=1.35, hjust=1, fill=NULL),
              size=3) +
    scale_fill_manual(values=c("white","white","gray","gray"))+
    scale_x_discrete(position = "top") +
    theme_void() + ylim(-1,1.35) +
    theme(axis.text.y = element_text(size=9, margin=margin(0,0,0,-3)),
          panel.grid=element_blank(),
          legend.position="none",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=9, margin=margin(0,0,1,0)),
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(expression(atop(italic("ECOLI_large")*" (70 edges),"*italic(" n = 5 000   "),
                            "   % Recall             % Precision            #E"))) +
   coord_flip(clip="off")
dev.off()
