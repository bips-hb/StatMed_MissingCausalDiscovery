library(graph) #ugraph
library(pcalg)
library(bnlearn) # Hamming distance
library(mice)
library(parallel)
library(ggplot2)

load("res_TilePlot.RData")


## prepare results for plotting

mean_vsparse_vweak <- apply(res_vsparse_vweak, 1, mean)
mean_vsparse_weak <- apply(res_vsparse_weak, 1, mean)
mean_vsparse_medium <- apply(res_vsparse_medium, 1, mean)
mean_vsparse_strong <- apply(res_vsparse_strong, 1, mean)
mean_vsparse_vstrong <- apply(res_vsparse_vstrong, 1, mean)

mean_sparse_vweak <- apply(res_sparse_vweak, 1, mean)
mean_sparse_weak <- apply(res_sparse_weak, 1, mean)
mean_sparse_medium <- apply(res_sparse_medium, 1, mean)
mean_sparse_strong <- apply(res_sparse_strong, 1, mean)
mean_sparse_vstrong <- apply(res_sparse_vstrong, 1, mean)

mean_medium_vweak <- apply(res_medium_vweak, 1, mean)
mean_medium_weak <- apply(res_medium_weak, 1, mean)
mean_medium_medium <- apply(res_medium_medium, 1, mean)
mean_medium_strong <- apply(res_medium_strong, 1, mean)
mean_medium_vstrong <- apply(res_medium_vstrong, 1, mean)

mean_dense_vweak <- apply(res_dense_vweak, 1, mean)
mean_dense_weak <- apply(res_dense_weak, 1, mean)
mean_dense_medium <- apply(res_dense_medium, 1, mean)
mean_dense_strong <- apply(res_dense_strong, 1, mean)
mean_dense_vstrong <- apply(res_dense_vstrong, 1, mean)

mean_vdense_vweak <- apply(res_vdense_vweak, 1, mean)
mean_vdense_weak <- apply(res_vdense_weak, 1, mean)
mean_vdense_medium <- apply(res_vdense_medium, 1, mean)
mean_vdense_strong <- apply(res_vdense_strong, 1, mean)
mean_vdense_vstrong <- apply(res_vdense_vstrong, 1, mean)


means <- rbind(mean_vsparse_vweak, mean_vsparse_weak, mean_vsparse_medium,
               mean_vsparse_strong, mean_vsparse_vstrong,
               mean_sparse_vweak, mean_sparse_weak, mean_sparse_medium,
               mean_sparse_strong, mean_sparse_vstrong,
               mean_medium_vweak, mean_medium_weak, mean_medium_medium,
               mean_medium_strong, mean_medium_vstrong,
               mean_dense_vweak, mean_dense_weak, mean_dense_medium,
               mean_dense_strong, mean_dense_vstrong,
               mean_vdense_vweak, mean_vdense_weak, mean_vdense_medium,
               mean_vdense_strong, mean_vdense_vstrong)

y <- apply(means, 1, function(x) { (x[2]-x[3]) / x[1]})
d <- factor(rep(c("very\nsparse","sparse","medium","dense","very\ndense"),
                each=5),
            levels=c("very\nsparse","sparse","medium","dense","very\ndense"))
w <- factor(rep(c("very\nweak","weak","medium","strong","very\nstrong"), 5),
            levels=c("very\nweak","weak","medium","strong","very\nstrong"))

ggdat <- data.frame(y, d, w)

levels(ggdat$d) <- c("very\nsparse","sparse","medium","dense","very\ndense")
levels(ggdat$w) <- c("very\nweak","weak","medium","strong","very\nstrong")

cairo_ps(file="TilePlot.eps",
    width=4.8, height=3.5)
  ggplot(data=ggdat, aes(x=d, y=w, fill=y)) + geom_tile() +
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) +
    xlab("Edge density") + ylab("Edge strength") +
    scale_fill_gradient2(name=NULL,
                         low=rgb(255-23, 255-99, 255-170, maxColorValue = 255),
                         mid="white",
                         high=rgb(23, 99, 170, maxColorValue = 255),
                         breaks=c(-0.5,0,0.5),
                         limits=c(-0.5,0.5),
                         labels=c("test-wise\ndeletion\nbetter","0",
                                  "multiple\nimputation\nbetter"),
                         guide=guide_colourbar(ticks.colour="transparent",
                                               barheight=7,
                                               barwidth=0.5)) +
    theme_bw() + 
    theme(panel.grid=element_blank(), axis.ticks=element_blank())
    
dev.off()



################################################################################
### This is how I chose the edge coefficients:

# the sample size was the same: n=500 for every trial such that multiple
# imputation is potentially helpful. The coefficients were chosen such that
# the power for detecting a marginal dependence was 10, 30, 50, 70 or 90%.

pow <- function(r) {
  A <- rnorm(500)
  B <- r*A + rnorm(500)
  dat <- data.frame(A,B)
  pval <- gaussCItest(x=1, y=2, S=NULL, suffStat=list(C=cor(dat), n=500))
  return(pval < 0.05)
}

res <- replicate(10000, pow(0.111))
mean(res)

## 0.03 for a power of 10%

## 0.065 for a power of 30%

## 0.088 for a power of 50%

## 0.111 for a power of 70%

## 0.145 for a power of 90%
