library(ggplot2)
d <- read.csv("calibration_summary.csv")
d$design <- ifelse(startsWith(d$scenario, "Complete"), "Complete", "Gaps + serial dependence")
theme_set(theme_minimal(base_size=14) + theme(
  legend.position="bottom", plot.title=element_text(face="bold",size=17),
  panel.grid.minor=element_blank(), plot.margin=margin(12,17,12,12)))
make <- function(file, labels, goal, span, title, x) {
  z <- d[match(d$metric,names(labels),nomatch=0)>0,]
  z$test <- factor(labels[z$metric], rev(unname(labels)))
  p <- ggplot(z,aes(rate,test,color=design,shape=design)) +
    geom_vline(xintercept=goal,color="#6B7D89",linetype=3) +
    geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.19,
      position=position_dodge(width=.58),linewidth=.75) +
    geom_point(position=position_dodge(width=.58),size=2.7) +
    scale_x_continuous(labels=scales::label_percent(accuracy=1),limits=span) +
    scale_color_manual(values=c("#245D80","#B75126")) +
    labs(title=title,x=x,y=NULL,color=NULL,shape=NULL)
  ggsave(file,p,width=8.4,height=if(length(labels)>3) 4.5 else 3.9,dpi=175,bg="white")
}
make("rejection.png",c(null_maxT="Null Max-|T|",null_tfce="Null TFCE",
  null_omnibus="Null omnibus",raw_epoch_family="Study omnibus, unadjusted",
  holm_epoch_family="Study omnibus, Holm"),.05,c(0,.23),
  "The study family changes false-positive rates","Empirical false-positive rate")
make("coverage.png",c(null_band="Raw band, null curve",active_band="Raw band, active curve",
  window_coverage="Planned window"),.95,c(.84,1),
  "Band and window coverage answer different questions","Empirical coverage")
