library(ggplot2)
d <- read.csv("calibration_summary.csv")
d$design <- ifelse(startsWith(d$scenario, "Complete"), "Complete", "Gaps + serial dependence")
theme_set(theme_minimal(base_size=13) + theme(
  legend.position="bottom",
  plot.title.position="plot",
  plot.title=element_text(face="bold", size=13.6, color="#173145"),
  plot.subtitle=element_text(size=10.5, color="#495E6E"),
  panel.grid.minor=element_blank(),
  plot.margin=margin(12,18,12,12)))
make <- function(file, labels, goal, span, title, subtitle, x) {
  z <- d[match(d$metric,names(labels),nomatch=0)>0,]
  z$test <- factor(labels[z$metric], rev(unname(labels)))
  p <- ggplot(z,aes(rate,test,color=design,shape=design)) +
    geom_vline(xintercept=goal,color="#6B7D89",linetype=3) +
    geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.19,
      position=position_dodge(width=.58),linewidth=.75) +
    geom_point(position=position_dodge(width=.58),size=2.7) +
    scale_x_continuous(labels=scales::label_percent(accuracy=1),limits=span) +
    scale_color_manual(values=c("#245D80","#B75126")) +
    labs(title=title,subtitle=subtitle,x=x,y=NULL,color=NULL,shape=NULL)
  ggsave(file,p,width=9.2,height=if(length(labels)>3) 4.6 else 4.0,dpi=180,bg="white")
}
make("rejection.png",c(null_maxT="Single trajectory: Max-T",null_tfce="Single trajectory: 1D TFCE",
  null_omnibus="Single trajectory: Omnibus",raw_epoch_family="3-parcel family: Unadjusted",
  holm_epoch_family="3-parcel family: Holm FWER"),.05,c(0,.23),
  "Empirical Type I Error: Single-Trajectory vs. Multi-Parcel Families",
  "180 synthetic cohorts (S = 20 participants); dotted line marks 5% nominal rate",
  "Empirical false-positive rate")
make("coverage.png",c(null_band="Simultaneous Max-T band (null)",active_band="Simultaneous Max-T band (active)",
  window_coverage="Planned window ([3,4] - [0,1] s)"),.95,c(.84,1),
  "Simultaneous 95% Trajectory Band and Window-Contrast Coverage",
  "180 synthetic cohorts (S = 20 participants); dotted line marks 95% nominal coverage",
  "Empirical 95% coverage")
