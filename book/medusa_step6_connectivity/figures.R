library(ggplot2)
d <- read.csv("ppi_draws.csv")
d$scenario <- factor(d$scenario, c("Genuine PPI","Shared-error null"))
d$method <- factor(d$method, c("Naive","M-MOMS","Known-error moments"))
theme_set(theme_minimal(base_size=14) + theme(
  legend.position="none", plot.title=element_text(face="bold",size=17),
  strip.text=element_text(face="bold"), panel.grid.minor=element_blank(),
  plot.margin=margin(12,16,12,12)))
ref <- unique(d[c("scenario","truth")])
p <- ggplot(d,aes(estimate, method, fill=method)) +
  geom_vline(data=ref,aes(xintercept=truth),linetype=3,color="#637887") +
  geom_boxplot(width=.5,outlier.alpha=.4) + facet_wrap(~scenario,ncol=1) +
  scale_fill_manual(values=c("#CAD8E2","#73A7C3","#82C0AD")) +
  labs(title="A predictor correction can amplify shared error",
       subtitle="Dotted lines mark the generating PPI slopes",
       x="Estimated population PPI slope",y=NULL)
ggsave("ppi_estimates.png",p,width=8.4,height=4.7,dpi=175,bg="white")
g <- data.frame(assumption=seq(0,.95,length.out=101))
g$estimate <- .7^2*(.65-g$assumption)
p <- ggplot(g,aes(assumption,estimate)) +
  geom_hline(yintercept=0,color="#637887",linetype=3) +
  geom_vline(xintercept=.65,color="#637887",linetype=3) +
  geom_line(color="#245D80",linewidth=1.2) +
  geom_point(data=data.frame(assumption=.65,estimate=0),color="#B75126",size=3) +
  labs(title="Shared-error assumptions can change the conclusion",
       subtitle="Population sensitivity calculation for the shared-error null",
       x="Assumed amplitude of task-dependent error correlation",
       y="Corrected PPI slope")
ggsave("ppi_sensitivity.png",p,width=8.4,height=3.7,dpi=175,bg="white")
