library(ggplot2)
set.seed(2022)
t <- seq(0, 32, .1)
dt <- diff(t)[1]
h <- dgamma(t, 6) - .15 * dgamma(t, 15)
h <- h / sum(h)
a <- data.frame(time=t, Decision=as.numeric(t>=3 & t<5 | t>=12 & t<16),
                Feedback=as.numeric(abs(t-5.5)<.055 | abs(t-16.5)<.055) * 8)
forward <- function(x) convolve(x, rev(h), type="open")[seq_along(t)]
d <- forward(a$Decision)
f <- forward(a$Feedback)
z <- rbind(data.frame(time=t, layer="Event sources", signal=c(a$Decision, a$Feedback/8), series=rep(c("Decision","Feedback"),each=length(t))),
           data.frame(time=t, layer="Predicted BOLD contributions", signal=c(d,f), series=rep(c("Decision","Feedback"),each=length(t))),
           data.frame(time=t, layer="Observed BOLD", signal=d+f+as.numeric(arima.sim(list(ar=.82),n=length(t),sd=.017)),series="Observed"))
z$layer <- factor(z$layer,c("Event sources","Predicted BOLD contributions","Observed BOLD"))
p <- ggplot(z,aes(time,signal,color=series)) + geom_line(linewidth=.95) +
  facet_wrap(~layer,ncol=1,scales="free_y") +
  scale_color_manual(values=c(Decision="#22757A",Feedback="#B35530",Observed="#223B52")) +
  labs(x="Seconds from run start",y=NULL,color=NULL) +
  theme_minimal(base_size=15) + theme(legend.position="bottom",panel.grid.minor=element_blank(),
     strip.text=element_text(face="bold"),plot.margin=margin(10,17,10,12))
ggsave("assets/forward-model.png",p,width=8.3,height=5.4,dpi=175,bg="white")
