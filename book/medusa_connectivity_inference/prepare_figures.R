args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
library(ggplot2)
s<-read.csv(file.path(root,'simulation_summary.csv'),check.names=FALSE)
labels<-c(Oracle='Oracle',Naive='Naive',Diagonal='Diagonal moments',
 'Zero shared covariance'='No shared error', 'Full moments, row SE'='Full + row SE',
 'Full moments, cluster SE'='Full + cluster SE',
 'Resolution inverse + full moments'='Undo K + full')
order<-unname(labels)
s$method_label<-factor(unname(labels[s$method]),levels=rev(order))
scenarios<-c('Serial + shared Gaussian','Serial + shared t5','Different resolution')
s$scenario_label<-factor(s$scenario,levels=scenarios,labels=c('Shared / serial Gaussian','Shared / serial t5','Different resolution'))
d<-s[s$scenario%in%scenarios,]
d$term<-factor(d$term,levels=c('Connectivity','PPI'))
d$bias_mcse<-d$empirical_sd/sqrt(d$usable)
ci<-t(vapply(seq_len(nrow(d)),function(i)binom.test(round(d$coverage[i]*d$usable[i]),d$usable[i])$conf.int,numeric(2)))
d$coverage_low<-ci[,1];d$coverage_high<-ci[,2]
theme<-theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
 strip.text=element_text(face='bold'),axis.title.y=element_blank(),
 plot.title=element_text(face='bold',size=14),plot.title.position='plot',
 plot.margin=margin(12,20,12,12))
p1<-ggplot(d,aes(x=bias,y=method_label))+geom_vline(xintercept=0,linetype='dashed',color='#64748B')+
 geom_errorbar(aes(xmin=bias-1.96*bias_mcse,xmax=bias+1.96*bias_mcse),orientation='y',width=.22,color='#2B6CB0')+
 geom_point(size=2.3,color='#2B6CB0')+facet_grid(term~scenario_label,scales='free_x')+
 labs(title='Correcting sampling variance does not correct a biased measurement mean',
 x='Mean coefficient bias (estimate minus truth)')+theme
p2<-ggplot(d,aes(x=coverage,y=method_label))+geom_vline(xintercept=.95,linetype='dashed',color='#64748B')+
 geom_errorbar(aes(xmin=coverage_low,xmax=coverage_high),orientation='y',width=.22,color='#2B6CB0')+
 geom_point(size=2.3,color='#2B6CB0')+facet_grid(term~scenario_label)+
 scale_x_continuous(limits=c(0,1),breaks=c(0,.5,.95),labels=scales::label_percent(accuracy=1))+
 labs(title='Interval coverage depends on the error model and the inference unit',x='Coverage of nominal 95% intervals')+theme
# Materialize plot layouts before opening graphics devices.
p1 <- ggplot_gtable(ggplot_build(p1))
p2 <- ggplot_gtable(ggplot_build(p2))
for(item in list(list(p=p1,file='bias.png'),list(p=p2,file='coverage.png')))
 ggsave(file.path(root,item$file),item$p,width=12,height=7.8,dpi=180,device=ragg::agg_png,bg='white')
write.csv(d,file.path(root,'figure_data.csv'),row.names=FALSE)
