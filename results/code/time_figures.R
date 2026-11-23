draw_time_figures <- function(fits,specs,figures) {
  d<-fits[fits$method %in% c("VB","VEM","Greedy"),]
  colors<-c(VB="#E41A1C",VEM="#377EB8",Greedy="#4DAF4A")
  groups<-split(d,interaction(d[c("setting","n","method")],drop=TRUE))
  boxes<-outliers<-list()
  for(x in groups) {
    v<-x$fit_sec;q<-quantile(v,c(.25,.5,.75),type=7,names=FALSE);iq<-q[3]-q[1]
    lo<-min(v[v>=q[1]-1.5*iq]);hi<-max(v[v<=q[3]+1.5*iq])
    boxes[[length(boxes)+1L]]<-data.frame(x[1,c("setting","n","method")],
      ymin=lo,lower=q[1],middle=q[2],upper=q[3],ymax=hi)
    ov<-v[v<lo | v>hi]
    if(length(ov))outliers[[length(outliers)+1L]]<-data.frame(setting=x$setting[1],n=x$n[1],method=x$method[1],seconds=ov)
  }
  b<-do.call(rbind,boxes);o<-do.call(rbind,outliers)
  # Precompute quartiles/whiskers on raw seconds BEFORE applying the log scale,
  # matching the original manuscript's matplotlib boxplots (type-7 quantiles).
  plot_setting<-function(i) {
    s<-specs[[i]];z<-b[b$setting==s$setting,];z$method<-factor(z$method,levels=names(colors))
    z<-z[order(z$method,z$n),];ov<-o[o$setting==s$setting,]
    parts<-strsplit(s$title,", ",fixed=TRUE)[[1]]
    fam<-switch(parts[1],"two community"="Two communities","bipartite"="Bipartite graph",
      "star"="Star graph",tools::toTitleCase(parts[1]))
    title<-paste0(s$setting,"  |  ",fam,"\n",paste(parts[-1],collapse=", "))
    ggplot(z,aes(n,middle,color=method,group=method))+
      geom_line(linewidth=.65)+
      geom_boxplot(aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,
        fill=method,group=interaction(method,n)),stat="identity",position="identity",width=11,
        linewidth=.5,outlier.shape=NA)+
      geom_point(data=ov,aes(n,seconds,color=method),inherit.aes=FALSE,size=1.2,alpha=.43)+
      scale_color_manual(values=colors,drop=FALSE)+
      scale_fill_manual(values=vapply(colors,function(c)adjustcolor(c,alpha.f=.16),character(1)),drop=FALSE)+
      scale_y_log10(labels=function(v)format(v,trim=TRUE,scientific=FALSE,big.mark=","))+
      scale_x_continuous(breaks=simulation_sizes)+coord_cartesian(xlim=c(23,427))+
      labs(x="Number of nodes, n",y="Fitting time (seconds)",title=title,color=NULL,fill=NULL)+
      theme_classic(base_size=12)+theme(legend.position="bottom",panel.grid.major.y=element_line(color="#dfe4e8",linewidth=.35),
        plot.title=element_text(size=13,lineheight=1.35,margin=margin(b=9)),
        plot.margin=margin(12,18,10,12),axis.line=element_line(color="#89939b"),legend.key.width=grid::unit(16,"pt"))
  }
  pdf(file.path(figures,"fitting_time_s1_s14.pdf"),width=8.45,height=5.4,useDingbats=FALSE)
  for(i in 1:14)print(plot_setting(i));invisible(dev.off())
  for(i in 1:14)ggsave(file.path(figures,sprintf("S%02d_fitting_time.png",i)),plot_setting(i),
    width=8.45,height=5.4,dpi=180,type="cairo",bg="white")
  invisible(b)
}
