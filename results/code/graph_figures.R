draw_graph_figures <- function(summary,specs,figures) {
  methods<-simulation_methods;alphas<-simulation_alphas
  colors<-setNames(c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),methods)
  types<-setNames(c("solid","longdash","dotdash","dotted"),methods)
  shapes<-setNames(c(16,17,15,18,4),sprintf("%.2f",alphas))
  plot_setting<-function(i,combined=FALSE) {
    s<-specs[[i]];d<-summary[summary$setting==s$setting,]
    d$method<-factor(d$method,levels=methods)
    d$alpha_level<-factor(sprintf("%.2f",d$alpha),levels=names(shapes))
    d$n_panel<-factor(d$n,levels=simulation_sizes,labels=paste("n =",simulation_sizes))
    d<-d[order(d$n,d$method,d$alpha),]
    parts<-strsplit(s$title,", ",fixed=TRUE)[[1]]
    family<-tools::toTitleCase(parts[1]);if(parts[1]=="two community")family<-"Two communities"
    title<-paste0(s$setting,". ",family," | ",paste(parts[-1],collapse=", "))
    ggplot(d,aes(FDR_mean,TDR_mean,color=method,group=method))+
      geom_vline(xintercept=alphas,linetype="dashed",color="grey80",linewidth=.28)+
      geom_path(aes(linetype=method),linewidth=if(combined).69 else .80)+
      geom_point(aes(shape=alpha_level),size=if(combined)1.85 else 2.4,stroke=.7)+
      facet_wrap(~n_panel,ncol=2)+scale_color_manual(values=colors,breaks=methods,drop=FALSE)+
      scale_linetype_manual(values=types,breaks=methods,drop=FALSE)+
      scale_shape_manual(values=shapes,limits=names(shapes),drop=FALSE)+
      scale_x_continuous(breaks=seq(0,.3,.1))+scale_y_continuous(breaks=c(0,.5,1))+
      coord_cartesian(xlim=c(-.009,.3),ylim=c(-.045,1.045),expand=FALSE)+
      labs(title=title,x=if(!combined)"Empirical FDR",y=if(!combined)"TDR",
        color=NULL,linetype=NULL,shape="Nominal alpha")+
      guides(color=guide_legend(order=1,override.aes=list(shape=NA)),linetype=guide_legend(order=1),
        shape=guide_legend(order=2,override.aes=list(color="#222222")))+
      theme_bw(base_size=if(combined)9 else 12)+
      theme(legend.position="bottom",legend.box="horizontal",legend.key.width=grid::unit(19,"pt"),
        legend.spacing.x=grid::unit(4,"pt"),legend.text=element_text(size=if(combined)8.5 else 10),
        panel.grid=element_blank(),panel.spacing.y=grid::unit(if(combined)10 else 17,"pt"),
        panel.spacing.x=grid::unit(if(combined)13 else 22,"pt"),
        strip.background=element_rect(fill="#efefef",color=NA),
        strip.text=element_text(size=if(combined)7.7 else 11,margin=margin(2,0,2,0)),
        plot.title=element_text(size=if(combined)9.3 else 13,face="bold",lineheight=1.13,
          margin=margin(b=if(combined)5 else 9)),axis.ticks.length=grid::unit(2,"pt"),
        axis.text=element_text(size=if(combined)7.5 else 10,color="#333333"),
        plot.margin=if(combined)margin(5,12,6,4)else margin(8,14,5,8))
  }
  for(i in 1:14)ggsave(file.path(figures,sprintf("S%02d_TDR_FDR.png",i)),plot_setting(i),
    width=9.2,height=7.4,dpi=220,bg="white",type="cairo")
  pdf(file.path(figures,"TDR_FDR_main_six_settings.pdf"),width=8.05,height=9.20,useDingbats=FALSE)
  grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(5,3,
    widths=grid::unit(c(.23,1,1),c("in","null","null")),
    heights=grid::unit(c(1,1,1,.23,.34),c("null","null","null","in","in")))))
  for(i in 1:6) print(plot_setting(i,TRUE)+theme(legend.position="none"),
    vp=grid::viewport(layout.pos.row=(i+1)%/%2,layout.pos.col=(i-1)%%2+2))
  grid::grid.text("TDR",rot=90,gp=grid::gpar(fontsize=9.6),vp=grid::viewport(layout.pos.row=1:3,layout.pos.col=1))
  grid::grid.text("Empirical FDR",gp=grid::gpar(fontsize=9.6),vp=grid::viewport(layout.pos.row=4,layout.pos.col=2:3))
  legend<-ggplotGrob(plot_setting(1,TRUE))
  key<-legend$grobs[[match("guide-box-bottom",legend$layout$name)]]
  grid::pushViewport(grid::viewport(layout.pos.row=5,layout.pos.col=1:3));grid::grid.draw(key)
  grid::popViewport(2);invisible(dev.off())
}
