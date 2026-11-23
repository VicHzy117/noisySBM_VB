# The sole public setting/order configuration: current main article + Supplement.
simulation_settings <- function() {
  ids<-paste0("S",1:14)
  family<-c("two_community","two_community","near_duplicate","bipartite","star","pa",
    "two_community","two_community","near_duplicate","bipartite","bipartite","bipartite","star","pa")
  titles<-c("two community, strong connectivity, strong signal","two community, strong connectivity, weak signal",
    "near-duplicate communities, weak signal","bipartite, strong connectivity, weak signal",
    "star, weak signal","preferential attachment, weak signal",
    "two community, weak connectivity, strong signal","two community, weak connectivity, weak signal",
    "near-duplicate communities, strong signal","bipartite, strong connectivity, strong signal",
    "bipartite, weak connectivity, strong signal","bipartite, weak connectivity, weak signal",
    "star, strong signal","preferential attachment, strong signal")
  mu<-c(2,1,NA,1,1,1,2,1,NA,2,2,1,2,2)
  qstrong<-matrix(c(.8,.1,.1,.8),2,2)
  qweak<-matrix(c(.55,.35,.35,.55),2,2)
  qnd<-matrix(c(.75,.12,.12,.12,.58,.57,.12,.57,.58),3,3,byrow=TRUE)
  out<-lapply(1:14,function(i) {
    K<-if(family[i]=="near_duplicate")3L else if(family[i] %in% c("two_community","bipartite"))2L else NA_integer_
    Q<-switch(family[i],two_community=if(i<=2)qstrong else qweak,near_duplicate=qnd,
      bipartite=matrix(c(0,if(i==11 || i==12).25 else .65,if(i==11 || i==12).25 else .65,0),2,2),NULL)
    M<-if(family[i]=="near_duplicate") if(i==3)matrix(c(1,1,1,1,1.5,1,1,1,.5),3,3,byrow=TRUE)
      else matrix(c(2,2,2,2,3,1,2,1,.5),3,3,byrow=TRUE) else if(!is.na(K))matrix(mu[i],K,K)else NULL
    list(setting=ids[i],number=i,title=titles[i],family=family[i],K_true=K,
      proportions=if(!is.na(K))rep(1/K,K)else NULL,Q=Q,mu=M,mu_scalar=mu[i],
      sigma2=if(!is.na(K))matrix(1,K,K)else NULL,null_mean=0,null_variance=1,signal_variance=1,
      K_max=if(i %in% c(1,2,7,8))3L else 4L,
      root_fraction=if(family[i]=="pa").2 else NA_real_,root_density=if(family[i]=="pa").5 else NA_real_,
      attachment_fraction=if(family[i]=="pa").75 else NA_real_)
  })
  names(out)<-ids;out
}
simulation_sizes <- c(50L,100L,200L,400L)
simulation_alphas <- c(.01,.05,.10,.15,.25)
simulation_methods <- c("VB","VEM","Greedy","SC")
