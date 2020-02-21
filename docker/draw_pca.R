Draw.PCA<-function(results=NA,
                   outfile=NA,
                   res=300,
                   width=12,
                   height=12,
                   group.colors=c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")){

  # Make PCA on normalized count
  
  
  # perform PCA and make PCA plots
  pca=prcomp(t(results$norm.mtx))  # compute PCA
  pca_df=merge(data.frame(pca$x), results$annotations, by=0)
  pca_plot<-ggplot(pca_df, aes(x=PC1, y=PC2, 
                               color=Group, label=Sample_ID))+
    geom_point(size=3)+ scale_color_manual(values=group.colors)
 
  
  if(!(is.na(outfile))){
    ggsave(outfile, plot = pca_plot,
         width = width, height = height, unit="in",
         dpi = res, limitsize = TRUE, device="png")
  }else{return(pca_plot)}
}
