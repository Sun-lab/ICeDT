#-------------------------------------------------------------------#
# PLOTTING ABERRANCE                                                #
#-------------------------------------------------------------------#
### Probability of Consistent Marker:
# No Weights:
summary(Testing_Orig$P_Consistent)
summary(c(Testing_Orig$PP_Consistent[which(SimData$Aberrant>0),]))
summary(c(Testing_Orig$PP_Consistent[which(SimData$Aberrant==0),]))

ICeDT_nw_PPAb1  = c(Testing_Orig$PP_Consistent[which(SimData$Aberrant==1),])
ICeDT_nw_PPAb2  = c(Testing_Orig$PP_Consistent[which(SimData$Aberrant==2),])
ICeDT_nw_PPAb3  = c(Testing_Orig$PP_Consistent[which(SimData$Aberrant==3),])
ICeDT_nw_PPnAb = c(Testing_Orig$PP_Consistent[which(SimData$Aberrant==0),])

# Weights
summary(Testing$P_Consistent)
summary(c(Testing$PP_Consistent[which(SimData$Aberrant>0),]))
summary(c(Testing$PP_Consistent[which(SimData$Aberrant==0),]))

ICeDT_w0_PPAb1  = c(Testing$PP_Consistent[which(SimData$Aberrant==1),])
ICeDT_w0_PPAb2  = c(Testing$PP_Consistent[which(SimData$Aberrant==2),])
ICeDT_w0_PPAb3  = c(Testing$PP_Consistent[which(SimData$Aberrant==3),])
ICeDT_w0_PPnAb = c(Testing$PP_Consistent[which(SimData$Aberrant==0),])

# Posterior Probability (Aberrant):
Ab1_Plot_Data = data.frame(Model=c(rep("ICeD-T (Weights)",length(ICeDT_w0_PPAb1))),
                           PostProb = c(ICeDT_w0_PPAb1))

Ab2_Plot_Data = data.frame(Model=c(rep("ICeD-T (Weights)",length(ICeDT_w0_PPAb2))),
                           PostProb = c(ICeDT_w0_PPAb2))

Ab3_Plot_Data = data.frame(Model=c(rep("ICeD-T (Weights)",length(ICeDT_w0_PPAb3))),
                           PostProb = c(ICeDT_w0_PPAb3))

# Plots:
Ab1plot <- ggplot(data = Ab1_Plot_Data,aes(PostProb))+geom_density(alpha=0.2)+
  labs(y="Density",x="Posterior Prob. Consistent",title="Down-Regulated")

Ab2plot <- ggplot(data = Ab2_Plot_Data,aes(PostProb))+geom_density(alpha=0.2)+
  labs(y="Density",x="Posterior Prob. Consistent",title="Up-Regulated")

Ab3plot <- ggplot(data = Ab3_Plot_Data,aes(PostProb))+geom_density(alpha=0.2)+
  labs(y="Density",x="Posterior Prob. Consistent",title="Tumor Expression")

lp <- get_legend(nAbPlot)

MultiPlot = plot_grid(Abplot+theme(legend.position="none"),nAbPlot+theme(legend.position="none"),PPlot+theme(legend.position="none"),NULL,lp,NULL,rel_widths = c(1,1,1),rel_heights = c(1,0.2),
                      nrow=2,ncol=3,labels = c("a","b","c","","",""))
