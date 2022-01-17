



AUCs_XGB_StrSampling <-  rbind(
  ci(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP, levels = c(0,1), direction = "<"))
)


AUCs_XGB_StrSampling[,c(1,2,3)] <- AUCs_XGB_StrSampling[,c(2,1,3)]
colnames(AUCs_XGB_StrSampling) <- c("AUC", "CI_low", "CI_High")



AUCs_XGB_StrSampling <- as.data.frame(AUCs_XGB_StrSampling)

AUCs_XGB_StrSampling$model_type <- c("Mechanistic", "Mechanistic", "Agnostic", "Agnostic")
AUCs_XGB_StrSampling$data_type <- c("Training", "Testing", "Training", "Testing")
AUCs_XGB_StrSampling$algorithm <- rep("XGB", 4)
AUCs_XGB_StrSampling$approach <- rep("Str.Sampling", 4)  


##  Plot the figure
#pd <- position_dodge(0.5)

#ggplot(AUCs_XGB_StrSampling, aes(x=AUC, y = model_type, group = data_type)) +
#  geom_point(position=pd) +
#  geom_errorbar(data=AUCs_XGB_StrSampling, aes(xmin=CI_low, xmax=CI_High, 
#                                color=model_type), width = 0.5, position=pd) 



###########
p = ggplot(data=AUCs_XGB_StrSampling,
           aes(x = data_type,y = AUC, ymin = CI_low, ymax = CI_High))+
  geom_pointrange(aes(col=data_type))+
  geom_hline(aes(fill=data_type),yintercept =0.5, linetype=2)+
  xlab('model_type')+ ylab("AUC (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_High,col=data_type),width=0.5,cex=1)+ 
  facet_wrap(~model_type,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p





