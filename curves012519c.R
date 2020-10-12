#%prediction versus confidence intervals
#https://stats.stackexchange.com/questions/16493/difference-between-confidence-intervals-and-prediction-intervals

#https://janhove.github.io/reporting/2017/05/12/visualising-models-2
#https://janhove.github.io/analysis/2018/07/04/bayesian-breakpoint-model
#nice plots
#https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
#http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

#read csv file
#plot FA Along tracts for group 1
#plot FA along tracts for group 2
# stats
library(ggplot2)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) # need to run every time you start R and want to use %>%
library(dplyr)

  #confidence intervals inspired from https://janhove.github.io/reporting/2017/05/12/visualising-models-2
  #Note that this confidence band is based solely on the fixed-effects. As such, it should be taken to reflect the uncertainty about the average participant’s regression curve;
  #the regression curves of individual participants can differ wildly from this average.

#for confidence intervals on the real data
#https://stackoverflow.com/questions/51993044/plot-time-series-with-ggplot-with-confidence-interval
mean_cl_quantile <- function(x, q = c(0.1, 0.9), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

g3file='/Users/alex/AlexBadea_MyPapers/E3E4_connectome_manuscript/51--57group3FA.csv'
g4file='/Users/alex/AlexBadea_MyPapers/E3E4_connectome_manuscript/51--57group4FA.csv'
outputpath='/Users/alex/AlexBadea_MyPapers/E3E4_connectome_manuscript/groupcomp/'

g3<-read.csv(g3file, header = TRUE)
g4<-read.csv(g4file, header = TRUE)

df3<-data.frame(g3)
df3<-cbind(df3,genotype='3')

df4<-data.frame(g4)
df4<-cbind(df4,genotype='4')

df<-rbind(df3,df4)

#these csv files only contain the top 6 clusters
#see how big are those clusters
clus3bundleid<-unique(g3$Bundle.ID)
clus4bundleid<-unique(g4$Bundle.ID)

length(clus3bundleid)
length(clus4bundleid)

vec <- vector()
bundlesize3 <- c(vec, length(clus3bundleid))
bundlesize4<-bundlesize3 


for (i in 1:length(clus3bundleid)){
  bundlesize3[i]<-length(which(g3$Bundle.ID==clus3bundleid[i]))
  bundlesize4[i]<-length(which(g4$Bundle.ID==clus4bundleid[i]))
 
}

#for all bundles

#https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
# %>%passes the previous dataframe into the next function as that functions df argument
muFA <- df %>% group_by(genotype) %>% summarize(grp.mean = mean(FA))
muFApointy <- df %>% group_by(genotype,Point.ID) %>% summarize(grp.meanpoint = mean(FA))
CIlowFApointy <- df %>% group_by(genotype,Point.ID) %>% summarize(grp.meanpoint = mean(FA))

ggplot(df, aes(x=FA, color=genotype, fill=genotype)) +
  geom_histogram( alpha=0.5, position="identity")+
  geom_density(alpha=0.6)+
  geom_vline(data = muFA, aes(xintercept = grp.mean, color = genotype), linetype = "dashed")+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())

ggsave(paste(outputpath,'51_57_allbundles_FAHistLines.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))

#https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
muLength <- df %>% group_by(genotype) %>% summarize(grp.mean = mean(Length))

ggplot(df, aes(x=Length, color=genotype, fill=genotype)) +
  geom_histogram( alpha=0.5, position="identity")+
  geom_density(alpha=0.6)+
  geom_vline(data = muLength, aes(xintercept = grp.mean, color = genotype), linetype = "dashed")+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())

ggsave(paste(outputpath,'51_57_allbundles_LengthHistLines.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

###################
####################
#####################



for (i in 1:length(clus3bundleid)){
dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
dftemp<-rbind(dftemp3,dftemp4)

ggplot(dftemp, aes(x=FA, color=genotype, fill=genotype)) +
  geom_histogram( alpha=0.5, position="identity")+
  geom_density(alpha=0.6)+
  ggtitle(paste('Bundle/Cluster', toString(i), sep= ' '))+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())

ggsave(paste(outputpath,'51_57_bundle', toString(i), '_FAHist.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
}



for (i in 1:length(clus3bundleid)){
  dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
  dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
  dftemp<-rbind(dftemp3,dftemp4)
  
  ggplot(dftemp, aes(x=Length, color=genotype, fill=genotype)) +
   geom_histogram( alpha=0.5, position="identity")+
   # geom_density(alpha=0.6)+
    geom_density(aes(color = genotype), size =1, alpha=0.6)+
    ggtitle(paste('Bundle/Cluster', toString(i), sep= ' '))+
    theme_classic()+
    scale_color_manual(values=c('chartreuse1','blueviolet'))+
    scale_fill_manual(values=c('chartreuse1','blueviolet'))+
    theme(legend.position="top")+
    theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
          axis.text.y = element_text(face="bold", size=14, angle=0),
          axis.line.x = element_line(colour = 'black', size=0.5),
          axis.line.y = element_line(colour = 'black', size=0.5),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank())
  
  ggsave(paste(outputpath,'51_57_bundle', toString(i), '_LengthHist.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
}

#write these results to tables
for (i in 1:length(clus3bundleid)){
  dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
  dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
  dftemp<-rbind(dftemp3,dftemp4)
  ksbundleFA<-ks.test(dftemp3$FA, dftemp4$FA, alternative = c("two.sided"),exact = NULL)
  ksbundleLength<-ks.test(dftemp3$Length, dftemp4$Length, alternative = c("two.sided"),exact = NULL)

#plot signal along fibers

dfstream3Bundlei<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
dfstream4Bundlei<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
dfstreamBundlei<-rbind(dfstream3Bundlei,dfstream4Bundlei)
dfstreamBundlei$genotypenorm <- ((-1)*as.numeric(dfstreamBundlei$genotype) - 3.5)

#this models does not work very well, but let's explore and make it better later
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#is-the-likelihood-ratio-test-reliable-for-mixed-models
fm1 <- lmer(
  formula = FA ~ poly(Point.ID,degree=4)+genotype+ (Point.ID|Streamlines.ID), data = dfstreamBundlei
  #+ , 
)
cV <- ranef(fm1, condVar = TRUE)   
ranvar <- attr(cV[[1]], "postVar")



FApred <- predict(fm1,data=dfstreamBundlei,re.form=NA)


ggplot(dfstreamBundlei, aes(y=FApred, x=Point.ID, color=genotype)) +
  geom_point()+
  #geom_line(aes(group = genotype), alpha = .3) +
  #geom_line(data = dfstreamBundlei, alpha = .8, size = 1) +
  stat_smooth() +
  #geom_histogram( alpha=0.5, position="identity")+
  #geom_density(alpha=0.6)+
  ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())
ggsave(paste(outputpath,'51_57_pred_allbundle', toString(i), '_FA.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)


dfstreamBundlei$lwr=dfstreamBundlei$FA-1.96*stderr(dfstreamBundlei$FA)
dfstreamBundlei$hi=dfstreamBundlei$FA+1.96*stderr(dfstreamBundlei$FA)

ggplot(dfstreamBundlei, aes(y=FA, x=Point.ID, color=genotype)) +
  #geom_smooth(stat = 'summary', group=genotype, fun.data = function(y) data.frame(ymin = quantile(y, .05), y = mean(y), ymax = quantile(y, .95))) +
  #geom_smooth(stat = 'summary', fun.data = function(y) data.frame(ymin = quantile(y, .025), y = mean(y), ymax = quantile(y, .975))) +
  #geom_smooth(stat = 'summary', fun.data = mean_cl_quantile)
  #stat_smooth() +
  geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
  ##bad crashes geom_smooth(aes(ymin = lwr, ymax = hi,color=genotype),stat = "identity")+
  ##bad crashes geom_ribbon(data=dfstreamBundlei,aes(ymin=lwr,ymax=hi),alpha=0.3)+
  #aes(ymin=lwr,ymax=hi,fill=genotype, colour=genotype),stat = "identity")+
  #geom_histogram( alpha=0.5, position="identity")+
  #geom_density(alpha=0.6)+
  ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())
ggsave(paste(outputpath,'51_57_alongbundleCI95', toString(i), '_FA.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)



#%start save stats per bundle: FA and Length
mytFA<-t.test(FA ~ genotype, data = dfstreamBundlei)
mytLength<-t.test(Length ~ genotype, data = dfstreamBundlei)

mytableFA<-as_data_frame(
  cbind(mytFA$data.name, i, mytFA$estimate[1], mytFA$estimate[2], mytFA$statistic, mytFA$p.value,  
        mytFA$conf.int[1], mytFA$conf.int[2] , mytFA$parameter)
  )
mytableLength<-as_data_frame(
  cbind(mytLength$data.name, i, mytLength$estimate[1], mytLength$estimate[2], mytLength$statistic, mytLength$p.value,  
        mytLength$conf.int[1], mytLength$conf.int[2] , mytLength$parameter)
)
                      # col.names=c('contrast', names(mytFA$estimate)[1],names(mytFA$estimate)[2],names(mytFA$statistic),'pvalue', 'CIlwr','CIhi', 'df'))
mycolnames<-c('contrast', 'bundle', names(mytFA$estimate)[1], names(mytFA$estimate)[2], names(mytFA$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df')

myfile<-paste(outputpath,'51_57_bundle', toString(i), '_stats.csv',sep='')
#add column names for teh first instance
write.table(mytableFA, file=myfile,col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#remove columnnames for subsequent rows
write.table(mytableLength, file=myfile, sep = "," , row.names = F,append=TRUE,col.names=FALSE)
#%end save stats per bundle: FA and Length


clus3streamid<-unique(dfstream3Bundlei$Streamlines.ID)
clus4streamid<-unique(dfstream4Bundlei$Streamlines.ID)
#by predict - this is your ?vbit Alex
m03 <- lm(FA~poly(Point.ID,4)*genotype, data=dfstreamBundlei)
newdf<-expand.grid(Point.ID <- seq(min(dfstream3Bundlei$ Point.ID), max(dfstream3Bundlei$Point.ID), length.out=100))

dfstreamBundlei$LoCI<-predict(m03,data=dfstreamBundlei,interval='confidence', level=0.95)[,2]
dfstreamBundlei$HiCI<-predict(m03,data=dfstreamBundlei,interval='confidence', level=0.95)[,3]


ggplot(dfstreamBundlei, aes(y=m03$fitted.values, x=Point.ID, color=genotype)) +
  geom_smooth(aes(ymin=LoCI,ymax=HiCI,fill=genotype, colour=genotype),stat = "identity")+
  ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
  theme_classic()+
  scale_color_manual(values=c('chartreuse1','blueviolet'))+
  scale_fill_manual(values=c('chartreuse1','blueviolet'))+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0), 
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())



ggsave(paste(outputpath,'51_57_bundle', toString(i), '_FAfitCI.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

}



# #start this bit not used
# m01 <- gls(FA~poly(Point.ID,4), data=dfstream3Bundlei)
# m02 <- lm(FA~poly(Point.ID,4), data=dfstream3Bundlei)
# newx <- seq(min(dfstream3Bundlei$Point.ID), max(dfstream3Bundlei$Point.ID), length.out=100)
# fit2<-predict(m02, newdata=data.frame(Point.ID=newx,interval='confidence'))
# #by hand
# V<-vcov(m02)
# X<-model.matrix(~poly(Point.ID,3),data=dfstream3Bundlei )
# se.fit<-sqrt(diag(X%*%V%*%t(X)))
# lwr<-fit2-1.96*se.fit
# upr<-fit2+1.96*se.fit
# #alpha<-0.95
# #z.val <- qnorm(1 - (1 - alpha)/2)
# #end this bit not used



#look into how good a model we fit

# fa1<-m03$fitted.values
# fa2<-(m03$model)$FA
# genotype<-(m03$model)$genotype
# dffit<-data_frame(fa1,fa2,genotype)
# ggplot(dffit, aes(y=fa1, x=fa2, color=genotype)) +
#   geom_line()+
#   geom_point()
# #end look at how good a model we fit