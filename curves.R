#read csv file
#plot FA Along tracts for group 1
#plot FA along tracts for group 2
# stats
library(ggplot2)
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

vec <- vector()
bundlesize3 <- c(vec, length(clus3bundleid))
bundlesize4<-bundlesize3 


for (i in 1:length(clus3bundleid)){
  bundlesize3[i]<-length(which(g3$Bundle.ID==clus3bundleid[i]))
  bundlesize4[i]<-length(which(g4$Bundle.ID==clus4bundleid[i]))
 
}

ggplot(df, aes(x=FA, color=genotype, fill=genotype)) +
  geom_histogram( alpha=0.5, position="identity")+
  geom_density(alpha=0.6)+
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

ggsave(paste(outputpath,'51_57_allbundles_FAHist.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))

ggplot(df, aes(x=Length, color=genotype, fill=genotype)) +
  geom_histogram( alpha=0.5, position="identity")+
  geom_density(alpha=0.6)+
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

ggsave(paste(outputpath,'51_57_allbundles_LengthHist.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)



#darkorchid2
#geom_histogram( fill="white", alpha=0.5, position="identity")+
#geom_histogram( alpha=0.5, position="identity")+

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
  
  ggsave(paste(outputpath,'51_57_bundle', toString(i), '_LengthHist.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
}



