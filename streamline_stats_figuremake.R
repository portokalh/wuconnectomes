#read csv file
#plot FA Along tracts for group 1
#plot FA along tracts for group 2
# stats
library(ggplot2)
library(readxl)

ratio = 1
if (ratio == 1) {ratio_str = 'all'} else {ratio_str = paste('ratio_',as.character(ratio),sep='')}
if (ratio == 1) {folder_str = ''} else {folder_str = paste('_',as.character(ratio))}

csvpath_folder = paste('/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Statistics_MDT_non_inclusive_symmetric',as.character(folder_str),'/',sep="")

#g3file = 'Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Statistics_MDT_non_inclusive_symmetric_100/APOE3_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_ratio_100_bundle_stats_fa.csv'
#g4file = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Statistics_MDT_non_inclusive_symmetric_100/APOE4_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_ratio_100_bundle_stats_fa.csv'

regions = list(22, 1)
regions_str = tolower(paste(my_atlas[regions[[1]],3],'_',my_atlas[regions[[1]],4],'_to_',my_atlas[regions[[2]],3],'_',my_atlas[regions[[2]],4],sep=""))

atlas_file='/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx'
my_atlas <- read_excel(atlas_file)

g3file = paste(csvpath_folder,'APOE3_',regions_str,'_',ratio_str,'_bundle_stats_fa.csv',sep="")
g4file = paste(csvpath_folder,'APOE4_',regions_str,'_',ratio_str,'_bundle_stats_fa.csv',sep="")

outputpath='/Users/jas/jacques/AD_decode_abstract/R_figures/'

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

ggplot(df, aes(x=fa, color=genotype, fill=genotype)) +
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

ggsave(paste(outputpath,regions_str,'_allbundles_FAHist','_',ratio_str,'.pdf',sep=''), plot = last_plot(), device = 'pdf',
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

ggsave(paste(outputpath,regions_str,'_allbundles_LengthHist','_',ratio_str,'.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)



#darkorchid2
#geom_histogram( fill="white", alpha=0.5, position="identity")+
#geom_histogram( alpha=0.5, position="identity")+

for (i in 1:length(clus3bundleid)){
  dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
  dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
  dftemp<-rbind(dftemp3,dftemp4)
  
  ggplot(dftemp, aes(x=fa, color=genotype, fill=genotype)) +
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
  
  ggsave(paste(outputpath,regions_str,'_bundle', toString(i), '_FAHist','_',ratio_str,'.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
}

for (i in 1:length(clus3bundleid)){
  dftemp3<-subset(df, (genotype=='3' & Bundle.ID==clus3bundleid[i]))
  dftemp4<-subset(df, (genotype=='4' & Bundle.ID==clus4bundleid[i]))
  dftemp<-rbind(dftemp3,dftemp4)
  #print(i)
  #if (i == 10) 
  #  {ratio_str = 'all'}
  
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
  
  ggsave(paste(outputpath,regions_str,'_bundle', toString(i), '_LengthHist','_',ratio_str,'.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
}



