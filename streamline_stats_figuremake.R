#read csv file
#plot FA Along tracts for group 1
#plot FA along tracts for group 2
# stats
library(ggplot2)
library(readxl)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) # need to run every time you start R and want to use %>%
library(dplyr)
library(Rmisc)

project = 'AMD'
ratio = 1

std_mean <- function(x) sd(x)/sqrt(length(x))

path_join <- function(path){
  if (Sys.info()['sysname']=='Windows') {
    newpath <- paste(unlist(path), collapse='\\')
  }
  else{
    newpath <- paste(unlist(path), collapse='/')
  }
  return(newpath)
}

my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}


if (Sys.info()[['user']]=='JacquesStout') {
  atlas_file='C:\\Users\\JacquesStout\\Documents\\Work\\atlases\\IITmean_RPI_index.xlsx'
  if (project == 'AMD') {
    mainpath = 'C:\\Users\\JacquesStout\\Documents\\Work\\AMD'
    outputpath= 'C:\\Users\\JacquesStout\\Documents\\Work\\AMD\\R_figures'
    
  }
}
if (Sys.info()[['user']]=='jas') {
  atlas_file='/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx'
  if (project=='AD_Decode') {
    mainpath = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/'
    outputpath='/Users/jas/jacques/AD_decode_abstract/R_figures/'
  }
  if (project == 'AMD') {
    mainpath = '/Volumes/Data/Badea/Lab/human/AMD/'
    outputpath= '/Users/jas/jacques/Whiston_article/R_figures/'
  }
}

if (ratio == 1) {ratio_str = 'all'} else {ratio_str = paste('ratio_',as.character(ratio),sep='')}
if (ratio == 1) {folder_str = ''} else {folder_str = paste('_',as.character(ratio))}

space_param = '_affinerigid'

my_atlas <- read_excel(atlas_file)
outputpath_singlebundle = 
  csvpath_folder = path_join(list(mainpath, paste('Statistics',space_param,'_non_inclusive_symmetric',as.character(folder_str),sep='')))

myYear = 'Initial'
if (myYear=='Initial'){
  groups =list('Paired Initial AMD','Paired Initial Control')
}
if (myYear=='2Year'){
  groups =list('Paired 2-YR AMD','Paired 2-YR Control')
}
#groups =list('Paired 2-YR AMD','Paired 2-YR Control')
#groups =list('Paired Initial AMD','Paired Initial Control')
#groups =list('1-Control','2-AMD')
group_colors <- list('blueviolet','chartreuse1')

regions_all <- vector(mode = 'list', length = 8)

#regions_all <- list(list(28, 9),list(62, 1),list(77, 43),list(61, 29))
regions_all <- list(list(62, 28),list(58, 45),list(77, 43),list(61, 29))
#references <- list(fa,md)

for (z in 3:4) {
  regions = regions_all[[z]]
  regions_str = tolower(paste(my_atlas[regions[[1]],3],'_',my_atlas[regions[[1]],4],'_to_',my_atlas[regions[[2]],3],'_',my_atlas[regions[[2]],4],sep=''))
  
  g3file = path_join(list(csvpath_folder,paste(groups[1],'_',regions_str,'_',ratio_str,'_bundle_stats.csv',sep='')))
  g4file = path_join(list(csvpath_folder,paste(groups[2],'_',regions_str,'_',ratio_str,'_bundle_stats.csv',sep='')))
  
  g3<-read.csv(g3file, header = TRUE)
  g4<-read.csv(g4file, header = TRUE)
  
  df3<-data.frame(g3)
  df3<-cbind(df3,genotype=groups[[1]])
  
  df4<-data.frame(g4)
  df4<-cbind(df4,genotype=groups[[2]])
  
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
  
  #write these results to tables
  for (i in 1:length(clus3bundleid)){
    dftemp3<-subset(df, (genotype==groups[[1]] & Bundle.ID==clus3bundleid[i]))
    dftemp4<-subset(df, (genotype==groups[[2]] & Bundle.ID==clus4bundleid[i]))
    dftemp<-rbind(dftemp3,dftemp4)
    ksbundlefa<-ks.test(dftemp3$fa, dftemp4$fa, alternative = c('two.sided'),exact = NULL)
    ksbundleLength<-ks.test(dftemp3$Length, dftemp4$Length, alternative = c('two.sided'),exact = NULL)
    
    clus3bundleid<-unique(g3$Bundle.ID)
    #plot signal along fibers
    
    dfstream3Bundlei<-subset(df, (genotype==groups[[1]] & Bundle.ID==clus3bundleid[i]))
    dfstream4Bundlei<-subset(df, (genotype==groups[[2]] & Bundle.ID==clus4bundleid[i]))
    dfstreamBundlei<-rbind(dfstream3Bundlei,dfstream4Bundlei)
    dfstreamBundlei$genotypenorm <- ((-1)*as.numeric(dfstreamBundlei$genotype) - 3.5)
    
    dfstreamBundlei$lwr=dfstreamBundlei$fa-1.96*STDERR(dfstreamBundlei$fa)
    dfstreamBundlei$hi=dfstreamBundlei$fa+1.96*STDERR(dfstreamBundlei$fa)
    
    #%start save stats per bundle: FA and Length
  
    mytLength<-t.test(Length ~ genotype, data = dfstreamBundlei)
    mytCoherence<-t.test(Streamline.coherence ~ genotype, data = dfstreamBundlei)
    mytFA<-t.test(fa ~ genotype, data = dfstreamBundlei)
    mytMD<-t.test(md ~ genotype, data = dfstreamBundlei)
    mytAD<-t.test(ad ~ genotype, data = dfstreamBundlei)
    mytRD<-t.test(rd ~ genotype, data = dfstreamBundlei)
  
    size_bundle <- list(length(unique(dftemp3$Streamlines.ID)),length(unique(dftemp4$Streamlines.ID)))
    #size_bundle_1 <- length(unique(dftemp3$Streamlines.ID))
    #size_bundle_2 <- length(unique(dftemp4$Streamlines.ID)) 
    groups_order <- list(2,1)
    
    mytableLength<-as_data_frame(
      cbind(mytLength$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytLength$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytLength$estimate[as.numeric(groups_order[2])], mytLength$statistic, mytLength$p.value,  
            mytLength$conf.int[1], mytLength$conf.int[2] , mytLength$parameter)
    )
    mytableCoherence<-as_data_frame(
      cbind(mytCoherence$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytCoherence$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytCoherence$estimate[as.numeric(groups_order[2])], mytCoherence$statistic, mytCoherence$p.value,  
            mytCoherence$conf.int[1], mytCoherence$conf.int[2] , mytCoherence$parameter)
    )
    mytableFA<-as_data_frame(
      cbind(mytFA$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytFA$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytFA$estimate[as.numeric(groups_order[2])], mytFA$statistic, mytFA$p.value,  
            mytFA$conf.int[1], mytFA$conf.int[2] , mytFA$parameter)
    )
    mytableMD<-as_data_frame(
      cbind(mytMD$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytMD$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytMD$estimate[as.numeric(groups_order[2])], mytMD$statistic, mytMD$p.value,  
            mytMD$conf.int[1], mytMD$conf.int[2] , mytMD$parameter)
    )
    mytableAD<-as_data_frame(
      cbind(mytAD$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytAD$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytAD$estimate[as.numeric(groups_order[2])], mytAD$statistic, mytAD$p.value,  
            mytAD$conf.int[1], mytAD$conf.int[2] , mytAD$parameter)
    )
    mytableRD<-as_data_frame(
      cbind(mytRD$data.name, i, as.numeric(size_bundle[as.numeric(groups_order[1])]), mytRD$estimate[as.numeric(groups_order[1])], as.numeric(size_bundle[as.numeric(groups_order[2])]), mytRD$estimate[as.numeric(groups_order[2])], mytRD$statistic, mytRD$p.value,  
            mytRD$conf.int[1], mytRD$conf.int[2] , mytRD$parameter)
    )
    # col.names=c('contrast', names(mytFA$estimate)[1],names(mytFA$estimate)[2],names(mytFA$statistic),'pvalue', 'CIlwr','CIhi', 'df'))
    mycolnames<-c('contrast', 'bundle',paste('num str ',groups[as.numeric(groups_order[1])],sep=''), names(mytFA$estimate)[as.numeric(groups_order[1])], paste('num str ',groups[as.numeric(groups_order[2])],sep=''), names(mytFA$estimate)[as.numeric(groups_order[2])], names(mytFA$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df')
    
    myfile<-paste(outputpath,regions_str,'_bundle', toString(i), '_stats.csv',sep='')
    if (!file.exists(myfile)){
    #add column names for the first instance
    write.table(mytableLength, file=myfile, col.names = mycolnames , sep = ',' , row.names = F,append=TRUE)
    #remove columnnames for subsequent rows
    write.table(mytableCoherence, file=myfile, sep = ',' , row.names = F,append=TRUE,col.names=FALSE)
    write.table(mytableFA, file=myfile, sep = ',' , row.names = F,append=TRUE,col.names=FALSE)
    write.table(mytableMD, file=myfile, sep = ',' , row.names = F,append=TRUE,col.names=FALSE)
    write.table(mytableAD, file=myfile, sep = ',' , row.names = F,append=TRUE,col.names=FALSE)
    write.table(mytableRD, file=myfile, sep = ',' , row.names = F,append=TRUE,col.names=FALSE)
    }
    
    clus3streamid<-unique(dfstream3Bundlei$Streamlines.ID)
    clus4streamid<-unique(dfstream4Bundlei$Streamlines.ID)
    
    #This is the part where you should probably loop it by references
    
    ggplot(dfstreamBundlei, aes(y=fa, x=Point.ID, color=genotype)) +
      geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
      ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
      theme_classic()+ xlab('Points')+ ylab('FA')+
      scale_color_manual(values=c(group_colors[1],group_colors[2]))+
      scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
      theme(legend.position='top')+
      theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
            axis.text.y = element_text(face='bold', size=14, angle=0),
            axis.line.x = element_line(colour = 'black', size=0.5),
            axis.line.y = element_line(colour = 'black', size=0.5),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),
            panel.background = element_blank())
    ggsave(path_join(list(outputpath,paste(myYear,'_',regions_str,'_alongbundle', toString(i), '_FA.pdf',sep=''))), plot = last_plot(), device = 'pdf',
           scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)
    
    ggplot(dfstreamBundlei, aes(y=Local.coherence, x=Point.ID, color=genotype)) +
      geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
      ggtitle(paste('Coherence along streamlines', toString(i), sep= ' '))+
      theme_classic()+ xlab('Points')+ ylab('Coherence')+
      scale_color_manual(values=c(group_colors[1],group_colors[2]))+
      scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
      theme(legend.position='top')+
      theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
            axis.text.y = element_text(face='bold', size=14, angle=0),
            axis.line.x = element_line(colour = 'black', size=0.5),
            axis.line.y = element_line(colour = 'black', size=0.5),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),
            panel.background = element_blank())
    ggsave(path_join(list(outputpath,paste(myYear,'_',regions_str,'_alongbundle', toString(i), '_coherence.pdf',sep=''))), plot = last_plot(), device = 'pdf',
           scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)
    
    ggplot(dfstreamBundlei, aes(y=(Local.coherence*fa), x=Point.ID, color=genotype)) +
      geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
      ggtitle(paste('FA weighted by coherence along streamlines', toString(i), sep= ' '))+
      theme_classic()+ xlab('Points')+ ylab('FAw')+
      scale_color_manual(values=c(group_colors[1],group_colors[2]))+
      scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
      theme(legend.position='top')+
      theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
            axis.text.y = element_text(face='bold', size=14, angle=0),
            axis.line.x = element_line(colour = 'black', size=0.5),
            axis.line.y = element_line(colour = 'black', size=0.5),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),
            panel.background = element_blank())
    ggsave(path_join(list(outputpath,paste(myYear,'_',regions_str,'_alongbundle', toString(i), '_FAwcoherence.pdf',sep=''))), plot = last_plot(), device = 'pdf',
           scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)
    
    m03 <- lm(fa~poly(Point.ID,4)*genotype, data=dfstreamBundlei)
    newdf<-expand.grid(Point.ID <- seq(min(dfstream3Bundlei$ Point.ID), max(dfstream3Bundlei$Point.ID), length.out=100))
    
    dfstreamBundlei$LoCI<-predict(m03,data=dfstreamBundlei,interval='confidence', level=0.95)[,2]
    dfstreamBundlei$HiCI<-predict(m03,data=dfstreamBundlei,interval='confidence', level=0.95)[,3]
    
    #geom_smooth is what adds confidence intervals
    
    ggplot(dfstreamBundlei, aes(y=m03$fitted.values, x=Point.ID, color=genotype)) +
      geom_smooth(aes(ymin=LoCI,ymax=HiCI,fill=genotype, colour=genotype),stat = 'identity')+
      ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
      theme_classic()+
      scale_color_manual(values=c(group_colors[1],group_colors[2]))+
      scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
      theme(legend.position='top')+
      theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
            axis.text.y = element_text(face='bold', size=14, angle=0), 
            axis.line.x = element_line(colour = 'black', size=0.5),
            axis.line.y = element_line(colour = 'black', size=0.5),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),
            panel.background = element_blank())
    
    ggsave(path_join(list(outputpath,paste(myYear,'_',regions_str,'_bundle', toString(i), '_FAfitCI.pdf',sep=''))), plot = last_plot(), device = 'pdf',
           scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)
    
    
    
  }
}


# ggplot(dfstreamBundlei, aes(y=Point.weight.to.centroid, x=Point.ID, color=genotype)) +
#   #geom_smooth(stat = 'summary', group=genotype, fun.data = function(y) data.frame(ymin = quantile(y, .05), y = mean(y), ymax = quantile(y, .95))) +
#   #geom_smooth(stat = 'summary', fun.data = function(y) data.frame(ymin = quantile(y, .025), y = mean(y), ymax = quantile(y, .975))) +
#   #geom_smooth(stat = 'summary', fun.data = mean_cl_quantile)
#   #stat_smooth() +
#   geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
#   ##bad crashes geom_smooth(aes(ymin = lwr, ymax = hi,color=genotype),stat = 'identity')+
#   ##bad crashes geom_ribbon(data=dfstreamBundlei,aes(ymin=lwr,ymax=hi),alpha=0.3)+
#   #aes(ymin=lwr,ymax=hi,fill=genotype, colour=genotype),stat = 'identity')+
#   #geom_histogram( alpha=0.5, position='identity')+
#   #geom_density(alpha=0.6)+
#   ggtitle(paste('Point weight along streamlines', toString(i), sep= ' '))+
#   theme_classic()+
#   scale_color_manual(values=c(group_colors[1],group_colors[2]))+
#   scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
#   theme(legend.position='top')+
#   theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
#         axis.text.y = element_text(face='bold', size=14, angle=0),
#         axis.line.x = element_line(colour = 'black', size=0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5),
#         # panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())
# ggsave(path_join(list(outputpath,paste(regions_str,'_alongbundle', toString(i), '_pointweight.pdf',sep=''))), plot = last_plot(), device = 'pdf',
#        scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)







# #this models does not work very well, but let's explore and make it better later
# #http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#is-the-likelihood-ratio-test-reliable-for-mixed-models
# fm1 <- lmer(
#   formula = fa ~ poly(Point.ID,degree=4)+genotype+ (Point.ID|Streamlines.ID), data = dfstreamBundlei
#   #+ , 
# )
# cV <- ranef(fm1, condVar = TRUE)   
# ranvar <- attr(cV[[1]], 'postVar')
# 
# fapred <- predict(fm1,data=dfstreamBundlei,re.form=NA)
# 
# ggplot(dfstreamBundlei, aes(y=fapred, x=Point.ID, color=genotype)) +
#   geom_point()+
#   #geom_smooth(data=dfstreamBundlei, se = TRUE, level=0.95, aes(color=genotype,fill=genotype))+
#   #geom_line(aes(group = genotype), alpha = .3) +
#   #geom_line(data = dfstreamBundlei, alpha = .8, size = 1) +
#   stat_smooth() +
#   #geom_histogram( alpha=0.5, position='identity')+
#   #geom_density(alpha=0.6)+
#   ggtitle(paste('FA along streamlines', toString(i), sep= ' '))+
#   theme_classic()+
#   scale_color_manual(values=c(group_colors[1],group_colors[2]))+
#   scale_fill_manual(values=c(group_colors[1],group_colors[2]))+
#   theme(legend.position='top')+
#   theme(axis.text.x = element_text(face='bold',  size=14, angle=0),
#         axis.text.y = element_text(face='bold', size=14, angle=0),
#         axis.line.x = element_line(colour = 'black', size=0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5),
#         # panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())
# ggsave(path_join(list(outputpath,paste(regions_str,'_pred_allbundle', toString(i), '_FA.pdf',sep=''))), plot = last_plot(), device = 'pdf',
#        scale = 1, width = 5, height = 5, units = c('in'),dpi = 300)

