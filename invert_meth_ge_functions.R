#invert_meth_ge_functions.R
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())



# SET UP BIOPROJECTS ------------------------------------------------------

bioprojectDirs = read.table('bioprojects/bioprojects.txt', stringsAsFactors = FALSE)$V1
names(bioprojectDirs) = sapply(bioprojectDirs, function(x) strsplit(x, '_')[[1]][1])
dat_names = names(bioprojectDirs)




# FUNCTIONS ---------------------------------------------------------------

#get inverse log2 labels for axes
log2_to_percent = function(x) {
  round(2^x, 2)*100
}

get_pearson_cor = function(dat, xcol, ycol){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pearsonCor=cor(x=dat[,xcol],
                 y=dat[,ycol])
  r=round(pearsonCor, digits=2)
  return(r)
}


#same as r2 version but for pearson correlation 
plot_scatter_pearsonCor_annotated = function(dat, xcol, ycol, xlab, ylab, ALPHA=0.1, ylim=FALSE){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pearsonCor=cor(x=dat[,xcol],
                 y=dat[,ycol])
  r=round(pearsonCor, digits=2)
  print(summary(lm1))
  plt = dat %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha = ALPHA) +
    labs(x=xlab,
         y=ylab)
  pbuild = ggplot_build(plt)
  yrange = pbuild$layout$panel_params[[1]]$y.range
  xrange = pbuild$layout$panel_params[[1]]$x.range
  if (length(ylim)>1){
    print('setting y limits')
    plt=plt + lims(y=ylim) +
      annotate("text", x = xrange[1], y = ylim[2],
               label = paste('italic(r) ==', r), parse=TRUE, color='black',
               hjust=0)
  } else {
    plt=plt +
      annotate("text", x = xrange[1], y = yrange[2],
               label = paste('italic(r) ==', r), parse=TRUE, color='black',
               hjust=0)
  }
  return(plt)
}


#
plot_shared_x = function(plotList, xlab){
  revList = lapply(plotList, function(x) return(x + theme(axis.title.x = element_blank())))
  top = plot_grid(plotlist=revList, nrow=1)
  xlab = ggdraw() + draw_label(xlab)
  plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1), nrow=nrow)
}



#
plot_shared_x_y = function(plotList, xlab, ylab, relXlab=1/15, relYlab=1/15){
  revList = lapply(plotList, function(x) return(x + theme(axis.title = element_blank())))
  top = plot_grid(plotlist=revList, nrow=1)
  pxlab = ggdraw() + draw_label(xlab)
  pylab = ggdraw() + draw_label(ylab, angle=90)
  wX = plot_grid(top, pxlab, nrow=2, rel_heights=c(1, relXlab))
  wXY = plot_grid(pylab, wX, nrow=1, rel_widths=c(relYlab,1))
  return(wXY)
}

#get inverse log2 labels for axes
log2_to_percent = function(x) {
  round(2^x, 2)*100
}

#get inverse log2 labels for axes
log10_to_percent = function(x) {
  round(10^x, 6)*100
}



plot_RNA_PCA <- function (df, coldat, intgroup = "treatInfo", ntop = 25000, returnData = F, pcs = 2, pc1 = 1, pc2 = 2, main = "\n", SIZE = 4, legendTitle=NULL, xInvert=1) 
{
  df=as.matrix(df)
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(coldat[, intgroup, 
                                      drop = FALSE])
  group <- intgroup.df[,intgroup]
  d <- data.frame(pca$x[,1:pcs], group = group, 
                  intgroup.df, name = colnames(df))
  
  d$my_letter = coldat$my_letter
  attr(d, "percentVar") <- percentVar[1:2]
  #Invert X if needed
  d[,paste('PC', pc1, sep = '')]=d[,paste('PC', pc1, sep = '')]*xInvert
  
  g = ggplot(data = d, 
             aes_string(x = paste('PC', pc1, sep = ''),
                        y = paste('PC', pc2, sep = ''), color = "group")) + 
    geom_point(size = SIZE) +
    xlab(paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")) + 
    ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")) + 
    coord_fixed() 
  g = g + labs(subtitle=main, color=legendTitle)
  g = g + theme(legend.position='right',
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank())
  print(g)
  if (returnData == T){
    return(d)
  }
  else{
    return(g)
  }
}



#function to load a methylKit_regionCounts.Rdata and its assocaited bisulfite_treat_table.txt (when running methylKitWindowsV2.R)
#returns a dataframe with methylation proportions for each region for each sample with labels from the id column in bisulfite_treat_table.txt
get_meth_proportion_from_reg_counts = function(rgFile, bs_treat_table){
  ll=load(rgFile)
  tt = read_tsv(bs_treat_table)
  rc = data.frame(my.reg.counts) %>% 
    as_tibble()
  pos = rc %>% 
    dplyr::select(c('chr', 'start', 'end'))
  ccols = colnames(rc)[grep('numCs', colnames(rc))]
  snums = as.numeric(sub('numCs', '', ccols))
  meth = rc %>% 
    dplyr::select(grep('numCs', colnames(rc))) %>% 
    data.frame()
  unmeth = rc %>% 
    dplyr::select(grep('numTs', colnames(rc))) %>% 
    data.frame()
  mdat = pos
  samples = ccols
  for (i in snums){
    sample = samples[i]
    methCol = paste('numCs', i, sep='')
    unmethCol = paste('numTs', i, sep='')
    mVec = meth[,methCol]
    umVec = unmeth[,unmethCol]
    rVec = mVec / (mVec + umVec)
    mdat[,sample] = rVec
  }
  colnames(mdat)[4:ncol(mdat)] = tt$id
  return(mdat)
}


#function to plot pca from normalized expression counts
build_pca = function(df,
                     coldata,
                     ntop = 25000,
                     pcs = 6){
  #get row varainces and select top
  df = as.matrix(df)
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  #build pca
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d = cbind(data.frame(pca$x[,1:pcs]),
            coldata)
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}


#function to plot PCA scatterplot from output from build_rld_pca
plot_rld_pca = function(df,
                        group_col = 'treatment',
                        pc1 = 1,
                        pc2 = 2,
                        subtitle = "",
                        size = 4,
                        legend_title=NULL,
                        x_invert=1,
                        legend_position = 'none',
                        fix_coords = TRUE){
  #select PCs to plot and nvert X if needed
  plt_df = tibble(x = df[,paste('PC', pc1, sep = '')]*x_invert,
                  y = df[,paste('PC', pc2, sep = '')],
                  col = factor(df[,group_col]))
  #pull out the percent variances
  percentVar = attr(df, "percentVar")[c(pc1,pc2)]
  #build axis labs
  xlab = paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")
  ylab = paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")
  g = plt_df %>% 
    ggplot(aes(x=x,
               y=y,
               color=col)) + 
    geom_point(size = size) +
    labs(x = xlab,
         y = ylab,
         subtitle=subtitle,
         color=legend_title) +
    theme(legend.position=legend_position,
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
  if (fix_coords){
    g = g + coord_fixed() 
  }
  return(g)
}



#function to plot bimodal distirbution
plot_bimodal = function(lvl_df, logo, logo_scale){
  Ngenes = length(unique(lvl_df$name))
  hist = lvl_df %>% 
    ggplot(aes(x=lfrac_meth)) +
    geom_histogram() +
    scale_x_continuous(labels = log2_to_percent) +
    labs(subtitle=paste(Ngenes, 'genes'), x='% methylation') +
    theme(axis.title = element_blank())
  ggdraw(hist) + draw_image(logo,
                            x=logo_corner - logo_scale,
                            y=logo_corner-logo_scale, 
                            width = logo_scale,
                            height = logo_scale)
}




