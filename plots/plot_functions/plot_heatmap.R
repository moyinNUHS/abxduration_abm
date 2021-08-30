library(parallel)
library(ggplot2)
library(ggpubr)
library(svMisc)

scaleFUN <- function(x) sprintf("%.3f", x) #3 decimals for axis labels 

data_heatmap <- function(model = model, params_ranges) {
  
  ##parameters to plot from the 3 models 
  if (model == 'simple') {
    params_list = c('repop.s', 'mu', 'pi_ssr', 'bif', 'abx.s', 'abx.r')
  } else if (model == 'binary') {
    params_list = c('repop.s', 'repop.r', 'mu', 'pi_ssr', 'bif', 'abx.s', 'abx.r')
  } else if (model == 'frequency') {
    params_list = c('total_prop', 'k', 's_growth', 'r_growth', 'r_thres', 'r_trans', 'pi_ssr', 'abx.s', 'abx.r')
  }
  
  #object to store all plots
  plots_data = vector("list", length(params_list))
  
  for (i in 1:length(params_list)) { # i refers to the rows of the plot matrix
    for (j in 1:(length(params_list)-1)) { # j refers to the columns of the plot matrix

      #show progress
      progress(i, progress.bar = TRUE)
      Sys.sleep(0.01)
      if (j == (length(params_list)-1)) cat(paste0(i, " in ", length(params_list), "th row is Done! \n"))
      
      if (i == j) { #blank plots 
        
        y_name = params_list[i] 
        y_seq = params_ranges$esbl[[i]]
        
        x_name = params_list[j] 
        x_seq = params_ranges$esbl[[j]]
        
        outcome.df = cbind.data.frame(x = rep(x_seq, each = pixels), 
                                      y = rep(y_seq, pixels),
                                      outcome = 0)
        
      } else { #filled plots 
        
        if (i > j) { # esbl plots 
          
          y_name = params_list[i] 
          y_seq = params_ranges$esbl[[i]]
          
          x_name = params_list[j] 
          x_seq = params_ranges$esbl[[j]]
          
          # Parameters not explored are considered fixed as defined in default_params
          para.feed = c(para[[1]]$esbl, para[[2]]$esbl)
          default.para.list = para.feed[which(names(para.feed) %in% parameters_diff_prevalence)]
          
        } else { # cpe plots 
          
          y_name = params_list[i] 
          y_seq = params_ranges$cpe[[i]]
          
          x_name = params_list[j] 
          x_seq = params_ranges$cpe[[j]]
          
          # Parameters not explored are considered fixed as defined in default_params
          para.feed = c(para[[1]]$cpe, para[[2]]$cpe)
          default.para.list = para.feed[which(names(para.feed) %in% parameters_diff_prevalence)]
          default.para.list$abx.r = 0
        }
        
        # format parameter input into lists to run model 
        default.para.df = as.data.frame(matrix(as.vector((rep(unlist(default.para.list), each = pixels**2))), ncol=pixels**2, byrow = T))
        rownames(default.para.df) = names(default.para.list)
        colnames(default.para.df) = 1:pixels**2
        default.para.df = default.para.df[match(parameters_diff_prevalence, rownames(default.para.df)),]
        default.para.df[grep(x_name,rownames(default.para.df)),] = rep(x_seq, each = pixels)
        default.para.df[grep(y_name,rownames(default.para.df)),] = rep(rep(y_seq, pixels), each = length(grep(y_name, rownames(default.para.df))))
        feed.list = as.list(default.para.df)
        
        ###run model 
        if (model == 'simple'){
          
          foo = function (x) {
            diff_prevalence(n.bed= x[1], max.los=x[2], prop_R=x[3], prop_S=x[4],
                            bif=x[5], pi_ssr=x[6], repop.s=x[7], 
                            mu=x[8], abx.s=x[9], abx.r=x[10], p.infect=x[11], 
                            cum.r.1=x[12], p.r.day1=x[13], 
                            short_dur=x[14], long_dur=x[15])[3]
          }
          
          save_runs = mclapply(feed.list, foo,  mc.cores = 11)
          
        } else if (model == 'binary'){
          
          foo=function (x) { 
            diff_prevalence(n.bed=x[1], max.los=x[2], prop_R=x[3], prop_r=x[4], 
                            prop_Sr=x[5], prop_S=x[6],
                            bif=x[7], pi_ssr=x[8], repop.s=x[9], repop.r=x[10], 
                            mu=x[11], abx.s=x[12], abx.r=x[13], 
                            p.infect=x[14], cum.r.1=x[15], p.r.day1=x[16], 
                            short_dur=x[17], long_dur=x[18])[3]
          }
          
          save_runs = mclapply(feed.list, foo, mc.cores = 11)
          
        } else if (model == 'frequency'){
          
          foo = function (x) { 
            diff_prevalence(n.bed=x[1], max.los=x[2], p.infect=x[3], 
                            cum.r.1=x[4], p.r.day1=x[5],
                            K=x[6], total_prop=x[7], prop_R=x[8], pi_ssr=x[9], 
                            r_trans=x[10], r_growth=x[11], r_thres=x[12], s_growth=x[13],
                            abx.s=x[14], abx.r=x[15], short_dur=x[16],long_dur=x[17])[3]
          }
          
          save_runs = mclapply(feed.list, foo,  mc.cores = 11)
        }
        
        outcome.df = cbind.data.frame(x = rep(x_seq, each=pixels), 
                                      y = rep(y_seq, pixels),
                                      outcome = unlist(save_runs))
      }
      
      plots_data[[i]][[j]] = outcome.df
    }
  }
  
  return(plots_data)
}

plot_heatmap <- function (model = model, plots_data) {
  
  font.size.axislab = 14
  font.size = 7
  
  ##parameters to plot from the 3 models 
  if (model == 'simple') {
    params_list = c('repop.s', 'mu', 'pi_ssr', 'bif', 'abx.s', 'abx.r')
  } else if (model == 'binary') {
    params_list = c('repop.s', 'repop.r', 'mu', 'pi_ssr', 'bif', 'abx.s', 'abx.r')
  } else if (model == 'frequency') {
    params_list = c('total_prop', 'k', 's_growth', 'r_growth', 'r_thres', 'r_trans', 'pi_ssr', 'abx.s', 'abx.r')
  }
  
  plots = vector("list", length(plots_data))
  
  for (i in 1:length(plots_data)) {
    for (j in 1:(length(plots_data)-1)) {
      
      outcome.df = plots_data[[i]][[j]]
      y_name = params_list[i]
      x_name = params_list[j]
      
      if (i == j & j == 1) {
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          ylab(y_name) + 
          scale_fill_gradient2(low="#388697", mid="white", high="#CC2936", 
                               midpoint = 0, guide = NULL) +
          theme(axis.line=element_blank(),axis.text.x=element_blank(),
                axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.y = element_text(size=font.size.axislab),
                axis.title.x=element_blank(),legend.position="none",
                panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank())
        
      } else if (i == j) { #diagonal plots to be filled with white spaces 
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          scale_fill_gradient2(low="#388697", mid="white", high="#CC2936", 
                               midpoint = 0, guide = NULL) +
          theme_void()
        
      } else if (j == 1 & i != length (plots_data)) { #1st column to keep y axis label on left 
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          scale_fill_gradient2(low="#388697", mid="grey95", high="#CC2936", 
                               midpoint = 0, limits=range(outcome.df$outcome),
                               breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3))+0.001,
                                          as.numeric(format(round(0),0), nsmall=0), 
                                          as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))-0.001),
                               name = "Difference in resistance carriers in wards given\nlong vs short antibiotic durations") +
          scale_y_continuous(labels=scaleFUN) +
          ylab(y_name)+
          theme_minimal()+
          theme(legend.position = 'none', 
                axis.title.y = element_text(size=font.size.axislab),
                axis.text.y = element_text(size=font.size), 
                axis.title.x = element_blank(),
                axis.text.x = element_blank())
        
      } else if (j == 1 & i == length(plots_data)){ #left bottom corner plot to keep all labels 
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          scale_fill_gradient2(low="#388697", mid="grey95", high="#CC2936", 
                               midpoint = 0, limits=range(outcome.df$outcome),
                               breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3))+0.001,
                                          as.numeric(format(round(0),0), nsmall=0), 
                                          as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))-0.001),
                               name = "Difference in resistance carriers in wards given\nlong vs short antibiotic durations") +
          ylab(y_name)+
          xlab(x_name)+
          scale_y_continuous(labels=scaleFUN) +
          scale_x_continuous(labels=scaleFUN) +
          theme_minimal()+
          theme(legend.position = 'none',
                axis.title.y = element_text(size=font.size.axislab),
                axis.text.y = element_text(size=font.size),
                axis.title.x = element_text(size=font.size.axislab),
                axis.text.x = element_text(size=font.size))
        
      } else if (j > 1 & i == length(plots_data)){ #plots on the bottom to keep x axis labels 
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          scale_fill_gradient2(low="#388697", mid="grey95", high="#CC2936", 
                               midpoint = 0, limits=range(outcome.df$outcome),
                               breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3))+0.001,
                                          as.numeric(format(round(0),0), nsmall=0), 
                                          as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))-0.001),
                               name = "Difference in resistance carriers in wards given\nlong vs short antibiotic durations") +
          xlab(x_name)+
          theme_minimal()+
          scale_x_continuous(labels=scaleFUN) +
          theme(legend.position = 'none', 
                axis.title.x = element_text(size=font.size.axislab),
                axis.text.x = element_text(size=font.size),
                axis.title.y=element_blank(),
                axis.text.y=element_blank())
        
      } else if (j > 1 &  i < length(plots_data) ) { #plots in the middle to remove all axis labels 
        
        plot = ggplot(outcome.df, aes(x, y)) +
          geom_tile(aes(fill = outcome)) +
          scale_fill_gradient2(low="#388697", mid="grey95", high="#CC2936", 
                               midpoint = 0, limits=range(outcome.df$outcome),
                               breaks = c(as.numeric(format(round(min(outcome.df$outcome),3),nsmall=3))+0.001,
                                          as.numeric(format(round(0),0), nsmall=0), 
                                          as.numeric(format(round(max(outcome.df$outcome),3),nsmall=3))-0.001),
                               name = "Difference in resistance carriers in wards given\nlong vs short antibiotic durations") +
          theme_minimal()+
          scale_y_continuous(labels=scaleFUN) +
          scale_x_continuous(labels=scaleFUN) +
          theme(legend.position = 'none', 
                text = element_text(size=font.size),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.text.x=element_blank())
      }
      plots[[i]][[j]] = ggplotGrob(plot)
    }
  }
  
  combinedplot = ggpubr::ggarrange(plotlist = unlist(plots, recursive=FALSE), ncol = (length(plots_data)-1), nrow = length(plots_data))
  
  return(combinedplot)
}

