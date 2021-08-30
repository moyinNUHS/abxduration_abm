###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################################## Plot cobweb ####################################
###################################################################################

require(plotrix) #load MASS package

plot_cobweb <- function(output) {
  
  x.samples=output[['data']]
  y.output=output$res[,3,]
  #x.samples<-x.samples[, match(rownames(prcc), names(x.samples))] 
  outcome.df<-as.data.frame(cbind(x.samples,y.output)) #matrix with parameter values in columns and outcome in last column
  for (i in 1:nrow(outcome.df)) {                      #label the rows of parameter values that produced top 20% of the outcomes
    if (outcome.df$y.output[i]<quantile(outcome.df$y.output,probs = 0.8)) { 
      outcome.df$top[i] <-0 } else {
        outcome.df$top[i] <-1
      }
  }
  
  for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced bottom 20% of the outcomes
    if (outcome.df$y.output[i]<quantile(outcome.df$y.output,probs = 0.2)) { 
      outcome.df$bottom[i] <-1 } else {
        outcome.df$bottom[i] <-0
      }
  }
  
  blue<-alpha("lightskyblue", alpha=0.5)
  red<-alpha("#E85D75", alpha=0.5)
  colors<- c(blue, red) #1 for parameters that produced top outcomes and one for the bottom
  outcome.df$top<- as.factor(outcome.df$top)
  outcome.df$bottom<- as.factor(outcome.df$bottom)
  outcome.df=outcome.df[,-grep('dur',colnames(outcome.df))]
  outcome.df<- outcome.df[which(outcome.df$top==1 | outcome.df$bottom==1),]
  
  parcoordlabel<-function (x, col = 1, lty = 1,  lblcol="black",...) 
  {
    df <- as.data.frame(x)
    pr <- lapply(df, pretty)
    rx <- lapply(pr, range, na.rm = TRUE)
    x <- mapply(function(x,r) {
      (x-r[1])/(r[2]-r[1])
    },
    df, rx)
    matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty, 
            xlab = "", ylab = "Sampled parameter values", cex.lab = 1.5, axes = FALSE,...)
    axis(1, at = 1L:ncol(x), labels = c(colnames(x)), las = 2, cex.axis = 1.2)
    for (i in 1L:ncol(x)) {
      lines(c(i, i), c(0, 1), col = "grey")
      text (c(i, i), seq(0,1,length.out=length(pr[[i]])), labels = pr[[i]], 
            xpd = NA, col=lblcol, cex = 1.2)
    }
    invisible()
  }
  parcoordlabel(outcome.df[,c(1:(ncol(x.samples)-2))], col = colors[outcome.df$top])
  
}