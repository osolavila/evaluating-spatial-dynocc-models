

##Colext number of parameters (assuming no factors)
get.nP <- function(psiformula = ~ 1, # First-year occupancy
                   gammaformula = ~ 1, # Colonization
                   epsilonformula = ~ 1, # Extinction
                   pformula = ~ 1) {
  
  formulas <- list(psiformula, gammaformula, epsilonformula,  pformula)
  ##Get formula right side as string
  formulas <- lapply(formulas, function(x){
    if (typeof(x)=="character") as.character(as.formula(x))[2]
    else as.character(x)[2]
  })
  
  ##Get vector of number of parameters
  num_param <- sapply(formulas, function(x)
    if (x=="1") 1
    else length(strsplit(x, "+", fixed=T)[[1]])+1)
  
  return(sum(num_param))
  
  
}


##Grouped boxplot
groupped.boxplot.helper <- function(df, group1.name, group2.name){
  legend <- unique(df[,group1.name])
  labels <- unique(df[,group2.name])
  ##length group1
  group_1 <- length(legend)
  ##length group2
  group_2 <- length(labels)
  position.boxplots <- c()
  position.axis <- c()
  for (i in 1:group_2){
    position.boxplots <- c(position.boxplots, 1:group_1+(group_1+1)*(i-1) )
    position.axis <- c(position.axis, group_1/2+ 0.5 + (group_1+1)*(i-1))
  }
  return(list(legend=legend,
              labels=labels,
              position.boxplots=position.boxplots,
              position.axis=position.axis))
}

# boxplot.utils <- groupped.boxplot.helper(model.performance, "variable",
#                                        "sp_group")
# 
# colors  <- c("#332288", "#F0E442", "#0072B2")
# boxplot(data= model.performance, value~ variable+sp_group,
#         at = boxplot.utils[["position.boxplots"]], col = colors,  xaxt = 'n',
#         ylim=c(0,50), xlab="", ylab="")
# axis(side = 1, at = boxplot.utils[["position.axis"]],
#      labels = boxplot.utils[["labels"]], lwd.ticks = FALSE,
#      las=1, cex.axis=0.7)
# legend("topleft", fill = colors, 
#        legend = boxplot.utils[["legend"]] , ncol = 2)