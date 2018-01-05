findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
 if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
 stop("vectors must be same length")
 arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


errorBarPlot <- function(vals,sampleType,col=standardColors(),legend=TRUE,elwd=2,ylim=NA,xlim=NA,length=0.1,...){
 if(is.null(dim(vals))) vals = cbind(vals,vals)
 yy <- t(findFromGroups(vals,sampleType))
 col = col[1:dim(yy)[1]]
 ee <- t(findFromGroups(vals,sampleType,sd))
 if(is.na(ylim[1])) ylim = c(0,max(ee+yy))
 if(is.na(xlim[1])) xlim = c(0,((dim(ee)[2]+1)*dim(ee)[1])*1.4)
 barx <- barplot(yy, beside=TRUE,col=col,legend.text=legend,ylim=ylim,xlim=xlim,...)
 error.bar(barx,yy,ee,lwd=elwd,length=length)
}

mouse2human2 <- function (mouse, m2h){
 # Convert mouse to human symbols
 rownames(m2h) = m2h$Mou
 noHumn = which(!(mouse%in%m2h$Mouse_Symbol))
 humn = which((mouse%in%m2h$Mouse_Symbol))
 mouse[humn] = as.character(m2h[mouse[humn],1])
 mouse[noHumn] = toupper(mouse[noHumn])
 return(mouse)
}

t.test.l <- function(x){
  l  = length(x)
  tt = t.test(x[1:(l/2)],x[(l/2+1):l],paired=FALSE)
  out = c(tt$est,tt$stat,tt$p.val)
  if(is.na(out[2])) out[2] = 0
  if(is.na(out[3])) out[3] = 1
  return(out)
}

getAnovaPvalforApply <- function(x,varLabels,varWeights=NULL){
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))  
  return(aov.out[[1]]$'Pr(>F)'[1])
}



meanEx <- function(x) {if(sum(x)==0) return(0); return(mean(x[x>0]));}


t.test.l.paired <- function(x){
  l  = length(x)
  tt = t.test(x[1:(l/2)],x[(l/2+1):l],paired=TRUE)
  out = c(tt$est,tt$stat,tt$p.val)
  if(is.na(out[2])) out[2] = 0
  if(is.na(out[3])) out[3] = 1
  return(out)
}

  
getSpecificityScore <- function(propExpr,returnScore=FALSE) {
  # GET THE SPECIFICITY SCORE TO DETERMINE GENES THAT WILL GO INTO DENDROGRAM
  # This function is very similar to the beta marker gene score and was what was used to build the trees
  # Marker "specificity" score is combination of specificity and sparsity
  # propExpr = proportions of cells in a given cluster with CPM/FPKM > 1 (or 0, HCT uses 1)
  keep.cl <- colnames(propExpr)
  max.scores <- sapply(1:(length(keep.cl) - 1), function(x) {
    y <- c(rep(1, x), rep(0, length(keep.cl) - x))
    d1 <- dist(y)
    sum(d1^2) * sd(d1) / mean(d1)
    })
  marker.score <- apply(propExpr, 1, function(x) {
    d1 <- dist(x)
    sum(d1^2) * sd(d1) / mean(d1) / max(max.scores)
  })
  marker.score[is.na(marker.score)] <- 0
  if(returnScore) return(marker.score)
  scoreRank = rank(-marker.score)
  return(scoreRank)
}


calc_beta <- function(y, spec.exp = 2) {
  d1 <- as.matrix(dist(y))
  eps1 <- 1e-10
  # Marker score is combination of specificity and sparsity
  score1 <- sum(d1^spec.exp) / (sum(d1) + eps1)
  return(score1)
}

#####################################################################
# FUNCTIONS FOR BUILDING AND PLOTTING THE TREE ARE BELOW

getDend <- function(input,distFun = function(x) return(as.dist(1-cor(x)))){
 distCor  = distFun(input) 
 avgClust = hclust(distCor,method="average")
 dend = as.dendrogram(avgClust)
 dend = labelDend(dend)[[1]]
}

labelDend <- function(dend,n=1)
  {  
    if(is.null(attr(dend,"label"))){
      attr(dend, "label") =paste0("n",n)
      n= n +1
    }
    if(length(dend)>1){
      for(i in 1:length(dend)){
        tmp = labelDend(dend[[i]], n)
        dend[[i]] = tmp[[1]]
        n = tmp[[2]]
      }
    }
    return(list(dend, n))
  }

reorder.dend <- function(dend, l.rank,verbose=FALSE)
  {
    tmp.dend = dend
    sc=sapply(1:length(dend), function(i){
      l = dend[[i]] %>% labels
      mean(l.rank[dend[[i]] %>% labels])
    })
    ord = order(sc)
	if(verbose){
      print(sc)
	  print(ord)
    }
	if(length(dend)>1){
      for(i in 1:length(dend)){
        if(ord[i]!=i){
          dend[[i]]= tmp.dend[[ord[i]]]
        }
        if(length(dend[[i]])>1){
          dend[[i]]=reorder.dend(dend[[i]],l.rank)
        }
      }
    }
    return(dend)
  }


# This function builds the tree and plots the dendrogram
buildAndPlotTree <- function(exprData,cl,l.rank,topNgenes = 1200,...){

  # Get median expression per cluster and the proportions
  normDat    = log2(exprData+1)
  medianExpr = do.call("cbind", tapply(names(cl), cl, function(x) rowMedians(normDat[,x]))) 
  rownames(medianExpr) = rownames(normDat) 
  medianExpr = medianExpr[apply(medianExpr,1,max)>0,]
  normDat    = normDat[rownames(medianExpr),]
  propExpr   = do.call("cbind", tapply(names(cl), cl, function(x) rowMeans(normDat[,x]>1))) 
  propExpr   = propExpr[,colnames(medianExpr)]
  rownames(propExpr) = rownames(medianExpr) 

  # Calculate the marker score and plot the dendrogram
  specificityScoreRank <- getSpecificityScore(propExpr)

  # Build and reorder the dendrogram
  dend <- getDend (medianExpr[specificityScoreRank<=topNgenes,])
  dend <- reorder.dend(dend,l.rank)
  dend <- collapse_branch(dend, 0.01)
  dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 2) 

  # Plot the dendrogram
  plot(dend,...)

}
  

