#' Single Cell Linkage by Distance Estimation is SLIDE
#'
#' @param infected A dataframe of protein expression levels in an infected subset of cells. All columns must be numeric expression levels.
#' @param uninfected A dataframe of protein expression levels in uninfected cells. All columns must be numeric expression levels. This can be interpreted as the population of cells, of which the infected cells are a subset.
#' @param cutoff This is the fold-change in distances being tested from an infected cell to its nearest uninfected and infected cell. The default value is 1.2. Determine an exact value using the \code{bootstrap_cutoff} function.
#' @return A list with containing five items. Non-parametric testing of up or downregulation in protein expression between the \code{infected} and \code{uninfected} cells, for both balanced and unbalanced comparisons. Three, it returns a boxplot mirroring these results. Four, the SLIDE Wilcox Rank Sum test results for evidence of cellular remodeling in the infected cells. Fifth, the mean ratio of the distances from each infected cell to its nearest uninfected cell and nearest infected cell.
#' @examples 
#' slide(I_sig, UN_sig)
#' @references
#' Sen, N., Mukherjee, G., and Arvin, A.M. (2015). Single cell mass cytometry reveals remodeling of human T cell phenotypes by varicella zoster virus. Methods 90, 85–94.

#' @importFrom grDevices dev.control dev.off pdf recordPlot rgb
#' @importFrom graphics boxplot par title
#' @importFrom stats p.adjust quantile wilcox.test
#' @export
slide<-function(infected ,uninfected, cutoff = 1.2)
{

  #Unncessecary now that input is pre-subsetted. ok for now
  V.sub<-infected
  UN.sub<-uninfected

  #obj skeleton to be iteratively filled
  UN.sub.1<-V.sub;

  #make sure there are equal columns in the comparisons
  if (ncol(V.sub) != ncol(UN.sub))
    stop("Infected and uninfected data frames must have equal number of columns")

  #make sure all columns are numeric inputs
  if (all(sapply(V.sub,class)=="numeric") == F | all(sapply(UN.sub,class)=="numeric") == F)
    stop("Dataframes must contain only numeric values")

  #check for numeric cols
  #try(apply())

  #returns zero distance values if there are less than 25 cells in subset
  if (dim(V.sub)[1] < 25)
  {
    li <- list('cell.dist' = rbind(0, 0),
               'protein.dist' = rbind(0, 0))
    return(list(li, warning = "Less than 25 cells in subset")) #double check this warning
  }


  #MY VERSION to test for speed USING SWEEP
  for (i in 1:nrow(V.sub))
  {
    temp<-V.sub[i,];
    temp1 <- apply(abs(sweep(UN.sub, 2, as.numeric(temp))),1,sum)
    UN.sub.1[i,]<-UN.sub[which.min(temp1),] ;
  }

  dist.v.protein.wise<-apply(abs(V.sub-UN.sub.1),1,mean) #distances to nearest uninf cell
  dist.v.cell.wise<-apply(abs(V.sub-UN.sub.1),2,mean)	#mean distance per protein


  ### Nested Function for non-parametric testing of protein expression
  test<-function(m1,m2)
  {
    t=matrix(0,ncol(m1),2);
    for (i in 1:ncol(m1)){
      x1=m1[,i]
      x2=m2[,i]
      t[i,]<-c(suppressWarnings(wilcox.test(x1,x2, alternative = "less")$p.value),
               suppressWarnings(wilcox.test(x1,x2, alternative = "greater")$p.value));
    }
    t1<-apply(t,1,min);
    t2<-apply(t,1,which.min);
    t3 <- ifelse(t2==1, "Down", "Up")
    return(data.frame(t1,t3))
  }

  ### Tests for sigDif Protein expression: Imbalanced comparison
  p.raw<-test(V.sub, UN.sub);
  p.adj<- p.adjust(p.raw[,1], 'BY')
  imbalPvals <- data.frame(names(UN.sub),round(p.adj, 16), round(p.raw[,1],16), p.raw[,2]);
  names(imbalPvals) <- c("Protein", "P.adj", "P.raw", "Infected Expression")
  ### #Tests for sigDif Protein expression: Balanced comparison
  p.raw<-test(V.sub, UN.sub.1);
  p.adj<- p.adjust(p.raw[,1], 'BY');
  balPvals <- data.frame(names(UN.sub.1),round(p.adj, 16), round(p.raw[,1],16), p.raw[,2]);
  names(balPvals) <- c("Protein", "P.adj", "P.raw", "Infected Expression")

  ### Second sweep
  UN.sub.2<-UN.sub.1

  for (i in 1:nrow(UN.sub.1))
  {
    temp<-UN.sub.1[i,];
    t2 <- apply(abs(sweep(UN.sub, 2, as.numeric(temp))),1,sum)
    ind.2<-which.min(t2);
    t3<-UN.sub[-c(ind.2),]
    ind.2<-which.min(t2[-c(ind.2)])
    UN.sub.2[i,]<-t3[ind.2,];
  }


  dist.u.protein.wise<-apply(abs(UN.sub.2-UN.sub.1),1,mean) #distances to nearest uninf cell
  dist.u.cell.wise<-apply(abs(UN.sub.2-UN.sub.1),2,mean) # mean distance per protein

  ### Plot differences in Protein Expression
  pdf(NULL)
  dev.control(displaylist="enable")
  par(bg="white",mfrow=c(2,1))

  boxplot(UN.sub,range=0,xlab='',ylab='',col='red',border='gray80',boxwex=0.8)
  boxplot(V.sub,range=0,xlab='',ylab='',add=T,col=rgb(0.5, 0.8, 0, 0.5),border='black',boxwex=0.25)
  title("Imbalanced Comparison (Green = inf, Red = uninf)")

  boxplot(UN.sub.1,range=0,xlab='',ylab='',col='red',border='gray80',boxwex=0.8)
  boxplot(V.sub,range=0,xlab='',ylab='',add=T,col=rgb(0.5, 0.8, 0, 0.5),border='black',boxwex=0.25)
  title("Balanced Comparison")

  plotExpress <- recordPlot()
  invisible(dev.off())

  #Distances
  li<-list('cell.dist'=rbind(dist.v.protein.wise,dist.u.protein.wise),'protein.dist'=rbind(dist.v.cell.wise,dist.u.cell.wise));

  #Run the rank sum test within function for now
  d1<-li$cell.dist[1,];
  d2<-li$cell.dist[2,];
  testRemodeling <- suppressWarnings(wilcox.test(d1/d2, cutoff, alternative = "greater"))
  meanRatio <- mean(d1/d2)

  #List with all relevant objects to return
  final <- list('imbalancedExpression'= imbalPvals,
                'balancedExpression' = balPvals,
                'plotExpression'=plotExpress,
                'testRemodeling'=testRemodeling,
                'meanRatio'=meanRatio)
  return(final)
}
