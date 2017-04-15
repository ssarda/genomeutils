#' @title Hclust in colour
#' @description Prints pretty hclust plots
#' @param hclust hclust object
#' @param labels a character vector of labels of the leaves of the tree
#' @param lab.col colour for the labels; NA=default device foreground colour
#' @param xlab title for x-axis
#' @param sub subtitle
#' @return Function
#' 
#' @export
color_hclust <- function( hclust, labels=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, xlab="", sub="")
{
    y <- rep(hclust$height,2)
    x <- as.numeric(hclust$merge)
    y <- y[which(x<0)]
    x <- x[which(x<0)]
    x <- abs(x)
    y <- y[order(x)]
    x <- x[order(x)]
    plot( hclust, labels=FALSE, hang=hang, xlab=xlab, sub=sub )
    text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=labels[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA )
}





#' @title MA plot
#'
#' @description Takes two vectors x and y and plots 
#' M=y-x versus A=(x+y)/2. 
#' A smooth curve is added to show trends.
#' 
#' @param x a numeric vector
#' @param y a numeric vector
#' @param n a numeric value to subsample data
#' @param subset index of the points to be plotted 
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param curve.add if TRUE a smooth curve is fit to the data and displayed. The function loess is used to fit the curve.
#' @param curve.col a numeric value that determines the color of the smooth curve
#' @param curve.span is passed on to loess
#' @param curve.lwd the line width for the smooth curve
#' @param curve.n a numeric value that determines the sample size used to fit the curve. This makes fitting the curve faster with large datasets
#' @param ... further arguments passed to plot
#' @return Function
#' @export
ma_plot <- function(x,y,n=10000,subset=NULL,xlab=NULL,ylab=NULL,
                   curve.add=TRUE,curve.col=2,curve.span=1/2,
                   curve.lwd=2,curve.n=2000,...)
{
    if(length(x)!=length(y)) stop("Length of 'x' and 'y' must be the same.")  
    if(is.null(xlab)) xlab="A"
    if(is.null(ylab)) ylab="M"
    if(!is.null(subset))
	{
        x=x[subset]
        y=y[subset]
    }
    if(length(x)>n)
	{
        ind=sample(length(x),10000)
        x=x[ind]
        y=y[ind]
    }
    m=y-x
    a=(x+y)/2
    plot(a,m,xlab=xlab,ylab=ylab,...)
    if(curve.add)
	{
        o=order(a)
        aa=a[o]
        mm=m[o]
        o=seq(1,length(aa),len=curve.n)
        aa=aa[o]
        mm=mm[o]
        fit1=loess(mm~aa,span=curve.span,degree=1)
        lines(aa,fit1$fitted,col=curve.col,lwd=curve.lwd)
    }
}


#' @title Fancy Heatmap
#' @description Plots a heatmap of clstering within a Gene x Sample Matrix of top n genes.
#' @param X Gene x Sample Matrix
#' @param top Pick top genes with highest variability
#' @param ... Additional arguments to pass to heatmap.2
#' @return Function
#' @export
fancy_heatmap <- function(X, top, ...)
{
	hmcol <- RColorBrewer::colorRampPalette(brewer.pal(9, "GnBu"))(100)
	rv <- genefilter::rowVars(X)
	idx <- order(-rv)[1:top]
	heatmap.2(X[idx,],trace="none", col=hmcol, ...)
}



#' @title Smooth Histogram
#' @description a smooth histogram with unit indicator. The advantage is interpretability. If input is a matrix, curves are plotted together.
#' @param z the data
#' @param unit the unit which determines the y axis scaling and is drawn
#' @param bw arguments to density
#' @param n arguments to density
#' @param from arguments to density
#' @param to arguments to density
#' @param plotHist a logical: should an actual histogram be drawn under curve?
#' @param add a logical: add should the curve be added to existing plot?
#' @param xlab x-axis title, defaults to no title
#' @param ylab y-axis title, defaults to no title
#' @param xlim range of the x-axis
#' @param ylim range of the y-axis
#' @param main an overall title for the plot
#' @return Function
#' @export
smooth_hist <- function(z, unit,
                  bw="nrd0", n, from, to,
                  plotHist = FALSE, add = FALSE,
                  xlab, ylab = "Frequency",
                  xlim, ylim, main)
{

  if (is.data.frame(z) && is.numeric(as.matrix(z)))
  {
    z <- as.matrix(z)
  }
  if (is.matrix(z) && ncol(z) == 1)
  {
    z <- as.vector(z)
  }
  stopifnot(is.vector(z) | is.matrix(z))
  
  if (missing(xlab))
  {
    xlab <- deparse(substitute(z))
  }
  if (missing(main))
  {
    main <- paste("Smooth histogram of", deparse(substitute(z)))
  }
  
  if (is.vector(z))
  {
    if (missing(unit))
	{
      unit <- .5 * sd(z)
    }
  }
  else if (is.matrix(z))
  {
    if (missing(unit))
	{
      unit <- .5 * mean(apply(z, 2, sd))
    }
  }

  smooth_hist.onesample <- function(z, unit, bw, n, from, to, plotHist, add, xlab, ylab,xlim, ylim, main, maxz, ...)
{
    lz <- length(z)
    if (missing(n))
	{
      n <- 512
    }
    d <- density(z, bw=bw, n=n, from=from, to=to)
    ymax <- max(d$y * lz * unit)
    if (missing(xlim))
	{
      xlim <- c(min(z) - unit, max(z) + unit)
    }
    if (missing(ylim))
	{
      ylim <- c(0, 1.2 * ymax)
    }
    if (plotHist)
	{
      h <- hist(z, breaks = seq(from = min(z) - unit, to = max(z) + unit, by = unit), plot = FALSE)
      ylim[2] <- max(ylim[2], 1 * 1 * h$count)
      plot(h, col = "grey", main = main, xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim)
    }
	else
	{
      if (!add) {
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
      }
    }

    lines(d$x, d$y * lz * unit)
    if (!add) {
      arrows(maxz - unit, 1.1*ymax, maxz, 1.1*ymax, angle = 90, code = 3, length = 0.1)
      text(maxz - unit, 1.1*ymax, paste("unit =", sprintf("%.2g", unit)), pos = 2)
    }
  }

  if (is.vector(z))
  {
    smooth_hist.onesample(z=z, unit=unit, bw=bw, n=n, from=from, to=to, plotHist=plotHist, add=add, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, maxz=max(z))
  }
  else if (is.matrix(z))
  {
    if (missing(xlim))
	{
      xlim <- c(min(z) - unit, max(z) + unit)
    }
    smooth_hist.onesample(z=z[,1], unit=unit, bw=bw, n=n, from=from, to=to,  plotHist=FALSE, add=FALSE, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, maxz=max(z))
    for (j in 2:ncol(z))
	{
      smooth_hist.onesample(z=z[,j], unit=unit, bw=bw, n=n, from=from, to=to, plotHist=FALSE, add=TRUE)
    }
  }
}
