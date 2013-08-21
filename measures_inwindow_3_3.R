# Перед запуском нужно установить эти 3 пакета запуском команд:
# install.packages("sna")
# install.packages("QuACN")
# install.packages("fNonlinear")
# Пакеты graph нужно скачать вручную из http://www.bioconductor.org/packages/release/bioc/html/graph.html
# http://www.bioconductor.org/packages/2.12/bioc/html/RBGL.html
# http://www.bioconductor.org/packages/2.12/bioc/html/BiocGenerics.html
# Архивы также сохранены в папке RequiredStrangePackages
# и установить командой:
# install.packages(file.choose(), repos=NULL)
# указав путь к зип-файлу пакета.

library(sna)
library(QuACN)
library(fNonlinear)

recurrencePlotNew =
function(x, m, d, end.time, eps, nt = 10, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a recurrence plot

    # Arguments
    #   x - time series
    #   m - embedding dimension
    #   d - time delay
    #   end.time - ending time (as no. of observations)
    #   eps - neighbourhood threshold
    #   nt - observations in each step
    #   ... - further parameters to be passed to plot

    # Value:
    #   Produces the recurrence plot, as proposed by Eckmann et al. (1987).
    #   To reduce the number of points plotted (especially with highly
    #   sampled data), each nt observations, one single point is plotted.

    # FUNCTION:

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    w = (0:(m-1))*d
    .dist = function(i, j) { sum((series[i+w]-series[j+w])^2) }

    .checkEmbParms(series, m, d)
    if (eps <= 0) stop("eps must be positive")
    nt = as.integer(nt)
    if (nt<=0) nt = 1
    n = length(series)-(m-1)*d
    if(end.time > n) end.time = n
    eps = eps^2
    xyz = .embeddPSR(series, m = m, d = d)[1:end.time, ]

    # Plot:
    if (doplot)  {
        plot(0, xlim = c(0, end.time), ylim = c(0, end.time), type = "n",
            main = "Recurrence Plot", xlab = "i", ylab = "j")}
        res=matrix(0,end.time,end.time)
        for(i in seq(1, end.time, by = nt))
            for(j in seq(i,end.time, by = nt))
                if(.dist(i,j) < eps) 
                 {
                  if (doplot) points(c(i, j), c(j, i), ...)
                  res[i,j]=1
                  res[j,i]=1}
                else
                 {res[i,j]=0
                  res[j,i]=0}
               # }
    #}
    # Return Value:
    # invisible()
    return (res)
}

graphEntropy <- function(adj, type="SoleValverde") {
  if (type == "SoleValverde") {
    return(graphEntropySoleValverde(adj))
  }
  else {
    return(graphEntropyWang(adj))
  }
}

LacasaAdjMatrix =
function(x, end.time, doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates an adjancency matrix for Lacasa visibility method

    # Arguments
    #   x - time series
    #   end.time - ending time (as no. of observations)
    #   ... - further parameters to be passed to plot

    # Settings:
    if (class(x) == "timeSeries") x = as.vector(x)
    series = as.ts(x)
    # w = (0:(m-1))*d
    .visib = function(i, j) {
     a1=(series[j]-series[i])/(j-i)
     a0=series[i]-a1*i
     if (i+1<j)
     for (ii in seq(i+1,j-1,by=1))
     { 
       if (series[ii]>a0+ii*a1)
        {
         return (-1)
        } 
     }
     return (1)
    }

    #.checkEmbParms(series, m, d)
    #if (eps <= 0) stop("eps must be positive")
    nt = 1 # as.integer(nt)
    # if (nt<=0) nt = 1
    # n = length(series)-(m-1)*d
    if(end.time > n) end.time = n
    # eps = eps^2
    # xyz = .embeddPSR(series, m = m, d = d)[1:end.time, ]

    # Plot:
    if (doplot)  {
        plot(0, xlim = c(0, end.time), ylim = c(0, end.time), type = "n",
            main = "Recurrence Plot", xlab = "i", ylab = "j")}
        res=matrix(0,end.time,end.time)
        for(i in seq(1, end.time, by = nt))
            for(j in seq(i,end.time, by = nt))
                if(.visib(i,j) > 0) 
                 {
                  if (doplot) points(c(i, j), c(j, i), ...)
                  res[i,j]=1
                  res[j,i]=1}
                else
                 {res[i,j]=0
                  res[j,i]=0}
               # }
    #}
    # Return Value:
    # invisible()
    return (res)
}


graphEntropySoleValverde <- function(adj) {
  # Calculate Sole & Valverde, 2004 graph entropy
  # Uses Equations 1 and 4
  # First we need the denominator of q(k)
  # To get it we need the probability of each degree
  # First get the number of nodes with each degree
  library(sna)
  existingDegrees = sna::degree(adj)/2
  maxDegree = nrow(adj) - 1
  allDegrees = 0:maxDegree
  degreeDist = matrix(0, 3, length(allDegrees)+1) # Need an extra zero prob degree for later calculations
  degreeDist[1,] = 0:(maxDegree+1)
  for(aDegree in allDegrees) {
    degreeDist[2,aDegree+1] = sum(existingDegrees == aDegree)
  }
  # Calculate probability of each degree
  for(aDegree in allDegrees) {
    degreeDist[3,aDegree+1] = degreeDist[2,aDegree+1]/sum(degreeDist[2,])
  }
  # Sum of all degrees mult by their probability
  sumkPk = 0
  for(aDegree in allDegrees) {
    sumkPk = sumkPk + degreeDist[2,aDegree+1] * degreeDist[3,aDegree+1]
  }
  # Equivalent is sum(degreeDist[2,] * degreeDist[3,])
  # Now we have all the pieces we need to calculate graph entropy
  graphEntropy = 0
  for(aDegree in 1:maxDegree) {
    q.of.k = ((aDegree + 1)*degreeDist[3,aDegree+2])/sumkPk
    # 0 log2(0) is defined as zero
    if (q.of.k != 0) {
      graphEntropy = graphEntropy + -1 * q.of.k * log2(q.of.k)
    }
  }
  return(graphEntropy)
}

graphEntropyWang <- function(adj) {
  # Calculate Wang, 2008 graph entropy
  # Uses Equation 14
  # bigN is simply the number of nodes
  # littleP is the link probability.  That is the same as graph density calculated by sna with gden().
  bigN = nrow(adj)
  littleP = gden(adj)
  graphEntropy = 0
  if (littleP != 1 && littleP != 0) {
    graphEntropy = -1 * .5 * bigN * (bigN - 1) * (littleP * log2(littleP) + (1-littleP) * log2(1-littleP))
  }
  return(graphEntropy)
}



y <- read.table("dax_04.txt")
y=t(y)
wind=500
m=1	#   m - embedding dimension
d=1 	#   d - time delay
end.time=wind    	#   end.time - ending time (as no. of observations)
eps=0.1    	#   eps - neighbourhood threshold
nt=1    	#   nt - observations in each step
graph_type=2 # 1 - crp, 2 - visibility
logfile="dax_04_logsR_v3_wind_500.txt"

n=length(y)
mas_ge=matrix(0,n-wind,2)
fileConn<-file(logfile, open="w")
#writeLines("t\tGraphEntropy\n", fileConn)
writeLines("t\ttopologicalEntr\tbertz\tradialCentric\tvertexDegree\tgraphVertexComplexity\tedgeEqualityMIC\tedgeMagnitudeMIC\tsymmetryIndex\tdistanceDegreeMIC\tdistanceDegreeEquality\n", fileConn)
close(fileConn)



for (i in seq(1,n-wind,by=10)) #n-wind)
{
 print (i)
 y_fragm=y[1,i:(i+wind)]
 meanvalue=mean(y_fragm)
 stand_dev=sd(y_fragm)
 y_fragm=(y_fragm-meanvalue)/stand_dev
 # plot(y_fragm)
 # graph_type= 2 # 1 - crp, 2 - visibility
 if (graph_type==1)
  { adj = recurrencePlotNew(y_fragm,m,d,end.time,eps,nt,doplot=FALSE)
  }
 else 
  {
  adj = LacasaAdjMatrix (y_fragm,end.time,doplot=FALSE)
  }
 
 graphnelg=as(adj,"graphNEL")
 topinfc = topologicalInfoContent(graphnelg)
 entrp=topinfc$entropy
 #bon1=bonchev1(graphnelg)
 #bon2=bonchev2(graphnelg)
 #bon3=bonchev3(graphnelg)
 berz=bertz(graphnelg)
 rc=radialCentric(graphnelg)
 vdi=vertexDegree(graphnelg)
 #bal1=balabanlike1(graphnelg)
 #bal2=balabanlike2(graphnelg)
 gvc=graphVertexComplexity(graphnelg)
 #gdc=graphDistanceComplexity(graphnelg)
 #ibi=informationBondIndex(graphnelg)
 edgEqMic=edgeEqualityMIC(graphnelg)
 edgMgMic=edgeMagnitudeMIC(graphnelg)
 symmInd=symmetryIndex(graphnelg)
 ddmic=distanceDegreeMIC(graphnelg)
 ddeq=distanceDegreeEquality(graphnelg)
 #ddcomp=distanceDegreeCompactness(graphnelg)
 #infLayerInd=informationLayerIndex(graphnelg)


 fileConn<-file(logfile, open="a")
 writeLines(paste(toString(i),toString(entrp),toString(berz),toString(rc),toString(vdi),toString(gvc),toString(edgEqMic),toString(edgMgMic),toString(symmInd),toString(ddmic),toString(ddeq),sep="\t"), fileConn)
 close(fileConn)

 mas_ge[i,1]=i
 mas_ge[i,2]=ge
}
plot(mas_ge[,1],mas_ge[,2])
