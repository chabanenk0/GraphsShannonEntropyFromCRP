library(sna)

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
    if (doplot) # {
        plot(0, xlim = c(0, end.time), ylim = c(0, end.time), type = "n",
            main = "Recurrence Plot", xlab = "i", ylab = "j")
        res=matrix(0,end.time,end.time)
        for(i in seq(1, end.time, by = nt))
            for(j in seq(i,end.time, by = nt))
                if(.dist(i,j) < eps) 
                 {
                  if (doplot) points(c(i, j), c(j, i), ...)
                  res[i,j]=1
                  res[j,i]=1}
                else
                 {res[i,j]=1
                  res[j,i]=1}
                }
    #}
    # Return Value:
    # invisible()
    return (res)
}



library(fNonlinear)
y <- read.table("/Users/Андрей/Dropbox/SolovievChabCommonFolder/sp_91.txt")
y=t(y)
wind=250
m=1	#   m - embedding dimension
d=1 	#   d - time delay
end.time=wind    	#   end.time - ending time (as no. of observations)
eps=0.1    	#   eps - neighbourhood threshold
nt=1    	#   nt - observations in each step
logfile="/Users/Андрей/Dropbox/SolovievChabCommonFolder/sp_91_logsR.txt"
n=length(y)
mas_ge=matrix(0,n-wind,2)
fileConn<-file(logfile, open="w")
writeLines("t\tGraphEntropy\n", fileConn)
close(fileConn)



for (i in seq(1,100,by=10)) #n-wind)
{
 print (i)
 y_fragm=y[1,i:(i+wind)]
 meanvalue=mean(y_fragm)
 stand_dev=sd(y_fragm)
 y_fragm=(y_fragm-meanvalue)/stand_dev
 # plot(y_fragm)
 adj = recurrencePlotNew(y_fragm,m,d,end.time,eps,nt,doplot=FALSE)
 ge = graphEntropy (adj)
 fileConn<-file(logfile, open="a")
 writeLines(paste(toString(i),toString(ge),sep="\t"), fileConn)
 close(fileConn)


 mas_ge[i,1]=i
 mas_ge[i,2]=ge
}
plot(mas_ge[,1],mas_ge[,2])
