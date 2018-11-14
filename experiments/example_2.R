## -----------------------------------------------------------------------------
## Example 2
##
##       Example illustrating the difference between CycleSampler
##       and a Maximum Entropy model.
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla example_2.R
##         
## -----------------------------------------------------------------------------

## --------------------------------------------------
## Load libraries and set the random seed
## --------------------------------------------------

library(cyclesampler)
set.seed(42)

## --------------------------------------------------
## The smallest network 1-2-3 with weights of 0.3 and 0.6
## --------------------------------------------------

data <- matrix(c(1,2,2,3,0.3,0.6),2,3)

# Allowed edge weights in [0,1]
a <- c(0,0)
b <- c(1,1)

## Allow node weights in [0.25,1]
A <- c(0.25,0.25,0.25)
B <- c(1.5,1.5,1.5)

## Append data with loops
aux <- addselfloops(data,A,B)
datal <- rbind(data,aux$data)
al <- c(a,aux$a)
bl <- c(b,aux$b)


## --------------------------------------------------
## Make a maxent model
## --------------------------------------------------
Y <- maxentsampler(data)
junk <- Y$optimlambda(tol=0)

## Obtain 10000 samples from weight of vertex 2
samples_me <- replicate(1000,Y$sample())

## --------------------------------------------------
## Make a cycle model with loops
## --------------------------------------------------
X <- cyclesampler(datal,a=al,b=bl)
samples_cy <- replicate(1000,{
    X$samplecycles2(1000)
    X$getstate()[1:2]
})

## --------------------------------------------------
## Plot figures
## --------------------------------------------------

pdf("fig1a.pdf")
plot(c(0,1),c(0,1),bty="n",type="n",
     main="",
     xlab=expression(paste(w^"*")(group("{",list(1,2),"}"))),
     ylab=expression(paste(w^"*")(group("{",list(2,3),"}"))))
points(t(samples_me),pch=4)
points(0.3,0.6,col="red",pch=19)
lines(c(0.25,1,1,0.5,0.25,0.25),c(0.25,0.25,0.5,1,1,0.25),col="red",lwd=2)
dev.off()

pdf("fig1b.pdf")
plot(c(0,1),c(0,1),bty="n",type="n",
     main="",
     xlab=expression(paste(w^"*")(group("{",list(1,2),"}"))),
     ylab=expression(paste(w^"*")(group("{",list(2,3),"}"))))
points(t(samples_cy),pch=4)
points(0.3,0.6,col="red",pch=19)
lines(c(0.25,1,1,0.5,0.25,0.25),c(0.25,0.25,0.5,1,1,0.25),col="red",lwd=2)
dev.off()

## --------------------------------------------------