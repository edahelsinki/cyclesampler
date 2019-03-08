###########################################################################
## Copyright (c) 2018 Kai Puolamaki <kai.puolamaki@iki.fi>
## Copyright (c) 2018 University of Helsinki
## 
## Permission to use, copy, modify, and distribute this software for any
## purpose with or without fee is hereby granted, provided that the above
## copyright notice and this permission notice appear in all copies.
## 
## THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
## WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
## ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
## WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
## ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
## OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
###########################################################################


#' find values of lambdas in maxent method
#'
#' @param n number of values
#' @param m number of constraints
#' @param nu number of constraints per value, vector of length n
#' @param u vector of constraints
#' @param s vector of expectations (e.g., row/column sums)
#' @param l initial values of lambdas (not used)
#' @param tol accuracy of lambdas in root finding
#' @param tole accuracy of expectations
#' @param tolz z-normalized accuracy of expectations
#' @param maxiter maximum number of iterations
#' @param check boolean, whether to do sanity checks for inputs
#' @return list containing new lambdas, tolerance values (tol, tole,
#'     tolz) as well as numer of iterations done.
#'
#' @export
findlambdaC <- function(n,m,nu,u,s,l=rep(0.0,m),
                        tol=0.001,tole=-1.0,tolz=-1.0,
                        maxiter=1000,initlambda=TRUE,
                        check=TRUE) {

    
    if(check) {
        if(any(nu<1)) {
            stop("findlambdaC: empty constraint")
        }
        if(length(nu)!=n || length(u)!=sum(nu) ||
           length(s)!=m || length(l)!=m) {
            stop(sprintf("findlambdaC: length of a vector is incorrect.",
                         length(nu),n,length(u),sum(nu),length(s),length(l),m))
        }
        uniqu <- unique(u)
        if(length(uniqu)!=m || any(sort(uniqu)!=0:(m-1))) {
            stop("findlambdaC: constraints should be in 0:(m-1) with none empty.")
        }
        starts <- c(1,1+cumsum(nu[-length(nu)]))
        ends <- c(starts[-1]-1,length(u))
        if(any(mapply(function(i,j) any(duplicated(u[i:j])),starts,ends))) {
            stop("findlambdaC: duplicate items in a constraint.")
        }
    }
        
    res <- .C("findlambdasC",
              as.integer(n),             # 1
              as.integer(m),             # 2
              as.integer(nu),           # 3
              as.integer(u),             # 4
              as.double(s),              # 5
              as.double(l),              # 6
              as.double(tol),          # 7
              as.double(tole),        # 8
              as.double(tolz),        # 9
              as.integer(maxiter), # 10
              if(initlambda) as.integer(1) else as.integer(0))
    
    list(l=res[[6]],tol=res[[7]],tole=res[[8]],tolz=res[[9]],
         iter=maxiter-res[[10]])
}

#' Scales the values of data to (0,1) interval
#' 
#' @param x data vector of values
#' @return Data vectorscaled to (0,1) interval.
#'
#' @export
scaledata <- function(x,minv=min(x),maxv=max(x),drange=maxv-minv+1) {
    0.5/drange+(x-minv)/drange
}


#' Descales the values of data from (0,1) back to original format
#'
#' 
#' @param data data matrix where columns indicate row indices, column
#'     indices and values, respectively
#' @param toint if TRUE values are rounded to next integers.
#' @return Data matrix with value (third column) scaled to (0,1)
#'     interval.
#'
#' @export
descaledata <- function(x,minv=1,maxv=5,drange=maxv-minv+1,toint=TRUE) {
    x <- drange*(x-0.5/drange)+minv
    if(toint) x <- pmax(minv,pmin(maxv,round(x)))
    x
}


#' Makes "normalized" version of the data matrix such that the rows are
#' indexed from 1:nrow, columns from 1:ncol, and that there are no empty
#' rows or columns.
#'
#' @param data data matrix where columns indicate row indices, column
#'     indices and values, respectively
#' @param scale_values Should the values be scaled to the (0, 1)
#'     interval. Default is \code{FALSE}.
#' @return Data matrix with rows and columns re-indexed so that there
#'     are no empty rows and columns.
#'
#' @export
normalizedata <- function(data,
                          minv=min(data[,3]),
                          maxv=max(data[,3]),
                          drange=maxv-minv+1,
                          scale_values = FALSE) {
    data[,1] <- nonzero(data[,1])$p[data[,1]]
    data[,2] <- nonzero(data[,2])$p[data[,2]]+max(data[,1])

    if (scale_values)
        data[,3] <- scaledata(data[,3],minv=minv,maxv=maxv,drange=drange)

    data
}


#' Makes training/test data split
#'
#' @param data Triplets where columns indicate row indices, column
#'     indices and values, respectively
#' @return Training set and test set.
#'
#' @export
trtesplit <- function(data,nte=round(0.014*dim(data)[1]),limit=1) {
    ## Sample test data set of size nte
    idx <- sample.int(dim(data)[1],nte)
    ## Find rows/columns with at least limit items in training data
    data[,1] <- nonzero(data[-idx,1],nbins=max(data[,1]),
                        limit=limit)$p[data[,1]]
    data[,2] <- nonzero(data[-idx,2],nbins=max(data[,2]),
                        limit=limit)$p[data[,2]]
    ## Training data
    tr <- data[-idx,]
    ## Clean out entries which have less than limit rows/columns
    tr <- tr[!is.na(tr[,1]) & !is.na(tr[,2]),]
    ## Test data
    te <- data[idx,]
    ## Clean out entries which are in empty rows/columns in training data
    te <- te[!is.na(te[,1]) & !is.na(te[,2]),]
    list(tr=tr,te=te)
}


#' Finds mapping from indices so that they are mapped to continuous range
#' where each index has some instances. E.g., if s==c(1,3,4) we obtain
#' mapping p <- c(1,NA,2,3), i.e., p[s]==1:3.
#'
#' @param s List of indices (positive integers)
#' @return Mapping as described above.
#'
#' @export
nonzero <- function(s,nbins=max(s),limit=1) {
    idx <- (1:nbins)[tabulate(s,nbins=nbins)>=limit]
    p <- rep(NA,nbins)
    p[idx] <- 1:length(idx)
    list(p=p,idx=idx)
}


#' Samples of exponential distribution parametrized by x
#'
#' @param x (A vector of) parameters to exponential distribution
#' @return A vector of length length(x) of samples in interval [0,1].
#'
#' @export
sampleU01 <- function(x,n=length(x),y=runif(n),eps=.Machine$double.eps^0.25) {
    ifelse(abs(x)<eps,y,log(1+(exp(x)-1)*y)/x)
}


#' expectation of exponential distribution defined in [0,1]
#'
#' @param x parameter of distribution in interval [-Inf,Inf]
#' @return The expectation in interval [0,1].
#'
#' @export
fexpU01 <- function(x,eps=.Machine$double.eps^0.25) {
    ifelse(abs(x)<eps,0.5,(exp(x)-exp(x)/x+1/x)/(exp(x)-1))
}


#' Efficient implementation of a standard FIFO queue using cyclic
#' array.
#'
#' @param size maximum size of the queue.
#' @return Returns methods $isempty(), $push(x), and $eject().
#'
#' @export
makequeue <- function(size) {
    queue <- rep(NA,size)
    first <- last <- NULL
    list(
        isempty=function() { is.null(first) },
        push=function(x) {
            if(is.null(last)) {
                first <<- 1
                last <<- 1
            } else {
                last <<- if(last<size) last+1 else 1
            }
            queue[last] <<- x
        },
        eject=function() {
            if(is.null(first)) {
                NULL
            } else {
                res <- queue[first]
                if(first==last) {
                    first <<- NULL
                    last <<- NULL
                } else {
                    first <<- if(first<size) first+1 else 1
                }
                res
            }
        })
}


#' Performs Breadth-first search for a tree and finds a spanning tree.
#'
#' @param graph a list whose elements correspond to nodes and contain
#'     a vector of nodes neighbours.
#' @param root index of the root.
#' @return A vector which contains the index of the parent node. The
#'     root node's parent node has a index of zero.
#'
#' @export
BFS <- function(graph,root=sample.int(length(graph),size=1),
                Q=makequeue(length(graph)),
                bft=rep(NA,length(graph)),
                depth=rep(NA,length(graph)),
                st=rep(NA,length(graph)),
                treecount=1) {
    ## Example in Wikipedia, see
    ## https://en.wikipedia.org/wiki/Breadth-first_search
    ## graph <- list(c(2,3,5),c(1,6),c(1,7,8),8,c(1,10),c(2,9),3,c(3,4,10),c(6,10),c(5,8,9))
    ## BFS(graph,root=1)
    ## Should output: 0 1 1 8 1 2 3 3 6 5

    while(root>0) {
        bft[root] <- 0
        depth[root] <- 0
        st[root] <- treecount
        Q$push(root)
        while(!Q$isempty()) {
            current <- Q$eject()
            for(node in graph[[current]]) {
                if(is.na(bft[node])) {
                    bft[node] <- current
                    depth[node] <- depth[current]+1
                    st[node] <- treecount
                    Q$push(node)
                }
            }
        }
        root <- if(any(is.na(bft))) sample(which(is.na(bft)),size=1) else 0
        treecount <- treecount+1
    }
    list(bft=bft,depth=depth,st=st)
}


#' Finds the triplet indices that match the links in the spanning tree.
#'
#' @param bft spanning tree, e.g., output by BFS
#' @param triplets nX2 matrix that contains the links.
#' @return Vector of length length(bft) that contains the triplet
#'     indices that match the links in the spanning tree.
#'
#' @export
bft2idx <- function(bft,triplets) {
    bftpairs <- matrix(c(bft,1:length(bft)),length(bft),2)
    p <- bftpairs[,1]>0

    aux <- merge(
        data.frame(i=pmin(bftpairs[p,1],bftpairs[p,2]),
                   j=pmax(bftpairs[p,1],bftpairs[p,2]),
                   j1=1:sum(p)),
        data.frame(i=pmin(triplets[,1],triplets[,2]),
                   j=pmax(triplets[,1],triplets[,2]),
                   j2=1:dim(triplets)[1]),
        by=c("i","j"))

    idx    <- rep(NA,length(bft))
    idx[p] <- aux[order(aux[,"j1"]),"j2"]
    idx
}


#' Auxiliary function that indexes a vector.
#'
#' @param rows Vector of integers (e.g., row indices related to
#'     values), containing values in [1,m] where m is the number of
#'     constraints
#' @return A list whose i'th element is a list containing the indices
#'     of values related to the i'th constraint.
#'
#' @export
findrows <- function(rows) {
    p <- order(rows)
    u <- c((1:length(rows))[!duplicated(rows[p])],length(rows)+1)
    lapply(1:(length(u)-1),function(i) p[u[i]:(u[i+1]-1)])
}


#' Makes a graph of triplets.
#'
#' @param triplets nX2 matrix containing the links.
#' @return List where the i'th element is the list of nodes node i
#'     links to.
#'
#' @export
makegraph <- function(triplets) {
    tmp <- c(triplets[, 2], triplets[, 1])
    lapply(findrows(c(triplets[, 1], triplets[, 2])), function(i) tmp[i])
}


#' Finds 1 cycle from the spanning tree.
#'
#' @param bft Spannig tree.
#' @param depth depth of the items in spanning tree.
#' @param idx Triplet indices of the links in spanning tree.
#' @param idx0 Cycle link that does not appear in the spanning tree.
#' @param left triplets[idx0,1]
#' @param right triplets[idx0,2]
#' @return Cycle extracted from the spanning tree.
#'
#' @export
find1cycle <- function(bft,depth,idx,idx0,left,right) {
    p <- idx0
    m <- NULL
    fleft <- fright <- FALSE
    while(left!=right) {
        if(depth[left]>depth[right]) {
            if(fleft) {
                p <- c(p,idx[left])
            } else {
                m <- c(m,idx[left])
            }
            fleft <- !fleft
            left <- bft[left]
            if(left<1) stop("find1cycle: left<1!")
        } else {
            if(fright) {
                p <- c(p,idx[right])
            } else {
                m <- c(m,idx[right])
            }
            fright <- !fright
            right <- bft[right]
            if(right<1) stop("find1cycle: right<1!")
        }
    }
    c(p,m)
}

#' Finds a random term to add to the cycle.
#'
#' @param x random states in the cycle where the first half are the
#'     "positive" and the second half "negative" entries.
#' @return A random vector to number to add to the cycle.
#'
#' @export
delta1cycle <- function(x) {
    n <- length(x)/2
    p <- 1:n
    x[p] <- 1-x[p]
    runif(1,min=-min(x),max=min(1-x))*c(rep(-1,n),rep(1,n))
}


#' Samples n cycles and alters the random state accordingly
#'
#' @param n Number of cycles to sample
#' @param bft Spannig tree.
#' @param depth depth of the items in spanning tree.
#' @param idx Triplet indices of the links in spanning tree.
#' @param idx0 Vector of cycle links that do not appear in the spanning tree.
#' @param triplets nX2 matrix of triplest (links)
#' @param randomstate Current random state.
#' @return New random state.
#'
#' @export
samplecycles <- function(n,bft,depth,idx,idx0,triplets,randomstate) {
    for(j in 1:n) {
        i <- sample(idx0,size=1)
        cycle <- find1cycle(bft,depth,idx,i,triplets[i,1],triplets[i,2])
        randomstate[cycle] <- (
            randomstate[cycle]+delta1cycle(randomstate[cycle]))
    }
    randomstate
}


#' C interface that is functionally equivalent to
#' samplecycles. Samples n cycles and alters the random state
#' accordingly
#'
#' @param n Number of cycles to sample
#' @param bft Spannig tree.
#' @param depth depth of the items in spanning tree.
#' @param idx Triplet indices of the links in spanning tree.
#' @param idx0 Vector of cycle links that do not appear in the spanning tree.
#' @param triplets nX2 matrix of triplest (links)
#' @param randomstate Current random state.
#' @return New random state.
#'
#' @export
samplecyclesC <- function(n,bft,depth,idx,idx0,
                          triplets,randomstate) {
    idx[is.na(idx)] <- 0
    res <- .C("samplecyclesC",
              as.integer(n),               # 1
              as.integer(length(bft)),     # 2
              as.integer(length(idx0)),    # 3
              as.integer(bft-1),           # 4
              as.integer(depth),           # 5
              as.integer(idx-1),           # 6
              as.integer(idx0-1),          # 7
              as.integer(triplets[,1]-1),  # 8
              as.integer(triplets[,2]-1),  # 9
              as.double(randomstate))      # 10
    
    res[[10]]
}


#' C interface that is functionally equivalent to
#' samplecycles. Samples n cycles and alters the random state
#' accordingly
#'
#' @param n Number of cycles to sample
#' @param bft Spannig tree.
#' @param depth depth of the items in spanning tree.
#' @param st subtree indices
#' @param idx Triplet indices of the links in spanning tree.
#' @param idx0 Vector of cycle links that do not appear in the
#'     spanning tree.
#' @param slength length of cycles corresponding to idx0
#' @param odds array with the number of odd cycles for each subgraph
#' @param triplets nX2 matrix of triplest (links)
#' @param randomstate Current random state.
#' @param a lower bound for the weight of each triplet
#' @param b upper bound for the weight of each triplet
#' @param what define what to return (string). One of \code{all} returning all
#'     parameters, \code{randomstate} returning the new random state,
#'     or \code{acceptance_rate} returning the fraction of accepted
#'     moves.
#' @return New random state.
#'
#' @export
samplecycles2C <- function(n,bft,depth,st,idx,idx0,slength,
                           odds,triplets,
                           randomstate,
                           a=rep(0,dim(triplets)[1]),
                           b=rep(1,dim(triplets)[1]), what = "randomstate")
{
    idx[is.na(idx)] <- 0

    tmp <- .C("samplecycles2C",
       as.integer(n),                      # 1
       as.integer(length(bft)),            # 2
       as.integer(dim(triplets)[1]),       # 3
       as.integer(length(idx0)),           # 4
       as.integer(bft-1),                  # 5
       as.integer(depth),                  # 6
       as.integer(st-1),                   # 7
       as.integer(idx-1),                  # 8
       as.integer(idx0-1),                 # 9
       as.integer(slength),                # 10
       as.integer(sapply(odds,length)),    # 11
       as.integer(do.call("c",odds)-1),    # 12
       as.integer(triplets[,1]-1),         # 13
       as.integer(triplets[,2]-1),         # 14
       as.double(a),                       # 15
       as.double(b),                       # 16
       as.integer(0),                      # 17 - accepted moves
       as.double(randomstate))             # 18 - random state

    switch(what,
           "all"             = tmp,
           "randomstate"     = tmp[[18]],
           "acceptance_rate" = tmp[[17]] / n)
}


#' Finds node at which the tree defined by bft merges if starting from
#' nodes left and right.
#'
#' @param bft Spannig tree.
#' @param depth depth of the items in spanning tree.
#' @param left node
#' @param right node that is connected to left
#' @return Node at which the trees with leafs at left and right merge.
#'
#' @export
findmerge <- function(bft,depth,left,right) {
    while(left!=right) {
        if(depth[left]<depth[right]) {
            right <- bft[right]
        } else {
            left <- bft[left]
        }
    }
    left
}


#' Cyclesampler object.
#'
#' @param data n X 3 matrix containing triplets (first 2 columns and
#'     the random state (data[, 3]). The data matrix should be
#'     normalized so that the row and column indices are continous,
#'     i.e., there is no node to which there are no links.
#' @return Returns two function $getstate() that gives the current
#'     random state and $samplecycles(n) that samples n cycles.
#'
#' @export
cyclesampler <- function(data,
                         a=rep(0, dim(data)[1]),
                         b=rep(1, dim(data)[1]),
                         useC=TRUE,
                         nukezeros=TRUE) {

    triplets <- as.matrix(data[,1:2])

    if(nukezeros)
        triplets <- matrix(nonzero(triplets)$p[triplets], dim(triplets)[1], 2)

    randomstate <- data[,3]

    if (any(randomstate < a) | any(randomstate > b))
        stop("cyclesampler: edge weights outside bounds.\n")

    graph <- makegraph(triplets)

    ## order the graph according to node strengths
    nw      <- get_node_weights(data)
    graph   <- lapply(graph, function(nl) nl[order(nw[nl], decreasing = TRUE)])
    aux     <- BFS(graph, root = which.max(nw))

    bft     <- aux$bft
    depth   <- aux$depth
    st      <- aux$st
    idx     <- bft2idx(bft,triplets)
    idx0    <- setdiff(1:dim(triplets)[1],idx)
    idx0st  <- st[triplets[idx0,1]]

    slength <- sapply(idx0,
                      function(i) 
                          (1+depth[triplets[i,1]]+depth[triplets[i,2]]
                              -2*depth[findmerge(bft,depth,
                                                 triplets[i,1],
                                                 triplets[i,2])]))
    idx0odd <- slength%%2==1
    odds    <- lapply(1:max(st),function(i) idx0[idx0odd & idx0st==i])
    f       <- if(useC) samplecyclesC else samplecycles

    list(
        getstate = function() { randomstate },
        setstate = function(rs) { randomstate <<- rs },
        samplecycles = function(n=1) {
            randomstate <<- f(n,bft,depth,idx,idx0,triplets,randomstate)
        },
        samplecycles2 = function(n=1, what = "all") {
        tmp <- samplecycles2C(
                n,bft,depth,st,idx,idx0,slength,odds,triplets,randomstate,
            a=a,b=b, what = what)
        randomstate <<- tmp[[18]]
        acceptance_rate <<- tmp[[17]] / n
    },

    getncycles = function() { length(idx0) },
    getacceptancerate = function() { acceptance_rate }
    )
}


#' Maxentsampler object.
#'
#' @param data nX3 matrix containing triplets (first 2 columns and the
#'     random state (data[,3]). The data matrix should be normalized
#'     so that the random states are in the interval [0,1] and that
#'     the row and column indices are continous, i.e., there is no
#'     node to which there are no links.
#' @return Returns functions $optimlambda() that finds the parameters
#'     lambda, $sample() that produces one sample, and $list() that gives
#'     a listing of internal parameters.
#' 
#' @export
maxentsampler_old <- function(data) {
    triplets <- data[,1:2]
    triplets[,2] <- triplets[,2]-min(triplets[,2])+1
    if(any(triplets<1) || any(data[,3]<0 | 1<data[,3]))
        stop("maxentsampler: invalid data.\n")
    n <- dim(triplets)[1]
    rows <- findrows(triplets[,1])
    cols <- findrows(triplets[,2])
    m <- length(rows)+length(cols)
    rsums <- sapply(rows,function(i) sum(data[i,3]))
    csums <- sapply(cols,function(i) sum(data[i,3]))
    lambdar <- rep(0,length(rows))
    lambdac <- rep(0,length(cols))

    list(
        optimlambda=function(tol=0.001,tole=-1.0,tolz=-1.0,maxiter=1000,
                             initlambda=FALSE) {
            res <- findlambdaC(n=n,
                               m=m,
                               nu=rep(2,n),
                               u=c(t(triplets[,1:2]
                                     +(rep(1,n) %o% c(-1,length(rows)-1)))),
                               s=c(rsums,csums),
                               l=c(lambdar,lambdac),
                               tol=tol,
                               tole=tole,
                               tolz=tolz,
                               maxiter=maxiter,
                               initlambda=initlambda)
            
            lambdar <<- (res$l)[1:length(rows)]
            lambdac <<- (res$l)[-(1:length(rows))]
            ## The probability distribution is not changed if we add a
            ## constant to the row lambdas and decrease the same
            ## constant from columns lambdas. This means that even if
            ## the resulting distribution is same we can have
            ## different sets of values for lambdas. To make lambdas
            ## unique we choose to have the mean of row lambdas zero,
            ## i.e., we add to row lambdas minus their mean and add
            ## the mean of row lambdas to column lambdas,
            ## respectively.
            mlambdar <- mean(lambdar)
            lambdar <<- lambdar-mlambdar
            lambdac <<- lambdac+mlambdar
            res
        },
        sample=function() {
            sampleU01(lambdar[triplets[,1]]+lambdac[triplets[,2]])
        },
        exp=function() {
            fexpU01(lambdar[triplets[,1]]+lambdac[triplets[,2]])
        },
        list=function() {
            list(lambdar=lambdar,lambdac=lambdac,rsums=rsums,csums=csums)
        }
    )
}


#' Maxentsampler2 object.
#'
#' @param data nX3 matrix containing triplets (first 2 columns and the
#'     random state (data[,3]). The data matrix should be normalized
#'     so that the random states are in the interval [0,1] and that
#'     the row and column indices are continous, i.e., there is no
#'     node to which there are no links.
#' @return Returns functions $optimlambda() that finds the parameters
#'     lambda, $sample() that produces one sample, and $list() that gives
#'     a listing of internal parameters.
#' 
#' @export
maxentsampler <- function(data) {
    triplets <- as.matrix(data[,1:2])
    rs <- data[,3]
    rs2 <- rep(rs,2)
    uu <- sort(unique(c(triplets)))

    if(any(rs<0) || any(1<rs) || any(uu!=1:length(uu)))
        stop("maxentsampler: invalid data.\n")

    n <- dim(triplets)[1]
    rows <- findrows(c(triplets[,1],triplets[,2]))
    m <- length(rows)
    rsums <- sapply(rows,function(i) sum(rs2[i]))
    lambdar <- rep(0,m)

    list(
        optimlambda=function(tol=0.001,tole=-1.0,tolz=-1.0,maxiter=1000,
                             initlambda=FALSE) {
            res <- findlambdaC(n=n,
                               m=m,
                               nu=rep(2,n),
                               u=c(t(triplets))-1,
                               s=rsums,
                               l=lambdar,
                               tol=tol,
                               tole=tole,
                               tolz=tolz,
                               maxiter=maxiter,
                               initlambda=initlambda)
            
            lambdar <<- res$l
            res
        },
        sample=function() {
            sampleU01(lambdar[triplets[,1]]+lambdar[triplets[,2]])
        },
        exp=function() {
            fexpU01(lambdar[triplets[,1]]+lambdar[triplets[,2]])
        },
        list=function() {
            list(lambdar=lambdar,rsums=rsums)
        }
    )
}


#' addselfloops
#'
#' @param data nX3 matrix containing triplets (first 2 columns and the
#'     random state (data[,3]). The data matrix should be normalized
#'     so that 
#'     the row and column indices are continous, i.e., there is no
#'     node to which there are no links. A and B contain the lower
#'     and upper bounds for the vertex weights, respectively.
#' @return Returns self-loops and lower and upper bounds for their
#'     weights.
#' 
#' @export
addselfloops <- function(data,A,B) {
  triplets <- as.matrix(data[,1:2])
  rs <- data[,3]
  rs2 <- rep(rs,2)
  uu <- sort(unique(c(triplets)))
  if(any(uu!=1:length(uu))) 
    stop("addselfloops: invalid data.\n")
  idx <- which(!is.na(A) & !is.na(B))
  rows <- findrows(c(triplets[,1],triplets[,2]))
  W <- sapply(rows[idx],function(i) sum(rs2[i]))
  a <- W-B[idx]
  b <- W-A[idx]
  if(any(0<a) || any(b<0))
    stop("addselfloops: inconsistent vertex constraint.\n")
  list(data=matrix(c(uu[idx],uu[idx],rep(0,length(idx))),
                    length(idx),3),
       a=a,b=b)
}

#' Get node weights
#'
#' @param data A 3-column matrix with the triplets.
#'
#' @return The Frobenius norm for each of the n_samples with respect to the starting state.
#'
#' @export
get_node_weights <- function(data) {
    rs2  <- rep(data[, 3], 2)
    rows <- findrows(c(data[, 1], data[, 2]))
    sapply(rows, function(i) sum(rs2[i]))
}
