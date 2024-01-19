
# first prepare for the data
require(metafor)
require( mixmeta )
require(cubature)
require(parallel)
require(tidyr)
require(dplyr)
require(pracma)
options(warn = -1)

data.x <- data.frame(
    fp = c(3, 3, 9, 7, 6, 4, 3, 2, 19, 2, 7, 1, 5, 11, 0, 1, 3, 1), 
    tp = c(0, 1, 2, 0, 5, 1, 1, 1, 1, 1, 0, 0, 3, 6, 0, 0, 1, 4),
    n0 = c(117,  35, 195, 136, 157, 139, 177,  39, 103, 122,  64,  58, 175,
           180, 105, 262, 362,  69),
    n1 = c(116,  44, 208, 130, 151,  98, 174,  74,  97, 113,  66,  70, 188,
           187, 118, 252, 345,  64)
) %>% dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>% dplyr::mutate( yall = fp + tp ) %>%
    dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ) | (fp == 0 ) | (tn == 0) )*0.5 ) %>%
    dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ), 
                   v = 1/( tp + cx ) + 1/( tn + cx ) + 
                       1/( fn + cx ) + 1/( fp + cx )) %>%
    dplyr::mutate( tmp = y/sqrt(v) ) %>%
    dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )

nstudy = nrow(data.x)

# obtain the estimate with the reml
# est.o <- mixmeta(y, S=v, data=data.x, method = "ml")
# est.o.vec <- c(est.o$coefficients, sqrt(est.o$Psi)) %>% `names<-`(NULL)

data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( tp = 0:data.x$yall[ind], yall = data.x$yall[ind], n1 = data.x$n1[ind], n0 = data.x$n0[ind]  ) %>%
        dplyr::mutate( fp = yall - tp ) %>% dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>%
        dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ) | (fp == 0 ) | (tn == 0) )*0.5 ) %>%
        dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ), 
                       v = 1/( tp + cx ) + 1/( tn + cx ) + 
                           1/( fn + cx ) + 1/( fp + cx ))  %>%
        dplyr::mutate( tmp = y/sqrt(v) ) %>% dplyr::mutate(index = ind) %>%
        dplyr::mutate( addict = log( choose(yall, tp) ) ) %>%
        dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )
})
data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
    x1 = which(data.opt.vec[['index']] == ind)
    c(min(x1), max(x1))
} )

yindex <- with( data.x, cumsum( c( fp + tp + 1 ) ) )
xindex <- with( data.x, tp + 1 + c(0, yindex[1:(nstudy-1)]  )  )

yall <- data.x$yall
tmp <- data.x$tmp

tp.all <- data.opt.vec$tp
yall.all <- data.opt.vec$yall
n1.all <- data.opt.vec$n1
n0.all <- data.opt.vec$n0
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp
addict.all <- data.opt.vec$addictx

# first, for BN
estimate.bn  <- function(p = 0.7){
    
    start.p <- c(-0.71, 0.28, 0.1)
    
    if (p==1){
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(log(n1/n0) + tau*tvec + mu)))
                dbinom( tp, yall, p1)*dnorm(tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral)
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            return(-l1 )
            
        }
        result <- try(nlminb(start.p[1:2], llk.o, 
                             lower=c(-Inf, 0.01), upper = c(Inf,10)), silent=F)
        
    }else{
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            beta <- par[3]
            
            ## ESTIMATE ALPHA
            
            f.a <- function(alpha) {
                pox <- pnorm(alpha + beta*tmp.all)
                prob.x <- cubature::hcubature( f = function(tvec){
                    
                    p1 = 1/(1 + exp(-(log(n1.all/n0.all) + tau*tvec + mu)))
                    vec = dbinom(tp.all, yall.all, p1)*pox
                    
                    sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
                }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral
                
                sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
                
            }
            
            alpha.opt <- uniroot(f.a, lower = -10, upper = 10, extendInt = 'yes')$root
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(log(n1/n0) + tau*tvec + mu)))
                dbinom( tp, yall, p1)*dnorm(tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral)
            
            pox <- pnorm(alpha.opt + beta*tmp.all)
            prob.denom <- cubature::hcubature( f = function(tvec){
                
                p1 = 1/(1 + exp(-(log(n1.all/n0.all) + tau*tvec + mu)))
                vec = dbinom(tp.all, yall.all, p1)*pox
                
                sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            l2 <- sum(pnorm(alpha.opt + beta*tmp, log.p = TRUE), na.rm = TRUE)
            l3 <- sum( log( prob.denom), na.rm = TRUE )
            
            ll <- l1 + l2 - l3
            
            return(-ll )
            
        }
        
        result <- try(nlminb(start.p, llk.o, 
                             lower=c(-Inf, 0.01, -10), upper = c(Inf,10,10)), silent=F)
        
    }
    
    # return( c(result$par, rep(NA, 4 - length(result$par)) ) )
    
    tau = result$par[2]
    if(p == 1){
        hess <- numDeriv::hessian( llk.o, x = result$par )[1:2, 1:2]
    }else{
        hess <- numDeriv::hessian( llk.o, x = result$par )
    }

    fisher <- diag( solve(hess)[1:2, 1:2] )

    lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )

    lower.bound[2] <- exp( log( tau ) - 1.96*sqrt(fisher)[2]/tau )
    upper.bound[2] <- exp( log( tau ) + 1.96*sqrt(fisher)[2]/tau )

    return(list(
        result = result,
        lower = lower.bound,
        mean = result$par[1:2],
        upper = upper.bound,
        beta.opt = ifelse(p==1, NA, result$par[3])
    ))
    
}

result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn, mc.cores = 10L)

result.bn.matrix <- sapply( result.bn, FUN = function(x) c(x$mean, x$beta.opt) ) %>% t() %>% `rownames<-`(NULL)
print(result.bn.matrix)

# second, for HN
estimate.hn  <- function(p = 0.7){
    
    start.p <- c(-0.71, 0.28, 0.1)
    
    if(p==1){
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            
            prob.prior <- cubature::hcubature(
                f = function(tvec){
                    weight.log <- addict.all + (tvec*tau + mu)*tp.all
                    weight <- exp( weight.log - unlist(
                        lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
                    ) )
                    
                    weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )
                    
                    (weight[xindex]/weight.sum.bypart)*dnorm( tvec )
                    
                }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
            )$integral
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            return(-l1 )
            
        }
        result <- try(nlminb(start.p[1:2], llk.o, 
                             lower=c(-Inf, 0.01), upper = c(Inf,10)), silent=F)
    }else{
        llk.o <- function(par){
            
            mu   <- par[1]
            tau <- par[2]
            beta <- par[3]
            
            f.a <- function(alpha) {
                pox = pnorm(alpha + beta*tmp.all)
                prob.x <- cubature::hcubature(
                    f = function(tvec){
                        weight.log <- addict.all + (tvec*tau + mu)*tp.all
                        
                        weight <- exp( weight.log - unlist(
                            lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
                        ) )
                        
                        denom <- weight*pox
                        sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/
                            sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )*dnorm(tvec)
                        
                    }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
                )$integral
                
                sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
                
            }
            
            alpha.opt <- uniroot(f.a, lower = -10, upper = 10, extendInt = 'yes')$root
            
            pox = pnorm(alpha.opt + beta*tmp.all)
            prob.x.mean <- cubature::hcubature(
                f = function(tvec){
                    weight.log <- addict.all + (tvec*tau + mu)*tp.all
                    weight <- exp( weight.log - unlist(
                        lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
                    ) )
                    
                    weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )
                    denom <- weight*pox
                    
                    c(weight[xindex]/weight.sum.bypart, 
                      sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/weight.sum.bypart)*dnorm( tvec )
                    
                }, lowerLimit = -5, upperLimit = 5, fDim = 2*nstudy, tol = 1e-4
            )$integral
            
            l1 <- sum( log( prob.x.mean[1:nstudy]), na.rm = TRUE )
            l2 <- sum(pnorm(alpha.opt + beta * tmp , log.p = TRUE), na.rm = TRUE)
            l3 <- sum( log( prob.x.mean[(nstudy+1):(2*nstudy)]), na.rm = TRUE )
            
            ll <- l1 + l2 - l3
            
            return(-ll )
            
        }
        
        result <- try(nlminb(start.p, llk.o, 
                             lower=c(-Inf, 0.01, -10), upper = c(Inf,10,10)), silent=F)
    }
    
    # return( c(result$par, rep(NA, 4 - length(result$par)) ) )
    
    tau = result$par[2]
    if(p == 1){
        hess <- numDeriv::hessian( llk.o, x = result$par )[1:2, 1:2]
    }else{
        hess <- numDeriv::hessian( llk.o, x = result$par )
    }

    fisher <- diag( solve(hess)[1:2, 1:2] )

    lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )

    lower.bound[2] <- exp( log( tau ) - 1.96*sqrt(fisher)[2]/tau )
    upper.bound[2] <- exp( log( tau ) + 1.96*sqrt(fisher)[2]/tau )

    return(list(
        result = result,
        lower = lower.bound,
        mean = result$par[1:2],
        upper = upper.bound,
        beta.opt = ifelse(p==1, NA, result$par[3])
    ))
    
}

result.hn <- mclapply((1:10)*0.1, FUN = estimate.hn, mc.cores = 10L)
# as.data.frame( result.hn ) %>% t() %>% `rownames<-`(NULL) %>% print()

result.hn.matrix <- sapply( result.hn, FUN = function(x) c(x$mean, x$beta.opt) ) %>% t() %>% `rownames<-`(NULL)
print(result.hn.matrix)


# second, for nn
estimate.nn  <- function(p = 0.7){
    
    start.p <- c(-0.71, 0.28, 0.1)
    eps = .Machine$double.eps^0.5
    
    if(p==1){
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            
            prob.prior <- with(data.x, dnorm( y, mean = mu, sd = sqrt(tau^2 + v) ) )
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            return(-l1 )
            
        }
        result <- try(nlminb(start.p[1:2], llk.o, 
                             lower=c(-Inf, eps), upper = c(Inf,10)), silent=F)
    }else{
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            beta <- par[3]
            
            f.a <- function(alpha) {
                
                prob.denom = with(data.x, pnorm( (alpha + beta*mu/sqrt(v))/sqrt( 1 + beta^2*(1 + tau^2/v) ) ))
                sum( 1/prob.denom, na.rm = TRUE ) - nstudy/p
             }
            
            alpha.opt <- uniroot(f.a, lower = -3, upper = 3, extendInt = 'downX')$root
            
            l1 <- with(data.x, dnorm( y, mean = mu, sd = sqrt(tau^2 + v) ) )
            l2 <- pnorm(alpha.opt + beta * tmp)
            l3 <- with(data.x, pnorm( (alpha.opt + beta*mu/sqrt(v))/sqrt( 1 + beta^2*(1 + tau^2/v) ) ))
            
            ll <- sum( log(l1), na.rm = TRUE ) + sum( log(l2), na.rm = TRUE ) - sum( log(l3), na.rm = TRUE )
            
            return(-ll )
            
        }
        
        result <- try(nlminb(start.p, llk.o, 
                             lower=c(-Inf, eps, -10), upper = c(Inf, 10, 10)), silent=F)
    }
    
    # return( c(result$par, rep(NA, 4 - length(result$par)) ) )
    
    tau = result$par[2]
    if(p == 1){
        hess <- numDeriv::hessian( llk.o, x = result$par )[1:2, 1:2]
    }else{
        hess <- numDeriv::hessian( llk.o, x = result$par )
    }

    fisher <- diag( solve(hess)[1:2, 1:2] )

    lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )

    lower.bound[2] <- exp( log( tau ) - 1.96*sqrt(fisher)[2]/tau )
    upper.bound[2] <- exp( log( tau ) + 1.96*sqrt(fisher)[2]/tau )

    return(list(
        result = result,
        lower = lower.bound,
        mean = result$par[1:2],
        upper = upper.bound,
        beta.opt = ifelse(p==1, NA, result$par[3])
    ))
    
}

result.nn <- mclapply((1:10)*0.1, FUN = estimate.nn, mc.cores = 10L)
result.nn.matrix <- sapply( result.nn, FUN = function(x) c(x$mean, x$beta.opt) ) %>% t() %>% `rownames<-`(NULL)
print(result.nn.matrix)

# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.mean <- round( sapply( result.hn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.lower <- round( sapply( result.hn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.upper <- round( sapply( result.hn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
nn.mean <- round( sapply( result.nn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
nn.lower <- round( sapply( result.nn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
nn.upper <- round( sapply( result.nn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

mu.table <- data.frame(
    # p = (10:1)*0.1,
    HN = paste0(hn.mean[, 1], '(', hn.lower[, 1], ', ', hn.upper[, 1], ')'), 
    BN = paste0(bn.mean[, 1], '(', bn.lower[, 1], ', ', bn.upper[, 1], ')'), 
    NN = paste0(nn.mean[, 1], '(', nn.lower[, 1], ', ', nn.upper[, 1], ')')
) %>% `rownames<-`((10:1)*0.1)

tau.table <- data.frame(
    # p = (10:1)*0.1,
    HN = paste0(hn.mean[, 2], '(', hn.lower[, 2], ', ', hn.upper[, 2], ')'), 
    BN = paste0(bn.mean[, 2], '(', bn.lower[, 2], ', ', bn.upper[, 2], ')'), 
    NN = c(rep('<10^(-8)', 5), '1.35times10^(-5)', rep('<10^(-8)', 4) )
) %>% `rownames<-`((10:1)*0.1)

beta.table <- data.frame(
    # p = (10:1)*0.1,
    HN = round( sapply( result.hn, FUN = function(x) x$beta.opt ), 3 )[10:1], 
    BN = round( sapply( result.bn, FUN = function(x) x$beta.opt ), 3 )[10:1], 
    NN = round( sapply( result.nn, FUN = function(x) x$beta.opt ), 3 )[10:1]
) %>% `rownames<-`((10:1)*0.1)

xtable::xtable(mu.table, digits = 3) %>% print()
xtable::xtable(tau.table, digits = 3) %>% print()
xtable::xtable(beta.table, digits = 3) %>% print()

