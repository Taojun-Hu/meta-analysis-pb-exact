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
    y1 = c(07, 08, 03, 06, 01, 01, 17, 03, 02), 
    y0 = c(11, 08, 14, 10, 01, 08, 15, 07, 03), 
    t1  = c(1491, 1638, 11484, 1446, 0370, 0786, 6840, 1156, 0400), 
    t0  = c(1988, 1460, 10962, 1502, 0482, 0913, 6840, 1015, 0440)
) %>% dplyr::mutate( yall = y1+y0 ) %>%
    dplyr::mutate( cx = (y0==0 | y1==0)*0.5 ) %>%
    dplyr::mutate( y = log( ( (y1+cx)/t1 )/( (y0+cx)/t0 ) ), 
                   v = 1/(y1+cx) + 1/(y0+cx) ) %>%
    dplyr::mutate( tmp = y/sqrt(v) )

nstudy = nrow(data.x)

data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( y1 = 0:data.x$yall[ind], yall = data.x$yall[ind], t1 = data.x$t1[ind], t0 = data.x$t0[ind]  ) %>%
        dplyr::mutate( y0 = yall - y1 ) %>%
        dplyr::mutate( cx = (y0==0 | y1==0)*0.5 ) %>%
        dplyr::mutate( y = log( ( (y1+cx)/t1 )/( (y0+cx)/t0 ) ), 
                       v = 1/(y1+cx) + 1/(y0+cx) ) %>%
        dplyr::mutate( tmp = y/sqrt(v) ) %>% 
        dplyr::mutate(index = ind) 
})
data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
    x1 = which(data.opt.vec[['index']] == ind)
    c(min(x1), max(x1))
} )

yall <- data.x$yall
tmp <- data.x$tmp

y1.all <- data.opt.vec$y1
yall.all <- data.opt.vec$yall
t1.all <- data.opt.vec$t1
t0.all <- data.opt.vec$t0
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp

estimate.bnicr  <- function(p = 0.7){
    
    start.p <- c(-0.71, 0.28, 0.1)
    # start.p <- c(est.o.vec, 0.1)
    
    if (p==1){
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(log(t1/t0) + tau*tvec + mu)))
                dbinom( y1, yall, p1)*dnorm(tvec)
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
                    
                    p1 = 1/(1 + exp(-(log(t1.all/t0.all) + tau*tvec + mu)))
                    vec = dbinom(y1.all, yall.all, p1)*pox
                    
                    sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
                }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral
                
                sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
                
            }
            
            alpha.opt <- uniroot(f.a, lower = -10, upper = 10, extendInt = 'yes')$root
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(log(t1/t0) + tau*tvec + mu)))
                dbinom( y1, yall, p1)*dnorm(tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral)
            
            pox <- pnorm(alpha.opt + beta*tmp.all)
            prob.denom <- cubature::hcubature( f = function(tvec){
                
                p1 = 1/(1 + exp(-(log(t1.all/t0.all) + tau*tvec + mu)))
                vec = dbinom(y1.all, yall.all, p1)*pox
                
                sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )$integral
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            l2 <- sum( log( pnorm(alpha.opt + beta*tmp) ) , na.rm = TRUE)
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

result.bnicr <- mclapply((1:10)*0.1, FUN = estimate.bnicr, mc.cores = 10L)

result.bnicr.matrix <- sapply( result.bnicr, FUN = function(x) c(x$mean, x$beta.opt) ) %>% t() %>% `rownames<-`(NULL)
print(result.bnicr.matrix)

bnicr.mean <- round( sapply( result.bnicr, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bnicr.lower <- round( sapply( result.bnicr, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bnicr.upper <- round( sapply( result.bnicr, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

final.table <- data.frame(
    # p = (10:1)*0.1,
    theta = paste0(bnicr.mean[, 1], '(', bnicr.lower[, 1], ', ', bnicr.upper[, 1], ')'), 
    tau   = paste0(bnicr.mean[, 2], '(', bnicr.lower[, 2], ', ', bnicr.upper[, 2], ')'), 
    beta  = round( sapply( result.bnicr, FUN = function(x) x$beta.opt ), 3 )[10:1]
) %>% `rownames<-`((10:1)*0.1)

xtable::xtable(final.table) %>% print()

# second, for nn
# estimate.nn  <- function(p = 0.7){
# 
#     start.p <- c(-0.71, 0.28, 0.1)
#     # start.p <- c(est.o.vec, 0.1)
#     eps = .Machine$double.eps^0.5
# 
#     if(p==1){
#         llk.o <- function(par){
# 
#             mu   <- par[1]
#             tau  <- par[2]
# 
#             prob.prior <- with(data.x, dnorm( y, mean = mu, sd = sqrt(tau^2 + v) ) )
# 
#             l1 <- sum( log( prob.prior), na.rm = TRUE )
#             return(-l1 )
# 
#         }
#         result <- try(nlminb(start.p[1:2], llk.o,
#                              lower=c(-Inf, eps), upper = c(Inf,10)), silent=F)
#     }else{
#         llk.o <- function(par){
# 
#             mu   <- par[1]
#             tau  <- par[2]
#             beta <- par[3]
# 
#             f.a <- function(alpha) {
# 
#                 prob.denom = with(data.x, pnorm( (alpha + beta*mu/sqrt(v))/sqrt( 1 + beta^2*(1 + tau^2/v) ) ))
#                 sum( 1/prob.denom, na.rm = TRUE ) - nstudy/p
#             }
# 
#             alpha.opt <- uniroot(f.a, lower = -10, upper = 10, extendInt = 'yes')$root
# 
#             l1 <- with(data.x, dnorm( y, mean = mu, sd = sqrt(tau^2 + v) ) )
#             l2 <- pnorm(alpha.opt + beta * tmp)
#             l3 <- with(data.x, pnorm( (alpha.opt + beta*mu/sqrt(v))/sqrt( 1 + beta^2*(1 + tau^2/v) ) ))
# 
#             ll <- sum( log(l1), na.rm = TRUE ) + sum( log(l2), na.rm = TRUE ) - sum( log(l3), na.rm = TRUE )
# 
#             return(-ll )
# 
#         }
# 
#         result <- try(nlminb(start.p, llk.o,
#                              lower=c(-Inf, eps, -10), upper = c(Inf, 10, 10)), silent=F)
#     }
# 
#     # return( c(result$par, rep(NA, 4 - length(result$par)) ) )
# 
#     tau = result$par[2]
#     if(p == 1){
#         hess <- numDeriv::hessian( llk.o, x = result$par )[1:2, 1:2]
#     }else{
#         hess <- numDeriv::hessian( llk.o, x = result$par )
#     }
# 
#     fisher <- diag( solve(hess)[1:2, 1:2] )
# 
#     lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
#     upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )
# 
#     lower.bound[2] <- exp( log( tau ) - 1.96*sqrt(fisher)[2]/tau )
#     upper.bound[2] <- exp( log( tau ) + 1.96*sqrt(fisher)[2]/tau )
# 
#     return(list(
#         result = result,
#         lower = lower.bound,
#         mean = result$par[1:2],
#         upper = upper.bound,
#         beta.opt = ifelse(p==1, NA, result$par[3])
#     ))
# 
# }
# 
# result.nn <- mclapply((1:10)*0.1, FUN = estimate.nn, mc.cores = 10L)
# result.nn.matrix <- sapply( result.nn, FUN = function(x) c(x$result$par, rep(NA, 4 - length(x$result$par)) ) ) %>% t() %>% `rownames<-`(NULL)
# print(result.nn.matrix)
