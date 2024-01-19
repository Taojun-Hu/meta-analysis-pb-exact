library( dplyr )
library( parallel )

args <- commandArgs(TRUE)
if ( length(args) == 0 ){
    s = 25
}else{
    arg_data <- do.call('rbind', strsplit(sub("^--", "", args), "=")) %>% as.data.frame()
    argsl <- as.list(as.character(arg_data$V2)) %>% `names<-`(arg_data$V1)
    
    if ( is.null(argsl$s)){
        s = 25
    }else{
        s = as.numeric( argsl$s )
    }
}

n.datasets <- 1000
t.theta <- -3
t.sigma <- 0.3
beta <- 1
set.seed(1234)

# global setting up!

data.generate <- lapply(
    1:n.datasets, FUN = function(ind){
        study.theta <- rnorm( s, mean = t.theta, sd = t.sigma)
        study.pi <- plogis(study.theta)
        n_all <- runif( s, min = 200, max = 400 ) %>% round()
        y_all <- rbinom(s, size = n_all, prob = study.pi)
        y.revise <- ifelse( ( y_all == 0 )| (y_all == n_all), y_all + 0.5, y_all )
        n.revise <- ifelse( ( y_all == 0 )| (y_all == n_all), n_all + 1, n_all )
        t.par <- log( y.revise/(n.revise - y.revise) )/sqrt( ( 1/y.revise + 1/(n.revise - y.revise) ) )

        expectation <- sqrt(n_all)*t.theta/(exp(-t.theta/2) + exp(t.theta/2)) -1/(8*sqrt(n_all))*(
            t.theta*(exp(t.theta/2) + exp(-t.theta/2))
        ) + 1/2*t.sigma**2*(
            sqrt(n_all)*(
                exp(-3*t.theta/2)*(1 + t.theta/4) - exp(3*t.theta/2)*(1 - t.theta/4) + 
                    exp(-t.theta/2)*(1 - 5*t.theta/4) - exp(t.theta/2)*(1 + 5*t.theta/4)
            )/(
                exp(-t.theta/2) + exp(t.theta/2)
            )**4 - 1/(8*sqrt(n_all))*(
                exp(t.theta/2)*(1+t.theta/4) + exp(-t.theta/2)*(t.theta/4 - 1)
            )
        ) + 1/24*3*t.sigma**4*(
            1/16*exp(-t.theta/2)*(
                sqrt(n_all)*exp(t.theta)*(
                    exp(3*t.theta)*(176 - 76*t.theta) + exp(4*t.theta)*(t.theta - 8) + 
                        230*exp(2*t.theta)*t.theta + t.theta - 4*exp(t.theta)*(19*t.theta + 44) + 8
                )/(1 + exp(t.theta))**5 - 1/(8*sqrt(n_all))*(t.theta + exp(t.theta)*(t.theta + 8) - 8)
            )
        )
        
        variance <- (
            sqrt(n_all)*(exp(t.theta/2)*(1 - t.theta/2) + exp(-t.theta/2)*(1+t.theta/2))/(
                exp(t.theta/2) + exp(-t.theta/2)
            )**2 -1/(8*sqrt(n_all))*(exp(t.theta/2)*(1 + t.theta/2) + exp(-t.theta/2)*(1-t.theta/2))
        )**2*t.sigma**2 +(
            sqrt(n_all)*(
                exp(-3*t.theta/2)*(1 + t.theta/4) - exp(3*t.theta/2)*(1 - t.theta/4) + 
                    exp(-t.theta/2)*(1 - 5*t.theta/4) - exp(t.theta/2)*(1 + 5*t.theta/4)
            )/(
                exp(-t.theta/2) + exp(t.theta/2)
            )**4 - 1/(8*sqrt(n_all))*(
                exp(t.theta/2)*(1+t.theta/4) + exp(-t.theta/2)*(t.theta/4 - 1)
            )
        )**2*1/4*3*t.sigma**4 + (
            1 + 1/2*t.theta*(1 - exp(t.theta))/(1 + exp(t.theta)) + 1/2*t.sigma**2*(
                exp(2*t.theta)*(t.theta-2) - exp(t.theta)*(t.theta + 2)
            )/(1 + exp(t.theta))**3
        )**2 + 1/4*(1 - exp(2*t.theta) - 2*t.theta*exp(t.theta))**2/(1 + exp(t.theta))**4*t.sigma**2
        
        pselect <- function(alpha){
            
            select.tmp = pnorm( (beta*expectation + alpha)/sqrt(1 + beta**2*variance) )
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/0.7 )
        }
        
        alpha.set <- uniroot(pselect, lower = -10, upper = 10, extendInt = 'yes')$root
        
        t.prop <- pnorm( alpha.set + beta*t.par )
        t.sample <- sapply( t.prop, FUN = function(x) {rbinom(1, 1, x)} ) %>% as.logical()
        
        data.org = data.frame( y = y_all, n = n_all )
        data.select = data.org[t.sample, ] %>% `rownames<-`(NULL)  
        
        return( list( data.org = data.org, data.select = data.select, alpha = alpha.set) )
    }
)


alpha.all <- sapply( data.generate, FUN = function(x) x[['alpha']], simplify = TRUE )
data.org.all <- lapply( data.generate, FUN = function(x) x[['data.org']] )
data.select.all <- lapply( data.generate, FUN = function(x) x[['data.select']] )

alpha.kt <- mean( alpha.all, na.rm = TRUE )
save( list =  c('alpha.all', 'data.org.all', 'data.select.all'), file =  paste0('data_', s, '.rda')  )

# printout 1
print('alpha')
print(mean( alpha.kt ))
print('Average events proportion (without select)')
op.org.all <- sapply( data.org.all, FUN = function(x){ c(sum(x[['y']]), sum(x[['n']]))  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.org.all[1]/op.org.all[2])
print('Average events proportion (after select)')
op.select.all <- sapply( data.select.all, FUN = function(x){ c(sum(x[['y']]), sum(x[['n']]))  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.select.all[1]/op.select.all[2])

write.csv( (
    list(alpha = alpha.kt, oporg = op.org.all[1]/op.org.all[2], 
         opselect = op.select.all[1]/op.select.all[2] ) %>% as.data.frame()
), file = paste0('summary_1_', s, '.csv' ) )

rm(list = c('n.datasets', 'beta'))

estimate.fun.all <- function(data = data.org, start.ini = c(-3.6, 0.28)){
    nstudy = nrow( data )
    y = data[['y']]
    n = data[['n']]
    
    start.p <- start.ini
    llk <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dbinom( y, n, prob = plogis( sigma.ini*theta + theta.ini ) )*dnorm(theta)
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )[['integral']]
        
        -sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
    }
    
    result <- try(nlminb(start.p, llk, 
                         lower=c(-Inf, 0.01), upper = c(Inf, 10)), silent=T)
    
    # hess <- try( numDeriv::hessian( llk, x = result$par ), silent = TRUE )
    # 
    # if (class(hess) != 'try-error'){
    #     
    #     theta.opt <- result$par[1]
    #     sigma.opt <- result$par[2]
    #     
    #     fisher <- diag( solve(hess)[1:2, 1:2] )
    #     
    #     lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    #     upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )
    #     
    #     lower.bound[2] <- exp( log( sigma.opt ) - 1.96*sqrt(fisher)[2]/sigma.opt )
    #     upper.bound[2] <- exp( log( sigma.opt ) + 1.96*sqrt(fisher)[2]/sigma.opt )
    #     
    #     coverage <- ( (lower.bound <= c(t.theta, t.sigma)) | (upper.bound >= c(t.theta, t.sigma)) )
    #     
    # }else{
    #     lower.bound <- upper.bound <- fisher <- c(NA, NA)
    #     coverage <- c(NA, NA)
    # }
    
    lower.bound <- upper.bound <- fisher <- c(NA, NA)
    coverage <- c(NA, NA)
    
    return(list( result = result, lower.bound = lower.bound, upper.bound = upper.bound,
                 fisher = fisher, coverage = coverage ))
    
}

# data.select <- data.select.all[[1]]
estimate.fun.nn <- function(data = data.select, start.ini = c(1.0, 0.8, 0.2), method = 'nlm',  p = 0.7){
    nstudy = nrow(data)
    y = data[['y']]
    n = data[['n']]
    
    g.tmp <- (data.frame( y = y, n = n ) %>%
                  dplyr::mutate( ytmp = ifelse( (y==0) | (y==n), y+0.5, y ),
                                 ntmp = ifelse( (y==0) | (y==n), n+1, n )) %>%
                  dplyr::mutate( ttmp = log( ytmp/(ntmp - ytmp) )/sqrt(
                      ( 1/ytmp + 1/(ntmp - ytmp) )
                  ) ))
    
    yp <- with(g.tmp, log(ytmp/(ntmp - ytmp)))
    v = with( g.tmp, 1/ytmp + 1/(ntmp - ytmp) )
    
    start.p = start.ini
    
    t.tmp = g.tmp[['ttmp']]
    
    llk <- function( par ){
        theta.ini = par[1]
        sigma.ini = par[2]
        beta.ini = par[3]
        
        pselect <- function(alpha, print = FALSE){
            
            select.tmp <- pnorm(
                ( beta.ini*theta.ini/sqrt(v) + alpha )/sqrt(
                    1 + beta.ini**2*(1 + sigma.ini**2/v)
                )
            )
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/p  )
        }
        
        alpha.opt = uniroot(pselect, lower = -10, upper = 10, extendInt = 'yes')$root
        
        prob.prior <- dnorm( yp, theta.ini, sd = sqrt( sigma.ini**2 + v )  )
        
        prob.denom <- pnorm(
            ( beta.ini*theta.ini/sqrt(v) + alpha.opt )/sqrt(
                1 + beta.ini**2*(1 + sigma.ini**2/v)
            )
        )
        
        l1 = sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
        l2 = sum( pnorm( alpha.opt + beta.ini*t.tmp, log.p = TRUE  ), na.rm = TRUE )
        l3 = sum( log( prob.denom ), na.rm = TRUE )
        
        return( -(l1 + l2 - l3) )
    }
    
    if( method == 'optim'){
        result <- try(
            optim( par = start.p, fn = llk, method = 'L-BFGS-B', lower = c(-Inf, 0.01, -Inf), upper = c(Inf, 10, Inf)),
            silent = TRUE
        )
    }else{
        
        result <- try(nlminb(start.p, llk,
                             lower=c(-Inf, 0.01, -10), upper = c(Inf, 10, 10)), silent=F)
    }
    
    # if(p == 1){
    #     hess <- try( numDeriv::hessian( llk, x = result$par )[1:2, 1:2], silent = TRUE )
    # }else{
    #     hess <- try( numDeriv::hessian( llk, x = result$par ), silent = TRUE )
    # }
    # 
    # if (class(hess) != 'try-error'){
    #     
    #     theta.opt <- result$par[1]
    #     sigma.opt <- result$par[2]
    #     
    #     fisher <- diag( solve(hess)[1:2, 1:2] )
    #     
    #     lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    #     upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )
    #     
    #     lower.bound[2] <- exp( log( sigma.opt ) - 1.96*sqrt(fisher)[2]/sigma.opt )
    #     upper.bound[2] <- exp( log( sigma.opt ) + 1.96*sqrt(fisher)[2]/sigma.opt )
    #     
    #     coverage <- ( (lower.bound <= c(t.theta, t.sigma)) | (upper.bound >= c(t.theta, t.sigma)) )
    #     
    # }else{
    #     lower.bound <- upper.bound <- fisher <- c(NA, NA)
    #     coverage <- c(NA, NA)
    # }
    
    lower.bound <- upper.bound <- fisher <- c(NA, NA)
    coverage <- c(NA, NA)
    
    return(list( result = result, lower.bound = lower.bound, upper.bound = upper.bound,
                 fisher = fisher, coverage = coverage ))
    
}

estimate.fun.asymp <- function(data = data.select, start.ini = c(-3.6, 0.28, 0.4), method = 'optim',  p = 0.7){
    nstudy = nrow(data)
    y = data[['y']]
    n = data[['n']]
    
    llk.com <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dbinom( y, n, prob = plogis( sigma.ini*theta + theta.ini ) )*dnorm(theta)
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )[['integral']]
        
        -sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
    }
    
    start.p <- try(nlminb(start.ini[1:2], llk.com, lower=c(-Inf, 0.01), upper = c(Inf, 10)), silent=T)[['par']]
    if(start.p[2] %in% c(0.01, 10) ){
        start.p[2] = start.ini[2]
    }
    start.p <- c(start.p, start.ini[3])
    
    g.tmp <- (data.frame( y = y, n = n ) %>%
                  dplyr::mutate( ytmp = ifelse( (y==0) | (y==n), y+0.5, y ), 
                                 ntmp = ifelse( (y==0) | (y==n), n+1, n )) %>%
                  dplyr::mutate( ttmp = log( ytmp/(ntmp - ytmp) )/sqrt(
                      ( 1/ytmp + 1/(ntmp - ytmp) )
                  ) ))
    
    # start.p = start.ini
    
    t.tmp = g.tmp[['ttmp']]
    
    llk <- function( par ){
        theta.ini = par[1]
        sigma.ini = par[2]
        beta.ini = par[3]
        
        expectation <- sqrt(n)*theta.ini/(exp(-theta.ini/2) + exp(theta.ini/2)) -1/(8*sqrt(n))*(
            theta.ini*(exp(theta.ini/2) + exp(-theta.ini/2))
        ) + 1/2*sigma.ini**2*(
            sqrt(n)*(
                exp(-3*theta.ini/2)*(1 + theta.ini/4) - exp(3*theta.ini/2)*(1 - theta.ini/4) + 
                    exp(-theta.ini/2)*(1 - 5*theta.ini/4) - exp(theta.ini/2)*(1 + 5*theta.ini/4)
            )/(
                exp(-theta.ini/2) + exp(theta.ini/2)
            )**4 - 1/(8*sqrt(n))*(
                exp(theta.ini/2)*(1+theta.ini/4) + exp(-theta.ini/2)*(theta.ini/4 - 1)
            )
        ) + 1/24*3*sigma.ini**4*(
            1/16*exp(-theta.ini/2)*(
                sqrt(n)*exp(theta.ini)*(
                    exp(3*theta.ini)*(176 - 76*theta.ini) + exp(4*theta.ini)*(theta.ini - 8) + 
                        230*exp(2*theta.ini)*theta.ini + theta.ini - 4*exp(theta.ini)*(19*theta.ini + 44) + 8
                )/(1 + exp(theta.ini))**5 - 1/(8*sqrt(n))*(theta.ini + exp(theta.ini)*(theta.ini + 8) - 8)
            )
        )
        
        variance <- (
            sqrt(n)*(exp(theta.ini/2)*(1 - theta.ini/2) + exp(-theta.ini/2)*(1+theta.ini/2))/(
                exp(theta.ini/2) + exp(-theta.ini/2)
            )**2 -1/(8*sqrt(n))*(exp(theta.ini/2)*(1 + theta.ini/2) + exp(-theta.ini/2)*(1-theta.ini/2))
        )**2*sigma.ini**2 +(
            sqrt(n)*(
                exp(-3*theta.ini/2)*(1 + theta.ini/4) - exp(3*theta.ini/2)*(1 - theta.ini/4) + 
                    exp(-theta.ini/2)*(1 - 5*theta.ini/4) - exp(theta.ini/2)*(1 + 5*theta.ini/4)
            )/(
                exp(-theta.ini/2) + exp(theta.ini/2)
            )**4 - 1/(8*sqrt(n))*(
                exp(theta.ini/2)*(1+theta.ini/4) + exp(-theta.ini/2)*(theta.ini/4 - 1)
            )
        )**2*1/4*3*sigma.ini**4 + (
            1 + 1/2*theta.ini*(1 - exp(theta.ini))/(1 + exp(theta.ini)) + 1/2*sigma.ini**2*(
                exp(2*theta.ini)*(theta.ini-2) - exp(theta.ini)*(theta.ini + 2)
            )/(1 + exp(theta.ini))**3
        )**2 + 1/4*(1 - exp(2*theta.ini) - 2*theta.ini*exp(theta.ini))**2/(1 + exp(theta.ini))**4*sigma.ini**2
        
        pselect <- function(alpha, print = FALSE){
            
            select.tmp <- pnorm( (alpha + beta.ini*expectation)/sqrt(1+beta.ini**2*variance ) )
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/p  )
        }
        
        alpha.opt = uniroot(pselect, lower = -10, upper = 10, extendInt = 'yes')$root
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dbinom( y, n, prob = plogis( sigma.ini*theta + theta.ini ) )*dnorm(theta)
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )[['integral']]
        
        prob.denom <- pnorm( (alpha.opt + beta.ini*expectation)/sqrt(1+beta.ini**2*variance ) )
        
        l1 = sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
        l2 = sum( pnorm( alpha.opt + beta.ini*t.tmp, log.p = TRUE  ), na.rm = TRUE )
        l3 = sum( log( prob.denom + 1e-30 ), na.rm = TRUE )
        
        return( -(l1 + l2 - l3) )
    }
    
    if( method == 'optim'){
        result <- try(
            optim( par = start.p, fn = llk, method = 'L-BFGS-B', lower = c(-Inf, 0.01, -Inf), upper = c(Inf, 10, Inf)), 
            silent = TRUE
        )
    }else{
        result <- try(nlminb(start.p, llk, 
                             lower=c(-Inf, 0.01, -10), upper = c(Inf, 10, 10)), silent=F)
    }
    
    # if(p == 1){
    #     hess <- try( numDeriv::hessian( llk, x = result$par )[1:2, 1:2], silent = TRUE )
    # }else{
    #     hess <- try( numDeriv::hessian( llk, x = result$par ), silent = TRUE )
    # }
    # 
    # if (class(hess) != 'try-error'){
    #     
    #     theta.opt <- result$par[1]
    #     sigma.opt <- result$par[2]
    #     
    #     fisher <- diag( solve(hess)[1:2, 1:2] )
    #     
    #     lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    #     upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )
    #     
    #     lower.bound[2] <- exp( log( sigma.opt ) - 1.96*sqrt(fisher)[2]/sigma.opt )
    #     upper.bound[2] <- exp( log( sigma.opt ) + 1.96*sqrt(fisher)[2]/sigma.opt )
    #     
    #     coverage <- ( (lower.bound <= c(t.theta, t.sigma)) | (upper.bound >= c(t.theta, t.sigma)) )
    #     
    # }else{
    #     lower.bound <- upper.bound <- fisher <- c(NA, NA)
    #     coverage <- c(NA, NA)
    # }
    
    lower.bound <- upper.bound <- fisher <- c(NA, NA)
    coverage <- c(NA, NA)
    
    return(list( result = result, lower.bound = lower.bound, upper.bound = upper.bound,
                 fisher = fisher, coverage = coverage ))
    
}

result.asymp <- mclapply( data.select.all, FUN = estimate.fun.asymp, start.ini = c(-3.6, 0.28, 0.4), method = 'nlm', p = 0.7, 
                          mc.cores = 100L)
result.com <- mclapply( data.org.all, FUN = estimate.fun.all, start.ini = c(-3.6, 0.28), mc.cores = 100L )
result.com.noadjust <-  mclapply( data.select.all, FUN = estimate.fun.all, start.ini = c(-3.6, 0.28), mc.cores = 100L )
result.nn <- mclapply( data.select.all, FUN = estimate.fun.nn, start.ini = c(-3.6, 0.28, 0.4), method = 'nlm', p = 0.7,
                          mc.cores = 100L)


tidy.asymp <- sapply(
    result.asymp, FUN = function(x){
        if ( class(x[['result']]) == 'try-error'){
            return(c(NA, NA))
        }else{
            if(x[['result']]$convergence == 0){ x[['result']]$par[1:2]}else{c(NA, NA)}
        }
    }, simplify = TRUE
)

# keepindex <- ((tidy.asymp[2, ] != 0.01) & (tidy.asymp[2, ] != 10) )
# tidy.asymp <- tidy.asymp[, keepindex]

rowMeans(tidy.asymp, na.rm = TRUE) %>% print()
print.asymp <- rbind( rowMeans(tidy.asymp, na.rm = TRUE),
                      apply(tidy.asymp-c(t.theta, t.sigma), 1, FUN = function(x) mean( (x - mean(x, na.rm = TRUE))**2, na.rm = TRUE ) ),
                      sapply( result.asymp, FUN = function(x) x[['coverage']] ) %>% rowMeans(na.rm = TRUE) )
# 
tidy.com <- sapply(
    result.com, FUN = function(x){
        if ( class(x[['result']]) == 'try-error'){
            return(c(NA, NA))
        }else{
            if(x[['result']]$convergence == 0){ x[['result']]$par[1:2]}else{c(NA, NA)}
        }
    }, simplify = TRUE
)
# keepindex <- ((tidy.com[2, ] != 0.01) & (tidy.com[2, ] != 10) )
# tidy.com <- tidy.com[, keepindex]

rowMeans(tidy.com, na.rm = TRUE) %>% print()
print.com <- rbind( rowMeans(tidy.com, na.rm = TRUE),
                      apply(tidy.com-c(t.theta, t.sigma), 1, FUN = function(x) mean( (x - mean(x, na.rm = TRUE))**2, na.rm = TRUE ) ),
                      sapply( result.com, FUN = function(x) x[['coverage']] ) %>% rowMeans(na.rm = TRUE) )

tidy.com.noadjust <- sapply(
    result.com.noadjust, FUN = function(x){
        if ( class(x[['result']]) == 'try-error'){
            return(c(NA, NA))
        }else{
            if(x[['result']]$convergence == 0){ x[['result']]$par[1:2]}else{c(NA, NA)}
        }
    }, simplify = TRUE
)
# keepindex <- ((tidy.com.noadjust[2, ] != 0.01) & (tidy.com.noadjust[2, ] != 10) )
# tidy.com.noadjust <- tidy.com.noadjust[, keepindex]

rowMeans(tidy.com.noadjust, na.rm = TRUE) %>% print()
print.com.noadjust <- rbind( rowMeans(tidy.com.noadjust, na.rm = TRUE),
                      apply(tidy.com.noadjust-c(t.theta, t.sigma), 1, FUN = function(x) mean( (x - mean(x, na.rm = TRUE))**2, na.rm = TRUE ) ),
                      sapply( result.com.noadjust, FUN = function(x) x[['coverage']] ) %>% rowMeans(na.rm = TRUE) )

tidy.nn <- sapply(
    result.nn, FUN = function(x){
        if ( class(x[['result']]) == 'try-error'){
            return(c(NA, NA))
        }else{
            # if(x[['result']]$convergence == 0){ x[['result']]$par[1:2]}else{c(NA, NA)}
            x[['result']]$par[1:2]
        }
    }, simplify = TRUE
)
# keepindex <- ((tidy.nn[2, ] != 0.01) & (tidy.nn[2, ] != 10) )
# tidy.nn <- tidy.nn[, keepindex]

rowMeans(tidy.nn, na.rm = TRUE) %>% print()
print.nn <- rbind( rowMeans(tidy.nn, na.rm = TRUE),
                      apply(tidy.nn-c(t.theta, t.sigma), 1, FUN = function(x) mean( (x - mean(x, na.rm = TRUE))**2, na.rm = TRUE ) ),
                      sapply( result.nn, FUN = function(x) x[['coverage']] ) %>% rowMeans(na.rm = TRUE) )

save(list =  c(paste0('result.', c('asymp', 'com', 'com.noadjust', 'nn') ), 
               paste0('tidy.', c('asymp', 'com', 'com.noadjust', 'nn') )), 
     file = paste0('allsave_', s, '.rda' ))

rio::export(list(asymp = (print.asymp %>% round(4) %>% as.data.frame()), 
                 com = (print.com %>% round(4) %>% as.data.frame()), 
                 com.noadjust = (print.com.noadjust %>% round(4) %>% as.data.frame()), 
                 nn = (print.nn %>% round(4) %>% as.data.frame())  ), 
            file = paste0('printout_', s, '.xlsx' ))
