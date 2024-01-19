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

data.generate <- lapply(
    1:n.datasets, FUN = function(ind){
        study.theta <- rnorm( s, mean = t.theta, sd = t.sigma)
        t <- runif( s, min = 250, max = 350 )
        y.give <- rpois( s, lambda = t*exp( study.theta ) )
        
        y.revise <- ifelse( ( y.give == 0 ), y.give + 0.5, y.give )
        t.par <- log( y.revise/t )/sqrt( ( 1/y.revise ) )
        
        expectation = sqrt(t)*exp(t.theta/2)*t.theta - 1/(8*sqrt(t))*t.theta*exp(-t.theta/2) + 1/2*(
            sqrt(t)*exp(t.theta/2)*(1 + t.theta/4) - 1/(8*sqrt(t))*exp(-t.theta/2)*(t.theta/4-1)
        )*t.sigma**2 + 1/24*3*t.sigma**4*(
            sqrt(t)*exp(t.theta/2)*(t.theta/16 + 1/2) - 1/(8*sqrt(t))*exp(-t.theta/2)*(t.theta/16 - 1/2)
        )

        variance = (
            sqrt(t)*exp(t.theta/2)*(1 + t.theta/2) - 1/(8*sqrt(t))*exp(-t.theta/2)*(1 - t.theta/2)
        )**2*t.sigma**2 + (
            sqrt(t)*exp(t.theta/2)*(1 + t.theta/4) - 1/(8*sqrt(t))*exp(-t.theta/2)*(t.theta/4-1)
        )**2*1/4*t.sigma**4 + 1/4*(t.sigma**2 + t.theta**2) + 1 + t.theta

        pselect <- function(alpha){
            
            select.tmp = pnorm(
                (beta*expectation + alpha)/(
                    sqrt( 1 + beta^2*variance )
                )
            )
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/0.7 )
        }
        
        alpha.set <- uniroot( pselect, lower = -10, upper = 10, extendInt = 'yes' )$root
        
        t.prop <- pnorm( alpha.set + beta*t.par )
        t.sample <- sapply( t.prop, FUN = function(x) {rbinom(1, 1, x)} ) %>% as.logical()
        
        data.org = data.frame( y = y.give, t = t )
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
op.org.all <- sapply( data.org.all, FUN = function(x){ c(sum(x[['y']]), sum(x[['t']]))  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.org.all[1]/op.org.all[2])
print('Average events proportion (after select)')
op.select.all <- sapply( data.select.all, FUN = function(x){ c(sum(x[['y']]), sum(x[['t']]))  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.select.all[1]/op.select.all[2])

write.csv( (
    list(alpha = alpha.kt, oporg = op.org.all[1]/op.org.all[2], 
         opselect = op.select.all[1]/op.select.all[2] ) %>% as.data.frame()
), file = paste0('summary_1_', s, '.csv' ) )

rm(list = c('n.datasets', 'beta'))

estimate.fun.all <- function(data = data.org, start.ini = c(-3.6, 0.28)){
    nstudy = nrow( data )
    y = data[['y']]
    t = data[['t']]
    
    start.p <- start.ini
    
    llk <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dpois( y, lambda = t*exp( theta*sigma.ini + theta.ini ) )*dnorm( theta )
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

estimate.fun.nn <- function(data = data.select, start.ini = c(-3.6, 0.28, 0.4), method = 'nlm',  p = 0.7){
    nstudy = nrow(data)
    y = data[['y']]
    t = data[['t']]
    
    g.tmp <- (data.frame( y = y, t = t ) %>%
                  dplyr::mutate( ytmp = ifelse( (y==0), y+0.5, y )) %>%
                  dplyr::mutate( ttmp = log(ytmp/t)/sqrt(
                      ( 1/ytmp )
                  ) ))
    yp = with(g.tmp, log(ytmp/t))
    v = with( g.tmp, 1/ytmp )
    
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
    t = data[['t']]
    
    g.tmp <- (data.frame( y = y, t = t ) %>%
                  dplyr::mutate( ytmp = ifelse( (y==0), y+0.5, y )) %>%
                  dplyr::mutate( ttmp = log(ytmp/t)/sqrt(
                      ( 1/ytmp )
                  ) ))
    
    llk.com <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dpois( y, lambda = t*exp( theta*sigma.ini + theta.ini ) )*dnorm( theta )
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )[['integral']]
        
        -sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
    }
    
    start.p <- try(nlminb(start.ini[1:2], llk.com, lower=c(-Inf, 0.01), upper = c(Inf, 10)), silent=T)[['par']]
    if(start.p[2] %in% c(0.01, 10) ){
        start.p[2] = start.ini[2]
    }
    start.p <- c(start.p, start.ini[3])
    
    t.tmp = g.tmp[['ttmp']]
    v = with(g.tmp, 1/ytmp)
    
    llk <- function( par ){
        theta.ini = par[1]
        sigma.ini = par[2]
        beta.ini = par[3]
        
        expectation = sqrt(t)*exp(theta.ini/2)*theta.ini - 1/(8*sqrt(t))*theta.ini*exp(-theta.ini/2) + 1/2*(
            sqrt(t)*exp(theta.ini/2)*(1 + theta.ini/4) - 1/(8*sqrt(t))*exp(-theta.ini/2)*(theta.ini/4-1)
        )*sigma.ini**2 + 1/24*3*sigma.ini**4*(
            sqrt(t)*exp(theta.ini/2)*(theta.ini/16 + 1/2) - 1/(8*sqrt(t))*exp(-theta.ini/2)*(theta.ini/16 - 1/2)
        )
        
        variance = (
            sqrt(t)*exp(theta.ini/2)*(1 + theta.ini/2) - 1/(8*sqrt(t))*exp(-theta.ini/2)*(1 - theta.ini/2)
        )**2*sigma.ini**2 + (
            sqrt(t)*exp(theta.ini/2)*(1 + theta.ini/4) - 1/(8*sqrt(t))*exp(-theta.ini/2)*(theta.ini/4-1)
        )**2*1/4*sigma.ini**4 + 1/4*(sigma.ini**2 + theta.ini**2) + 1 + theta.ini
        
        pselect <- function(alpha, print = FALSE){
            
            select.tmp = pnorm(
                (beta.ini*expectation + alpha)/(
                    sqrt( 1 + beta.ini^2*variance )
                )
            )
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/p  )
        }
        
        alpha.opt = try( uniroot( pselect, lower = -10, upper = 10, extendInt = 'yes' )$root )
        
        prob.prior <- cubature::hcubature( f = function(theta){
            dpois( y, lambda = t*exp( theta*sigma.ini + theta.ini ) )*dnorm( theta )
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4 )[['integral']]
        
        prob.denom <- pnorm(
            (beta.ini*expectation + alpha.opt)/(
                sqrt( 1 + beta.ini^2*variance )
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
rowMeans(tidy.com.noadjust, na.rm = TRUE) %>% print()
print.com.noadjust <- rbind( rowMeans(tidy.com.noadjust, na.rm = TRUE),
                             apply(tidy.com.noadjust-c(t.theta, t.sigma), 1, FUN = function(x) mean( (x - mean(x, na.rm = TRUE))**2, na.rm = TRUE ) ),
                             sapply( result.com.noadjust, FUN = function(x) x[['coverage']] ) %>% rowMeans(na.rm = TRUE) )

tidy.nn <- sapply(
    result.nn, FUN = function(x){
        if ( class(x[['result']]) == 'try-error'){
            return(c(NA, NA))
        }else{
            if(x[['result']]$convergence == 0){ x[['result']]$par[1:2]}else{c(NA, NA)}
        }
    }, simplify = TRUE
)
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
