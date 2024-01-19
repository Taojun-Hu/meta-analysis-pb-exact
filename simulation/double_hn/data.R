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
beta <- 0.5
set.seed(1234)

data.generate <- lapply(
    1:n.datasets, FUN = function(ind){
        study.theta <- rnorm( s, mean = t.theta, sd = t.sigma)
        n0_all <- runif( s, min = 200, max = 400 ) %>% round()
        n1_all <- runif( s, min = 200, max = 400 ) %>% round()
        y_all <- runif( s, min = 10, max = 20) %>% round()
        
        data.opt <- lapply( data.frame( y_all = y_all, n1 = n1_all, n0 = n0_all ) %>% t()  %>% as.data.frame() %>% `colnames<-`(NULL) %>% 
                                as.list(), FUN = function(x){
                                    data.frame( tp = 0:x[1], n1 = x[2], n0 = x[3] ) %>%
                                        dplyr::mutate( fp = x[1] - tp  ) %>%
                                        dplyr::mutate( fn = n1 - tp, tn = n0 - fp )
                                }  )
        
        prob.sample <- mapply(FUN = function(data.x, theta.x){
            weight.log <- with( data.x, log( choose( fn + tp, tp ) ) + log( choose( fp + tn, fp )) + 
                                    theta.x*tp)
            exp( weight.log - max(weight.log) )/sum( exp( weight.log - max(weight.log) ) )
            
        }, data.opt, study.theta, SIMPLIFY = FALSE)
        
        y_i <- mapply( FUN = function(x, y){
            sample( 0:x, size = 1, prob = y  )
        }, y_all, prob.sample, SIMPLIFY = TRUE )
        
        data.tmp <- data.frame(
            tp = y_i, fp = y_all - y_i
        ) %>% dplyr::mutate( fn = n1_all - tp, tn = n0_all - fp ) %>%
            dplyr::mutate( correction = ( tp ==0 | fn ==0 | fp ==0 | tn == 0 )*0.5 ) %>%
            dplyr::mutate( y = log( ( tp + correction )*( tn + correction )/( fp + correction )/( fn + correction ) ) ) %>%
            dplyr::mutate( v = 1/( tp + correction ) + 1/( tn + correction ) + 1/( fp + correction ) + 1/( fn + correction ) ) %>%
            dplyr::select( -correction ) %>%
            dplyr::mutate( tmp = y/sqrt(v) )
        
        t.par = with( data.tmp, y/sqrt(v) )
        
        # data.tmp is the generated whole dataset
        
        # GENERATE data.opt for imputing ALPHA
        
        expectation <- (sqrt(y_all)*t.theta/(exp(-(t.theta + log(n1_all/n0_all))/2) + exp((t.theta + log(n1_all/n0_all))/2)) -1/(8*sqrt(y_all))*(
            t.theta*(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))
        ) + 1/2*t.sigma**2*(
            sqrt(y_all)*(
                exp(-3*(t.theta + log(n1_all/n0_all))/2)*(1 + t.theta/4) - exp(3*(t.theta + log(n1_all/n0_all))/2)*(1 - t.theta/4) + 
                    exp(-(t.theta + log(n1_all/n0_all))/2)*(1 - 5*t.theta/4) - exp((t.theta + log(n1_all/n0_all))/2)*(1 + 5*t.theta/4)
            )/(
                exp(-(t.theta + log(n1_all/n0_all))/2) + exp((t.theta + log(n1_all/n0_all))/2)
            )**4 - 1/(8*sqrt(y_all))*(
                exp((t.theta + log(n1_all/n0_all))/2)*(1+t.theta/4) + exp(-(t.theta + log(n1_all/n0_all))/2)*(t.theta/4 - 1)
            )
        ) + 1/24*3*t.sigma**4*(
            4*sqrt(y_all)*(
                -6*(1/2*exp((t.theta + log(n1_all/n0_all))/2) -1/2*exp(-(t.theta + log(n1_all/n0_all))/2) )**3/( exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2) )**4 + 
                    6*(1/4*exp(-(t.theta + log(n1_all/n0_all))/2) + 1/4*exp((t.theta + log(n1_all/n0_all))/2))*(1/2*exp((t.theta + log(n1_all/n0_all))/2) - 1/2*exp(-(t.theta + log(n1_all/n0_all))/2))/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**3 - 
                    (1/8*exp((t.theta + log(n1_all/n0_all))/2) - 1/8*exp(-(t.theta + log(n1_all/n0_all))/2))/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**2
            ) + sqrt(y_all)*t.theta*(
                24*(1/2*exp((t.theta + log(n1_all/n0_all))/2) - 1/2*exp(-(t.theta + log(n1_all/n0_all))/2))**4/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**5 - 
                    36*(1/4*exp(-(t.theta + log(n1_all/n0_all))/2) + 1/4*exp((t.theta + log(n1_all/n0_all))/2))*(1/2*exp((t.theta + log(n1_all/n0_all))/2) - 1/2*exp(-(t.theta + log(n1_all/n0_all))/2))**2/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**4 +
                    8*(1/8*exp((t.theta + log(n1_all/n0_all))/2) - 1/8*exp(-(t.theta + log(n1_all/n0_all))/2))*(1/2*exp((t.theta + log(n1_all/n0_all))/2) - 1/2*exp(-(t.theta + log(n1_all/n0_all))/2))/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**3 - 
                    (1/16*exp((t.theta + log(n1_all/n0_all))/2) +1/16*exp(-(t.theta + log(n1_all/n0_all))/2))/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**2 + 
                    6*(1/4*exp(-(t.theta + log(n1_all/n0_all))/2) + 1/4*exp((t.theta + log(n1_all/n0_all))/2))**2/(exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2))**3
            ) - 1/(8*sqrt(y_all))*(exp((t.theta + log(n1_all/n0_all))/2)*(1/2 + (t.theta + log(n1_all/n0_all))/16) + exp(-(t.theta + log(n1_all/n0_all))/2)*((t.theta + log(n1_all/n0_all))/16 - 1/2) )
        ) )
        
        variance <- ((
            sqrt(y_all)*(
                sqrt(n1_all/n0_all)*exp(t.theta/2)*(1 - t.theta/2) + sqrt(n0_all/n1_all)*exp(-t.theta/2)*(1 + t.theta/2)
            )/( exp((t.theta + log(n1_all/n0_all))/2) + exp(-(t.theta + log(n1_all/n0_all))/2) )**2 - 
                1/(8*sqrt(y_all))*(
                    sqrt(n1_all/n0_all)*exp(t.theta/2)*(1 + t.theta/2) + sqrt(n0_all/n1_all)*exp(-t.theta/2)*(1 - t.theta/2)
                )
        )**2*t.sigma**2 + (
            sqrt(y_all)*(
                exp(-3*(t.theta + log(n1_all/n0_all))/2)*(1 + t.theta/4) - exp(3*(t.theta + log(n1_all/n0_all))/2)*(1 - t.theta/4) + 
                    exp(-(t.theta + log(n1_all/n0_all))/2)*(1 - 5*t.theta/4) - exp((t.theta + log(n1_all/n0_all))/2)*(1 + 5*t.theta/4)
            )/(
                exp(-(t.theta + log(n1_all/n0_all))/2) + exp((t.theta + log(n1_all/n0_all))/2)
            )**4 - 1/(8*sqrt(y_all))*(
                exp((t.theta + log(n1_all/n0_all))/2)*(1+t.theta/4) + exp(-(t.theta + log(n1_all/n0_all))/2)*(t.theta/4 - 1)
            )
        )**2*1/4*t.sigma**4 + cubature::hcubature( f = function(theta){
            ((1/2 - 1/(1 + exp(theta*t.sigma + t.theta + log(n1_all/n0_all))))*(
                theta*t.sigma + t.theta
            ) + 1)**2*dnorm(theta)
        }, lowerLimit = -10, upperLimit = 10, tol=1e-4, fDim = s )[['integral']] )

        pselect <- function(alpha){
            
            select.tmp = pnorm(
                ( beta*expectation + alpha )/(
                        1 + beta**2*variance)) 
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/0.7 )
        }
        
        alpha.set <- uniroot(pselect, lower = -10, upper = 10, extendInt = 'yes')$root
        
        t.prop <- pnorm( alpha.set + beta*t.par )
        t.sample <- sapply( t.prop, FUN = function(x) {rbinom(1, 1, x)} ) %>% as.logical()
        
        data.org = data.tmp
        data.select = data.org[t.sample, ] %>% `rownames<-`(NULL)  
        
        return( list( data.org = data.org, data.select = data.select, alpha = alpha.set) )
    }
)


alpha.all <- sapply( data.generate, FUN = function(x) x[['alpha']], simplify = TRUE )
data.org.all <- lapply( data.generate, FUN = function(x) x[['data.org']] )
data.select.all <- lapply( data.generate, FUN = function(x) x[['data.select']] )


# alpha.kt <- mean( alpha.all, na.rm = TRUE )
# data.select.all <- lapply( data.org.all, FUN = function(x){
#     
#     
#     t.par <- x[['tmp']]
#     
#     t.prop <- pnorm( alpha.kt + beta*t.par )
#     t.sample <- sapply( t.prop, FUN = function(x) {rbinom(1, 1, x)} ) %>% as.logical()
#     
#     data.select = x[t.sample, ] %>% `rownames<-`(NULL)  
#     
#     data.select
#     
# }  )
# data.select.all <- lapply( data.generate, FUN = function(x) x[['data.select']] )

alpha.kt <- mean( alpha.all, na.rm = TRUE )

save( list =  c('alpha.all', 'data.org.all', 'data.select.all'), file =  paste0('data_', s, '.rda')  )

print('alpha')
print(mean( alpha.kt ))
print('Average events proportion (without select)')
op.org.all <- sapply( data.org.all, FUN = function(x){ c(sum(x[['tp']] + x[['fp']] ), 
                                                         sum(x[['tp']] + x[['fp']] + x[['tn']] + x[['fn']] ))  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.org.all[1]/op.org.all[2])
print('Average events proportion (after select)')
op.select.all <- sapply( data.select.all, FUN = function(x){ c(sum(x[['tp']] + x[['fp']] ), 
                                                               sum(x[['tp']] + x[['fp']] + x[['tn']] + x[['fn']] ) )  }  ) %>% rowSums() %>% `names<-`(NULL)
print( op.select.all[1]/op.select.all[2])

write.csv( (
    list(alpha = alpha.kt, oporg = op.org.all[1]/op.org.all[2], 
         opselect = op.select.all[1]/op.select.all[2] ) %>% as.data.frame()
), file = paste0('summary_1_', s, '.csv' ) )

rm(list = c('n.datasets', 'beta'))

estimate.fun.all <- function(data = data.org, start.ini = c(-3.6, 0.28)){
    nstudy = nrow( data )
    tp = data[['tp']]
    y_all <- with( data, tp + fp )
    n1 <- with( data, tp + fn )
    n0 <- with( data, fp + tn )
    
    start.p <- start.ini
    
    # prepare the data with an imputation
    data.opt <- lapply(
        data %>% dplyr::select(c('tp', 'tn', 'fp', 'fn'))  %>% t() %>% as.data.frame() %>% as.list(), 
        FUN = function(x){
            y_all <- x[1] + x[3]
            n1 <- x[1] + x[2]
            n0 <- x[3] + x[4]
            
            data.frame( tp = 0:y_all, y_all = y_all, n1 = n1, n0 = n0 ) %>%
                dplyr::mutate( fp = y_all - tp ) %>% dplyr::mutate( tn = n1 - tp, fn = n0 - fp ) %>% dplyr::relocate(
                    c(y_all, n1, n0), .after = everything()
                )
        }
    ) %>% do.call(what = 'rbind') %>%
        dplyr::mutate(
            cx = ( (tp == 0 ) | (tn == 0) | (fp == 0) | ( fn ==0) ) *0.5
        ) %>%
        dplyr::mutate( y = log( ( tp + cx )*( fn + cx )/( tn + cx ) / ( fp + cx ) ), 
                       v = 1/( tp + cx ) + 1/( fn + cx) + 1/(tn + cx) + 1/( fp + cx )) %>% dplyr::select( -cx ) %>%
        `rownames<-`(NULL)
    
    tp.all <- data.opt[['tp']]
    tn.all <- data.opt[['tn']]
    fp.all <- data.opt[['fp']]
    fn.all <- data.opt[['fn']]
    n1.all <- data.opt[['n1']]
    n0.all <- data.opt[['n0']]
    y_all.all <- data.opt[['y_all']]
    y.all <- data.opt[['y']]
    v.all <- data.opt[['v']]
    
    index.all <- mapply( FUN = function(x, y) rep(x, y + 1), seq_along( y_all ), y_all, SIMPLIFY = FALSE ) %>% 
        do.call(what = 'c')
    index.all <- lapply( 1:nstudy, FUN = function(k) which(index.all == k)   )
    
    xindex.all <- with(data, c(0, cumsum( fp + tp + 1 )[1:(nstudy - 1  )]) + tp + 1 )
    yindex.all <- cumsum(y_all + 1)
    
    llk <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature(
            f = function(theta){
                weight.log <- log( choose(n1.all, tp.all) ) + log( choose( n0.all, fp.all ) ) + (theta*sigma.ini + theta.ini)*tp.all
                
                weight <- exp( weight.log - unlist(
                    lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.all[[k]]]), y_all[k] + 1 )   )
                ) )
                
                weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.all[[k]]]), simplify = TRUE )
                
                weight[xindex.all]/weight.sum.bypart*dnorm( theta)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
        )[['integral']]
        
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
    tp = data[['tp']]
    y_all <- with( data, tp + fp )
    n1 <- with( data, tp + fn )
    n0 <- with( data, fp + tn )
    t.tmp = data[['tmp']]
    yp = with(data, y)
    v = with( data, v)
    
    start.p = start.ini
    
    llk <- function( par ){
        print(par)
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
    
    nstudy = nrow( data )
    tp = data[['tp']]
    y_all <- with( data, tp + fp )
    n1 <- with( data, tp + fn )
    n0 <- with( data, fp + tn )
    t.tmp = data[['tmp']]
    
    # prepare the data with an imputation
    data.opt <- lapply(
        data %>% dplyr::select(c('tp', 'tn', 'fp', 'fn'))  %>% t() %>% as.data.frame() %>% as.list(), 
        FUN = function(x){
            y_all <- x[1] + x[3]
            n1 <- x[1] + x[2]
            n0 <- x[3] + x[4]
            
            data.frame( tp = 0:y_all, y_all = y_all, n1 = n1, n0 = n0 ) %>%
                dplyr::mutate( fp = y_all - tp ) %>% dplyr::mutate( tn = n1 - tp, fn = n0 - fp ) %>% dplyr::relocate(
                    c(y_all, n1, n0), .after = everything()
                )
        }
    ) %>% do.call(what = 'rbind') %>%
        dplyr::mutate(
            cx = ( (tp == 0 ) | (tn == 0) | (fp == 0) | ( fn ==0) ) *0.5
        ) %>%
        dplyr::mutate( y = log( ( tp + cx )*( fn + cx )/( tn + cx ) / ( fp + cx ) ), 
                       v = 1/( tp + cx ) + 1/( fn + cx) + 1/(tn + cx) + 1/( fp + cx )) %>% dplyr::select( -cx ) %>%
        `rownames<-`(NULL)
    
    tp.all <- data.opt[['tp']]
    tn.all <- data.opt[['tn']]
    fp.all <- data.opt[['fp']]
    fn.all <- data.opt[['fn']]
    n1.all <- data.opt[['n1']]
    n0.all <- data.opt[['n0']]
    y_all.all <- data.opt[['y_all']]
    y.all <- data.opt[['y']]
    v.all <- data.opt[['v']]
    
    index.all <- mapply( FUN = function(x, y) rep(x, y + 1), seq_along( y_all ), y_all, SIMPLIFY = FALSE ) %>% 
        do.call(what = 'c')
    index.all <- lapply( 1:nstudy, FUN = function(k) which(index.all == k)   )
    
    xindex.all <- with(data, c(0, cumsum( fp + tp + 1 )[1:(nstudy - 1  )]) + tp + 1 )
    yindex.all <- cumsum(y_all + 1)
    
    # estimate the initial value
    llk.com <- function(par){
        theta.ini = par[1]
        sigma.ini = par[2]
        
        prob.prior <- cubature::hcubature(
            f = function(theta){
                weight.log <- log( choose(n1.all, tp.all) ) + log( choose( n0.all, fp.all ) ) + (theta*sigma.ini + theta.ini)*tp.all
                
                weight <- exp( weight.log - unlist(
                    lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.all[[k]]]), y_all[k] + 1 )   )
                ) )
                
                weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.all[[k]]]), simplify = TRUE )
                
                weight[xindex.all]/weight.sum.bypart*dnorm( theta)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
        )[['integral']]
        
        -sum( log( prob.prior + 1e-30 ), na.rm = TRUE )
    }
    
    start.p <- try(nlminb(start.ini[1:2], llk.com, lower=c(-Inf, 0.01), upper = c(Inf, 10)), silent=T)[['par']]
    if(start.p[2] %in% c(0.01, 10) ){
        start.p[2] = start.ini[2]
    }
    start.p <- c(start.p, start.ini[3])
    
    llk <- function( par ){
        theta.ini = par[1]
        sigma.ini = par[2]
        beta.ini = par[3]
        
        expectation <- (sqrt(y_all)*theta.ini/(exp(-(theta.ini + log(n1/n0))/2) + exp((theta.ini + log(n1/n0))/2)) -1/(8*sqrt(y_all))*(
            theta.ini*(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))
        ) + 1/2*sigma.ini**2*(
            sqrt(y_all)*(
                exp(-3*(theta.ini + log(n1/n0))/2)*(1 + theta.ini/4) - exp(3*(theta.ini + log(n1/n0))/2)*(1 - theta.ini/4) + 
                    exp(-(theta.ini + log(n1/n0))/2)*(1 - 5*theta.ini/4) - exp((theta.ini + log(n1/n0))/2)*(1 + 5*theta.ini/4)
            )/(
                exp(-(theta.ini + log(n1/n0))/2) + exp((theta.ini + log(n1/n0))/2)
            )**4 - 1/(8*sqrt(y_all))*(
                exp((theta.ini + log(n1/n0))/2)*(1+theta.ini/4) + exp(-(theta.ini + log(n1/n0))/2)*(theta.ini/4 - 1)
            )
        ) + 1/24*3*sigma.ini**4*(
            4*sqrt(y_all)*(
                -6*(1/2*exp((theta.ini + log(n1/n0))/2) -1/2*exp(-(theta.ini + log(n1/n0))/2) )**3/( exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2) )**4 + 
                    6*(1/4*exp(-(theta.ini + log(n1/n0))/2) + 1/4*exp((theta.ini + log(n1/n0))/2))*(1/2*exp((theta.ini + log(n1/n0))/2) - 1/2*exp(-(theta.ini + log(n1/n0))/2))/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**3 - 
                    (1/8*exp((theta.ini + log(n1/n0))/2) - 1/8*exp(-(theta.ini + log(n1/n0))/2))/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**2
            ) + sqrt(y_all)*theta.ini*(
                24*(1/2*exp((theta.ini + log(n1/n0))/2) - 1/2*exp(-(theta.ini + log(n1/n0))/2))**4/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**5 - 
                    36*(1/4*exp(-(theta.ini + log(n1/n0))/2) + 1/4*exp((theta.ini + log(n1/n0))/2))*(1/2*exp((theta.ini + log(n1/n0))/2) - 1/2*exp(-(theta.ini + log(n1/n0))/2))**2/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**4 +
                    8*(1/8*exp((theta.ini + log(n1/n0))/2) - 1/8*exp(-(theta.ini + log(n1/n0))/2))*(1/2*exp((theta.ini + log(n1/n0))/2) - 1/2*exp(-(theta.ini + log(n1/n0))/2))/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**3 - 
                    (1/16*exp((theta.ini + log(n1/n0))/2) +1/16*exp(-(theta.ini + log(n1/n0))/2))/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**2 + 
                    6*(1/4*exp(-(theta.ini + log(n1/n0))/2) + 1/4*exp((theta.ini + log(n1/n0))/2))**2/(exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2))**3
            ) - 1/(8*sqrt(y_all))*(exp((theta.ini + log(n1/n0))/2)*(1/2 + (theta.ini + log(n1/n0))/16) + exp(-(theta.ini + log(n1/n0))/2)*((theta.ini + log(n1/n0))/16 - 1/2) )
        ) )
        
        variance <- ((
            sqrt(y_all)*(
                sqrt(n1/n0)*exp(theta.ini/2)*(1 - theta.ini/2) + sqrt(n0/n1)*exp(-theta.ini/2)*(1 + theta.ini/2)
            )/( exp((theta.ini + log(n1/n0))/2) + exp(-(theta.ini + log(n1/n0))/2) )**2 - 
                1/(8*sqrt(y_all))*(
                    sqrt(n1/n0)*exp(theta.ini/2)*(1 + theta.ini/2) + sqrt(n0/n1)*exp(-theta.ini/2)*(1 - theta.ini/2)
                )
        )**2*sigma.ini**2 + (
            sqrt(y_all)*(
                exp(-3*(theta.ini + log(n1/n0))/2)*(1 + theta.ini/4) - exp(3*(theta.ini + log(n1/n0))/2)*(1 - theta.ini/4) + 
                    exp(-(theta.ini + log(n1/n0))/2)*(1 - 5*theta.ini/4) - exp((theta.ini + log(n1/n0))/2)*(1 + 5*theta.ini/4)
            )/(
                exp(-(theta.ini + log(n1/n0))/2) + exp((theta.ini + log(n1/n0))/2)
            )**4 - 1/(8*sqrt(y_all))*(
                exp((theta.ini + log(n1/n0))/2)*(1+theta.ini/4) + exp(-(theta.ini + log(n1/n0))/2)*(theta.ini/4 - 1)
            )
        )**2*1/4*sigma.ini**4 + cubature::hcubature( f = function(theta){
            ((1/2 - 1/(1 + exp(theta*sigma.ini + theta.ini + log(n1/n0))))*(
                theta*sigma.ini + theta.ini
            ) + 1)**2*dnorm(theta)
        }, lowerLimit = -10, upperLimit = 10, tol=1e-4, fDim = nstudy )[['integral']] )
        
        pselect <- function(alpha, print = FALSE){
            
            select.tmp <- pnorm(
                ( beta.ini*expectation+ alpha )/(
                    1 + beta.ini**2*variance) ) 
            
            return( mean( 1/select.tmp, na.rm = TRUE ) - 1/p  )
        }
        
        alpha.opt = uniroot(pselect, lower = -10, upper = 10, extendInt = 'yes')$root
        
        prob.prior <- cubature::hcubature(
            f = function(theta){
                weight.log <- log( choose(n1.all, tp.all) ) + log( choose( n0.all, fp.all ) ) + (theta*sigma.ini + theta.ini)*tp.all
                
                weight <- exp( weight.log - unlist(
                    lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.all[[k]]]), y_all[k] + 1 )   )
                ) )
                
                weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.all[[k]]]), simplify = TRUE )
                
                weight[xindex.all]/weight.sum.bypart*dnorm( theta)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
        )[['integral']]
        
        prob.denom <- pnorm(
            ( beta.ini*expectation+ alpha.opt )/(
                1 + beta.ini**2*variance )
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
