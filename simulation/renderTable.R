library( tidyr )
library( rio )

theta_xl <- -3 + c(-1.96, 0, 1.96)*0.3
print(plogis( theta_xl ))

table_all <- summary_all <- vector(mode = 'list', length = 5)

ind = 1
for (dir_x in c('double_bn', 'double_pn', 'double_hn', 'single_bn', 'single_pn')){
    s_all <- c(25)
    list_summary <- list_theta_mean <- list_theta_mse <- list_theta_sd <- list_sigma_sd <-
        list_sigma_mean <- list_sigma_mse <- vector(mode = 'list', length = length(s_all)) 
    
    for (i in seq_along(s_all)){
        s <- s_all[i]
        x <- read.csv(paste0(dir_x, '/summary_1_', s, '.csv'))
        list_summary[[i]] = x[, 2:4] %>% unlist() %>% `names<-`(NULL)
        
        y <- rio::import_list( paste0(dir_x, '/printout_', s, '.xlsx') )
        
        load( paste0(dir_x, '/allsave_', s, '.rda') )
        
        # list_theta_mean[[i]] <- lapply( y, FUN = function(t){t[11,  ] %>% unlist()} ) %>% unlist()
        list_theta_mse[[i]] <- lapply( y, FUN = function(t){t[21,  ] %>% unlist()} ) %>% unlist()
        # list_sigma_mean[[i]] <- lapply( y, FUN = function(t){t[1, 2] %>% unlist()} ) %>% unlist()
        list_sigma_mse[[i]] <- lapply( y, FUN = function(t){t[2, 2] %>% unlist()} ) %>% unlist()
        
        list_theta_mean[[i]] <- list(
            asymp = mean(tidy.asymp[1,  ], na.rm = TRUE),
            com = mean(tidy.com[1,  ], na.rm = TRUE),
            com.asymp = mean(tidy.com.noadjust[1,  ], na.rm = TRUE),
            nn = mean(tidy.nn[1,  ], na.rm = TRUE)
        ) %>% unlist()
        
        list_sigma_mean[[i]] <- list(
            asymp = mean(tidy.asymp[2,  ], na.rm = TRUE),
            com = mean(tidy.com[2,  ], na.rm = TRUE),
            com.asymp = mean(tidy.com.noadjust[2,  ], na.rm = TRUE),
            nn = mean(tidy.nn[2,  ], na.rm = TRUE)
        ) %>% unlist()
        
        
        list_theta_sd[[i]] <- list(
            asymp = sd(tidy.asymp[1,  ], na.rm = TRUE),
            com = sd(tidy.com[1,  ], na.rm = TRUE),
            com.asymp = sd(tidy.com.noadjust[1,  ], na.rm = TRUE),
            nn = sd(tidy.nn[1,  ], na.rm = TRUE)
        ) %>% unlist()
        
        list_sigma_sd[[i]] <- list(
            asymp = sd(tidy.asymp[2,  ], na.rm = TRUE),
            com = sd(tidy.com[2,  ], na.rm = TRUE),
            com.asymp = sd(tidy.com.noadjust[2,  ], na.rm = TRUE),
            nn = sd(tidy.nn[2,  ], na.rm = TRUE)
        ) %>% unlist()
        
    }
    
    table_mean <- rbind(
        list_theta_mean[[1]], 
        list_sigma_mean[[1]]
        # list_theta_mean[[2]], 
        # list_sigma_mean[[2]],
        # list_theta_mean[[3]], 
        # list_sigma_mean[[3]]
        # list_sigma_mean[[3]],
        # list_theta_mean[[4]], 
        # list_sigma_mean[[4]]
    ) %>% as.matrix() %>% round(3) %>% as.data.frame() %>% `rownames<-`(NULL) %>%
        `colnames<-`(c('asymp', 'com', 'com.noadjust', 'nn'))
    
    table_sd <- rbind(
        list_theta_sd[[1]], 
        list_sigma_sd[[1]]
        # list_theta_sd[[2]], 
        # list_sigma_sd[[2]],
        # list_theta_sd[[3]], 
        # list_sigma_sd[[3]]
        # list_sigma_sd[[3]],
        # list_theta_sd[[4]], 
        # list_sigma_sd[[4]]
    ) %>% as.matrix() %>% round(3) %>% as.data.frame() %>% `rownames<-`(NULL) %>%
        `colnames<-`(c('asymp', 'com', 'com.noadjust', 'nn'))
    
    table_final <- data.frame(
        asymp = paste0(table_mean[['asymp']], '(', table_sd[['asymp']], ')'), 
        com = paste0(table_mean[['com']], '(', table_sd[['com']], ')'), 
        com.noadjust = paste0(table_mean[['com.noadjust']], '(', table_sd[['com.noadjust']], ')'), 
        nn = paste0(table_mean[['nn']], '(', table_sd[['nn']], ')')
    )
    
    table_all[[ind]] = table_final
    summary_all[[ind]] = list_summary[[1]]
    
    write.csv(table_final, paste0(dir_x, '.csv'))
    
    ind = ind+1
}

table_xlsx <- do.call('rbind', table_all)
summary_xlsx <- do.call('rbind', summary_all) %>% as.data.frame() %>%
    `colnames<-`(paste0('V', 1:3)) %>%
    dplyr::mutate( V2 = V2*100, V3 = V3*100 ) %>%
    as.matrix() %>% round(2) %>% as.data.frame()



