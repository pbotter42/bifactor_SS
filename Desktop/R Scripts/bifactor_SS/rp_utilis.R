# for some reason this onle works with 17 quad points

## norm_dist_2d
norm_dist_2d <- function(theta_gen, theta_spec) { 
  dist_m <- matrix(0,length(theta_gen),length(theta_gen)) 
  phi <- 0#toDo 
  for (i in 1:length(theta_gen)) { 
    for (j in 1:length(theta_spec)) { 
      dist_m[i,j] <- exp(-0.5*(theta_gen[i]^2+theta_spec[j]^2-2*phi*theta_gen[i]*theta_spec[j])/(1-phi*phi)) 
    } 
  } 
  dist_m <- dist_m/sum(dist_m) 
  return(dist_m) 
} 

## marg_dist_2d
marg_dist_2d <- function(dist_2d) {
  n_quad <- nrow(dist_2d)
  dist_m <- rep(0, n_quad)
  for (i in 1:n_quad) {
    for (j in 1:n_quad) {
      dist_m[i] <- dist_m[i] + dist_2d[i,j]
    }
  }
  return(dist_m)
}

# ts_list
ts_list <- function(a_gen, a_spec, b, rp, rp_n, ic_ref, theta) { 
  ic_n <- length(unique(ic_ref)) # calculating the number of unique item clusters 
  gaussian_1d <- marg_dist_2d(norm_dist_2d(theta, theta)) 
  ts_list <- list() # list containing  trace surfaces for all rp_n 
  
  # compute trace surfaces and storeing them in a list indexed by ic 
  ts_list <- list() 
  jth_item <- 0 # this references the item number irredardless of ic 
  for (ic_counter in 1:ic_n) { 
    ts_list[[ic_counter]] <- list() 
    ic_nitem <- length(ic_ref[ic_ref == ic_counter]) # finding the number of items in that particular ic 
    for(j in 1:ic_nitem) { 
      jth_item <- jth_item + 1 
      ts_list[[ic_counter]][[j]] <- trace_surface(a_gen = a_gen[jth_item], 
                                                  a_spec = a_spec[jth_item], 
                                                  b = b[jth_item], 
                                                  theta_gen = theta, 
                                                  theta_spec = theta) 
      
    } 
  }   
  return(ts_list) 
} 

##trace_surface
trace_surface <- function(a_gen, a_spec, b, theta_gen, theta_spec) { 
  ts <- matrix(nrow = length(theta), ncol = length(theta)) 
  for(q1 in 1:length(theta_gen)) {    
    for(q2 in 1:length(theta_spec)) {    
      ts[q1,q2] <- 1/(1 + exp(-(a_gen * theta_gen[q1] + a_spec * theta_spec[q2] + b))) 
    }    
  }  
  return(ts) 
} 

## ic_item_ref_cal2
# creates a list that contains the responses by rp -> ic -> item 
ic_item_ref_cal2 <- function(rp, ic_ref) { 
  rp_n <- nrow(rp) 
  ic_n <- length(unique(ic_ref)) 
  ic_item_ref <- list() # list referenced by ic that contains all responses for items in the ic 
  for(i in 1:rp_n) { 
    jth_item <- 0 
    ic_item_ref[[i]] <- list() 
    for(ic_counter in 1:ic_n) { 
      ic_nitem <- length(ic_ref[ic_ref == ic_counter]) # finding the number of items in that particular ic 
      ic_item_ref[[i]][[ic_counter]] <- numeric(length = ic_nitem) 
      for(j in 1:ic_nitem) { 
        jth_item <- jth_item + 1 
        ic_item_ref[[i]][[ic_counter]][[j]] <- rp[i,jth_item] 
      } 
    } 
  } 
  return(ic_item_ref) 
} 
## ic_item_ref_cal3
# creates a list that contains the responses by item -> rp 
ic_item_ref_cal3 <- function(rp) { 
  ic_item_ref <- list() 
  rp_n <- nrow(rp) 
  item_n <- ncol(rp) 
  for(j in 1:item_n) { 
    ic_item_ref[[j]] <- list() 
    for(i in 1:rp_n) { 
      ic_item_ref[[j]][[i]] <- unname(rp[i,j]) 
    } 
  } 
  return(ic_item_ref) 
} 
## LE_iis
#function that calculates:Li, Lis, Ei, and Eis 
LE_iis <- function(ts_list, a_gen, a_spec, b, rp, rp_n, rp_count, ic_ref, theta) { 
  ic_n <- length(unique(ic_ref)) # calculating the number of unique item clusters 
  gaussian_1d <- marg_dist_2d(norm_dist_2d(theta, theta)) 
  ic_item_ref <- ic_item_ref_cal2(rp, ic_ref) 
  L_is_list <- list() 
  E_is_list <- list() 
  E_i_list <- list() 
  L_i_list <- list()
  logLik <- 0
  
  for(i in 1:rp_n) { 
    L_is_list[[i]] <- list() 
    E_is_list[[i]] <- list() 
    E_i_list[[i]] <- list() 
    E_i_list[[i]][[1]] <- rep(1, length(theta)) 
    L_i_list[[i]] <- numeric() 
    for (ic_counter in 1:ic_n) { 
      L_is_list[[i]][[ic_counter]] <- list() 
      L_is_list[[i]][[ic_counter]][[1]] <- matrix(rep(1, length(theta)^2), 
                                                  nrow = length(theta), 
                                                  ncol = length(theta)) 
      E_is_list[[i]][[ic_counter]] <- numeric(length = length(theta)) 
      ic_nitem <- length(ic_item_ref[[i]][[ic_counter]]) # finding the number of items in that particular ic 
      for(j in 1:ic_nitem) { # the jth item in each ic 
        if (ic_item_ref[[i]][[ic_counter]][[j]] == 1) { 
          L_is_list[[i]][[ic_counter]][[1]] <- L_is_list[[i]][[ic_counter]][[1]] *  
            ts_list[[ic_counter]][[j]] 
        } else { 
          L_is_list[[i]][[ic_counter]][[1]] <- L_is_list[[i]][[ic_counter]][[1]] *  
            (1 - ts_list[[ic_counter]][[j]]) 
        } 
      } 
      for (q in 1:length(theta)) {  
        temp <- L_is_list[[i]][[ic_counter]][[1]][q,]*gaussian_1d # temp vector containing the product  
        E_is_list[[i]][[ic_counter]][q] <- sum(temp) # integrating out the specific dim 
      } 
    } 
    for(ic_counter2 in 1:ic_n) { 
      E_i_list[[i]][[1]] <- E_is_list[[i]][[ic_counter2]]*E_i_list[[i]][[1]] # for rp i, taking the product of each ic joint lik with the specific dims integrated out 
    } 
    temp2 <- E_i_list[[i]][[1]]*gaussian_1d # temp vector containing the product of each ic joint lik with the specific dim integrated out 
    L_i_list[[i]] <- sum(temp2)
    logLik <- rp_count[i]*log(L_i_list[[i]])
  } 
  out <- list()   
  out[["L_is_list"]] <- L_is_list 
  out[["E_is_list"]] <- E_is_list 
  out[["E_i_list"]] <- E_i_list 
  out[["L_i_list"]] <- L_i_list
  out[["logLik"]] <- logLik
  return(out) 
}
# r_ij
r_ij <- function(LE_iis, a_gen, a_spec, b, rp, rp_n, ic_ref, theta) { 
  
  L_is_list <- LE_iis[["L_is_list"]] 
  E_is_list <- LE_iis[["E_is_list"]] 
  E_i_list <- LE_iis[["E_i_list"]] 
  L_i_list <- LE_iis[["L_i_list"]] 
  
  ic_n <- length(unique(ic_ref)) # calculating the number of unique item clusters 
  ic_item_ref <- ic_item_ref_cal2(rp, ic_ref) 
  r_ij_list <- list() 
  for(i in 1:rp_n) { 
    jth_item <- 0 
    r_ij_list[[i]] <- list() 
    for(ic_counter in 1:ic_n) { 
      ic_nitem <- length(ic_item_ref[[ic_counter]]) # finding the number of items in that particular ic 
      for(j in 1:ic_nitem) { 
        jth_item <- jth_item + 1 
        r_ij_list[[i]][[j]] <- as.matrix((E_i_list[[i]][[1]]/E_is_list[[i]][[ic_counter]])*(L_is_list[[i]][[ic_counter]][[1]]/L_i_list[[i]])) 
      } 
    } 
  } 
  return(r_ij_list) 
} 
## r_ij_transform
#change the format of r_ij, so that instead of being testee -> ic -> item w/in ic. it is item -> testee, or actually it is r_ij (of which there is only one per ts) -> testee 
r_ij_transform <- function(r_ij_object, ic_ref) { 
  out <- list() 
  for(i in 1:length(r_ij_object)) { 
    counter <- 0 
    out[[i]] <- list() 
    for(ic_counter in 1:length(r_ij_object[[i]])) { 
      counter <- counter + 1 
      out[[i]][[counter]] <- as.matrix(r_ij_object[[i]][[ic_counter]]) 
    } 
  } 
  
  out_t <- purrr::transpose(out) 
  out_t2 <- list() 
  jth_item <- 0 
  for(ic_counter in 1:length(unique(ic_ref))) { # so for "ic_counter" in the number of ics 
    #jth_item <- 0 
    ic_nitem <- length(ic_ref[ic_ref == ic_counter]) # the number of intems in each ic 
    for(j in 1:ic_nitem) { 
      jth_item <- jth_item + 1 
      out_t2[[jth_item]] <- list() 
      out_t2[[jth_item]] <- out_t[[ic_counter]] 
    } 
  } 
  return(out_t2) 
} 

##r_jk_allInOne
r_jk_allInOne <- function(a_gen, 
                          a_spec, 
                          b, 
                          rp, 
                          rp_n, # the tota number of rp
                          rp_count, # number of testees with each rp
                          ic_ref, 
                          theta) { 
  # out is a list that will contain an item wise list, that also have sub lists by correct or incorrect responses. 
  # put differently this list contains the E tables for each response 
  out <- list() 
  
  ts_list_object <- ts_list(a_gen = a_gen, a_spec = rep(.4, 6), b = d, rp = rp, rp_n = rp_n, ic_ref = ic_ref, theta = theta)  
  
  LE_iis_object <- LE_iis(ts_list = ts_list_object, 
                          a_gen = a_gen, a_spec = rep(.4, 6), b = d, rp = rp, rp_n = rp_n,
                          rp_count = rp_count, ic_ref = ic_ref, theta = theta)
  
  r_ij_object <- r_ij(LE_iis = LE_iis_object,
                      a_gen = a_gen, a_spec = rep(.4, 6), b = d, rp = rp, rp_n = rp_n, ic_ref = ic_ref, theta = theta)
  
  ic_item_ref <- ic_item_ref_cal3(rp) 
  r_ij_transformed <- r_ij_transform(r_ij_object = r_ij_object, ic_ref = ic_ref) 
  E_tables <- list() 
  E_tables[["p"]] <- list() 
  E_tables[["q"]] <- list()   
  for(j in 1:length(r_ij_transformed)) { # for j in 1:the number of items 
    zeroTraceSurface <- matrix(rep(0,length(theta)*length(theta)), 
                               nrow = length(theta), 
                               ncol = length(theta))   
    E_tables[["p"]][[j]] <- zeroTraceSurface # initializing the E table for item j with a matrix of zeros 
    E_tables[["q"]][[j]] <- zeroTraceSurface # initializing the E table for item j with a matrix of zeros 
    for(i in 1:length(r_ij_transformed[[j]])) {# for i in 1:number of testees 
      if (ic_item_ref[[j]][[i]] == 1) { 
        E_tables[["p"]][[j]] <- E_tables[["p"]][[j]] + rp_count[i]*r_ij_transformed[[j]][[i]] 
      } else { 
        E_tables[["q"]][[j]] <- E_tables[["q"]][[j]] + rp_count[i]*r_ij_transformed[[j]][[i]] 
      } 
    } 
  } 
  out[["E_tables"]] <- E_tables
  out[["logLik"]] <- LE_iis_object[["logLik"]]
  return(out) 
}

## bifactorll
bifactorll <- function(p, r1, r0, theta) { 
  a_gen <- p[1] 
  a_spec <- p[2] 
  b <- p[3] 
  
  ts <- trace_surface(a_gen, a_spec, b, theta, theta) 
  # doing a single sum across the the trace surface may be a problem. but i am sticking with it for now
  l <- (-1)*(sum(r1*(log(ts))) + (sum(r0*(log(1.0-ts))))) 
  #l <- (sum(r1*(log(ts))) + (sum(r0*(log(1.0-ts))))) 
  return(l)
}

bifactorEM <- function(a_gen_start, 
                       a_spec_start, 
                       b_start, 
                       rp, 
                       rp_n,
                       rp_count,
                       ic_ref, 
                       theta,
                       niter = 25) {
  
  out <- list()
  out[[1]] <- list()
  paramIterList <- list()
  logLik_vector <- numeric()
  paramIterList[[1]] <- list()
  paramIterList[[1]][[1]] <- numeric()
  paramIterList[[1]][[2]] <- numeric()
  paramIterList[[1]][[3]] <- numeric()
  r_jk_allInOne_out <- r_jk_allInOne(a_gen = a_gen_start, a_spec = a_spec_start,b = b_start, 
                                     rp = rp, rp_n = nrow(rp), rp_count = rp_count, ic_ref = ic_ref,
                                     theta = theta)
  E_tables <- r_jk_allInOne_out[["E_tables"]]
  r1 <- E_tables[["p"]]
  r0 <- E_tables[["q"]]
  logLik_vector[1] <- r_jk_allInOne_out[["logLik"]]
  
  for(j in 1:length(a_spec_start)) {
    out[[1]][[j]] <- nlm(f = bifactorll, p = c(a_gen_start[j],a_spec_start[j], b_start[j]),r1 = r1[[j]], r0 = r0[[j]], theta = theta)
    paramIterList[[1]][[1]][j] <- out[[1]][[j]]$estimate[1]
    paramIterList[[1]][[2]][j] <- out[[1]][[j]]$estimate[2]
    paramIterList[[1]][[3]][j] <- out[[1]][[j]]$estimate[3]
  }
  
  
  for(iter in 2:niter) {
    out[[iter]] <- list()
    paramIterList[[iter]] <- list()
    paramIterList[[iter]][[1]] <- numeric()
    paramIterList[[iter]][[2]] <- numeric()
    paramIterList[[iter]][[3]] <- numeric()
    
    r_jk_allInOne_out <- r_jk_allInOne(a_gen = unlist(paramIterList[[iter-1]][[1]]),
                                       a_spec = unlist(paramIterList[[iter-1]][[2]]),
                                       b = unlist(paramIterList[[iter-1]][[3]]), 
                                       rp = rp, rp_n = nrow(rp), rp_count = rp_count, ic_ref = ic_ref,
                                       theta = theta)
    E_tables <- r_jk_allInOne_out[["E_tables"]]
    r1 <- E_tables[["p"]]
    r0 <- E_tables[["q"]]
    logLik_vector[iter] <- r_jk_allInOne_out[["logLik"]]
    
    for(j in 1:length(a_spec_start)) {
      out[[iter]][[j]] <- nlm(f = bifactorll, p = c(paramIterList[[iter-1]][[1]][j],
                                                    paramIterList[[iter-1]][[2]][j],
                                                    paramIterList[[iter-1]][[3]][j]),
                              r1 = r1[[j]], r0 = r0[[j]], theta = theta)
      paramIterList[[iter]][[1]][j] <- out[[iter]][[j]]$estimate[1]
      paramIterList[[iter]][[2]][j] <- out[[iter]][[j]]$estimate[2]
      paramIterList[[iter]][[3]][j] <- out[[iter]][[j]]$estimate[3]
    }
  }
  out2 <- list()
  out2[["paramIterList"]] <- paramIterList
  out2[["logLik"]] <- logLik_vector
  return(out2)
}

## test

theta <- seq(-4, 4, .5)
a <- matrix(c( 
  .8,.4,NA, 
  .4,.4,NA, 
  .7,.4,NA, 
  .8,NA,.4, 
  .4,NA,.4, 
  .7,NA,.4),ncol=3,byrow=TRUE) 

d <- c(-1,1.5,0.0,-1.5,1,2) 
items <- rep("2PL", 6) 
dataset <- simdata(a,d,2000,items) 
dataset2 <- matrix(nrow = nrow(dataset), 
                   ncol = ncol(dataset)) 
for(i in 1:nrow(dataset2)) { 
  for(j in 1:ncol(dataset2)) { 
    if (dataset[i,j] == 1) { 
      dataset2[i,j] = TRUE 
    } else { 
      dataset2[i,j] = FALSE 
    } 
  } 
} 

dataset.1 <- as.data.frame(dataset) %>%  
  group_by(Item_1,Item_2,Item_3,Item_4,Item_5,Item_6) %>%  
  summarise(count = n()) 
dataset2 <- dataset[,-ncol(dataset.1)] 
rp_n <- dataset.1[,ncol(dataset.1)] 

# collapsed data set
ds_short <- as.data.frame(dataset) %>%  
  group_by(Item_1,Item_2,Item_3,Item_4,Item_5,Item_6) %>%  
  summarise(count = n())
ds_short_rp <- ds_short[,-ncol(ds_short)] %>% 
  as.data.frame()
ds_short_rpn <- ds_short[,ncol(ds_short)] %>% 
  as.data.frame()
ds_short_rpn2 <- as.numeric(ds_short_rpn[,1])

bifactorEM_foo <- bifactorEM(a_gen_start = c(.8,.4,.7,.8,.4,.7), 
                             a_spec_start = rep(.4, 6), 
                             b_start = d, 
                             rp = ds_short_rp, 
                             rp_n = nrow(ds_short_rp), 
                             rp_count = ds_short_rpn2,
                             ic_ref = c(1, 1, 1, 
                                        2, 2, 2), 
                             theta = theta,
                             niter = 100)

