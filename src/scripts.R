library(tidyverse)
theme_set(theme_light())

# library(GGally)
# library(grafify)

# library(rosalia)

# library(corpcor)
library(arm)
# library(doParallel)
# library(ramify)

# library(Rlab)
library(tensor)
# library(nortest)
library(modelr)
# library(XNomial)
# library(AMR)
library(rstatix)

library(entropy)
library(rsample)

getOption("warn")
options(warn = -1)

method = 'shrink'

scaling = 2
plain_reg = 1
bayes_reg = 1

compute_n_params = function(n_unique, order, n_spp=8){
    n_params = 1 + n_spp
    if (order >=2){n_params = n_params + n_spp * (n_spp-1) /2}
    if (order >=3){n_params = n_params + n_spp * (n_spp-1)/3 * (n_spp-2)/2}
    return(n_params)
}                  

compare_preds = function(pred, obs, min_species_for_comp = 1, only_dark_patches = FALSE){
    
    pred = pred %>%
        as_tibble() %>%
        filter(min_species_for_comp <= rowSums(.[1:ncol(pred)]))  %>%
        rename_with(~colnames(obs)) %>%
        group_by_all() %>%
        summarize(pred = n()) %>%
        ungroup()

    obs = obs %>%
        as_tibble() %>%
        filter(min_species_for_comp <= rowSums(.[1:ncol(obs)])) %>%
        group_by_all() %>%
        summarize(obs = n()) %>%
        ungroup()    
        
    if (only_dark_patches == T){
        results = obs %>%
            left_join(pred)
        } else {
        results = obs %>%
            full_join(pred)
        }
    
    results = results %>%
        replace_na(list(pred = 0, obs = 0)) %>% 
        ungroup() %>%
        mutate(
           # n_species = rowSums(.[1:ncol(obs)]),
           n_obs = sum(obs),
           n_pred = sum(pred),
           pred_prop = pred / sum(pred),
           obs_prop = obs / sum(obs),
           diff_prop = obs_prop - pred_prop,
           pred_prop_err_2 = pred_prop * (1 - pred_prop) / sum(pred),
           obs_prop_err_2 = obs_prop * (1 - obs_prop) / sum(obs),
           resid_prop = diff_prop / sqrt( pred_prop_err_2 + obs_prop_err_2 ),

           pred_bayes = (pred + bayes_reg) / (sum(pred) + 2 * bayes_reg),
           obs_bayes = (obs + bayes_reg) / (sum(obs) + 2 * bayes_reg),
           diff_bayes = obs_bayes - pred_bayes,
           obs_bayes_err_2 = obs_bayes * (1 - obs_bayes) / sum(obs_bayes),
           resid_bayes = diff_bayes / sqrt(obs_bayes_err_2),

           pred_scaled = pred * sum(obs) / sum(pred),
           diff_count = obs - pred_scaled,
           resid_count = diff_count / sqrt((pred_scaled + obs)/scaling),

           resid_count_b = diff_count / sqrt(pred_scaled),

           pred_reg=replace(pred, pred==0, plain_reg),
           obs_reg=replace(obs, obs==0, plain_reg),
           pred_reg = pred_reg * sum(obs_reg) / sum(pred_reg),
           diff_reg = obs_reg - pred_reg,
           resid_reg = diff_reg / sqrt(obs_reg),
            
            obs_shrink = n_obs * freqs(obs, method=method, verbose=F),
            pred_shrink = freqs(pred, method=method, verbose=F)

        ) %>%
        ungroup()
    
    summarized = results %>%
        summarize(
        # order = order,
        # x2_props = sum(resid_prop^2) / n_dofs,
        n_unique = n(),# - compute_n_params(n_unique, order, n_spp = ncol(x_nat)),
        x2_counts = sum(resid_count^2),#/ n_unique,
        # G_value = g.test(obs, pred) %>% pluck(1,1),
        pearson = cor.test(obs, pred, method='pearson', na.rm=T)$estimate,
        kendall = cor.test(obs, pred, method='kendall', na.rm=T)$estimate,

        Gstat_shrink = Gstat(obs_shrink, pred_shrink) %>% pluck('stat'),
        # pval_G_s = Gstat(obs_shrink, pred_shrink) %>% pluck('pval'),
        Gstat_obs = Gstat(obs, pred_shrink) %>% pluck('stat'),
        # pval_G_o = Gstat(obs, pred_shrink) %>% pluck('pval'),
        # pval_Chisq = chi2stat(obs_shrink, pred_shrink) %>% pluck('pval'),
        Chisq = chi2stat(obs_shrink, pred_shrink) %>% pluck('stat'),

        # x2_red = x2_counts / n_dofs,
        # x2_count_b = sum(resid_count_b^2) / n_dofs,
        # x2_reg = sum(resid_reg^2) / n_dofs,
        # x2_bayes = sum(resid_bayes^2) / n_dofs,
        # Shapiro_Wilk = shapiro.test(resid_count)[[2]],
        # Anderson_Darling = ad.test(resid_count)[[2]],
        # Kolmogorov_Smirnov = ks.test("pnorm",
        #                              mean=mean(resid_count),
        #                              sd=sd(resid_count)),
        # Lilliefors = lillie.test(resid_count)[[2]],
        # Shapiro_Francia = sf.test(resid_count)[[2]],
        # Pearson = pearson.test(resid_count)[[2]]
             ) %>%
        suppressMessages()

    return(list('results' = results, 'summarized' = summarized))
    # return(summarized)
}


compare_vae_preds = function(all_outcomes, all_targets, min_species_for_comp = 1, only_dark_patches = FALSE){

    all_targets = all_targets %>%
        as_tibble() %>%
        group_by(sampling_method, data_type, sample_id) %>%
        nest()

    all_outcomes = all_outcomes %>%
        as_tibble() %>%
        group_by(sampling_method, data_type, sample_id) %>%
        nest()

    results = all_outcomes %>%
        inner_join(all_targets, by = join_by(sampling_method, data_type, sample_id)) %>%
        mutate(results = map2(data.x, data.y,
                      \(outcome, target)
                          compare_preds(outcome, target,
                                                  min_species_for_comp, only_dark_patches) %>%
                              pluck('summarized'))) %>%
        suppressMessages() %>%
        dplyr::select(c(data_type, sample_id, sampling_method, results)) %>%
        unnest()
    
    return(results)
}

compare_with_self = function(data, pluck='summarized', min_species_for_comp=1, n_replicas=100){

    comp = tibble()

    for (i in 1:n_replicas){

        data_splits = data %>% as_tibble() %>%
                        resample_partition(c(train=train_frac, test=test_frac))

        comp = data_splits %>% pluck('train') %>% as.data.frame() %>%
                        compare_preds(data_splits %>% pluck('test') %>%
                                          as.data.frame(),
                                      min_species_for_comp) %>%
                        pluck(pluck) %>%
                        mutate(train_frac = train_frac,
                               test_frac = test_frac) %>%
                        bind_rows(comp) %>%
                        suppressMessages()
    }
    
    return(comp)
}

melt_beta = function(beta){
        beta %>%
        as_tibble() %>%
        mutate(row = colnames(beta)) %>%
        pivot_longer(cols=!row)
}

heatmap_betas = function(beta){
    beta %>%
    melt_beta() %>%
    ggplot(aes(y=row, x=name, fill=value, label=value)) +
        geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
    #     # geom_text(aes(label = value), color = "black", size = 4) +
        scale_fill_gradient2(low = "blue",
                           mid = "white",
                           high = "red") +
        # scale_y_discrete(limits=rev) +
        # scale_y_reverse() +
        coord_fixed()
}

extract_coefs = function(x, order = 3, min_species_for_fit = 0, symmetrize_betas = T){
    
    species_names = colnames(x)
    
    x = as.matrix(x)

    alpha = matrix(0, ncol(x))
    beta_fit = matrix(0, ncol(x), ncol(x))
    gamma_tensor_fit = array(0, dim=c(ncol(x), ncol(x), ncol(x)))
    gamma_matrix_fit = array(0, dim=c(ncol(x), (ncol(x) - 1) * (ncol(x) - 2) / 2)) #array(0, c(ncol(x), ncol(x), ncol(x)))
    gamma_tensor_compact = array(0, dim=c(ncol(x), ncol(x), ncol(x)))
    gamma_tensor_sym = array(0, dim=c(ncol(x), ncol(x), ncol(x)))
    # coef_se = matrix(0, ncol(x), ncol(x))
    
    results = tibble()
    
    for(i in 1:ncol(x)){
        if(var(x[,i]) > 0){
            
            X = x %>% as_tibble() %>%
            # if (include_empty_patches == FALSE){
            # X = X %>% 
                filter(min_species_for_fit <= rowSums(pick(-i))) %>% #if min_species_for_fit = == 1, we only fit on patches that have at least one other species present, as we don't know the number of "empty" patches
                as.matrix() 
            # }
            # X = X %>% as.matrix()
            
            if (order == 1){mod = bayesglm(X[,i] ~ 1, family = binomial)}
            if (order == 2){mod = bayesglm(X[,i] ~ X[ ,-i], family = binomial)}
            if (order == 3){
                X_ = X[ , -i]
                XX = array(0, dim=c(nrow(X), (ncol(X) - 1) * (ncol(X) - 2) / 2))

                for (s in 1:nrow(X)){
                    jk = 0
                    for (j in 1:(ncol(X_) - 1)){ # assume j < k, symmetrize later
                        for (k in (j+1):ncol(X_)){
                            jk = jk + 1
                            XX[s, jk] = X_[s, j] * X_[s, k]
                        }
                    }
                }

                mod = bayesglm(X[,i] ~ X[, -i] + XX[ , ], family = binomial)#, prior.df = Inf)
            }

            # mod = glm(X[,i] ~ 1, family = binomial)
            
            alpha[i] = coef(mod)[1]
            
            if (order >= 2){
                estimcoef = coef(mod)[2:ncol(x)]
                                
                beta_fit[i, -i] = estimcoef   
            }
            
            if (order == 3){
                
                estimcoef = coef(mod)[-(1:ncol(x))]#, dim=(nrow(x), ncol(x)))

                gamma_matrix_fit[i, ] = estimcoef

                jk = 0
                for (j in 1:(ncol(X_) - 1)){
                    for (k in (j+1):ncol(X_)){ # careful for later: we only picked the case j < k, so half of the cases
                        jk = jk + 1
                        J = (1:(ncol(X)))[-i][j]
                        K = (1:(ncol(X)))[-i][k]

                        gamma_tensor_fit[i, J, K] = gamma_matrix_fit[i, jk]
                    }
                }                
        }
    }
        }
    
    # if (symmetrize_betas == T){
    beta_sym = (beta_fit + t(beta_fit)) / 2 #symmetrize because the gammas are not very symmetric if min_species_for_fit = 1
    beta_compact = pull_upper_triangle(beta_sym) %>% as.matrix()
        # }
    
    gamma_tibble_compact = tibble(Var1=character(),
                                 Var2=character(),
                                 Var3=character(),
                                 value=numeric())

    for (I in 1:(ncol(x) - 2)){
        for (J in (I+1):(ncol(x) - 1)){
            for (K in (J+1):ncol(x)){ #we have enforced I < J < K
                gamma_tensor_compact[I, J, K] = (gamma_tensor_fit[I, J, K] + gamma_tensor_fit[J, I, K] + gamma_tensor_fit[K, I, J]) /3 # + gamma_tensor_fit[J, K, I]) /3 instead of symmetrizing, we take the average the 3 cases with J < K (out of 6 possible orderings in total, and only consider I < J < K
                gamma_tensor_sym[I, J, K] = gamma_tensor_compact[I, J, K] * 1/6
                gamma_tensor_sym[I, K, J] = gamma_tensor_compact[I, J, K] * 1/6
                gamma_tensor_sym[J, I, K] = gamma_tensor_compact[I, J, K] * 1/6
                gamma_tensor_sym[J, K, I] = gamma_tensor_compact[I, J, K] * 1/6
                gamma_tensor_sym[K, I, J] = gamma_tensor_compact[I, J, K] * 1/6
                gamma_tensor_sym[K, J, I] = gamma_tensor_compact[I, J, K] * 1/6
                
                gamma_tibble_compact = gamma_tibble_compact %>% add_row(Var1 = species_names[I], Var2 = species_names[J], Var3=species_names[K], value=gamma_tensor_compact[I,J,K])
                
                
            }
        }
    }
    
    species_names = colnames(x)
    rownames(alpha) = species_names
    dimnames(beta_fit) <- list(species_names, species_names)
    dimnames(beta_sym) <- list(species_names, species_names)
    dimnames(beta_compact) <- list(species_names, species_names)
    dimnames(gamma_tensor_fit) <- list(species_names, species_names, species_names)
    dimnames(gamma_tensor_compact) <- list(species_names, species_names, species_names)
    dimnames(gamma_tensor_sym) <- list(species_names, species_names, species_names)

    return(list('alpha' = alpha, 'beta_fit' = beta_fit, 'beta_sym' = beta_sym, 'beta_compact' = beta_compact, 'gamma_matrix_fit' = gamma_matrix_fit, 'gamma_tensor_fit' = gamma_tensor_fit, 'gamma_tensor_compact' = gamma_tensor_compact, 'gamma_tensor_sym' = gamma_tensor_sym, 'gamma_tibble_compact'=gamma_tibble_compact))
    }

extend_frac_empties = function(x, frac){ #frac is the fraction of non-empty sites, i.e. grows with coverage
    x = x %>% filter(0.5 < rowSums(.[1:ncol(x)]))
    if ((frac > 0) & (frac < 1)){
        n_extra_sites = ((1-frac)/frac  * nrow(x)) %>% as.integer()
        x = matrix(0, n_extra_sites, ncol(x)) %>% as_tibble %>% rename_with(~colnames(x)) %>% bind_rows(x)# %>% filter(0.5 < rowSums(.[1:ncol(x)])))
        }
    return(x)
}


rbern = function(n, prob){
    rbinom(n = n, size = 1, prob = prob)
    }

generate_straight_mc = function(results, n_spp, n_sites, n_gibbs){
    
    alpha = results$alpha
    beta = results$beta_compact
    gamma_tensor_sym = results$gamma_tensor_sym
    
    # pred = matrix(plogis(rep(1, n_sites) %*% t(alpha)), nrow = n_sites, ncol = n_spp)
    pred = matrix(rbern(n = n_sites * n_spp, plogis(rep(1, n_sites) %*% t(alpha))), n_sites, n_spp)

    for(i in 1:n_gibbs){
        cubic_contribution = matrix(0, n_sites, n_spp)

        for(i in 1:n_spp){ 
            build_product = array(0, dim=c(n_sites, n_spp, n_sites))
            build_product[ , i, ] = tensor(tensor(pred, results$gamma_tensor_sym[i, , ], 2, 1), pred, 2, 2) * 3#/ 2
            
            for (s in 1:n_sites){
                cubic_contribution[s, i] = build_product[s, i, s]
                }

            pred[,i] = rbern(nrow(pred),
                             plogis(
                                 cubic_contribution[, i] +
                                 pred %*% results$beta[ , i] +
                                 results$alpha[i]# +
                                 # env %*% alpha_env[ ,j] +
                             )
                            )
        }
    }

    pred = pred %>%
        as_tibble()
    #     rename_with(~colnames(x))

    # pred=pred>0

    return(pred)
}


apply_conditionals = function(unique_pred, count, sp, results){
        
    unique_pred = unique_pred %>% as.matrix()
    colnames(unique_pred) = NULL
    n_spp = ncol(unique_pred)
    cubic_contribution = matrix(0, ncol(unique_pred))
    out_pred = t(matrix(unique_pred, n_spp, count))
    # out_pred = 
    
    cubic_contribution[sp] = tensor(tensor(results$gamma_tensor_sym[sp, , ], unique_pred, 1, 2), unique_pred, 1, 2) * 3#[1]

    out_pred[, sp] = rbern(count,
                         plogis(
                             cubic_contribution[sp] +
                             unique_pred %*% results$beta_sym[, sp] + #since beta is already symmetrized
                             # results$beta_fit[sp, ] %*% t(unique_pred)  + #if not symmetrized
                             results$alpha[sp]# +
                             # env %*% alpha_env[ ,sp] +
                        )
                         ) #%>% as.matrix()
    # }

    return(out_pred %>% as_tibble())
}


rbern = function(n, prob){
    rbinom(n = n, size = 1, prob = prob)
    }

generate_shortcut_mc = function(results, n_spp, n_sites, n_gibbs){

    pred = matrix(rbern(n = n_sites * n_spp, plogis(rep(1, n_sites) %*% t(results$alpha))), n_sites, n_spp) %>% as_tibble()

    for(i in 1:n_gibbs){
        for (sp in 1:n_spp){
        
            pred = pred %>%
                group_by_all() %>%
                summarize(count = n(), .groups = "drop") %>%
                group_by(across(-count)) %>%
                ungroup() %>%
                mutate(id = row_number()) %>%
                nest(data = colnames(pred)) %>%
                # mutate(result = map2(data, count, apply_conditionals(sp))) %>%
                mutate(result = map2(data, count, \(data, count) apply_conditionals(data, count, sp, results))) %>%
                dplyr::select(result) %>%
                unnest() %>%
                # as.matrix() %>%
                # as_tibble() %>%
                rename_with(~colnames(pred)) #%>%
                # as.matrix()
        }
    }

    # pred=pred>0

    return(pred %>% as_tibble())
}

append_errors = function(x_distr, x_target, x2s, min_species_for_comp, min_species_for_fit, order, frac, n_sites, n_gibbs, datat_type, extra_tag){
    print(frac)
    
    x2s = x_distr %>%
    # resample_bootstrap() %>%
    # as.data.frame() %>%
    extend_frac_empties(frac) %>%
    extract_coefs(order = order, min_species_for_fit = min_species_for_fit) %>%
    generate_shortcut_mc(n_spp=ncol(x_distr), n_sites=n_sites, n_gibbs=n_gibbs) %>%
    compare_preds(x_target, min_species_for_comp) %>%
    pluck('summarized') %>%
    mutate(min_species_for_comp = min_species_for_comp,
           order = order,
           frac = frac,
           min_species_for_fit = min_species_for_fit,
           n_sites = n_sites,
           n_gibbs = n_gibbs,
           data_type = data_type,
           extra_tag = extra_tag
          ) %>%
    bind_rows(x2s) %>%
    suppressMessages()
    
    return(x2s)
}

append_results = function(x_distr, x_target, results_collected,
                          data_type,
                          min_species_for_comp, min_species_for_fit,
                          order, frac, n_sites, n_gibbs){
    
    results = x_distr %>%
    # resample_bootstrap() %>%
    # as.data.frame() %>%
    extend_frac_empties(frac) %>%
    extract_coefs(order = order, min_species_for_fit = min_species_for_fit) %>%
    generate_shortcut_mc(n_spp=ncol(x_distr), n_sites=n_sites, n_gibbs=n_gibbs) %>%
    compare_preds(x_target, min_species_for_comp) %>%
    pluck('results') %>%
    mutate(min_species_for_comp = min_species_for_comp,
           order = order,
           frac = frac,
           n_sites = n_sites,
           n_gibbs = n_gibbs,
           data_type = data_type,
           date_time = Sys.time()
          ) %>%
    bind_rows(x2s) %>%
    suppressMessages()
    
    return(x2s)
}

get_coefs = function(splits, frac = 0.2, order = 2, min_species_for_fit = 1, symmetrize_betas = T){
    splits %>%
    as_tibble() %>%
    extend_frac_empties(frac) %>%
    extract_coefs(order, min_species_for_fit, symmetrize_betas)
}

tidyup_alphas = function(coefs){
    coefs %>%
        pluck('alpha') %>%
        # as_tibble() %>%
        as.data.frame.table() %>%
        mutate(Freq = as.numeric(Freq),
              value = Freq) %>%
        as_tibble() %>%
        # filter(value != 0) %>%
        drop_na(value)
        # filter(Freq != 0)
        # rename(value = V1) %>%
        # rowid_to_column()
}

never_together = -99

tidyup_betas = function(coefs, symmetrize = TRUE){

    if (symmetrize==TRUE){
        
        final = coefs %>%
            pluck('beta_compact') %>%
            as.data.frame.table() %>%
            as_tibble() %>%
            mutate(Freq = as.numeric(Freq),
                  value = Freq) %>%
            # filter(value != 0) %>%
            drop_na(value)

        } else {
    
        temp = coefs %>%
            pluck('beta_fit') %>%
            as.data.frame.table() %>%
            as_tibble() %>%
            mutate(Freq = as.numeric(Freq),
                  value = Freq,
                  transposed = FALSE) %>%
            # filter(value != 0) %>%
            drop_na(value)

        final = coefs %>%
            pluck('beta_fit') %>%
            t() %>%
            as.data.frame.table() %>%
            as_tibble() %>%
            mutate(Freq = as.numeric(Freq),
                  value = Freq,
                  transposed = TRUE) %>%
            # filter(value != 0) %>%
            drop_na(value) %>%
            bind_rows(temp)
        }
    
    return(final)
        
}

tidyup_gammas = function(coefs){
    
    gammas = coefs %>%
        pluck('gamma_tibble_compact') %>%
        drop_na(value)
    
    return(gammas)
    }

append_coefs = function(data, alphas, betas, gammas, frac, order, min_species_for_fit, symmetrize_betas = TRUE, data_type, tag){
    
    coefs = get_coefs(data, frac, order, min_species_for_fit, symmetrize_betas)

    alphas = tidyup_alphas(coefs) %>%
        mutate(min_species_for_fit = factor(min_species_for_fit),
               order = factor(order),
               frac = frac,
               data_type = data_type,
               tag = tag) %>%
        bind_rows(alphas)

    betas = tidyup_betas(coefs) %>%
        mutate(min_species_for_fit = factor(min_species_for_fit),
               order = factor(order),
               frac = frac,
               data_type = data_type,
              tag = tag) %>%
        bind_rows(betas)

    gammas = tidyup_gammas(coefs) %>%
        mutate(min_species_for_fit = factor(min_species_for_fit),
               order = factor(order),
               frac = frac,
               data_type = data_type,
              tag = tag) %>%
        bind_rows(gammas)

    return(list('alphas' = alphas, 'betas' = betas, 'gammas' = gammas))
}