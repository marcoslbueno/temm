library(bnlearn)
source("sd_estimation_dbn.R")

theo_factors = function (data, vars, debug = F) {
  theoretical_factors = list()
  for (i in 1:length(vars)) {
    # i=4
    j = which(colnames(data) == vars[i])
    l = unique(data[[j]])
    if (is.factor(l)) { l = levels(l)[l] }
    theoretical_factors[[i]] = factor(l)
  }
  names(theoretical_factors) = vars
  theoretical_factors
}

selectdata = function (data, va_select) {
  data = do.call(rbind, data)
  data = data[,va_select,drop=FALSE]
  data
}

remove.factors = function(t, attrs_dynamic) {
  for (j in 1:length(t)) {
    # j=1
    case = t[[j]]
    for (i in 1:ncol(case)) {
      # i = 1
      if (is.factor(case[,i]) ) { case[,i] = levels(case[,i])[case[,i]] 
      } else { case[,i] = as.character(case[,i]) }
    }
    t[[j]] = case[,attrs_dynamic,drop=FALSE]
  } 
  t
}

subgroup_dbn = function (appls, t, t_compl, new_descriptor, size, theoretical_factors, attrs_dynamic, debug = F) {
    set.seed(0)
    dbn_s = learn.dbn.bpic(datalist = t, theoretical_factors = theoretical_factors,
                           va_select = attrs_dynamic, text = "sd-searchgreedy", F); if (debug) cat("ok\n")
    
    set.seed(0)
    dbn_compl = learn.dbn.bpic(datalist = t_compl, theoretical_factors = theoretical_factors,
                               va_select = attrs_dynamic, text = "sd-searchgreedy", F)    ; if (debug) cat("ok\n")
    
    ll_datasub_sub = logLik.dbnv2(dt = t, dtpaired = dbn_s$train_paired[[1]], dbn = dbn_s, F)
    ll_datasub_compl = logLik.dbnv2(dt = t, dtpaired = dbn_s$train_paired[[1]],  dbn = dbn_compl, F); if (debug) cat("ok\n")
    
    ll_datacompl_compl = logLik.dbnv2(dt = t_compl, dtpaired = dbn_compl$train_paired[[1]], dbn = dbn_compl, F)
    ll_datacompl_sub = logLik.dbnv2(dt = t_compl, dtpaired = dbn_compl$train_paired[[1]], dbn = dbn_s, F); if (debug) cat("ok\n")
    
    datasize = sum(sapply(t, nrow))
    bic_datasub_sub = bic.dbn(dbn = dbn_s, ll = ll_datasub_sub, datasize = datasize)
    bic_datasub_compl = bic.dbn(dbn = dbn_compl, ll = ll_datasub_compl, datasize = datasize)
    
    datasize_compl = sum(sapply(t_compl, nrow))
    bic_datacompl_compl = bic.dbn(dbn = dbn_compl, ll = ll_datacompl_compl, datasize = datasize_compl)
    bic_datacompl_sub = bic.dbn(dbn = dbn_s, ll = ll_datacompl_sub, datasize = datasize_compl)
    
    llright = ll_datasub_sub + ll_datacompl_compl
    llwrong = ll_datasub_compl + ll_datacompl_sub
    
    bicright = bic_datasub_sub + bic_datacompl_compl
    bicwrong = bic_datasub_compl + bic_datacompl_sub
    
    diff = llright - llwrong
    bic_diff = bicwrong - bicright
    
    r = list(ll_datasub_sub = ll_datasub_sub, ll_datasub_compl = ll_datasub_compl, 
             ll_datacompl_compl = ll_datacompl_compl, ll_datacompl_sub = ll_datacompl_sub, 
             bic_datasub_sub = bic_datasub_sub, bic_datasub_compl = bic_datasub_compl, 
             bic_datacompl_compl = bic_datacompl_compl, bic_datacompl_sub = bic_datacompl_sub, 
             diff = diff, bic_diff = bic_diff,
             appls = appls, descriptor = new_descriptor, 
             size =  size, dbn_s = dbn_s, dbn_compl = dbn_compl)
  }
  
  subgroup_sampling = function (appls, appls_compl, t, new_descriptor, new_descriptor_flat, 
                                theoretical_factors, attrs_dynamic, size, debug) {
    boot_rs = list()
    for (boot in 1:nboot) {
    # boot = 1
    # You can use the line below to replace the "for" loop to enable parallel processing
    #   with the doParallel pkg
    # boot_rs = foreach(boot = 1:nboot) %do% {
      set.seed(boot)
      
      # d_previous = d_train[node$subgroup_lbinded[,"application"]]
      if (any(appls %in% appls_compl)) { 
        warning("Attention: there is common sequence data between the two subgroups to be compared!") 
      }
      d_previous = d_train[c(appls, appls_compl)]
      
      sample = sample(x = 1:length(d_previous), size = length(d_previous), replace = F)
      Ns = as.numeric(sapply(d_previous,nrow))
      goal = sum(sapply(t, nrow))
      cur = 0
      subsample = c()
      for (k in 1:length(d_previous)) {
        # k = 1
        if (cur >= goal) { break }
        cur = cur + Ns[sample[k]]
        subsample = c(subsample, sample[k])
      }
      
      # boot_subgroup = node$subgroup_lbinded[subsample, ]
      # boot_appls = boot_subgroup[, "application"]
      # boot_t = d_previous[boot_appls]
      # boot_t_compl = d_previous[!(names(d_previous) %in% boot_appls)]
      
      boot_appls = c(appls, appls_compl)[subsample]
      boot_t = d_previous[boot_appls]
      boot_t_compl = d_previous[!(names(d_previous) %in% boot_appls)]
      
      boot_t = remove.factors(boot_t, attrs_dynamic)
      boot_t_compl = remove.factors(boot_t_compl, attrs_dynamic)
      #
      boot_dbn_s = learn.dbn.bpic(datalist = boot_t, theoretical_factors = theoretical_factors,
                                  va_select = attrs_dynamic, text = "sd-searchgreedy", F)#; cat("ok\n")
      
      boot_dbn_compl = learn.dbn.bpic(datalist = boot_t_compl, theoretical_factors = theoretical_factors,
                                      va_select = attrs_dynamic, text = "sd-searchgreedy", F)#    ; cat("ok\n")
      #
      boot_ll_datasub_sub = logLik.dbnv2(dt = boot_t, dtpaired = boot_dbn_s$train_paired[[1]], dbn = boot_dbn_s, F)
      boot_ll_datasub_compl = logLik.dbnv2(dt = boot_t, dtpaired = boot_dbn_s$train_paired[[1]],  dbn = boot_dbn_compl, F)#; cat("ok\n")
      
      boot_ll_datacompl_compl = logLik.dbnv2(dt = boot_t_compl, dtpaired = boot_dbn_compl$train_paired[[1]], dbn = boot_dbn_compl, F)
      boot_ll_datacompl_sub = logLik.dbnv2(dt = boot_t_compl, dtpaired = boot_dbn_compl$train_paired[[1]], dbn = boot_dbn_s, F)#; cat("ok\n")
      #
      datasize = sum(sapply(boot_t, nrow))
      bootbic_datasub_sub = bic.dbn(dbn = boot_dbn_s, ll = boot_ll_datasub_sub, datasize = datasize)
      bootbic_datasub_compl = bic.dbn(dbn = boot_dbn_compl, ll = boot_ll_datasub_compl, datasize = datasize)
      
      datasize_compl = sum(sapply(boot_t_compl, nrow))
      bootbic_datacompl_compl = bic.dbn(dbn = boot_dbn_compl, ll = boot_ll_datacompl_compl, datasize = datasize_compl)
      bootbic_datacompl_sub = bic.dbn(dbn = boot_dbn_s, ll = boot_ll_datacompl_sub, datasize = datasize_compl)
      # 
      bootbicright = bootbic_datasub_sub + bootbic_datacompl_compl
      bootbicwrong = bootbic_datasub_compl + bootbic_datacompl_sub
      # 
      bootbic_diff = bootbicwrong - bootbicright
      
      if (debug) cat("Subgroup",new_descriptor_flat,"; Sample",boot,"; BIC diff =",bootbic_diff,"; Bic right =",bootbic_datasub_sub,
          bootbic_datacompl_compl,"; Bic wrong =",bootbic_datasub_compl, bootbic_datacompl_sub,"\n")
      
      boot_r = list(bootbic_datasub_sub = bootbic_datasub_sub, bootbic_datasub_compl = bootbic_datasub_compl, 
                    bootbic_datacompl_compl = bootbic_datacompl_compl, bootbic_datacompl_sub = bootbic_datacompl_sub, 
                    bootbic_diff = bootbic_diff,
                    boot_appls = boot_appls, descriptor = new_descriptor, 
                    size =  size)
      
      boot_r
      boot_rs[[length(boot_rs)+1]] = boot_r
    }
    if (debug) cat("\n")
    boot_rs
  }
  
expand = function(node, d_train, lbinded,
                  attrs_dynamic, theoretical_factors, ls_sizes, query_ls_sizes, debug = FALSE) {
  
  if (length(node$attrs_static_cur) == 0) {
    return (list(list_NewNodes = NULL, ls_sizes = ls_sizes, query_ls_sizes = query_ls_sizes)) }
  
  Max_Cov = 1 - Min_Cov
  attrs_static_feasible = c()
  for (i in 1:length(node$attrs_static_cur)) {
    # i=1
    var_chosen = node$attrs_static_cur[i]
    dom = node$subgroup_lbinded[[var_chosen]]; dom = unique(dom)
    sizes = c(); c = sum(as.numeric(sapply(d_train,nrow)))
    
    if (length(dom) == 1) { next }
    
    feasible = F
    for (v in dom) {
      # v = dom[1]
      subgroup = node$subgroup_lbinded[node$subgroup_lbinded[[var_chosen]] == v,]
      
      a = subgroup[,"application"]
      b = sapply(d_train[a], nrow)
      size = sum(as.numeric(b)) / c
      
      if (size >= Min_Cov & size <= Max_Cov) { feasible = T; break } 
    }
    if (feasible) { attrs_static_feasible = c(attrs_static_feasible, var_chosen)}
  }
  
  if (length(attrs_static_feasible) == 0) {
    return (list(list_NewNodes = NULL, ls_sizes = ls_sizes, query_ls_sizes = query_ls_sizes)) }
  
  rs = boots_rs = c()
  list_NewNodes = c()
  for (i in 1:length(attrs_static_feasible)) {
    # i=1
    var_chosen = attrs_static_feasible[i]
    dom = node$subgroup_lbinded[[var_chosen]]; dom = unique(dom)
    sizes = c(); c = sum(as.numeric(sapply(d_train,nrow)))
    
    for (v in dom) {
      # v = dom[1]
      subgroup = node$subgroup_lbinded[node$subgroup_lbinded[[var_chosen]] == v,]
      
      a = subgroup[,"application"]
      b = sapply(d_train[a], nrow)
      size = sum(as.numeric(b)) / c
      
      sizes = c(sizes, size)
    }
    
    # Subgroup and Enough size
    aux = order(sizes, decreasing = T)
    sizes = sizes[aux]
    dom = dom[aux]
    
    for (j in 1:length(dom)) {
      # j=1
      if ((sizes[j] < Min_Cov) | (sizes[j] > 1-Min_Cov) ) { next }
      
      newNode = list()
      newNode$attrs_static_cur = attrs_static_feasible[which(!(attrs_static_feasible == var_chosen))]
      
      value_chosen = dom[j]
      subgroup = node$subgroup_lbinded[node$subgroup_lbinded[[var_chosen]] == value_chosen,]
      new_descriptor = rbind(node$descriptor, c(var_chosen,value_chosen))
      new_descriptor_flat = paste0(apply(new_descriptor, 1, paste0, collapse="="), collapse = ",")
      
      newNode$descriptor = new_descriptor
      newNode$descriptor_flat = new_descriptor_flat
      newNode$subgroup_lbinded = subgroup
      
      if (debug) cat("Subgroup:", new_descriptor_flat,"; Size =",sizes[j],"\n")
      
      appls = subgroup[,"application"]
      t = d_train[appls]
      appls_compl = lbinded$application[!(lbinded$application %in% appls)]
      t_compl = d_train[appls_compl]
      
      if (debug) cat("* Subgroup Length appls = ",sizes[j],"; Length appls_compl =",1-sizes[j],"\n")
      
      t = remove.factors(t, attrs_dynamic)
      t_compl = remove.factors(t_compl, attrs_dynamic)
      
      r = subgroup_dbn(appls = appls, t = t, t_compl = t_compl, size = sizes[j], new_descriptor = new_descriptor, 
                       theoretical_factors = theoretical_factors, attrs_dynamic = attrs_dynamic, debug)
      
      newNode$r = r
      
      if (debug) cat("Subgroup:","BIC diff =",newNode$r$bic_diff,"; Bic right =",newNode$r$bic_datasub_sub,
          newNode$r$bic_datacompl_compl,"; Bic wrong =",newNode$r$bic_datasub_compl, newNode$r$bic_datacompl_sub,"\n")

      sizes[j]
      if (length(ls_sizes) == 0) { 
        found = 0
      } else {
        if (is.null(ls_sizes[[as.character(sizes[j]) ]] )) { found = 0
        } else { 
          if (debug) cat("* Found subgroup size in the list of sizes.\n\n")
          found = 1 
          boot_rs = ls_sizes[[as.character(sizes[j]) ]]
          query_ls_sizes = query_ls_sizes + 1
        }
      }
      
      if (!found) {
        if (debug) cat("* Didn't find subgroup size in the list of sizes. Will calculate the DFD. \n")
        
        boot_rs = subgroup_sampling(appls = appls, appls_compl = appls_compl, t = t, new_descriptor = new_descriptor,
                                    new_descriptor_flat = new_descriptor_flat,  theoretical_factors = theoretical_factors,
                                    attrs_dynamic = attrs_dynamic, size = sizes[j], debug = debug)

        ls_sizes[[length(ls_sizes) + 1]] = boot_rs
        names(ls_sizes)[length(ls_sizes)] = as.character(sizes[j])
      }
      
      newNode$boot_rs = boot_rs
      
      list_NewNodes[[length(list_NewNodes) + 1]] = newNode
    }
    
  }
  # save(list = "rs", file = paste( node$descriptor_flat,Dname,"-rs.RData"))
  # save(list = "boots_rs", file = paste( node$descriptor_flat, Dname,"-boots-rs.RData"))
  return (list(list_NewNodes = list_NewNodes, ls_sizes = ls_sizes, query_ls_sizes = query_ls_sizes))
}
