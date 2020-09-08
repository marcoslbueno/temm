
# setwd("~/Google Drive/marcos/phd/papers/ecml2020/github/")

source("sd_search_functions.R")
source("sd_randmodels.R")

# Nb of simulations to be run, each with its own ground truth models
nb_seeds = 10

for (seed in 1:nb_seeds){
  # Random state
  # seed = 1
  
  # Nb temporal targets in the temporal model (DBN/Markov chain)
  ndynamic = 17
  # Nb of sequences to sample
  numseq = 10
  # Nb time points in each sequence
  seqlen = 10
  
  # Set to TRUE to enforce the temporal model to be a Markov chain
  # This will collapse the temporal targets into a single temporal variable 
  flag_markov_chain = F
  
  # Control similarity of ground-truth models
  # If set to TRUE, add_noise must be set to a number between 0 and 1
  flag_similar_dbns = F
  add_noise = 0.10
  
  # Nb of noisy binary variables (subgroups formed by random assignment of DBN sequences)
  num_noisy_vars = 5
  
  # Subgroup size threshold
  Min_Cov = 0.05
  # Significance level for determining exceptional/non-exceptional
  SigLevel = 0.05
  
  # Nb of random subgroups to generate for a Distribution of False Discoveries
  nboot = 100
  
  # Verbose mode
  verbose = T
  
  #######################################
  # Generate ground truth models and data
  
  rand.dbn = function(ndynamic, seqlen, numseq, debug = FALSE) {
    # Transitions and initial models
    
    prob.node = 0.1
    dbni_structure = rand.structure.bn(numvar = ndynamic, prob.node = prob.node)
    nodes = nodes(dbni_structure)
    nodes(dbni_structure) = paste0(nodes(dbni_structure),"_0")
    
    dbni = buildrand.bn(domsizes = rep(2, ndynamic), g = dbni_structure)
    dbnt = buildrand.dbn(domsizes = rep(2,ndynamic))
    
    if (debug) graphviz.plot(dbnt$fit)
    if (debug) graphviz.plot(dbni$graph)
    
    # Unroll the DBN
    graph = structure.unrolled.tbn(nodes = nodes, finals = seqlen, 
                                   dbn_init = dbni$graph, dbns_paired = list(dbnt$graph))
    if (debug) graphviz.plot(graph)
    
    fit = fit.unrolled.tbn(bn = graph, nodes = nodes, finals = seqlen, 
                           fit_dbn_init = dbni$fit, fits_paired =  list(dbnt$fit))
    
    dbn_unroll = list(graph = graph, fit = fit)
    
    list(dbn0 = dbni, dbnt = dbnt, dbn_unroll = dbn_unroll, nodes = nodes)
  }
  
  ind = ground_dbn = ground_data = c()
  if (flag_similar_dbns) {
    
    for (i in 1:2) {
      # i=2
      if (verbose) cat("Generated ground truth model", i, "with difference control of",add_noise,"\n")
      
      set.seed(i + 2*seed)
      # Random DBN
      if (i==1) {
        dbn_simul = rand.dbn(ndynamic, seqlen, numseq, debug = FALSE)
        ground_dbn[[length(ground_dbn) + 1]] = dbn_simul  
        # Sample data
        d_simul = dbn.simulate(f = dbn_simul$dbn_unroll$fit, nodes = dbn_simul$nodes, n = rep(seqlen,numseq))  
      } else {
        
        dbn_simul = ground_dbn[[1]]
        
        aa = fit.unrolled.tbn(bn = dbn_simul$dbn_unroll$graph, 
                              nodes = dbn_simul$nodes, 
                              finals = seqlen, 
                              fit_dbn_init = dbn_simul$dbn0$fit, 
                              fits_paired = list(dbn_simul$dbnt$fit), 
                              add_noise = add_noise)
        
        node = "N1_1"
        c(aa[[node]]$prob)
        c(dbn_simul$dbn_unroll$fit[[node]]$prob)
        
        d_simul = dbn.simulate(f = aa, nodes = dbn_simul$nodes, n = rep(seqlen,numseq))  
      }
      
      ground_data[[length(ground_data) + 1]] = sequences.fromdataset(d_simul)
      
      ind[[length(ind) + 1]] = list(val = i, amount = numseq)
    }
    
  } else {
    for (i in 1:2) { #i=1
      if (verbose) cat("Generated ground truth model", i, "\n")
      set.seed(i + 2*seed)
      # Random DBN
      dbn_simul = rand.dbn(ndynamic, seqlen, numseq, debug = FALSE)
      ground_dbn[[length(ground_dbn) + 1]] = dbn_simul
      
      class(dbn_simul$dbnt$fit$N1_t2)
      
      # Sample data
      d_simul = dbn.simulate(f = dbn_simul$dbn_unroll$fit, nodes = dbn_simul$nodes, n = rep(seqlen,numseq))
      ground_data[[length(ground_data) + 1]] = sequences.fromdataset(d_simul)
      
      ind[[length(ind) + 1]] = list(val = i, amount = numseq)
    }
  }
  
  d_simul_mc = list()
  col = colnames(ground_data[[1]][[1]])
  ground_data_mc = c()#descriptor = ground_data[[i]]
  for (i in 1:length(ground_data)) {
    # i=1
    dataset = ground_data[[i]]
    
    dataset_mc = c()
    for (j in 1:length(dataset) ) {
      # j=1
      
      seq = c()
      for (k in 1:nrow(dataset[[j]]) ) {
        
        row = paste0(paste0(col,"="), (dataset[[j]])[k,], collapse = ",")  
        seq = rbind(seq, row)
      }
      colnames(seq) = "N1"
      dataset_mc[[length(dataset_mc) + 1]] = seq
    }
    ground_data_mc[[length(ground_data_mc) + 1]] = dataset_mc
  }
  
  if (flag_markov_chain) {
    ground_data = ground_data_mc
  }
  
  set.seed(seed)
  id = 1; d_train = names = c()
  for (i in 1:length(ind)) {
    # i=2
    for (j in 1:ind[[i]]$amount) {
      # j=1
      val = ind[[i]]$val
      rows = ground_data[[val]][[j]]
      
      # Noise
      for (k in 1:num_noisy_vars) {
        R = sample(x = c("1","2"), size = 1, prob = c(0.5, 0.5))
        rows = cbind(R, rows)  
      }
      # R1 = sample(x = c("1","2"), size = 1, prob = c(0.5, 0.5))
      
      # Id + Descriptors
      rows = cbind(val, rows)
      rows = cbind(id, rows)
      
      # rownames(rows) = rep(id, nrow(rows))
      colnames(rows)[1] = "application"
      colnames(rows)[2] = "A1"
      colnames(rows)[3:(3+num_noisy_vars-1)] = paste0("R", 1:num_noisy_vars)
      names = c(names, id)
      
      d_train[[length(d_train) + 1]] =  rows
      id = id+1
    }
  }
  names(d_train) = names
  
  flag_write_datasets = F
  if (flag_write_datasets) {
    d_train_csv = do.call(rbind, d_train)
    colnames(d_train_csv)[1] = "seq ID"
    
    if (flag_similar_dbns) { similardbns_str =   add_noise 
    } else { similardbns_str =  flag_similar_dbns }
    write.csv(d_train_csv, row.names = F, file = paste0("simuldata_n=",ndynamic,"_nsequence=",numseq
                                                        ,"_combine=",flag_markov_chain,"_controlsimilarity=",similardbns_str , "_randomstate=",seed, ".csv") )
  }
  
  attrs = colnames(d_train[[1]])
  
  d_flat = do.call(rbind, d_train)
  d_flat = as.data.frame(d_flat)
  theoretical_factors = theo_factors(d_flat, attrs)
  
  attrs_static = c("A1",paste0("R", 1:num_noisy_vars))
  attrs_dynamic = colnames(ground_data[[1]][[1]])# paste0("N", 1:ndynamic) 
  
  linetoget = 1
  l = (lapply(d_train, function(x) { x[linetoget, ] }))
  lbinded = do.call(rbind, l)
  rownames(lbinded) = lbinded[,"application"]
  lbinded = as.data.frame(lbinded)
  
  # lbinded will *not* be used as factors
  for (i in 1:ncol(lbinded)) {
    if (!is.factor(lbinded[,i])) {
      lbinded[,i] = as.character(lbinded[,i])
    } else {
      lbinded[,i] = levels(lbinded[,i])[lbinded[,i]]
    }
  }
  
  ##################
  # Search algorithm
  attrs_dynamic_0 = paste0(attrs_dynamic,"_0")
  attrs_dynamic_t1 = paste0(attrs_dynamic,"_t1")
  attrs_dynamic_t2 = paste0(attrs_dynamic,"_t2")
  
  # First call: no subgroup generated so far
  list_ToExpand = list()
  list_Except = list()
  list_NoExcept = list()
  list_NewNodes = list()
  list_sizes = list()
  query_list_sizes = 0
  
  node = list(attrs_static_cur = attrs_static, descriptor = c(), descriptor_flat = "root",
              subgroup_lbinded = lbinded, r = NULL, boot_rs = NULL)
  list_ToExpand[[1]] = node
  
  # Uncomment to use parallel processing with pkg doParallel, itertools, foreach
  # nb_cores = 1
  # registerDoParallel(nb_cores)
  
  while (length(list_ToExpand) > 0) {
    node = list_ToExpand[[1]]
    list_ToExpand = list_ToExpand[-1]
    
    # list_ToExpand = NULL
    if (flag_markov_chain) { Dmc = paste0("-markovchain")} else { Dmc = "" }
    if (flag_similar_dbns) { Dn = paste0("-similar",add_noise)} else { Dn = "" }
    Dx = paste0("-ndynamic",ndynamic)
    
    Dname = paste0("dsimul",seed,Dn,Dmc,Dx,"-numseq",numseq)
    result = expand(node = node, d_train = d_train, lbinded = lbinded, attrs_dynamic = attrs_dynamic, 
                    theoretical_factors = theoretical_factors, ls_sizes = list_sizes, query_ls_sizes = query_list_sizes, debug=verbose)
    
    list_NewNodes = result$list_NewNodes
    list_sizes = result$ls_sizes
    query_list_sizes = result$query_ls_sizes
    
    # if (verbose) cat("* Size of list os sizes =",length(list_sizes))
    
    if (length(list_NewNodes) == 0) { next }
    
    # Determine if new nodes are exceptional or non-exceptional
    except = noexcept = c()
    for (i in 1:length(list_NewNodes)) {
      # i=1
      
      node = list_NewNodes[[i]]
      
      boot_diffs = sapply(node$boot_rs, function(i) { i$bootbic_diff})
      if (!all(is.finite(boot_diffs))) { 
        warning("There are ", length(which(is.finite(boot_diffs) == FALSE)), 
                " random subgroup distances that not finite. Removing such distances...\n" ) 
        boot_diffs = boot_diffs[which(is.finite(boot_diffs) == TRUE)]
      }
      diff = node$r$bic_diff
      
      if ( (!is.finite(diff)) |  (!any(is.finite(boot_diffs))) ) {
        warning("Subgroup ",node$descriptor_flat,": Either subgroup distance is not finite Or all random distances are not finite. Subgroup will be considered as non-exceptional.\n")
        
        list_NewNodes[[i]]$pval = -Inf
        list_NewNodes[[i]]$zscore = -Inf
        
        noexcept[[length(noexcept)+1]] = c(i, -Inf, -Inf)
        
        next
      }
      
      # variable 'den' acts as a small epsilon to prevent division by zero
      if (length(unique(boot_diffs)) == 1) { den = 0.001 } else { den = sd(boot_diffs) }
      zscore_Subgroup = (diff - mean(boot_diffs)) / den
      pval_Subgroup = 2*pnorm(-abs(zscore_Subgroup))
      
      list_NewNodes[[i]]$pval = pval_Subgroup
      list_NewNodes[[i]]$zscore = zscore_Subgroup
      
      if (verbose) cat("Subgroup",node$descriptor_flat,"\n")
      
      if (pval_Subgroup < SigLevel) { except[[length(except)+1]] = c(i, zscore_Subgroup, pval_Subgroup) 
      } else { noexcept[[length(noexcept)+1]] = c(i, zscore_Subgroup, pval_Subgroup) } 
      
      if (verbose) cat("z score Subgroup =",zscore_Subgroup,"; p-value Subgroup =",pval_Subgroup,"\n\n")
    }
    
    # Update lists of exceptionals and non-exceptionals
    if (length(except) > 0) {
      j = sapply(except, function(x) {x[1]})
      list_Except = c(list_Except, list_NewNodes[j])
    }
    if (length(noexcept) > 0) {
      j = sapply(noexcept, function(x) {x[1]})
      list_NoExcept = c(list_NoExcept, list_NewNodes[j])
    }
    
    # Avoid expanding subgroups already expanded (with a different order of descriptors)
    j = sapply(except, function(x) {x[1]})
    if (length(j) > 0 ) {
      new_Except = list_NewNodes[j]  
      
      strings_Except = c()
      for (i in 1:length(new_Except)){ # i = 1
        t = new_Except[[i]]$descriptor
        desc_Except = matrix(t[order(t[,1],t[,2]),], ncol = 2)
        string = paste0(apply(desc_Except,1, function(x) { paste0(x,collapse = "=") }), collapse = ";")
        strings_Except = c(strings_Except, string)
      }
      strings_ToExpand =  c()
      if (length(list_ToExpand) > 0){
        for (j in 1:length(list_ToExpand)) {
          t = list_ToExpand[[j]]$descriptor
          desc_ToExpand = matrix(t[order(t[,1],t[,2]),], ncol = 2)
          string = paste0(apply(desc_ToExpand,1, function(x) { paste0(x,collapse = "=") }), collapse = ";")
          strings_ToExpand = c(strings_ToExpand, string)
        }
      }
      if (length(strings_Except) > 0){
        toAdd = new_Except[!(strings_Except %in% strings_ToExpand)]
        list_ToExpand = c(list_ToExpand,toAdd)
      }
    }
  }
  
  write_subgroups = function() {
    save(list = "list_Except", file = paste0(Dname,"-list_Except.RData"))
    save(list = "list_NoExcept", file = paste0(Dname,"-list_NoExcept.RData"))
    save(list = "list_NewNodes", file = paste0(Dname,"-list_NewNodes.RData"))
    save(list = "list_ToExpand", file = paste0(Dname,"-list_ToExpand.RData"))
    save(list = "list_sizes", file = paste0(Dname,"-list_sizes.RData"))
    save(list = "query_list_sizes", file = paste0(Dname,"-query_list_sizes.RData"))  
  }
  
  # Set to TRUE to write the exceptional/non-exceptional to disk
  # Required for post-processing classification measures
  flag_write_subgroups = T
  if (flag_write_subgroups) write_subgroups()
  
  print_results = function (nb_subgroup_show) {
    cat("\nFinal results\n")
    cat("#################################\n\n")
    
    cat("Total exceptional subgroups found:",length(list_Except),"\n")
    cat("Total non-exceptional subgroups found:",length(list_NoExcept),"\n\n")
    
    cat("Selection of Exceptional subgroups found:\n")
    cat("#################################\n")
    l1 = min(nb_subgroup_show, length(list_Except))
    for (i in 1:l1) {
      cat(i,"\n")
      cat("Subgroup Description:",list_Except[[i]]$descriptor_flat,"\n")
      cat("Size:",list_Except[[i]]$r$size,"\n")
      cat("Mismatch distance:",list_Except[[i]]$r$bic_diff,"\n")
      cat("p-value:",list_Except[[i]]$pval,"\n\n")
    }
    
    cat("Selection of Non-exceptional subgroups found:\n")
    cat("#################################\n")
    l2 = min(nb_subgroup_show, length(list_NoExcept))
    for (i in 1:l2) {
      cat(i,"\n")
      cat("Subgroup Description:",list_NoExcept[[i]]$descriptor_flat,"\n")
      cat("Size:",list_NoExcept[[i]]$r$size,"\n")
      cat("Mismatch distance:",list_NoExcept[[i]]$r$bic_diff,"\n")
      cat("p-value:",list_NoExcept[[i]]$pval,"\n\n")
    }
  }
  
  if (verbose) {
    cat("Simulation settings:\n\n")
    cat("Random state:",seed,"\n\n")
    
    cat("DESCRIPTOR VARIABLES\n")
    cat("Number of true variables: 1 (binary)\n")
    cat("Number of noisy variables:",num_noisy_vars,"\n\n")
    
    cat("TARGET VARIABLES\n")
    cat("Number of temporal targets in the temporal model:",ndynamic,"\n")
    
    cat("\nDATA\n")
    cat("Number of sequences sampled (per subgroup):", numseq,"\n")
    cat("Number of data points (time points) in each sequence:", seqlen,"\n")
    cat("Temporal targets collapsed, i.e. Markov chain representation?",flag_markov_chain,"\n")
    cat("Controlling the difference in ground truth models?",flag_similar_dbns,"\n")
    if (flag_similar_dbns) { cat("Maximal probability difference (delta):",add_noise,"\n") }
    
    cat("\nSEARCH PARAMETERS\n")
    cat("Subgroup size thresold:",Min_Cov,"\n")
    cat("Significance level for exceptionality test:",SigLevel,"\n")
  }
  
  # Nb. of subgroups to be shown
  nb_subgroup_show = 10 #length(list_Except)
  print_results(nb_subgroup_show)  
}
