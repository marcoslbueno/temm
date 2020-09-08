library(PRROC)

# Nb of simulations to be run, each with its own ground truth models
nb_seeds = 10

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

# Nb of noisy binary variables (subgroups formed by random assignment of temporal sequences)
num_noisy_vars = 5

# Significance level for determining exceptional/non-exceptional
SigLevel = 0.05

if (flag_markov_chain) { Dmc = paste0("-markovchain")} else { Dmc = "" }
if (flag_similar_dbns) { Dn = paste0("-similar",add_noise)} else { Dn = "" }
Dx = paste0("-ndynamic",ndynamic)

label_subgroups = function () {
  for (seed in 1:nb_seeds) {
    #seed=1
    print(seed)
    
    Dname = paste0("dsimul",seed,Dn,Dmc,Dx,"-numseq",numseq)
    # Dname = "bpic18"
    print(Dname)
    load(file = paste0(Dname,"-list_Except.RData"), verbose = T)
    load(file = paste0(Dname,"-list_NoExcept.RData"), verbose = T)
    
    # Return indexes of NON-redundant subgroups
    redundant_descriptors = function (list_Subgroups, debug = FALSE) {
      if (length(list_Subgroups) == 0) { return(c())}
      strings = c()
      for (i in 1:length(list_Subgroups)) {
        # i = 1
        t = list_Subgroups[[i]]$descriptor
        desc_Except = matrix(t[order(t[,1],t[,2]),], ncol = 2)
        string = paste0(apply(desc_Except,1, function(x) { paste0(x,collapse = "=") }), collapse = ";")
        strings = c(strings, string)  
      }
      
      i_NoRedund = c()
      for (i in 1:(length(strings) - 1) ) {
        # i=1
        a = strings[i] == strings[(i+1):length(strings)] #strings[i] %in% strings[(i+1):length(strings)]
        if (any(a)>0) {
          if (debug) { 
            print(i)
            cat("Subgroup",strings[i],"is redundant in",which(a)+i,"\n")
          }
        } else {
          i_NoRedund = c(i_NoRedund, i)
        }
      }
      return(i_NoRedund)
    }
    
    noRedund_Except = redundant_descriptors(list_Subgroups = list_Except)
    noRedund_NoExcept = redundant_descriptors(list_Subgroups = list_NoExcept)
    
    list_Except = list_Except[noRedund_Except]
    list_NoExcept = list_NoExcept[noRedund_NoExcept]
    
    rs = boots_rs = diffs = rows = rows_Except = rows_NoExcept = c()
    for (i in 1:length(list_Except)){
      # i = 70
      cat("Collecting exceptional subgroup",i,"\n")
      node = list_Except[[i]]
      
      row = c(node$descriptor_flat, nrow(node$descriptor), node$r$size, node$pval, node$zscore)
      rows_Except = rbind(rows_Except, row) 
    }
    
    for (i in 1:length(list_NoExcept)){
      cat("Collecting non-exceptional subgroup",i,"\n")
      node = list_NoExcept[[i]]
      
      row = c(node$descriptor_flat, nrow(node$descriptor), node$r$size, node$pval, node$zscore)
      rows_NoExcept = rbind(rows_NoExcept, row)  
    }
    
    colnames = c("Subgroup", "Nb Descriptors", "Size", "p-value", "z-score")
    colnames(rows_Except) = colnames(rows_NoExcept) = colnames
    rows = rbind(rows_Except, rows_NoExcept)
    
    aux = rbind (c("Exceptional subgroups", rep(x = "",length(colnames)-1) ), rows_Except)
    aux = rbind(aux, c("Non-Exceptional subgroups", rep(x = "",length(colnames)-1) ), rows_NoExcept)
    write.table(x = aux, file = paste0(Dname,"-pvaluesubgroups.csv"), quote = F, sep = ";", dec = ",", row.names = F, col.names = T)
    
    # Labeling subgroups
    tab = as.data.frame(rows[,1] )
    tab[,2] = 1 - as.numeric(rows[,4])
    colnames(tab) = c("subgroup", "prob.exceptional")
    
    # Ground truth labels:
    # 0 = false unitary subgroup
    # 1 = true unitary subgroup 
    # 2 = false specialized subgroup
    # 3 = true specialized subgroup
    gtruths = c()
    for (i in 1:nrow(tab)) { # i=100
      long_descriptor = regexpr(pattern = ",", text = tab[i,1], fixed = T)[1]
      if (long_descriptor > -1) { # the subgroup is specialized...
        contains_true = regexpr(pattern = "A", text = tab[i,1], fixed = T)[1]
        if (contains_true > -1) { # true + noisy
          gtruth = 3
        } else { 
          gtruth = 2
        }
      } else { # the subgroup is unitary...
        contains_true = regexpr(pattern = "A", text = tab[i,1], fixed = T)[1]
        if (contains_true > -1) { 
          gtruth = 1
        } else {
          gtruth = 0
        }
      }
      gtruths = c(gtruths, gtruth)
    }
    tab = cbind(tab, gtruths)
    
    save(list = "tab", file = paste0(Dname,"-labels.RData"))
  }
}

classification_measures = function() {
  rows = num_subgroups = aucs = rows_specialized = c()
  for (seed in 1:nb_seeds) {
    # seed = 1
    cat("Collecting labeling results, random state",seed,"\n")
    
    Dname = paste0("dsimul",seed,Dn,Dmc,Dx,"-numseq",numseq)
    load(file = paste0(Dname,"-labels.RData"), verbose = T)
    
    tab = tab[apply(tab, 1, function(i) { is.finite(as.numeric(i[2])) }), ]
    tab_specialized = tab[!( (tab[,3] == 0) | (tab[,3] == 1)) , ]
    tab = tab[(tab[,3] == 0 | tab[,3] == 1), ]
    
    if (nrow(tab) == 0) { next }
    
    # AUROC
    # plot(roc.curve(scores.class0 = tab$prob.exceptional, weights.class0 = tab$gtruths, curve = T))
    auc = (roc.curve(tab$prob.exceptional, weights.class0 = tab$gtruths))$auc
    if (!is.finite(auc)) { auc = 0 }
    aucs = c(aucs, auc)
    
    # Labeling the discovered subgroups (exceptional / non-exceptional) based on the search's threshold
    num_subgroups = c(num_subgroups, nrow(tab))
    pred = as.numeric(tab$prob.exceptional >= 1-SigLevel)
    
    pred_specialized = as.numeric(tab_specialized$prob.exceptional >= 1-SigLevel)
    tab_specialized[,1] = levels(tab_specialized[,1])[tab_specialized[,1]]
    tab_specialized = cbind(tab_specialized[,1:2], pred_specialized, tab_specialized[,3])
    colnames(tab_specialized)[4] = "gtruths"
    
    tab[,1] = levels(tab[,1])[tab[,1]]
    tab = cbind(tab[,1:2], pred, tab[,3])
    colnames(tab)[4] = "gtruths"
    
    # Confusion matrix measures
    pos = tab[tab$subgroup == "A1=1", ]
    pos = rbind(pos, tab[tab$subgroup == "A1=2", ])
    tp = pos[pos$pred == 1,]
    tp = nrow(tp)
    fn = 2 - tp
    
    neg = c()
    for (i in 1:num_noisy_vars) { 
      # i = 1
      str1 = paste0(c("R",i,"=","1"), collapse = "")
      str2 = paste0(c("R",i,"=","2"), collapse = "")
      neg = rbind(neg, tab[tab$subgroup==str1, ])
      neg = rbind(neg, tab[tab$subgroup==str2, ])
    }
    fp =  neg[neg$pred == 1, ]
    fp = nrow(fp)
    tn = num_noisy_vars*2 - fp
    
    recall = tp / (tp  + fn)
    if (!is.finite(recall)) { recall = 0 }
    
    prec = tp / (tp + fp)
    if (!is.finite(prec)) { prec = 0 }
    
    x = tp+ fn + fp + tn
    row = c(seed, round(auc,2), round(prec,2), round(recall,2), 
            tp, fn, fp, tn, x, nrow(tab_specialized))
    rows = rbind(rows, row)
    
    spec_a1 = tab_specialized[tab_specialized$gtruths == 3, ]
    if (nrow(spec_a1) > 0) {
      perc_a1correct = nrow(spec_a1[spec_a1$pred_specialized == 1,])/nrow(spec_a1)
    } else {
      perc_a1correct = NA
    }
    
    row_specialized = c(nrow(spec_a1), perc_a1correct)
    
    spec_noisy = tab_specialized[tab_specialized$gtruths == 2, ]
    if (nrow(spec_noisy) > 0) {
      perc_noisycorrect = nrow(spec_noisy[spec_noisy$pred_specialized == 0,])/nrow(spec_noisy)
    } else {
      perc_noisycorrect = NA
    }
    
    row_specialized = c(row_specialized, c(nrow(spec_noisy), perc_noisycorrect))
    rows_specialized = rbind(rows_specialized, row_specialized)
    
    # data = sort(1-spec_a1$prob.exceptional, decreasing = T)
    # boxplot(data, notch = T)
    # plot(data, type = "l",ylim = c(0,1), xlim = c(0,260))
    # stripchart(x = data, method =  "jitter", pch = 1, vertical = T)
    
  }
  colnames(rows) = c("Iteration", "AUC-ROC",  "Precision","Recall", "TP","FN","FP","TN","Total Unitary Subgroups", "Total Specialized Subgroups")
  iter = nrow(rows)
  
  rows = rbind(rows, apply(rows, 2, function(i) { mean(as.numeric( i[1:iter] )) } )); 
  rows = rbind(rows, apply(rows, 2, function(i) { sd(as.numeric( i[1:iter] )) } )); 
  
  rows[nrow(rows)-1, 1] = "mean"
  rows[nrow(rows), 1] = "sd"
  
  rows_specialized = cbind(1:nrow(rows_specialized), rows_specialized)
  colnames(rows_specialized) = c("Iteration", "#Subgroups A1+noisy", "Labeled Excep", "#Subgroups noisy", "Labeled Non-Excep")
  
  nona1 = sapply(1:nrow(rows_specialized), function(i) { !is.na(rows_specialized[i,3]) })
  nona2 = sapply(1:nrow(rows_specialized), function(i) { !is.na(rows_specialized[i,5]) })
  meanrow = c(0, mean(rows_specialized[nona1,2]), mean(rows_specialized[nona1,3]), 
              mean(rows_specialized[nona2,4]), mean(rows_specialized[nona2,5]))
  sdrow = c(0, sd(rows_specialized[nona1,2]), sd(rows_specialized[nona1,3]), 
            sd(rows_specialized[nona2,4]), sd(rows_specialized[nona2,5]))
  rows_specialized = rbind(rows_specialized, meanrow ) 
  rows_specialized = rbind(rows_specialized, sdrow ) 
  
  rows_specialized[nrow(rows_specialized)-1, 1] = "mean"
  rows_specialized[nrow(rows_specialized), 1] = "sd"
  
  write.table(rows, file = paste0("results-dsimul-unitary",Dn,Dmc,Dx,"-numseq",numseq,".csv"), sep = " ; ", row.names = F)
  write.table(rows_specialized, file = paste0("results-dsimul-specialized",Dn,Dmc,Dx,"-numseq",numseq,".csv"), sep = " ; ", row.names = F)
}

label_subgroups()

classification_measures()
