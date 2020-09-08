library(bnlearn)

#->DBN: Initial BN
dbn.initial <- function (sequences, doms, score = "aic", debug = F) {
  ncol = length((sequences[[1]] )[1,])
  nrow = length(sequences)
  dimnames = dimnames(sequences[[1]])
  train_initial = matrix(ncol = ncol, nrow = nrow)
  dimnames(train_initial) = dimnames

  j = 1
  for(seq in sequences) { 
    #seq = sequences[[1]]
    if (is.vector(seq)) {
      row = seq
    } else {
      row = seq[1,]
    }
    
    train_initial[j,] = row
    j = j+1
  }
  
  dimnames(train_initial)[[2]] = paste0(dimnames(train_initial)[[2]], "_0")
  
  if (!is.null(doms) & nrow(train_initial) > 1) {
    train_initial = convertfactors(data = train_initial, theoretical_factors = doms, debug = debug)
  } 
    
  if (nrow(train_initial) > 1) {
    # do not convert again, else we lose the original factors on each column of 'train_initial'
    if (!is.data.frame(train_initial)) {
      train_initial = as.data.frame(train_initial)  
    }
    
    dbn_init = tabu(train_initial, score = score)
    fit_dbn_init = bn.fit(dbn_init, train_initial, iss = 1, method = "bayes")  
  } else {
    dbn_init = empty.graph(nodes = dimnames[[2]])    
    
    cpts = vector(mode = "list", length = ncol(train_initial))
    for (i in 1:length(cpts)) { # i = 1
      pr = rep(0, length(doms[[i]]))
      pr[which(doms[[i]] == train_initial[1,i])] = 1
      
      cpts[[i]] = matrix(pr, ncol = length(pr), dimnames = list(NULL, doms[[i]]))
    }
    names(cpts) = dimnames[[2]]
    
    fit_dbn_init = custom.fit(dbn_init, dist = cpts)      
  }
  
  list (net = dbn_init, fit = fit_dbn_init)
}

#->DBN: Transitions, given by one or more conditional BNs
dbn.blacklist <- function (nodes1, nodes2) {
  bl = NULL
  for (i in 1:length(nodes1)) { #i=1
    #arcs between nodes of time t
    if (length(nodes1) > 1) {
      rows1 = tiers2blacklist(list (nodes1[i], nodes1[-i]))
      bl = rbind (bl, rows1)
    }
    
    #arcs from (t+1) to (t)
    rows2 = tiers2blacklist(list (nodes1, nodes2[i]))
    bl = rbind (bl, rows2)
  }
  bl
}

#->DBN: 3-Simulating data 
structure.unrolled.tbn <- function (nodes, finals, dbn_init, dbns_paired) {
  allnames = c(sapply(0:(finals[length(finals)] - 1), function(x) { paste0(nodes,"_",x)} ))
  bn = empty.graph(allnames)
  
  #getting arcs
  arcs = dbn_init$arcs
  if (length(arcs) != 0) {
    for (i in 1:length(arcs[,1])) { #i=1
      arc = arcs[i,]
      bn = set.arc(bn, from = arc[1], to = arc[2])
    }
  }
  
  begin = 1
  end = finals[1]
  for (i in 1:length(dbns_paired)) { #i=1
    bn = arcs.dbn (bn, dbns_paired[[i]], begin, end)
    begin = end
    end = finals[i+1]
  }
  
  bn
}


arcs.dbn <- function (bn, dbn, begin, end) {
  arcs = dbn$arcs
  if (length(arcs) != 0) {
    for (i in 1:length(arcs[,1])) { #i=6
      for (j in begin:(end-1)) { #j=1
        arc = arcs[i, ]
        arc = gsub(x=arc, pattern = "_t2", replacement = paste0("_",j), fixed = T)
        arc = gsub(x=arc, pattern = "_t1", replacement = paste0("_",j-1), fixed = T)
        bn = set.arc(bn, from = arc[1], to = arc[2])
      }
    }
  }
  
  bn
}

fit.unrolled.tbn <- function (bn, nodes, finals, fit_dbn_init, fits_paired, add_noise = NULL) {
  numpoints = finals[length(finals)]
  fit = list()#vector(mode = "list", length = (numpoints*length(nodes)) )
  
  for (i in 1:length(fit_dbn_init)) { # i=8 10
    n = fit_dbn_init[[i]]
    
    #no parents
    if (length(n$parents) == 0) {
      data = n$prob
      
      #This works only for binary variables!
      if (!is.null(add_noise)) {
        maxdiff = abs(1-data[1])
        maxdiff = min(add_noise, maxdiff)
        change = runif(n=1, min=0, max=maxdiff)
        data[1] = data[1] + change
        data[2] = data[2] - change
      }
      
      m = matrix(data, ncol = length(data), dimnames = list(NULL, dimnames(data)[[1]]))
      dimnames(m)
    } else {
      data = c(n$prob)
      
      #This works only for binary variables!
      if (!is.null(add_noise)) {
        for (j in seq(from=1, to = length(data), by = 2) ) {
          # j = 3
            maxdiff = abs(1-data[j])
            maxdiff = min(add_noise, maxdiff)
            change = runif(n=1, min=0, max=maxdiff)
            data[j] = data[j] + change
            data[j+1] = data[j+1] - change
        }
      }
      #Uniform replacement!
      data = as.numeric(gsub(pattern = NaN, replacement = 0.5, x = data))
      
      m = matrix(data, ncol = 1)
      dim(m) = dim(n$prob)
      dimnames(m) = dimnames(n$prob) 
    }
    
    fit[[i]] = m
  }
  
  begin = 1
  end = finals[1]
  for (i in 1:length(fits_paired)) { #i=1
    fit = cpt.tbn (fit, fits_paired[[i]], nodes, begin, end, add_noise)
    begin = end
    end = finals[i+1]
  }
  
  allnames = paste0(nodes, "_0")
  
  begin = 1
  end = finals[1]
  aux = c()
  for (i in 1:length(fits_paired)) {    #i=1
    aux = c(sapply (nodes, function(x) { paste0(x, "_", begin:(end-1)) } ))
    allnames = c(allnames, aux)
    
    begin = end
    end = finals[i+1]
  }
  
  names(fit) = allnames
  
  f = custom.fit(bn, dist = fit)
  f
}

cpt.tbn <- function (fit, fit_dbn, nodes, begin, end, add_noise = NULL) {
  
  base = length(nodes)
  
  for (i in 1:length(nodes)) { #i=1
    for (j in begin:(end-1)) { #j=1
      n = fit_dbn[[base + i]]
      
      #no parents
      if (length(n$parents) == 0) {
        data = n$prob
      
        if (!is.null(add_noise)){
          maxdiff = abs(1-data[1])
          maxdiff = min(add_noise, maxdiff)
          change = runif(n=1, min=0, max=maxdiff)
          data[1] = data[1] + change
          data[2] = data[2] - change
        }
        
        m = matrix(data, ncol = length(data), dimnames = list(NULL, dimnames(data)[[1]]))
      } else {
        #renaming
        names(dimnames(n$prob)) = gsub (x=names(dimnames(n$prob)), pattern = "_t1", 
                                        replacement = paste0("_",j-1), fixed = T)
        names(dimnames(n$prob)) = gsub (x=names(dimnames(n$prob)), pattern = "_t2", 
                                        replacement = paste0("_",j), fixed = T)
        
        data = c(n$prob)
        #Uniform replacement!
        data = as.numeric(gsub(pattern = NaN, replacement = 0.5, x = data))
        
        #! Works only for binary variables!
        if (!is.null(add_noise)) {
          for (j in seq(from=1, to = length(data), by = 2) ) {
            # j = 3
            maxdiff = abs(1-data[j])
            maxdiff = min(add_noise, maxdiff)
            change = runif(n=1, min=0, max=maxdiff)
            data[j] = data[j] + change
            data[j+1] = data[j+1] - change
          }
        }
        
        m = matrix(data, ncol = 1)
        dim(m) = dim(n$prob)
        dimnames(m) = dimnames(n$prob) 
      }
      
      fit[[length(fit) + 1]] = m
    }
  }
  
  fit
}

dbn.simulate <- function (f, nodes, n) {
  nslices = length(f)/length(nodes)
  if (any(n > nslices)) { 
    warning("Error: requested sample is larger than the number of slices in the DBN.") 
    return(NULL)
  }
  train_simul = vector()
  shift = length(nodes) #sempre de 3 em 3
  
  for (k in 1:length(n)) { #k=1
    data = rbn(f, 1)
    
    seqsize = n[k]
    begin = 1
    for (i in 1:seqsize) { #i=1
      row = as.vector(t( data[begin:(begin + shift - 1)] ) )
      train_simul = rbind (train_simul, row)
      
      begin = begin+shift
    }
  }
  
  dimnames(train_simul) = list(NULL, nodes)
  
  train = list(x = train_simul, N = n)
  train
}

bic.dbn = function(dbn, ll, datasize) {
  nparams = nparams(dbn$dbn_init$fit) + nparams(dbn$fits_paired[[1]])
  bic = log(datasize) * nparams - 2*ll
  bic
}

learn.dbn.bpic = function (datalist, theoretical_factors, va_select, text = "", debug = F) {
  N = unname(sapply(datalist, nrow))
  d = selectdata(data = datalist, va_select = va_select)
  d = convertfactors(d, theoretical_factors = theoretical_factors[va_select], debug = debug)
  if (!is.null(va_select)) {
    d = d[,va_select,drop=F]
  }
  d = as.matrix(d) #make it from a matrix of factors to a matrix of chars 
  
  if (debug) { print(sum(N) == nrow(d)) }
  
  datamhsmm = list(x = d, N = N)
  rownames(datamhsmm$x) = NULL
  class(datamhsmm) = "mhsmm"
  
  #->DBN: initial distr
  train = datamhsmm
  score = "bic"
  doms = theoretical_factors[va_select]
  
  sequences = sequences.fromdataset (train)
  dbn_init = dbn.initial (sequences, doms, score = score, debug = debug)
  
  #->DBN: transitions
  nodes = dimnames(train$x)[[2]]
  nodes1 = paste0(nodes,"_t1")
  nodes2 = paste0(nodes, "_t2")
  
  bl = dbn.blacklist(nodes1, nodes2)
  
  train_paired = list()
  for (i in 1:length(sequences)) {
    # i=1
    seq = sequences[[i]]
    if (length(seq) == 1) { next }
    rows = matrix(cbind(seq[1:(nrow(seq)-1),], seq[2:nrow(seq),]), nrow = nrow(seq)-1)
    train_paired = rbind(train_paired, rows)
  }
  colnames(train_paired) = c(nodes1, nodes2)
  train_paired = as.data.frame(train_paired)
  train_paired = list (train_paired)
  
  train_paired[[1]] = convertfactors(data = train_paired[[1]], theoretical_factors = c(doms,doms), debug = debug)
  
  dbns_paired = tabu (train_paired[[1]], blacklist = bl, score = score)
  dbns_paired = list(dbns_paired)
  fits_paired = bn.fit(dbns_paired[[1]], train_paired[[1]], iss = 1, method = "bayes")
  fits_paired = list(fits_paired)
  
  bn = f = NULL # unfolded objects
  
  dbn = list(net = bn, fit = f, dbn_init = dbn_init, nets_paired = dbns_paired, fits_paired = fits_paired, train_paired = train_paired)
  dbn
}

logLik.dbnv2 = function (dt, dtpaired, dbn, debug = F) {
  t0 = sapply(dt, function(i) { (i[1, ]) }, simplify = F); t0 = (do.call(rbind, t0)); colnames(t0) = nodes(dbn$dbn_init$fit)
  t0 = convertfactors(data = t0, theoretical_factors = theoretical_factors[attrs_dynamic], debug = F)
  llinit = logLik(object = dbn$dbn_init$fit, data = t0, by.sample = F)
  
  #
  nodes = paste0(colnames(dt[[1]]), "_t2")
  lltrans = logLik(object = dbn$fits_paired[[1]], data = dtpaired, nodes = nodes, by.sample = F)
  
  if (debug)  { cat("logLik init =",sum(llinit[is.finite(llinit)]), "logLik transitions =",sum(lltrans),"\n")
  }
  llinit + lltrans
}