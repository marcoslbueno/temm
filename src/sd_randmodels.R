library(bnlearn)
library (gtools)

# Generates a random BN graph
rand.structure.bn <-  function (numvar, method.g = "melancon", prob.node = NULL, max.degree = NULL) {
  if (method.g == "empty") {
    g <- empty.graph(paste("N",as.character(1:numvar), sep=""))
  } else {
    if (!is.null(prob.node)) {
      g <- random.graph(paste("N",as.character(1:numvar), sep=""), prob = prob.node)      
    } else {
      if (!is.null(max.degree)) {
        g <- random.graph(paste("N",as.character(1:numvar), sep=""), method = method.g, max.degree = max.degree)
      } else {
        g <- random.graph(paste("N",as.character(1:numvar), sep=""), method = method.g)
      }
    }
  }
  g
}

# Generates a random BN graph + probabilities
buildrand.bn <- function (domsizes, method.g = "melancon", g = NULL, prob.node = NULL) { 
  numvar <- length(domsizes)
  #   1 - Defined graph, Random graph, or Empty graph
  if (is.null(g)) {
    if (is.null(prob.node)) {
      g = rand.structure.bn (numvar = numvar, method.g = method.g)  
    } else {
      g = rand.structure.bn (numvar = numvar, method.g = method.g, prob.node = prob.node)  
    }
  } 
  
  #   2 - CPTs
  m1 <- matrix (domsizes)
  dimnames (m1) <- list(nodes(g))
  arcs <- arcs (g)
  dist <- list ()
  for (i in 1:length(nodes(g)))   {   
    node <- nodes(g)[i]  # "nodes(g)" returns in alphabetical order
    #   Get the parents
    pa <- arcs[arcs[,2]==node, ]
    if (!is.matrix(pa)) { 
      pa <- c(pa[1])  # 0 or 1 parent
    } else { 
      pa <- c(pa[,1])
    }       
    
    #     Generate rows for the CPT
    varsize <- domsizes[i]
    if (length(pa)==0) {
      pasizes <- c()
      prob <- prob.rows (varsize, pasizes)      
      cpt <- matrix (prob, ncol=varsize, dimnames=list(NULL, c(1:varsize)))
    } else {
      pasizes <- c()
      for (j in 1:length(pa)) { pasizes <- c(pasizes, m1[pa[j],1]) }
      prob <- prob.rows (varsize, pasizes) 
      
      cpt <- prob
      di <- c(varsize)
      li <- list(node=c(1:varsize))
      names(li)[1] <- node
      for (j in 1:length(pa))  {  
        li[[j+1]] <- c(1:pasizes[[j]])    
        names(li)[j+1] <- pa[j]
        di <- c(di, pasizes[[j]])
      }
      dim(cpt) <- di
      dimnames(cpt) <- li
    }
    dist[[i]] <- cpt
  }
  
  #    Construct a fit object
  names(dist) <- nodes(g)
  g.fit <- custom.fit (g, dist=dist)
  list("graph"=g, "fit"=g.fit)
}

prob.rows <- function (varsize, pasizes) {
  if (length(pasizes) == 0) {
    ncombs <- 1
  } else {
    ncombs <- prod (pasizes) 
  }
  prob <- rdirichlet(ncombs, rep(1,varsize))
  trans <- t(prob)
  allrows <- c(trans)
  allrows
}

# Generates a random Dynamic Bayesian Network 
rand.tstructure.dbn <- function (numvar, pr_trans_same = NULL, pr_trans_other = NULL, pr_node = NULL) {
  # 1 - random structure for nodes at t+1 
  # 2 - make a structure with nodes t and t+1, copying arcs from struct. t+1
  if (is.null(pr_node)) { pr_node = 0.1 }
  
  str_right = rand.structure.bn (numvar = numvar, prob.node = pr_node)
  nodes = nodes(str_right)
  nodes(str_right) = paste0(nodes,"_t2")
  
  str_full = rand.structure.bn (2*numvar, method.g = "empty")
  #graphviz.plot(str_full)
  nodes(str_full) = c(paste0(nodes,"_t1"), paste0(nodes,"_t2"))
  
  arcs(str_full) = arcs(str_right)
  #graphviz.plot(str_full)
  
  # 3 - make transition edges
  if (is.null(pr_trans_same)) { pr_trans_same = 0.8 } 
  if (is.null(pr_trans_other)) { pr_trans_other = 0.1 }
  
  for (nodel in nodes) { # nodel = nodes[1]
    for (noder in nodes) { # noder = nodes[1]
      if (nodel == noder) { 
        # transition to the same variable 
        if (runif(n = 1) < pr_trans_same) {
          row = c(paste0(nodel,"_t1"), paste0(noder,"_t2"))
          arcs(str_full) = rbind(arcs(str_full), row)
        }
      } else {
        # transition to other variable
        if (runif(n = 1) < pr_trans_other) {
          row = c(paste0(nodel,"_t1"), paste0(noder,"_t2"))
          arcs(str_full) = rbind(arcs(str_full), row)
        }
      }
    }
  }
  
  str_full
}


buildrand.dbn <- function (domsizes, pr_trans_same = NULL, pr_trans_other = NULL, pr_node = NULL){
  tstr = rand.tstructure.dbn (numvar = length(domsizes), pr_trans_same, pr_trans_other, pr_node)
  # nodes(tstr)
  # graphviz.plot(tstr)
  dbn = buildrand.bn (domsizes = c(domsizes,domsizes), g = tstr)
  dbn
}

sequences.fromdataset <- function (train) { 
  sequences = list()
  begin = 1
  for (i in 1:length(train$N)) { 
    # i=3; print(i)
    
    end = begin + train$N[i] - 1

    if (end == begin) {
      # ! Won't work unless rownames(train$x) is NULL
      sequences[[i]] = t(train$x[begin:end, ,drop=F])
      
    } else {
      sequences[[i]] = train$x[begin:end, ,drop=F]
    }
    
    begin = end + 1
  }
  
  sequences
}

convertfactors = function (data, theoretical_factors, debug = T) {
  data = as.data.frame(data)
  for (i in 1:ncol(data)) {
    # i=1
    
    if (is.factor(data[,i])){
      col = levels(data[,i])[data[,i]]
    } else {
      col = data[,i]
    }
    
    l = theoretical_factors[[i]]
    col = factor(x = col, levels = l)
    data[,i] = col
    
    if (!is.null(colnames(data))) {
      name=colnames(data)[i]
    } else {
      name = i
    }
    
    if (debug) {
      cat(name,length(unique(data[,i])), " Theoretical: ", length(theoretical_factors[[i]]),"\n")
    }
  }
  data
}
