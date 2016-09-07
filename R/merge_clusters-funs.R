start.k = 50 # choose initial k for merging
n.iter = 48 # choose number of merging iterations

alloc.test = cutree(cover.flex, start.k)
for (i in 1:n.iter) {
  print(paste0("Starting merging iteration ",i,":"))
  pairs = GeneratePairwise(alloc.test)
  #PairwiseCombinedSpeciesCount(pa, alloc.test, pairs)
  alloc.test.pairwise = PairwiseManyglm(pa, alloc.test, pairs)
  #alloc.test.pairwise.sorted = alloc.test.pairwise[with(alloc.test.pairwise, order(deltaSumAIC)),]
  # calculate summary and merge clusters using this information
  PairwiseSummary(alloc.test, alloc.test.pairwise)
  merge.pair = c(as.numeric(alloc.test.pairwise$target[alloc.test.pairwise$deltaSumAIC==min(alloc.test.pairwise$deltaSumAIC)]),
                 as.numeric(alloc.test.pairwise$test[alloc.test.pairwise$deltaSumAIC==min(alloc.test.pairwise$deltaSumAIC)]))
  alloc.test[alloc.test==merge.pair[1] | alloc.test==merge.pair[2]] = 1000 + i
  print("#############")
  print(paste0("merged cluster: ",merge.pair[1]," and ",merge.pair[2]))
  print("#############")
  assign(x=paste0("alloc.test.",start.k-i), value=alloc.test)
}

# fit and get AIC

alloc.list = list()
for (i in 1:length(2:49)) {
  alloc.list[[i]] = get(paste0("alloc.test.",c(2:49)[i]))
}
save(alloc.list, file="alloc.list.RData")

data = mvabund(pa)
alloc.sumAIC = numeric(length(alloc.list))
for (i in 1:length(alloc.list)) {
  alloc.sumAIC[i] = manyglm(formula=data~as.factor(alloc.list[[i]]), family="binomial")$AICsum
}

alloc.pairwise.AIC = data.frame(AICsum=alloc.sumAIC, groups=c(2:49))
load("pa.AICsum.RData")
pa.AICsum.pairwise = rbind(pa.AICsum, alloc.pairwise.AIC)
save(pa.AICsum.pairwise, file="pa.AICsum.pairwise.RData")
load("pa.AICsum.pairwise.RData")


# plot pairwise modelling merges
AIC.hclust = pa.AICsum.pairwise[1:49,]
AIC.pairwise = pa.AICsum.pairwise[50:97,]




# FUNCTIONS ---------------------------------------------------------------


GeneratePairwise = function(communities) {
  # generate pairwise tests
  communities = as.character(unique(communities))
  target.com = character(round(length(communities)^2/2))
  test.com = character(round(length(communities)^2/2))
  # loop across communities to create pairwise comparisons, without reverse tests
  ticker=0
  for (target in 1:(length(communities)-1)){
    for (test in target:(length(communities)-1)){
      ticker = ticker+1
      target.com[ticker] = communities[target]
      test.com[ticker] = communities[test+1]
    }
  }
  target.com=target.com[!target.com==""]
  test.com=test.com[!test.com==""]
  return(data.frame(target=target.com, test=test.com, stringsAsFactors=F))
}

## function that calculates the species count after removing species with 0 occurance in a pariwise combination
PairwiseCombinedSpeciesCount = function(multivar, predictor, pairs) {
  target = as.character(pairs$target)
  test = as.character(pairs$test)
  for (i in 1:length(target)){
    # subset to target sites
    target.rows = which(predictor==target[i])
    target.multivar = multivar[target.rows,]
    target.pred = predictor[target.rows]
    # subset to test sites
    test.rows = which(predictor==test[i])
    test.multivar = multivar[test.rows,]
    test.pred = predictor[test.rows]
    # create data objects for the pariwise test
    combined.multivar = rbind(target.multivar, test.multivar)
    combined.multivar = combined.multivar[,colSums(combined.multivar)>0] # remove species that don't occur in either community
    print(ncol(combined.multivar))
  }
}

PairwiseManyglm = function(multivar, predictor, pairs) {
  target = as.character(pairs$target)
  test = as.character(pairs$test)
  rank = numeric(nrow(pairs))
  nspecies = numeric(nrow(pairs))
  deltaSumAIC = numeric(nrow(pairs))
  for (i in 1:length(target)){
    # subset to target sites
    target.rows = which(predictor==target[i])
    target.multivar = multivar[target.rows,]
    target.pred = predictor[target.rows]
    # subset to test sites
    test.rows = which(predictor==test[i])
    test.multivar = multivar[test.rows,]
    test.pred = predictor[test.rows]
    # create data objects for the pariwise test
    combined.multivar = rbind(target.multivar, test.multivar)
    combined.multivar = mvabund(data.matrix(combined.multivar))
    #combined.multivar = mvabund(data.matrix(combined.multivar[,colSums(combined.multivar)>0])) # remove double absences
    combined.pred = as.factor(c(as.character(target.pred), as.character(test.pred)))
    # fit/test the model
    fit = manyglm(combined.multivar ~ combined.pred, family="negative.binomial")
    fit.null = manyglm(combined.multivar ~ 1, family="negative.binomial")
    # calculate % of species for which dAIC (from null) is >n
    dAIC = fit.null$aic-fit$aic
    rank[i] = sum(dAIC>4)/length(dAIC)
    nspecies[i] = ncol(combined.multivar)
    deltaSumAIC[i] = fit.null$AICsum - fit$AICsum
  }
  return(data.frame(pairs, rank, nspecies, deltaSumAIC, stringsAsFactors=F))
}

## function that calculates lowest ranked pairwise comparison and lowest ranked cluster
# alloc==allocation used, pairwisemanyglm=object from PairwiseManyglm() using alloc
PairwiseSummary = function(alloc, pairwisemanyglm) {
  # print cluster numbers
  print(table(alloc))
  # find lowest ranked pairwise comparison
  print(pairwisemanyglm[pairwisemanyglm$deltaSumAIC==min(pairwisemanyglm$deltaSumAIC),])
  # find cluster that has lowest mean ranks
  #   for (i in unique(alloc)) {
  #     print(paste0("Cluster ",i," mean delta sum-of-AIC:"))
  #     print(mean(pairwisemanyglm$deltaSumAIC[pairwisemanyglm$target==i | pairwisemanyglm$test==i]))
  #   }
}
