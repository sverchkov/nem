# ######################## #
# MCMC Sampling: Framework #
# ######################## #

# Method description: this method provides the framework for the MCMC sampling
# parameters: see runMCMC.m

mcmc = function(ratio,nrRuns,theta_init,changeList,logitprior.theta,prior.hidden,maxsteps_eminem,sd_val,signals,effects,changeHfreq,probVal,ep,lowMemFootprint=FALSE){

  # initialize the lists where the results of each sampling step are stored
  theta_list = list()
  graph2_list = list()
  ll_list = list()
  pP_list = list()
  nrSteps_list = list()
  acc_list = list()
  hidden_list = list()

  # Ensure lowMemFootprint is a logical
  lowMemFootprint <- !is.null(lowMemFootprint) && lowMemFootprint
  
  # the initial signals graph corresponds to the sampled graph, not to the local maximum -> it has to be optimized, too + the corresponding values have to be calculated
  graph2_i = theta_init
  res_EMiNEM = EMforMCMC(ratio,graph2_i,logitprior.theta,prior.hidden,maxsteps_eminem)
  theta_i = res_EMiNEM[["res"]]
  nrSteps_i = res_EMiNEM[["nrSteps"]]
  # calculate the log-likelihood for the local maximum
  ll_i = calcLikelihood(ratio,theta_i,prior.hidden)
  # calculate the prior for the local maximum -> see the paper for a derivation of this formula
  pP_i = sum(theta_i[changeList]*logitprior.theta[changeList])
  # set the initial acceptance value to "accepted"
  acc_i = 1
  
  # add the initial effects graph prior to the corresponding result-list
  hidden_list[[1]] = prior.hidden
  
  # initialize count variable for sampling
  i = 1
 
  # start the sampling process
  repeat{
  
    # add the current values to the corresponding results-lists before they are changed (also the last ones BEFORE the sampling is stopped)
    # Except if we're keeping a low memory footprint, then only keep the last value in the result list.
    if ( lowMemFootprint )
      result_index = 1
    else
      result_intex = i

    theta_list[[i]] = theta_i # Theta list needs to be maintained for EBayes
    graph2_list[[result_index]] = graph2_i
    ll_list[[result_index]] = ll_i
    pP_list[[result_index]] = pP_i
    nrSteps_list[[result_index]] = nrSteps_i
    acc_list[[result_index]] = acc_i
    
    # if maximal number of runs is reached -> stop the sampling process
    if (i==nrRuns) break()
    
    # if the number of steps after which the Empirical Bayes step should be applied is reached -> update the effects graph prior and add it to its result-list
    if(i%%changeHfreq==0){
      prior.hidden = updateHidden(ratio,prior.hidden,theta_list,signals,effects,i,changeHfreq)
      hidden_list[[i/changeHfreq+1]] = prior.hidden
    }

    # conduct the sampling step -> updated values (theta_i, ll_i, ...) are returned -> update the corresponding variables
    res_i = onestep_mcmc(ratio,theta_i,graph2_i,ll_i,pP_i,nrSteps_i,changeList,logitprior.theta,prior.hidden,maxsteps_eminem,sd_val,probVal,ep)
    theta_i = res_i[["theta"]]
    graph2_i = res_i[["graph2"]]
    ll_i = res_i[["ll"]]
    pP_i = res_i[["pP"]]
    nrSteps_i = res_i[["nrSteps"]]
    acc_i = res_i[["accepted"]]

    # reduce memory?
    rm("res_i")
    
    # update count variable
    i = i + 1
  
  }
  
  # return result-lists to the main function
  return(list(theta_list=theta_list,graph2_list=graph2_list,ll_list=ll_list,pP_list=pP_list,nrSteps_list=nrSteps_list,acc_list=acc_list,hidden_list=hidden_list))

}







