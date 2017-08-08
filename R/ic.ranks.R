#' Confidence intervals for ranks
#' 
#' This function calculates simultaneous confidence (sets) intervals (CIs) at a pre-specified level (1-alpha) for the ranks of centers mu_1,...,mu_n which are observed through a sample y using multiple testing techniques. Several possibilities are presented through a "Method" variable. There are bascially two main choices; one which uses the partitioing principle and the likelihood ratio test and the the other is based on Tukey's pairwise comparison procedure. See choices below, and for more details see the references.
#' @param y a real vector of observed data.
#' @param sigma a vector of standard deviations. If sigma is a single value, then we consider that all centers have the same standard deviation.
#' @param Method a character indicating the method used to produce the confidence intervals. The "Exact" produces confidence intervals using the partitioning principle and the likelihood ratio test. The "Bound" choice produces lower- or upper-bound confidence intervals (according to the "BoundChoice") for the ranks using a fast algorithm. The "Tukey" choice produces simultaneous confidence intervals for the ranks using Tukey's HSD. The "SeqTukey" produces simultaneous confidence intervals for the ranks using a sequential-rejective algorithm.
#' @param BoundChoice a character entry which is only relevant if the "Bound" choice is picked in the Method parameter. The default value is "Upper" which results in the upper-bound CIs for the ranks. If "Lower" is chosen, then the lower-bound CIs are generated.
#' @param ExactAlgo a character entry default to "Block". This parameter is only relevant when "Exact" choice in the Method parameter is chosen. This parameter gives the choice between two C++ implemented algorithms for an exact partitioning result. See below for more details.
#' @param alpha the significance level of the internal tests we perform (which corresponds to the FWER control of the corresponding multiple testing procedure). CIs have simultaneous significance level of 1-alpha.
#' @param control is a list of control parameters.
#' @param crit is the critical value for Tukey's HSD. If it is kept NULL, then it is calculated internally. The use of this parameter becomes handful in case the user wishes to make several simulations. By providing it, we avoid repeating a Monte-Carlo estimation of the quantile and thus we gain in execution time. In some cases (espcially when all centers have the same standard deviation), the critical value for Tukey's HSD can be found in some statistical tables.
#' @param trace a logical parameter which supresses the printing of the details of the method which was chosen. The default is TRUE (shows details).
#' @param adjustL a logical variable (default to FALSE) indicating if an adjustment on the lower bound according to the data must be considered (if possible). This choice is only relevenat if Method is chosen as "Bound" and BoundChoice is chosen as "Lower".
#' @param adjustU a logical variable (default to FALSE) which gives the user the choice to adjust the upper bound CIs through the parameter "n_adjust". This choice is only relevenat if Method is chosen as "Bound" and BoundChoice is chosen as "Upper".
#' @param n_adjust an integer-valued entry for advanced control over the lower- or upper-bound algorithms. When the "adjustL" parameter is TRUE, the new value of n_adjust is chosen automatically as the best adjustment on the lower affine bound of the chi-square quantiles according to the data. If adjustU is TRUE, then n_adjust contains the point at which the upper affine bound is tangent on the chi-square quantiles. Possible values {1,...,n-1}. If both adjustL and adjustU variables are left FALSE, then the default choice is that the lower affine bound passes between the chi-square quantiles at 1 and n-1 degrees of freedom, and the upper affine bound is tangent on n-1.
#' @param N the number of iterations used in order to calculate the Studentized range quantile for Tukey's algorithms.
#' @details The vector of observations need to be sorted. Otherwise, it is done internally. The observations are supposed to be realizations of independent Guassian distributions with unknown centers mu_1,...,mu_n and known standard deviations sigma = (sigma_1,...,sigma_n). 
#' @details The exact-partitioning confidence intervals are calculated using two algorithms; one which corresponds to the choice "Level" and another which corresponds to the choice "Block". Both choices use an algorithm with complexity 2^n, but the "Block" algorithm is generally the faster one. In the way it is constructed (the "Block" algorithm), the execution time depends on the data and not always on the size of the data. Thus, if you are lucky, you will get the result for 50 centers in a few seconds, and if you are not, then it might take up to 2 weeks. Both algorithms are written with a C++ code in order to optimize the performance and produce results in a reasonable time for a number of centers below 50. The "Block" algorithm requires lower and upper bounds for the confidence intervals. This is automatically calculated here using the option "Bound" from the Method parameter. Hypotheses in the "Level" algorithm are represented using the combinatorial number system whereas in the "Block" algorithm we use a sparse binary representation which runs faster but not convenient to the "Level" algorithm.
#' @details The lower- and upper-bound CIs are calculated with an algorithm whose complexity is of order n^3. The bracketing obtained from the lower and upper bounds is generally very narrow with a maximum gap of 1. Moreover, in regular situations, the lower and upper bounds coincide on at least 50 percent of the centers yielding the exact-partitioning result. Thus, the bracketing is an alternative for an exact-partitioning algorithm for medium- and large-size samples (n>50). When a calculus of the lower- and upper-bound CIs is required, the default choice is when no adjustment on neither the lower nor the upper bounds is taken into account. Thus, the lower affine bound of the chi-square is a line passing by the quantiles at 1 and n-1 degrees of freedom, whereas the upper affine bound is a line tangent on the chi-square quantiles at n-1 degrees of freedom. The adjustment on the lower bound CIs can in some contexts improve on the CIs and increase the number of centers where the lower and upper bounds coincide. The best option is to adjust for both the lower and upper bounds (separately) in case a complexity of n^3 is not considered high for the problem the you solve.
#' @details Both "Tukey" and "SeqTukey" are based on multiple comparison testing and are superior to the LR-based CIs if the centers are far apart from each others and if the standard deviations are not significantly different from each others among the centers. The sequential rejective variant of Tukey's HSD rejects at least as much as Tukey's HSD and thus produces generally shorter confidence intervals for the ranks. 
#' @author Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet. Correspondence to d.al_mohamad@@lumc.nl
#' @return a list of two vectors containing the lower and upper bounds of the confidence intervals for the sorted observed centers.
#' @references Simultaneous confidence sets for ranks using the partitioning principle - Technical report (Arxiv).
#' @references An improvement of Tukey's HSD with application to ranking institutions (Arxiv).
#' @examples
#' TrueCenters = 1:50
#' n = 10; alpha = 0.05; sigma = runif(n,min=0.5,max=1.5)
#' y = as.numeric(sapply(1:n, function(ll) rnorm(1,TrueCenters[ll],sd=sigma[ll])))
#' ind = sort.int(y, index.return = TRUE)$ix
#' y = y[ind]
#' sigma = sigma[ind] # The sigmas need to follow the order of the y's
#' res = ic.ranks(y, sigma, Method = "Exact", ExactAlgo = "Level",
#'    alpha = 0.05, control = list(trace = TRUE))
#' LowerExact = res$Lower; UpperExact = res$Upper
#' res = ic.ranks(y, sigma, Method = "Bound", BoundChoice = "Lower",
#'    control = list(adjustL = FALSE, adjustU = FALSE))
#' LowerL = res$Lower; UpperL = res$Upper
#' res = ic.ranks(y, sigma, Method = "Bound", BoundChoice = "Upper",
#'    control = list(adjustL = FALSE, adjustU = FALSE))
#' LowerU = res$Lower; UpperU = res$Upper
#' res = ic.ranks(y, sigma, Method = "Tukey")
#' LowerTuk = res$Lower; UpperTuk = res$Upper
#' res = ic.ranks(y, sigma, Method = "SeqTukey")
#' LowerTukSeq = res$Lower; UpperTukSeq = res$Upper
#' LowerExact 
#' LowerTuk 

#' @export
ic.ranks = function(y, sigma = rep(1,length(y)), Method = c("Exact","Bound","Tukey","SeqTukey"), BoundChoice = c("Upper", "Lower"), ExactAlgo = c("Level","Block"), alpha = 0.05, control = list(crit = NULL, trace = TRUE, adjustL = FALSE, adjustU = FALSE, n_adjust = n-1, N = 10^4))
{
if(is.null(Method)) Method = "SeqTukey"
trace = control$trace
if(is.null(trace)) trace = TRUE
if(!(Method %in% c("Exact","Bound","Tukey","SeqTukey"))) {print("Error! Method not supported."); return(0)}
n = length(y)
if(length(sigma) == 1) sigma = rep(sigma,n)
if(length(sigma) != n) {print("Error: sigma and y must have the same length!"); return(0)}
if(n == 1) return(1)
ind = sort.int(y, index.return = T)$ix
y = y[ind]
sigma = sigma[ind]
if(sum(ind != 1:n)) print("The sample had to be sorted in ascending way. Results are shown for the sorted sample.")
ranks = NULL
if(n <= 2 & Method == "Bound") {print("Upper- and Lower-bound CIs require at least three centers"); return(0)}

if(Method == "Exact"){
  crit = qchisq(1-alpha,1:(n-1))
  if(is.null(ExactAlgo)) ExactAlgo = "Block"
  if(!(ExactAlgo %in% c("Level","Block"))){print("Error! Could not recognize your choice for the type of the algorithm to be used."); return(0)}
  cat("\n Performing an exact partitioing procedure.\n")
  if(n>=45) cat("Execution time on regular computers might take (if you are not lucky) more than a week in worst case scenario. Check details in the help!")
  if(ExactAlgo == "Block")
  {
	res = ApproximatePartition(y, sigma, BoundChoice = "Lower", alpha,1)
	ind = which(res$Upper == n)[1]
	n_adjust = min(res$BlockMax[1:ind])+1
	res = ApproximatePartition(y,sigma,"Lower", alpha, n_adjust)
	Lower = res$Lower-1; Upper = res$Upper-1; MinBlock = res$BlockMax
	res = ApproximatePartition(y,sigma,"Upper", alpha, floor(n/2))
	MaxBlock = res$BlockMax
	ranks = PartitioningRankingBlock(y, sigma, crit, MinBlock, MaxBlock, Lower, Upper, n, trace)
	ranks = list(Lower = ranks[,1], Upper = ranks[,2])
  }
  if(ExactAlgo == "Level")
  {
	ranks = PartitioningRankingLevel(y, sigma, crit, n, trace)
	ranks = list(Lower = ranks[,1], Upper = ranks[,2])

  }
}
if(Method == "Bound")
{
  if(is.null(BoundChoice)) BoundChoice = "Upper"
  if(!(BoundChoice %in% c("Upper", "Lower"))){print("Error! Could not recognize your choice whether it is upper of lower bound."); return(0)}
  adjustL = control$adjustL
  if(is.null(adjustL)) adjustL = FALSE
  adjustU = control$adjustU
  if(is.null(adjustU)) adjustU = FALSE
  n_adjust = control$n_adjust
  if(is.null(n_adjust)) n_adjust = n-1

  if(BoundChoice == "Lower")
  {
	ranks = ApproximatePartition(y,sigma,"Lower", alpha, 1)
	if(adjustL == TRUE)
	{
	  ind = which(ranks$Upper == n)[1]
	  n_adjust = min(ranks$BlockMax[1:ind])+1
	  ranks = ApproximatePartition(y,sigma,"Lower", alpha, n_adjust)
	  if(trace == TRUE) cat(paste("\n Adjustment on the lower bound. Intersection with the chi-square quantile at n_adjust = ",n_adjust))
	}
	if(trace == TRUE)
	{
	  ranksU = ApproximatePartition(y,sigma,"Upper", alpha, n-1)
	  NbCorrectRanks = mean(ranksU$Upper == ranks$Upper & ranksU$Lower == ranks$Lower)
	  MaxErrorPerCI = max(ranksU$Upper - ranks$Upper, ranksU$Lower - ranks$Lower) 
	  cat(paste("\n Lower-bound confidence intervals for the ranks are calculated at simultaneous confidence level ",1-alpha))
	  cat(paste(paste("\n Percentage of correct confidence intervals is ", NbCorrectRanks),"\n"))
	  cat(paste("\n Maximum CI-length error per center is ",MaxErrorPerCI))
	}
	
  }else{
	# The upper choice is left as the alternative anyway.
	if(adjustU == FALSE) n_adjust = n-1 # keep regular bounds if no adjustment is required
	if(adjustU == TRUE & trace == TRUE) cat(paste("Adjustment on the upper bound by the user. Tangent on the chi-square quantile at n_adjust = ", n_adjust))
	ranks = ApproximatePartition(y,sigma,"Upper", alpha, n_adjust)
	if(trace == TRUE)
	{
	  ranksL = ApproximatePartition(y,sigma,"Lower", alpha, 1)
	  NbCorrectRanks = mean(ranks$Upper == ranksL$Upper & ranks$Lower == ranksL$Lower)
	  MaxErrorPerCI = max(ranks$Upper - ranksL$Upper, ranks$Lower - ranksL$Lower) 
	  cat(paste("\n Upper-bound confidence intervals for the ranks are calculated at simultaneous confidence level ",1-alpha))
	  cat(paste(paste("\n Percentage of correct confidence intervals is ", NbCorrectRanks),"\n"))
	  cat(paste("\n Maximum CI-length error per center is ",MaxErrorPerCI)); cat("\n")
	}

  }
}

if(Method == "Tukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 crit = control$crit
 if(is.null(crit)){
  x=t(mapply(rnorm,N,0,sigma))
  if(n<100){
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-alpha)   
	rm(d); rm(std); rm(S); rm(Cp) 
  }else
  {
	# Quantile calculus for larger sample size
	Diff = numeric(N)
	for(i in 1:N)
	{
	  for(j in 1:(n-1))
	  {
		Diff[i] = max(Diff[i],abs(x[j,i]-x[(j+1):n,i])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
	  }
	}
  
	crit = quantile(Diff,1-alpha)
	rm(Diff)
  }
 }
 rm(x)
 ranks = tukey(y,sigma,alpha, crit)
 cat(paste("\n Confidence intervals for ranks calculated using Tukey's HSD procedure at simultaneous level", 1-alpha))
}

if(Method == "SeqTukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 ranks = StepDownTukeySeqRej(y,sigma,alpha, N)
 if(trace == TRUE){
   cat(paste("\n Confidence intervals for ranks calculated using a sequential-rejective variant of Tukey's HSD procedure at simultaneous level", 1-alpha))
   cat(paste("\n Number of iterations = ",ranks$NbSteps))
   cat("\n")
 }
}
cat(paste(paste("\n Number of compared centers is ",n),"\n"))
return(list(Lower = ranks$Lower, Upper = ranks$Upper))
}

############################################
# Auxiliary functions. Internal usage only

ApproximatePartition = function(y, sigma, BoundChoice = c("Upper", "Lower"), alpha = 0.05, n_adjust)
{
#library(numDeriv)
critFun = function(x)
{
 if(x<=0) return(0)
 
 slop*x
}

n = length(y)
z = qchisq(1-alpha,1:(n-1))
slop = z[n_adjust] - z[n_adjust-1]; Intercept = z[n_adjust] - slop*n_adjust

if(BoundChoice == "Lower") {slop = (z[n-1] - z[n_adjust])/(n-n_adjust-1); Intercept = z[n_adjust] - slop*n_adjust}

# Calculate the matrix of contributions
LogL = matrix(0, nrow = n, ncol = n)
IndividContribBlock = matrix(0, nrow = n, ncol = n)
for(j in 2:n)
{
 for(i in (j-1):1)
 {
	LogL[i,j] = sum((y[i:j] - sum(y[i:j]/sigma[i:j]^2)/sum(1/sigma[i:j]^2))^2/sigma[i:j]^2)
	
	IndividContribBlock[i,j] = min(IndividContribBlock[i,i:(j-1)]+IndividContribBlock[(i+1):j,j], LogL[i,j] - critFun(j-i))
	# The intercept needs to be substracted so that we do not count it several times, because the sum of contributions is not the contribution of the sum
 }
}
# The minimum contribution must be either negative or zero

# Test the blocks by adding hypotheses to the left and to the right of it
Lower = 1:n; Upper = 1:n
BlockMax = numeric(n)
# Treat the case of mu_1=...=mu_j + hypothesis
for(j in (n-1):2)
{
	if(LogL[1,j] - critFun(j-1) + IndividContribBlock[j+1,n] - Intercept < 0) # Block [1,j] is accepted. Update the CIs.
	{
	  Lower[1:j] = 1
	  Upper[1:j] = pmax(Upper[1:j], j)
	  BlockMax[1] = j - 1
	  break
	}
}
# Treat the case of hypothesis + mu_i=...=mu_n
for(i in 2:(n-1))
{
	if(LogL[i,n] - critFun(n-i) + IndividContribBlock[1,i-1] - Intercept < 0) # Block [i,n] is accepted. Update the CIs.
	{
	  Lower[i:n] = pmin(Lower[i:n],i)
	  Upper[i:n] = n
	  break
	}
}
# Treat the case of hypothesis + mu_i=...=mu_j + hypothesis
for(i in 2:(n-1))
{#print(i)
 for(j in (n-1):(i+1))
 {
	if(LogL[i,j] - critFun(j-i) + IndividContribBlock[1,i-1] + IndividContribBlock[j+1,n] - Intercept < 0) # Block [i,j] is accepted. Update the CIs.
	{
	  Lower[i:j] = pmin(Lower[i:j],i)
	  Upper[i:j] = pmax(Upper[i:j], j)
	  BlockMax[i] = j-i
	  break
	}
 }
}
return(list(Lower = Lower, Upper = Upper, BlockMax = BlockMax))
}

####################################################
# The simple HSD
#library(multcomp)
tukey = function(y,sigma,alpha = 0.05, qq) {
  n=length(y)
  ranks=matrix(0,n,2)
  for(j in 1:n)
  {
	stat = (y[j]-y)/sqrt(sigma[j]^2+sigma^2)
	ranks[j,1]=1+sum(stat>qq)
	ranks[j,2]=n-sum(stat<(-qq))
  }

  return(list(Lower = ranks[,1], Upper = ranks[,2]))
}

###########################################################
# Step-down Tukey with sequential rejection
StepDownTukeySeqRej = function(y,sigma,alpha=0.05,N = 10^4)
{
  n = length(y)
 # Calculate the critical value
 Diff = NULL
 PosPairs = matrix(0,ncol=2,nrow = n*(n-1)/2)
 NegPairs = matrix(0,ncol=2,nrow = n*(n-1)/2)
 for(j in 1:(n-1))
 {
  PosPairs[(j-1)*n-j*(j-1)/2+1:(n-j),1] = rep(j,n-j)
  PosPairs[(j-1)*n-j*(j-1)/2+1:(n-j),2] = (j+1):n  
 }
 NegPairs[,2] = PosPairs[,1] # Those pairs stay untouched because they are never rejected
 NegPairs[,1] = PosPairs[,2] 

 x=t(mapply(rnorm,N,0,sigma))  
 Diff = numeric(N)
 for(k in 1:N)
 {
   for(j in 1:(n-1))
  {
	Diff[k] = max(Diff[k],abs(x[j,k]-x[(j+1):n,k])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
  }
 }

 qDown = quantile(Diff,1-alpha)
 rm(Diff)

 NbNRejOld = n*(n-1)/2; NbNRejNew = 0
 NBSteps = 0
 while(NbNRejOld>NbNRejNew)
 {
  NBSteps = NBSteps +1
  NewPairs = NULL

  for(i in 1:length(PosPairs[,1]))
  {
	T = abs(y[PosPairs[i,1]]-y[PosPairs[i,2]])/sqrt(sigma[PosPairs[i,1]]^2+sigma[PosPairs[i,2]]^2)
	if(T<qDown)
	{
	  NewPairs = rbind(NewPairs,PosPairs[i,])
	  #Diff = c(Diff,T)
	}
  }
 NbNRejOld = length(PosPairs[,1])
 PosPairs = NewPairs
 NbNRejNew = length(PosPairs[,1])
 
 # Calculate again the quantile  
 AllPairs = rbind(PosPairs, NegPairs)
 Diff = numeric(N)
 for(i in 1:N)
 {
	Diff[i] = max((x[AllPairs[,2],i]-x[AllPairs[,1],i])/sqrt(sigma[AllPairs[,1]]^2+sigma[AllPairs[,2]]^2))
 }
 
 qDown = quantile(Diff,1-alpha)
 rm(Diff)

 }
rm(x)
# Extract the ranks
ranks = matrix(0,nrow = n,ncol = 2)
ranks[,1] = 1:n
ranks[,2] = 1:n
indU = which(PosPairs[,1]==1)
if(length(indU) > 0) ranks[1,2] = PosPairs[indU[length(indU)],2]

for(i in 2:(n-1))
{
	indU = which(PosPairs[,1]==i)
	indL = which(PosPairs[,2]==i)
	if(length(indU) > 0) ranks[i,2] = PosPairs[indU[length(indU)],2]
	if(length(indL) > 0 ) ranks[i,1] = PosPairs[indL[1],1]	
}
indL = which(PosPairs[,2]==n)
if(length(indL) > 0) ranks[n,1] = PosPairs[indL[1],1]

return(list(Lower = ranks[,1], Upper = ranks[,2], NbSteps = NBSteps))
}



