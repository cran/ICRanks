#' Confidence intervals for ranks
#' 
#' This function calculates simultaneous confidence (sets) intervals (CIs) at a pre-specified level (1-alpha) for the ranks of centers mu_1,...,mu_n which are observed through a sample y using multiple testing techniques. Several possibilities are presented through a "Method" variable. There are bascially two main choices; one which uses the partitioing principle and the likelihood ratio test and the the other is based on Tukey's pairwise comparison procedure. See choices below, and for more details see the references.
#' @param y a real vector of observed data.
#' @param sigma a vector of standard deviations. If sigma is a single value, then we consider that all centers have the same standard deviation.
#' @param Method a character indicating the method used to produce the confidence intervals. The "Exact" produces confidence intervals using the partitioning principle and the likelihood ratio test. The "Bound" choice produces lower- or upper-bound confidence intervals (according to the "BoundChoice") for the ranks using a fast algorithm. The "Tukey" choice produces simultaneous confidence intervals for the ranks using Tukey's HSD. The "SeqTukey" produces simultaneous confidence intervals for the ranks using a sequential-rejective algorithm. The "Approximate" choice provides approximate confidence intervals which are shorter than exact ones by considering a subset of the partitions (the correctly ordered ones, see refs). The "TukeyNoTies" choice calculates a readustement for Tukey's method under the assumption that there are no ties and then use Tukey's method again with adjusted level.
#' @param BoundChoice a character entry which is only relevant if the "Bound" choice is picked in the Method parameter. The default value is "Upper" which results in the upper-bound CIs for the ranks. If "Lower" is chosen, then the lower-bound CIs are generated.
#' @param ApproxAlgo a character entry ("Upper" by default). This parameter controls which approximation is to be used. 
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
#' @details The "TukeyWithoutTies" method assumes that the true vector of parameters has no ties and therefore, instead of calculating a quantile q corresponding to mu=0 with set rank [1,n] for mu_i, we calculate a quantile corresponding to mu=0 with rank {i} for mu_i. The method provides shorter SCI for the ranks but is still conservative.
#' @details When the standard deviations are not the same for all the means, the methods based on the partitioning principle are no guaranteed to be the same. The "Block" algorith, however, is always compatible with the lower and upper CIs provided by option "Bound".
#' @details When the standard deviations are not the same the approximate methods based on the LRT are not guaranteed to cover and if the standard deviations are very different, the resulting SCIs are anticonservative.
#' @author Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet. Correspondence to d.al_mohamad@@lumc.nl (or diaa.almohamad@@gmail.com)
#' @return a list of two vectors containing the lower and upper bounds of the confidence intervals for the sorted observed centers.
#' @references Diaa Al Mohamad and Erik W. van Zwet and Jelle J. Goeman and Aldo Solari, Simultaneous confidence sets for ranks using the partitioning principle - Technical report (2017). https://arxiv.org/abs/1708.02729
#' @references Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet, An improvement of Tukey's HSD with application to ranking institutions (2017). https://arxiv.org/abs/1708.02428
#' @examples
#' n = 10; TrueCenters = 1:n
#' alpha = 0.05; sigma = runif(n,min=0.5,max=1.5)
#' y = as.numeric(sapply(1:n, function(ll) rnorm(1,TrueCenters[ll],sd=sigma[ll])))
#' ind = sort.int(y, index.return = TRUE)$ix
#' y = y[ind]
#' sigma = sigma[ind] # The sigmas need to follow the order of the y's
#' res = ic.ranks(y, sigma, Method = "ExactLR",alpha = 0.05, control = list(trace = TRUE))
#' LowerExact = res$Lower; UpperExact = res$Upper
#' res = ic.ranks(y, sigma, Method = "BoundLR", BoundChoice = "Lower",
#'    control = list(adjustL = FALSE, adjustU = FALSE))
#' LowerL = res$Lower; UpperL = res$Upper
#' res = ic.ranks(y, sigma, Method = "BoundLR", BoundChoice = "Upper",
#'    control = list(adjustL = FALSE, adjustU = FALSE, trace=FALSE))
#' LowerU = res$Lower; UpperU = res$Upper
#' res = ic.ranks(y, sigma, Method = "Tukey")
#' LowerTuk = res$Lower; UpperTuk = res$Upper
#' res = ic.ranks(y, sigma, Method = "SeqTukey")
#' LowerTukSeq = res$Lower; UpperTukSeq = res$Upper
#' res = ic.ranks(y, sigma, Method = "TukeyNoTies")
#' LowerTukNoTies = res$Lower; UpperTukNoTies = res$Upper
#' LowerExact 
#' LowerL
#' LowerU
#' LowerTuk 
#' @export
ic.ranks = function(y, sigma = rep(1,length(y)), Method = c("ExactLR","BoundLR","Tukey","SeqTukey","ApproximateLR", "TukeyNoTies"), BoundChoice = c("Upper", "Lower"), ApproxAlgo = c("Exact","Upper"), alpha = 0.05, control = list(crit = NULL, trace = TRUE, adjustL = FALSE, adjustU = FALSE, n_adjust = length(y)-1, N = 10^4))
{
# .................................................
# --------------------------
# Some checking before we start...
# --------------------------
# .................................................
if(length(Method) != 1) Method = "SeqTukey"
trace = control$trace
if(is.null(trace)) trace = TRUE
if(!(Method %in% c("ExactLR","BoundLR","Tukey","SeqTukey","ApproximateLR", "TukeyNoTies"))) {print("Error! Method not supported."); return(0)}
n = length(y)
if(length(sigma) == 1) sigma = rep(sigma,n)
if(length(sigma) != n) {print("Error: sigma and y must have the same length!"); return(0)}
if(n == 1) return(1)
ind = sort.int(y, index.return = T)$ix
y = y[ind]
sigma = sigma[ind]
if(sum(ind != 1:n)) print("The sample had to be sorted in ascending way. Results are shown for the sorted sample.")
ranks = NULL
if(n <= 2 & Method == "BoundLR") {print("Upper- and Lower-bound CIs require at least three centers"); return(0)}

# .................................................
# --------------------------
# The Full Partitioning Principle
# --------------------------
# .................................................
if(Method == "ExactLR")
{
  crit = qchisq(1-alpha,(n-1):1)
  ranks = PartitioningRankingLevel(y, sigma, crit, n, trace)
  ranks = list(Lower = ranks[,1], Upper = ranks[,2])
}
# .................................................
# --------------------------
# The Bracketing algorithm for the FULL Partitioning Principle
# --------------------------
# .................................................
if(Method == "BoundLR")
{
  if(length(BoundChoice)!= 1) BoundChoice = "Upper"
  if(!(BoundChoice %in% c("Upper", "Lower"))){print("Error! Could not recognize your choice whether it is upper of lower bound."); return(0)}
  adjustL = control$adjustL
  if(is.null(adjustL)) adjustL = FALSE
  adjustU = control$adjustU
  if(is.null(adjustU)) adjustU = FALSE
  n_adjust = control$n_adjust
  if(is.null(n_adjust)) n_adjust = n-1
  n_adjust = floor(n_adjust)
  if(n_adjust > n-1 | n_adjust<1){n_adjust = n-1; cat(paste("n_adjust can take values only between 1 and ", n-1)); cat(". Default value is considered.")}
 
  if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
  {
	if(trace==TRUE) cat("Process ended with trivial confidence intervals.\n")
	return(list(Lower = rep(1,n), Upper = rep(n,n)))
  }
  if(BoundChoice == "Lower")
  {
	if(trace == TRUE) cat('\n Calculate lower bounds for simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
	ranks = ApproximatePartition(y,sigma,"Lower", alpha, 1)
	if(adjustL == TRUE)
	{
	  ind = which(ranks$Upper == n)[1]
	  n_adjust = min(ranks$BlockMax[1:ind])+1
	  ranks = ApproximatePartition(y,sigma,"Lower", alpha, n_adjust, trace)
	  if(trace == TRUE) cat(paste("\n Adjustment on the lower bound. Intersection with the chi-square quantile curve at n_adjust = ",n_adjust))
	}	
  }else{
	# The upper choice is left as the alternative anyway.
	if(trace == TRUE) cat('\n Calculate upper (conservative) bounds for simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
	if(adjustU == FALSE) n_adjust = n-1 # keep regular bounds if no adjustment is required
	if(adjustU == TRUE & trace == TRUE) cat(paste("Adjustment on the upper bound by the user. Tangent on the chi-square quantile at n_adjust = ", n_adjust))
	ranks = ApproximatePartition(y,sigma,"Upper", alpha, n_adjust, trace)
  }
}
# .................................................
# --------------------------
# The Partitioning Principle when only the correctly ordered hypotheses are considered
# --------------------------
# .................................................
if(Method == "ApproximateLR")
{
 if(length(ApproxAlgo)!= 1) ApproxAlgo = "Upper"
 if(!(ApproxAlgo %in% c("Exact","Upper"))) {print("Error! Approximate algorithm not supported."); return(0)}
 if(trace == TRUE) cat('\n Calculate approximate simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
 if(ApproxAlgo == "Upper")
 {
   if(trace == TRUE) cat('\n A fast (cubic complex) algorithm is being used.\n')
   ranks = ApproximatePartitionCorrectOrder(y,sigma,"Upper",alpha,n-1)
 }else
 {
   crit = qchisq(1-alpha,(n-1):1)
   if(trace == TRUE) cat('\n A slow (exponentially complex) algorithm is being used.\n')
   res = ApproximatePartitionCorrectOrder(y, sigma, BoundChoice = "Lower", alpha,1)
   ind = which(res$Upper == n)[1]
   n_adjust = min(res$BlockMax[1:ind])+1
   res = ApproximatePartitionCorrectOrder(y,sigma,"Lower", alpha, n_adjust)
   Lower = res$Lower-1; Upper = res$Upper-1; MinBlock = res$BlockMax
   res = ApproximatePartitionCorrectOrder(y,sigma,"Upper", alpha, floor(n/2))
   MaxBlock = res$BlockMax
   ranks = PartitioningRankingBlockCorrectOrder(y, sigma, crit, MinBlock, MaxBlock, Lower, Upper, n, trace)
   ranks = list(Lower = ranks[,1], Upper = ranks[,2])
 }

}

# .................................................
# --------------------------
# General Tukey (ties are allowed)
# --------------------------
# .................................................
if(Method == "Tukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 crit = control$crit
 if(is.null(crit))
 {
 if(length(unique(sigma)) == 1 & alpha<0.3) # a quick calculation
 {
   crit = qtukey(1-alpha,n,Inf)/sqrt(2)
 }else
 {
 x=t(mapply(rnorm,N,0,sigma))
 if(n<100){
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-alpha)   
	rm(d); rm(std); rm(S); rm(Cp); rm(x)
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
	rm(x)
  }
 }
 }
 
 ranks = tukey(y,sigma,crit)
 if(trace == TRUE) cat(paste("\n Confidence intervals for ranks calculated using Tukey's HSD procedure at simultaneous level", 1-alpha))
}

# .................................................
# --------------------------
# Sequential Tukey
# --------------------------
# .................................................
if(Method == "SeqTukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 ranks = StepDownTukeySeqRej(y,sigma,alpha, N)
 if(trace == TRUE)
 {
   cat(paste("\n Confidence intervals for ranks calculated using a sequential-rejective variant of Tukey's HSD procedure at simultaneous level", 1-alpha))
   cat(paste("\n Number of iterations = ",ranks$NbSteps))
   cat("\n")
 }
}
# .................................................
# --------------------------
# Tukey under the assumption that there are no ties
# --------------------------
# .................................................
if(Method == "TukeyNoTies")
{
N = control$N
if(is.null(N)) N = 10^4
d = numeric(N) # a vector of differences that might be useful later
if(length(unique(sigma)) > 1 | alpha<0.3)
{
 x=t(mapply(rnorm,N,0,sigma))
 if(n<100)
 {
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	d = apply(abs(d),2,max)
	rm(std); rm(S); rm(Cp); rm(x)
  }else
  {
	# Quantile calculus for larger sample size
	d = numeric(N)
	for(i in 1:N)
	{
	  for(j in 1:(n-1))
	  {
		d[i] = max(d[i],abs(x[j,i]-x[(j+1):n,i])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
	  }
	}
  	rm(x)
  }
}

# Estimate the coverage of Tukey at alpha when mu=0
x=t(mapply(rnorm,N,0,sigma))
TukeyCoverage = function(a)
{
 q = 1
 if(length(unique(sigma)) > 1 | alpha<0.3)
 {
	q=quantile(d,1-a)
 }else
 {
	q = qtukey(1-a,n,Inf)/sqrt(2)
 }
 TrueLowerRank = 1:n; TrueUpperRank = 1:n
 coverageTuk = N
 for(i in 1:N)
 {
  y = x[,i]
  ind = sort.int(y, index.return = T)$ix
  y = y[ind]
  resTukey = tukey(y,sigma[ind], q)
  if(sum(TrueLowerRank[ind]<resTukey$Lower | TrueUpperRank[ind]>resTukey$Upper)>0) coverageTuk = coverageTuk - 1
 }
 coverageTuk/N
}

# Calculate a modified significance level
alphaTuk = uniroot(function(a)TukeyCoverage(a)-(1-alpha), c(alpha,0.9), maxiter=15)$root
rm(x)

# Use the new alpha to calculate the SCI for ranks
crit = 1
if(length(unique(sigma)) == 1 & alphaTuk<0.3) # a quick calculation
{
   crit = qtukey(1-alphaTuk,n,Inf)/sqrt(2)
}else
{
 x=t(mapply(rnorm,N,0,sigma))
 if(n<100)
 {
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-alphaTuk)   
	rm(d); rm(std); rm(S); rm(Cp); rm(x)
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
  
	crit = quantile(Diff,1-alphaTuk)
	rm(Diff)
	rm(x)
 }
}
 
 ranks = tukey(y,sigma, crit)
 if(trace == TRUE) 
 {
	cat(paste("\n Confidence intervals for ranks calculated using Tukey's HSD procedure at simultaneous level", 1-alpha))
	cat(paste("\n Rescaled significance level is ",alphaTuk)); cat(".")
 }
 
  
}

 if(trace == TRUE) cat(paste(paste("\n Number of compared centers is ",n),"\n"))
 return(list(Lower = ranks$Lower, Upper = ranks$Upper))
}

############################################
# Auxiliary functions. Internal usage only
ApproximatePartition = function(y, sigma, BoundChoice = c("Upper", "Lower"), alpha = 0.05, n_adjust, trace = FALSE)
{
critFun = function(x)
{
 if(x<=0) return(0)
 
 slop*x
}

n = length(y)
z = qchisq(1-alpha,1:(n-1))
slop = z[n_adjust] - z[n_adjust-1]; Intercept = z[n_adjust] - slop*n_adjust

if(BoundChoice == "Lower") {slop = (z[n-1] - z[n_adjust])/(n-n_adjust-1); Intercept = z[n_adjust] - slop*n_adjust}

### Part 0: calculate the SCI when no permutation is applied
# Calculate the matrix of contributions
if(trace == TRUE) cat('\n Caclulate simultaneous confidence intervals using the correctly ordered hypotheses.\n')
EmpOrder = 1:n
res = ApproximatePartitionCorrectOrder(y, sigma, BoundChoice = BoundChoice, alpha = alpha, n_adjust=n_adjust)
Lower = res$Lower; Upper = res$Upper
if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
minY = min(y)
maxY = max(y)
if(trace == TRUE) {cat('\n Forward permutations (');cat(n-1); cat(' permutations ).\n')}
### part1: permute (i,i+1,i+2,...)
for(I in 1:(n-1))#row index
{
 if(trace == TRUE) {cat(I); cat('.')}
 for(J in 1:(n-I))#column index
 {
	y_temp = y; sigma_temp = sigma
	EmpOrder_temp = 1:n
	# Permute
	for(s in J:(J+I-1))
	{
	  y_temp[s+1] = y[s]
	  sigma_temp[s+1] = sigma[s]
	  EmpOrder_temp[s+1] = EmpOrder[s] 
	}
	y_temp[J] = y[J+I]
	sigma_temp[J] = sigma[J+I]
	EmpOrder_temp[J] = EmpOrder[J+I]
	
	# Calculate the SCI for ranks
	res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)
	# Permute the ranks
	for(s in J:(J+I-1))
	{
	  Lower[s] = min(Lower[s], res$Lower[s+1])	
	  Upper[s] = max(Upper[s], res$Upper[s+1])
	}
	Lower[J+I] = min(Lower[J+I], res$Lower[J])
	Upper[J+I] = max(Upper[J+I], res$Upper[J])
	if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
 }

}

if(trace == TRUE) {cat('\n Applying backward permutations.\n')}
### Part 2: permute (i+4,i+3,..)
for(J in 1:(n-1))#column index
{
 if(trace == TRUE) {cat(J); cat('.')}
 I = 2
 while(I<=(n-J))#row index
 {
	y_temp = y; sigma_temp = sigma; EmpOrder_temp = 1:n
	# Permute
	for(s in J:(J+I-1))
	{
	  y_temp[s] = y[s+1]
	  sigma_temp[s] = sigma[s+1]
	  EmpOrder_temp[s] = EmpOrder[s+1]
	}
	y_temp[J+I] = y[J]
	sigma_temp[J+I] = sigma[J]
	EmpOrder_temp[J+I] = EmpOrder[J]
	
	# Calculate the SCI for ranks
	res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)

	# permute the ranks
	for(s in J:(J+I-1))
	{
	  Lower[s+1] = min(Lower[s+1], res$Lower[s])
	  Upper[s+1] = max(Upper[s+1], res$Upper[s])
	}
	Lower[J] = min(Lower[J], res$Lower[J+I])
	Upper[J] = max(Upper[J], res$Upper[J+I])
 	I = I+1
	if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))

 }

}
return(list(Lower = Lower, Upper = Upper))
}

# The function without permutations and the centers are assumed ordered.
ApproximatePartitionCorrectOrder = function(y, sigma, BoundChoice = c("Upper", "Lower"), alpha = 0.05, n_adjust)
{
#library(numDeriv)
n = length(y)
if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
{
	return(list(Lower = rep(1,n), Upper = rep(n,n), BlockMax = rep(n-1,n)))
}

critFun = function(x)
{
 if(x<=0) return(0)
 
 slop*x
}

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
for(i in 2:(n-2))
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

# The function with possibly wrong orderings. Wrong orderings which requires pooling are not accepted and are ignored.
ApproximatePartitionWrongOrder = function(y_temp, sigma_temp, EmpOrder_temp, critFun, Intercept)
{
# Make a few calculations before we start
maxY = max(y_temp) # could be calculated above
minY = min(y_temp)# could be calculated above
n = length(y_temp)
# Calculate the matrix of contributions
LogL = matrix(0, nrow = n, ncol = n)
#### If the Block has several sublocks, I need to keep track of the smallest average and the largest one !
for(j in 2:n)
{
 for(i in (j-1):1)
 {
	LogL[i,j] = sum((y_temp[i:j] - sum(y_temp[i:j]/sigma_temp[i:j]^2)/sum(1/sigma_temp[i:j]^2))^2/sigma_temp[i:j]^2)
 }
}
res = IndividContribs(y_temp,sigma_temp,LogL, critFun,1,n,minY,maxY,TRUE)
IndividContribBlock = res$IndividContribBlock
AverageBlock = res$AverageBlock

ContribRBlock = 0
ContribLBlock = 0
# Test the blocks by adding hypotheses to the left and to the right of it
Lower = EmpOrder_temp; Upper = EmpOrder_temp
BlockMax = numeric(n)
# Treat the case of mu_1=...=mu_j + hypothesis
for(j in (n-1):2)
{ # What if the minimum contributed block does not fulfill the ordering of the PAVA whereas another configuration which slightly higher contribution does ?!
  AvCurrBlock = sum(y_temp[1:j]/sigma_temp[1:j]^2)/sum(1/sigma_temp[1:j]^2)
  if(AvCurrBlock<=AverageBlock[j+1,n,1])
  {
	if(LogL[1,j] - critFun(j-1) + IndividContribBlock[j+1,n] - Intercept < 0) # Block [1,j] is accepted. Update the CIs.
	{
	  Lower[1:j] = 1
	  Upper[1:j] = pmax(Upper[1:j], j)
	  BlockMax[1] = j - 1
	  break
	}
  }else # The minimally contributing block does not respect the ordering. We need to find another one with slightly higher contrib but respects the order
  {
    ContribRBlock = IndividContribs(y_temp,sigma_temp,LogL,critFun,j+1,n,AvCurrBlock,maxY)
    if(AvCurrBlock<ContribRBlock$LowAv)
    {	  
	  if(LogL[1,j] - critFun(j-1) + ContribRBlock$Contrib - Intercept < 0) # Block [1,j] is accepted. Update the CIs.
	  {
	    Lower[1:j] = 1
	    Upper[1:j] = pmax(Upper[1:j], j)
	    BlockMax[1] = j - 1
	    break
	  }
    }
  }
}
# Treat the case of hypothesis + mu_i=...=mu_n
for(i in 2:(n-1))
{
  AvCurrBlock = sum(y_temp[i:n]/sigma_temp[i:n]^2)/sum(1/sigma_temp[i:n]^2)
  if(AverageBlock[1,i-1,2]<=AvCurrBlock)
  {
	if(LogL[i,n] - critFun(n-i) + IndividContribBlock[1,i-1] - Intercept < 0) # Block [i,n] is accepted. Update the CIs.
	{
	  Lower[i:n] = pmin(Lower[i:n],i)
	  Upper[i:n] = n
	  break
	}
  }else
  {
    ContribLBlock = IndividContribs(y_temp,sigma_temp,LogL,critFun,1,i-1,minY,AvCurrBlock)
    if(AvCurrBlock>ContribLBlock$UpAv)
    {
	  if(LogL[i,n] - critFun(n-i) + ContribLBlock$Contrib - Intercept < 0) # Block [i,n] is accepted. Update the CIs.
	  {
	    Lower[i:n] = pmin(Lower[i:n],i)
	    Upper[i:n] = n
	    break
	  }
	
    }
  }
}
# Treat the case of hypothesis + mu_i=...=mu_j + hypothesis
for(i in 2:(n-2))
{
 for(j in (n-1):(i+1))
 {
   AvCurrBlock = sum(y_temp[i:j]/sigma_temp[i:j]^2)/sum(1/sigma_temp[i:j]^2)
   if(AverageBlock[1,i-1,2]<=AvCurrBlock & AvCurrBlock<=AverageBlock[j+1,n,1])
   {
	if(LogL[i,j] - critFun(j-i) + IndividContribBlock[1,i-1] + IndividContribBlock[j+1,n] - Intercept < 0) # Block [i,j] is accepted. Update the CIs.
	{
	  Lower[i:j] = pmin(Lower[i:j],i)
	  Upper[i:j] = pmax(Upper[i:j], j)
	  BlockMax[i] = j-i
	  break
	}
   }else
  {
    ContribLBlock = IndividContribs(y_temp,sigma_temp,LogL,critFun,1,i-1,minY,AvCurrBlock)
    ContribRBlock = IndividContribs(y_temp,sigma_temp,LogL,critFun,j+1,n,AvCurrBlock,maxY)
    if(AvCurrBlock>ContribLBlock$UpAv & AvCurrBlock<ContribRBlock$LowAv)# THE CONDITION MUST BE ADAPTED TO THE MINIMUM CONTRIB BLOCK if it is not possible to bind with the left part, then move on !
    {
	if(LogL[i,j] - critFun(j-i) + ContribRBlock$Contrib  + ContribLBlock$Contrib - Intercept < 0) # Block [i,j] is accepted. Update the CIs.
	{
	  	Lower[i:j] = pmin(Lower[i:j],i)
	  	Upper[i:j] = pmax(Upper[i:j], j)
	  	BlockMax[i] = j-i
	  	break
	}
    }

  }

 }
}

return(list(Lower = Lower, Upper = Upper, BlockMax = BlockMax))
}

IndividContribs = function(y_temp,sigma_temp,LogL,critFun,K=1,L=length(y_temp),Binf=-Inf,Bsup=Inf, MatReturn = FALSE)
{
IndividContribBlock = matrix(0, nrow = L-K+1, ncol = L-K+1)
AverageBlock = array(dim = c(L-K+1,L-K+1,2)) # a block with a single observation has a mean equal to it.
for(i in 1:(L-K+1))
{
	AverageBlock[i,i,1:2] = y_temp[i+K-1]# I need to keep track of the minimum and maximum average (only these two values are necessary for the PAVA check)
}
### If the Block has several sublocks, I need to keep track of the smallest average and the largest one !
if(L-K+1>=2)
{
	for(j in 2:(L-K+1))
	{
	 for(i in (j-1):1)
	 {
		#LogL[i,j] = sum((y_temp[i:j] - sum(y_temp[i:j]/sigma_temp[i:j]^2)/sum(1/sigma_temp[i:j]^2))^2/sigma_temp[i:j]^2)
		# Minimum contribution is either attained on the whole block (i,j) or at the sum of the minimum of two adjacent sub-blocks of (i,j).
		IndividContribBlock[i,j] = min(0,LogL[i+K-1,j+K-1] - critFun(j-i))
		indMinContrib = 0; MinContribBlock = IndividContribBlock[i,j]
		for(s in 1:(j-i))
		{
		  if(AverageBlock[i,i+s-1,2] <= AverageBlock[i+s,j,1] & AverageBlock[i,i+s-1,1]>=Binf & AverageBlock[i+s,j,2]<=Bsup) # order is respected and no PAVA
		  {
			MinContribBlock = min(c(IndividContribBlock[i,i+s-1]+IndividContribBlock[i+s,j], LogL[i+K-1,j+K-1] - critFun(j-i)))
			if(IndividContribBlock[i,j]>MinContribBlock)
			{
			  IndividContribBlock[i,j] = MinContribBlock
			  indMinContrib = s
			}
		  }
		}
		if(indMinContrib == 0)
		{
		  if(abs(IndividContribBlock[i,j]) < 1e-10) # block contributes as mu_i<...<mu_j
		  {
			minYBlock_ind = which.min(y_temp[(K-1+i):(K-1+j)])
			maxYBlock_ind = which.max(y_temp[(K-1+i):(K-1+j)])
			if(minYBlock_ind==i & maxYBlock_ind==j)
			{
			  AverageBlock[i,j,1:2] = c(y_temp[K-1+minYBlock_ind],y_temp[K-1+maxYBlock_ind])
			}else
			{# The current block must be ignored so that its contribution should never imply nonrejection
			  AverageBlock[i,j,1:2] = c(y_temp[K-1+minYBlock_ind],y_temp[K-1+maxYBlock_ind])
			  IndividContribBlock[i,j] = Inf
			}
		  }else # block contributes as mu_i=...=mu_j
		  {
		  	AverageBlock[i,j,1:2] = sum(y_temp[(K-1+i):(K-1+j)]/sigma_temp[(K-1+i):(K-1+j)]^2)/sum(1/sigma_temp[(K-1+i):(K-1+j)]^2)
		  }
		}else
		{
		  # Update the averages in the minimum contributing block
		  AverageBlock[i,j,1] = AverageBlock[i,i+indMinContrib-1,1]
		  AverageBlock[i,j,2] = AverageBlock[i+indMinContrib,j,2]
		}
	 }
	}
}
if(MatReturn == TRUE)
{
	return(list(IndividContribBlock=IndividContribBlock, AverageBlock = AverageBlock))
}else
{
	return(list(Contrib = IndividContribBlock[1,L-K+1], UpAv = AverageBlock[1,L-K+1,2], LowAv = AverageBlock[1,L-K+1,1]))
}

}

####################################################
# The simple HSD
#library(multcomp)
tukey = function(y,sigma, qq) {
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
# set.seed(16021988)
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
ranks[1,2] = 1+sum(PosPairs[,1]==1)
for(i in 2:(n-1))
{
	ranks[i,1] = i - sum(PosPairs[,2] == i)
	ranks[i,2] = i + sum(PosPairs[,1]==i)
}
ranks[n,1] = n - sum(PosPairs[,2] == n)

return(list(Lower = ranks[,1], Upper = ranks[,2], NbSteps = NBSteps))
}




