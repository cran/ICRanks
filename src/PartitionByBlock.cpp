#include<cstdlib>
#include<cmath>
#include<string>
#include<ctime>
#include<Rcpp.h>
using namespace Rcpp;

double Summation(NumericVector x, NumericVector sigma, int LowI, int UppI)
{
	// LowI should be by default equal to 0
	// UppI should be by default equal to SampleSize - 1
	double SumRes = x[LowI]/(sigma[LowI]*sigma[LowI]);
	double SumSigma = 1/(sigma[LowI]*sigma[LowI]);
	for (int i = LowI + 1; i <= UppI; i++)
	{
		SumRes += x[i]/(sigma[i]*sigma[i]);
		SumSigma += 1/(sigma[i]*sigma[i]);
	}
	return SumRes/SumSigma;
}

double LogLikelihood(NumericVector y, NumericVector sigma, int LowI, int UppI) // Calculate the LogLiklihood for points with indices between LowI and UppI
{
	// Calculate the average
	double MeanY = Summation(y, sigma, LowI, UppI); // division by the number of elements is included
	double LogLikhood = 0;
	for (int i = LowI; i <= UppI; i++)
	{
		LogLikhood += (y[i] - MeanY)*(y[i] - MeanY) / (sigma[i] * sigma[i]);
	}
	return LogLikhood;
}

void BinaryConfig(unsigned long long int c, int* Config, int& l, const int& Shift)
{
	unsigned long long int residu = c;
	int counter = 0; l = 0;
	while (residu > 1)
	{
		if (residu % 2 == 1)
		{
			Config[l] = counter+Shift;
			l++;
		}

		residu /= 2;
		counter++;
	}
	if (residu == 1)
	{
		Config[l] = counter+Shift;
		l++;
	}
}
void RankUpdate(NumericVector& Lower, NumericVector& Upper, const int* InqPosi, const int& l, const int& n)
{
	for (int i = 0; i <= InqPosi[0]; i++)
	{
		Lower[i] = 0;
		if (InqPosi[0] > Upper[i]) {
			Upper[i] = InqPosi[0];
		}
	}
	int j = 0;
	
	while (j <= (l - 2))
	{
		for (int i = InqPosi[j] + 1; i <= InqPosi[j + 1]; i++)
		{
			if (InqPosi[j] + 1 < Lower[i]) {
				Lower[i] = InqPosi[j] + 1;
			}
			if (InqPosi[j + 1] > Upper[i]) {
				Upper[i] = InqPosi[j + 1];
			}
		}
		j++;
	}
	for (int i = InqPosi[l - 1] + 1; i < n; i++)
	{
		if (InqPosi[l - 1] + 1 < Lower[i]) {
			Lower[i] = InqPosi[l - 1] + 1;
		}
		Upper[i] = n - 1;
	}
}
unsigned long long int binomialCoeff(int n, int k)
{
	if (k>n) return 0;
	unsigned long long int res = 1;

	// Since C(n, k) = C(n, n-k)
	if (k > n - k){
		k = n - k;
	}

	// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res;
}
bool PAVACheck(const NumericVector y_temp, const NumericVector sigma_temp, const int l, const int* InqPosi, const int n)
{
	double MeanBlock1 = Summation(y_temp,sigma_temp,0,InqPosi[0]);
	double MeanBlock2;
	bool CheckPAVA = false;
	int k = 0;
	while (k <= (l - 2))
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
		else
		{
			MeanBlock1 = MeanBlock2;
		}
		k++;
	}
	if(CheckPAVA == false)
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
	}
	return CheckPAVA;
}
void CorrectPermutationsF(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector& Lower, NumericVector& Upper, int* InqPosi, const int l, const int n, const bool EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	NumericVector Lower_temp(n);
	NumericVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s] = y[s-1];
	  			sigma_temp[s] = sigma[s-1];
	  			Lower_temp[s] = s-1;
	  			Upper_temp[s] = s-1;
			}
			y_temp[j-1] = y[j+i-1];
			sigma_temp[j-1] = sigma[j+i-1];
			Lower_temp[j-1] = j+i-1;
			Upper_temp[j-1] = j+i-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck(y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);	
		  			Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
				}
				Lower[j+i-1] = fmin(Lower[j+i-1], Lower_temp[j-1]);
				Upper[j+i-1] = fmax(Upper[j+i-1], Upper_temp[j-1]);
			}
			else
			{
				if(EqSigma) break;
			}
			
		}
		if(j<n-i && EqSigma) break;
	}
}

void CorrectPermutationsB(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector& Lower, NumericVector& Upper, int* InqPosi, const int l, const int n, const bool EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	NumericVector Lower_temp(n);
	NumericVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s-1] = y[s];
	  			sigma_temp[s-1] = sigma[s];
	  			Lower_temp[s-1] = s;
	  			Upper_temp[s-1] = s;
			}
			y_temp[j+i-1] = y[j-1];
			sigma_temp[j+i-1] = sigma[j-1];
			Lower_temp[j+i-1] = j-1;
			Upper_temp[j+i-1] = j-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck( y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s] = fmin(Lower[s], Lower_temp[s-1]);	
		  			Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
				}
				Lower[j-1] = fmin(Lower[j-1], Lower_temp[j+i-1]);
				Upper[j-1] = fmax(Upper[j-1], Upper_temp[j+i-1]);
			}
			else
			{
				if(EqSigma) break;
			}
		}
		if(j<n-i && EqSigma) break;
	}
}
// This function is no longer used
NumericMatrix PartitioningRankingBlock(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector MinBlock, NumericVector MaxBlock, NumericVector Lower, NumericVector Upper, int n, bool trace)
{
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the vector of powers to 2.
	unsigned long long int* PowToN = new unsigned long long int[n];
	PowToN[0] = 1;
	for (int i = 1; i < n; i++)
	{
		PowToN[i] = 2 * PowToN[i-1];
	}

	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true){
			Rcout << "Process ended with trivial confidence intervals.\n";
		} 
	}
	else
	{
		int ConfigBase[2]; // It points out to the extremeties of the tested block.
		int* InqPosi = new int[n - 1];
		int* ConfigCompRight = new int[n - 1]; 
		int ConfigBaseLen = 2;
		int BlockConfig[2];
		for (int i = 0; i < n - 1; i++)
		{
			if(trace == true) {
				Rcout << i << ".";
			}
			if (MaxBlock[i] >= (MinBlock[i] + 1))
			{
				int k = MaxBlock[i];
				while (k >= MinBlock[i] + 1)
				{
					if (i > 0)
					{
						ConfigBase[0] = i - 1;
						ConfigBase[1] = k + i;
						ConfigBaseLen = 2;
					}
					if (i == 0)
					{
						ConfigBase[0] = k;
						ConfigBaseLen = 1;
					}
					if (k + i >= n-1)
					{
						ConfigBase[0] = i - 1;
						ConfigBaseLen = 1;
					}
					if (i == 0 && k + i >= n-1) {
						ConfigBaseLen = 0;
					}
					ConfigCompRight[0] = ConfigBase[0];
					ConfigCompRight[1] = ConfigBase[1];
					BlockConfig[0] = ConfigBase[0];
					BlockConfig[1] = ConfigBase[1];
					if (BlockConfig[1] > n - 1) {
						BlockConfig[1] = n - 1;
					}

					/********************************************************/
					/************ First Case: Left configuration ***********/
					/********************************************************/

					if (k + i >= n - 2) // substract 1 for subscribts diff
					{
						BlockConfig[1] = n - 1;
						// In this case, there is only one side where Partitioing must be done; the left side.
						// The partitions must cover mu_0 to mu_{i-2}
						unsigned long long int m = 2;
						unsigned long long int c = 1;
						if (i > 2)
						{
							int j = 0, l;
							m = PowToN[i-1];
							// Test the top hypothesis corresponding to mu_1=...=mu_{i-1}<mu_i=..=mu_n
							InqPosi[0] = ConfigBase[0];
							l = 1;
							if (ConfigBaseLen == 2) 
							{
								InqPosi[1] = ConfigBase[1];
								l++;
							}
							// Check if PAVA is required
			
							if(PAVACheck(y,sigma,l,InqPosi,n))
							{
								k--;
								continue;
							}
							// Inside each configuration, calculate each group's share in the likelihood
							Likelihood0 = LikelihoodMat[0][InqPosi[0]];
							j = 0;
							while (j <= (l - 2))
							{
								Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
								j++;
							}
							Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
							if (Likelihood0 < crit[l])
							{
								RankUpdate(Lower, Upper, InqPosi, l, n);
								break;
							}
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							// Check the block itself without any additions to the left or the right
							if(PAVACheck(y,sigma,l,InqPosi,n))
							{
								if(LikelihoodMat[InqPosi[0]+1][InqPosi[1]] < crit[n-InqPosi[1]+InqPosi[0]+1])
								{
									RankUpdate(Lower, Upper, InqPosi, l, n);
									break;
								}
							}
							for (c = 0; c<m ; c++)// When c = m-1, we get a binary representation of all ones. This is already tested with the initialization
							{

								BinaryConfig(c, InqPosi, l, 0);
								// Add ConfigBase to InqPosi. Keep in mind that we only use the l+ConfigBaseLen first elements.
								InqPosi[l] = ConfigBase[0];
								InqPosi[l + 1] = ConfigBase[1];
								// Update the length of InqPosi
								l += ConfigBaseLen; // If ConfigBase is empty, the previous attributions will not have effect on the configuration since its length is controled by l.


								/*********************************/
								if(PAVACheck(y,sigma,l,InqPosi,n))
								{
									continue;
								}
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][InqPosi[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									RankUpdate(Lower, Upper, InqPosi, l, n);
									break;// The block was accepted once, and no need to check others.
								}
							}
							if (c < m) break;
						}
					}



					/***************************************************************************/
					/************ Second Case: Right and possibly Left configuration ***********/
					/***************************************************************************/

					else
					{
						unsigned long long int m = 2;
						unsigned long long int c = 1; // I need the counter to check if I was able to find a hypothesis that was not rejected.
						//unsigned long long int c1 = 1;
						if (n - k - i > 2)
						{
							
							// Test the top level in this sub-partition
							int j = 0, l;
							m = PowToN[n - k - i - 2]; // The right part contains n-i-k-1 centers.
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							for (c = 0; c < m ; c++)
							{

								BinaryConfig(c, ConfigCompRight + ConfigBaseLen, l, k + i + 1);// ConfigBase if existed, is already included in ConfigCompRight.
								// Update the length of ConfigCompRight.
								l += ConfigBaseLen;

								/*********************************/
								if(PAVACheck(y,sigma,l,ConfigCompRight,n))
								{
									continue;
								}
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][ConfigCompRight[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[ConfigCompRight[j] + 1][ConfigCompRight[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[ConfigCompRight[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									
									RankUpdate(Lower, Upper, ConfigCompRight, l, n);
									break;// The block was accepted once, and no need to check others.
								}


								// Add a configuration to the left of ConfigBase in case there is suffcient place.
								/***************Third case: a configuration to the left and a configuration to the right*******************/

								// The top hypothesis in this sub-partitioning was just tested.
								unsigned long long int cc = 1;
								unsigned long long int mm = 2; // Initial values here are necessary to cancel the break hereafter in case, we do not enter the loop over cc.
								if (i > 1)
								{
									int ll;
									mm = PowToN[i - 1];
									// Notice that the upper hypothesis in this sub-partitioning corresponds to testing the block itself. This is already done in the initial CIs.
									for (cc = 0; cc<mm ; cc++) // If the left side does not contain any elements, this loop wont start.
									{

										BinaryConfig(cc, InqPosi, ll, 0);
										// Add now the configuration in the right side
										for (int iter = 0; iter < l; iter++) InqPosi[iter + ll] = ConfigCompRight[iter];
										// Update the length of InqPosi
										ll += l;

										/*********************************/
										if(PAVACheck(y,sigma,ll,InqPosi,n))
										{
											continue;
										}
										// Inside each configuration, calculate each group's share in the likelihood
										Likelihood0 = LikelihoodMat[0][InqPosi[0]];
										j = 0;
										while (j <= (ll - 2))
										{
											Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
											j++;
										}
										Likelihood0 += LikelihoodMat[InqPosi[ll - 1] + 1][n - 1];

										// Update the ranking
										if (Likelihood0 < crit[ll])
										{
											RankUpdate(Lower, Upper, InqPosi, ll, n);
											break;// The block was accepted once, and no need to check others.
										}
									}
									if (cc < mm) break;
								}

							}
							if (c < m) break;
						}
					}
					k--;
				}

			}

		}

		

		delete[] InqPosi;
		delete[] ConfigCompRight;

	}

	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	delete[] PowToN;
	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	return CIs;
}
void UnrankCombin(int* S, unsigned long long int m, int k, unsigned long long int** CnkMat)
{
	int i = k - 1;
	while (i >= 0)
	{
		int l = i;
		while (CnkMat[l][i + 1] <= m)
		{
			l++;
		}
		S[i]= l - 1;
		m = m - CnkMat[l - 1][i + 1];
		i--;
	}

}

// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevel(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace)
{
	// Detect equal-sigma for faster implementation
	bool EqSigma = true;
	for(int i=0;i<n-1;i++)
	{
		if(sigma[i] != sigma[i+1])
		{
			EqSigma = false;
			break;
		}
	}
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	NumericVector Lower(n);
	NumericVector Upper(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
	}
	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true) {
			Rcout << "Process ended with trivial confidence intervals.\n";
		}
	}
	else
	{
		int* InqPosi = new int[n];
		if(trace == true) {
				Rcout<<"Processed levels:";
		}
		for (int l = 1; l <= n - 2; l++)
		{
			//			int l = 3; int i = 8;
			if(trace == true) Rcout<<l<<".";
			unsigned long long int m = CnkMat[n - 1][l];
			for (unsigned long long int c = 0; c<m; c++)
			{
				UnrankCombin(InqPosi, c, l, CnkMat);
				// Inside each configuration, calculate each group's share in the likelihood
				Likelihood0 = LikelihoodMat[0][InqPosi[0]];
				int j = 0;
				while (j <= (l - 2))
				{
					Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
					j++;
				}
				Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
				// Update the ranking
				if (Likelihood0<crit[l])
				{
					RankUpdate(Lower, Upper, InqPosi, l, n);
					CorrectPermutationsF(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma);
					CorrectPermutationsB(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma);
					
				}


			}
			
		}
		delete[] InqPosi;
	}
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	return CIs;
}
// [[Rcpp::export]]
NumericMatrix PartitioningRankingBlockCorrectOrder(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector MinBlock, NumericVector MaxBlock, NumericVector Lower, NumericVector Upper, int n, bool trace)
{
	// Calculate the Likelihood matrix of the blocks
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the vector of powers to 2.
	unsigned long long int* PowToN = new unsigned long long int[n];
	PowToN[0] = 1;
	for (int i = 1; i < n; i++)
	{
		PowToN[i] = 2 * PowToN[i-1];
	}

	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true){
			Rcout << "Process ended with trivial confidence intervals.\n";
		} 
	}
	else
	{
		int ConfigBase[2]; // It points out to the extremeties of the tested block.
		int* InqPosi = new int[n - 1];
		int* ConfigCompRight = new int[n - 1]; 
		int ConfigBaseLen = 2;
		int BlockConfig[2];
		for (int i = 0; i < n - 1; i++)
		{
			if(trace == true) {
				Rcout << i << ".";
			}
			if (MaxBlock[i] >= (MinBlock[i] + 1))
			{
				int k = MaxBlock[i];
				while (k >= MinBlock[i] + 1)
				{
					if (i > 0)
					{
						ConfigBase[0] = i - 1;
						ConfigBase[1] = k + i;
						ConfigBaseLen = 2;
					}
					if (i == 0)
					{
						ConfigBase[0] = k;
						ConfigBaseLen = 1;
					}
					if (k + i >= n-1)
					{
						ConfigBase[0] = i - 1;
						ConfigBaseLen = 1;
					}
					if (i == 0 && k + i >= n-1) {
						ConfigBaseLen = 0;
					}
					ConfigCompRight[0] = ConfigBase[0];
					ConfigCompRight[1] = ConfigBase[1];
					BlockConfig[0] = ConfigBase[0];
					BlockConfig[1] = ConfigBase[1];
					if (BlockConfig[1] > n - 1) {
						BlockConfig[1] = n - 1;
					}

					/********************************************************/
					/************ First Case: Left configuration ***********/
					/********************************************************/

					if (k + i >= n - 2) // substract 1 for subscribts diff
					{
						BlockConfig[1] = n - 1;
						// In this case, there is only one side where Partitioing must be done; the left side.
						// The partitions must cover mu_0 to mu_{i-2}
						unsigned long long int m = 2;
						unsigned long long int c = 1;
						if (i > 2)
						{
							int j = 0, l;
							m = PowToN[i-1];
							// Test the top hypothesis corresponding to mu_1=...=mu_{i-1}<mu_i=..=mu_n
							InqPosi[0] = ConfigBase[0];
							l = 1;
							if (ConfigBaseLen == 2) 
							{
								InqPosi[1] = ConfigBase[1];
								l++;
							}
							// Inside each configuration, calculate each group's share in the likelihood
							Likelihood0 = LikelihoodMat[0][InqPosi[0]];
							j = 0;
							while (j <= (l - 2))
							{
								Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
								j++;
							}
							Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
							if (Likelihood0 < crit[l])
							{
								
								RankUpdate(Lower, Upper, InqPosi, l, n);
								
								break;
							}
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							// Check the block itself without any additions to the left or the right
							for (c = 0; c<m ; c++)
							{
								BinaryConfig(c, InqPosi, l, 0);
								// Add ConfigBase to InqPosi. Keep in mind that we only use the l+ConfigBaseLen first elements.
								InqPosi[l] = ConfigBase[0];
								InqPosi[l + 1] = ConfigBase[1];
								// Update the length of InqPosi
								l += ConfigBaseLen; // If ConfigBase is empty, the previous attributions will not have effect on the configuration since its length is controled by l.


								/*********************************/
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][InqPosi[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									RankUpdate(Lower, Upper, InqPosi, l, n);
									break;// The block was accepted once, and no need to check others.
								}
							}
							if (c < m) break;
						}
					}



					/***************************************************************************/
					/************ Second Case: Right and possibly Left configuration ***********/
					/***************************************************************************/

					else
					{
						unsigned long long int m = 2;
						unsigned long long int c = 1; // I need the counter to check if I was able to find a hypothesis that was not rejected.
						//unsigned long long int c1 = 1;
						if (n - k - i > 2)
						{
							
							// Test the top level in this sub-partition
							int j = 0, l;
							m = PowToN[n - k - i - 2]; // The right part contains n-i-k-1 centers.
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							for (c = 0; c < m ; c++)
							{
								BinaryConfig(c, ConfigCompRight + ConfigBaseLen, l, k + i + 1);// ConfigBase if existed, is already included in ConfigCompRight.
								// Update the length of ConfigCompRight.
								l += ConfigBaseLen;
								/*********************************/
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][ConfigCompRight[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[ConfigCompRight[j] + 1][ConfigCompRight[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[ConfigCompRight[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									RankUpdate(Lower, Upper, ConfigCompRight, l, n);
									break;// The block was accepted once, and no need to check others.
								}


								// Add a configuration to the left of ConfigBase in case there is suffcient place.
								/***************Third case: a configuration to the left and a configuration to the right*******************/

								// The top hypothesis in this sub-partitioning was just tested.
								unsigned long long int cc = 1;
								//unsigned long long int cc1 = 1;
								unsigned long long int mm = 2; // Initial values here are necessary to cancel the break hereafter in case, we do not enter the loop over cc.
								if (i > 1)
								{
									int ll;
									mm = PowToN[i - 1];
									// Notice that the upper hypothesis in this sub-partitioning corresponds to testing the block itself. This is already done in the initial CIs.
									for (cc = 0; cc<mm ; cc++) // If the left side does not contain any elements, this loop wont start.
									{

										BinaryConfig(cc, InqPosi, ll, 0);
										// Add now the configuration in the right side
										for (int iter = 0; iter < l; iter++) InqPosi[iter + ll] = ConfigCompRight[iter];
										// Update the length of InqPosi
										ll += l;

										/*********************************/
										// Inside each configuration, calculate each group's share in the likelihood
										Likelihood0 = LikelihoodMat[0][InqPosi[0]];
										j = 0;
										while (j <= (ll - 2))
										{
											Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
											j++;
										}
										Likelihood0 += LikelihoodMat[InqPosi[ll - 1] + 1][n - 1];

										// Update the ranking
										if (Likelihood0 < crit[ll])
										{
											RankUpdate(Lower, Upper, InqPosi, ll, n);
											break;// The block was accepted once, and no need to check others.
										}
									}
									if (cc < mm) break;
								}

							}
							if (c < m) break;
						}
					}
					k--;
				}

			}

		}

		

		delete[] InqPosi;
		delete[] ConfigCompRight;

	}

	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	delete[] PowToN;
	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	return CIs;
}
