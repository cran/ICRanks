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
	//double MeanY = Summation(y, sigma, LowI, UppI) / (UppI - LowI + 1);
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
void RankUpdate(NumericVector& Lower, NumericVector& Upper, int* InqPosi, const int& l, const int& n)
{
	for (int i = 0; i <= InqPosi[0]; i++)
	{
		Lower[i] = 0;
		//Upper[i] = fmax(Upper[i], InqPosi[0]);
		if (InqPosi[0] > Upper[i]) {
			Upper[i] = InqPosi[0];
		}
	}
	int j = 0;
	while (j <= (l - 2))
	{
		//NewLower[(Config[j]+1):Config[j+1]] = pmin(OldLower[(Config[j]+1):Config[j+1]], Config[j]+1)
		//NewUpper[(Config[j]+1):Config[j+1]] = pmax(OldUpper[(Config[j]+1):Config[j+1]], Config[j+1])
		for (int i = InqPosi[j] + 1; i <= InqPosi[j + 1]; i++)
		{
			//Lower[i] = fmin(Lower[i], InqPosi[j] + 1);
			if (InqPosi[j] + 1 < Lower[i]) {
				Lower[i] = InqPosi[j] + 1;
			}
			//Upper[i] = fmax(Upper[i], InqPosi[j + 1]);
			if (InqPosi[j + 1] > Upper[i]) {
				Upper[i] = InqPosi[j + 1];
			}
		}
		j++;
	}
	//NewLower[(Config[(l-1)]+1):n] = pmin(OldLower[(Config[(l-1)]+1):n], Config[l-1]+1)
	//NewUpper[(Config[(l-1)]+1):n] = n
	for (int i = InqPosi[l - 1] + 1; i < n; i++)
	{
		//Lower[i] = fmin(Lower[i], InqPosi[l - 1] + 1);
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
// [[Rcpp::export]]
NumericMatrix PartitioningRankingBlock(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector MinBlock, NumericVector MaxBlock, NumericVector Lower, NumericVector Upper, int n, bool trace)
{

	// Calculate the Likelihood matrix of the blocks.
	unsigned long long int NbEffectTests = 1;
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
			//cout<<LikelihoodMat[i][j]<<" ";
		}
		//cout<<endl;
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
		//bool check_significance = false;
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
						ConfigBase[0] = k + i;
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
						//unsigned long long int c1 = 1;
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
							NbEffectTests++;
							if (Likelihood0 < crit[l])
							{
								RankUpdate(Lower, Upper, InqPosi, l, n);
								break;
							}
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;

							for (c = 0; c<m-1 ; c++)// When c = m-1, we get a binary representation of all ones. This is already tested with the initialization
							{
								//c = (c1 + m/2) % m;
								//c = m-c1-2;
								NbEffectTests++;
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
						}
					
						if (c < m - 1) break;
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
							if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							for (c = 0; c < m-1 ; c++)
							{
								//c = (c1 + m/2) % m;
								//c = m - c1 - 2;
								NbEffectTests++;
								BinaryConfig(c, ConfigCompRight + ConfigBaseLen, l, k + i + 1);// ConfigBase if existed, is already included in ConfigCompRight.

								// Update the length of ConfigCompRight.
								l += ConfigBaseLen;

								// Check if the hypothesis is significant to the actual ranking

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
								if (i > 2)
								{
									int ll;
									mm = PowToN[i - 1];
									// Notice that the upper hypothesis in this sub-partitioning corresponds to testing the block itself. This is already done in the initial CIs.

									NbEffectTests++;
									for (cc = 0; cc<mm-1 ; cc++) // If the left side does not contain any elements, this loop wont start.
									{
										//cc = (cc1 + mm/2) % mm;
										//cc = mm - cc1 - 2;
										//if (c % (m / n) == 0) cout << c << " / " << m << endl;

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
								}

								if (cc < mm - 1) break;

							}
						}
						
						if (c < m - 1) break;
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
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
			//cout<<LikelihoodMat[i][j]<<" ";
		}
		//cout<<endl;
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
	//Rcout<<"Crit 1 = "<<crit[1]<<",Lik = "<<Likelihood0<<"\n";
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
				// Add a test of significance for the current hypothesis

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
					//Lower[1:Config[1]] = 1
					//Upper[1:Config[1]] = pmax(Upper[1:Config[1]], Config[1])
					for (int i = 0; i <= InqPosi[0]; i++)
					{
						Lower[i] = 0;
						Upper[i] = fmax(Upper[i], InqPosi[0]);
					}
					j = 0;
					while (j <= (l - 2))
					{
						//NewLower[(Config[j]+1):Config[j+1]] = pmin(OldLower[(Config[j]+1):Config[j+1]], Config[j]+1)
						//NewUpper[(Config[j]+1):Config[j+1]] = pmax(OldUpper[(Config[j]+1):Config[j+1]], Config[j+1])
						for (int i = InqPosi[j] + 1; i <= InqPosi[j + 1]; i++)
						{
							Lower[i] = fmin(Lower[i], InqPosi[j] + 1);
							Upper[i] = fmax(Upper[i], InqPosi[j + 1]);
						}
						j++;
					}
					//NewLower[(Config[(l-1)]+1):n] = pmin(OldLower[(Config[(l-1)]+1):n], Config[l-1]+1)
					//NewUpper[(Config[(l-1)]+1):n] = n
					for (int i = InqPosi[l - 1] + 1; i<n; i++)
					{
						Lower[i] = fmin(Lower[i], InqPosi[l - 1] + 1);
						Upper[i] = n - 1;
					}
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

