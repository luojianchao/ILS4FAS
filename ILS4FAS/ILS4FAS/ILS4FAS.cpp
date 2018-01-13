#include <iostream>
#include <fstream>
#include <time.h>
#include <cassert>
#include<string>
#include "data.h"
#define INF 100000
using namespace::std;

int max(int a, int b)
{
	return(a > b ? a : b);
}

int fnComputeT(int **nArrIniPop1, int nl, int nc, int *nJs);
void fnTranFir(int *M, int nT);
void fnComputeTime(int nC, int *pmatrix, int* tmatrix, int &finaltime);
void fnCompleteSequence(int **nArrIniPop1, int **nArrIniFs1, int nL, int nC);
void fnCompleteSequence1(int* nArrIniPop1, int *nArrIniFs1, int nC);
void fnRandomInitialize(int **nArrIniPop, int nL, int nC);
void fnLocalSearch(int **nArrIniPop1, int **nArrIniFs1, const int nMaxKeep, int *nArrFinial1, int nL);
void fnSearchLocal(int *nPmatrix, int nJss, int *nPmatrix1, int nC);
void fnEqual21L(int l, int m, int **arrL, int *arrR);
void fnEqual12(int l, int m, int **arrR, int *arrL);
void fnEDAPerturbation(int **nArrIniPop1, int nL, int nC, int **nArrIniPop2);
void fnBigDisturbance(int **nArrIniPop1, int nL, int nC, int **nArrIniPop2);
void fnRoundBet(int **nArrIniPop1, int **nArrIniFs1, int *nArrFinial1, int **nArrIniPop2, int **nArrIniFs2, int *nArrFinial2, int nIniPopNum);
void fnParetoDominanceAcceptanceCriterion(int **nArrIniPop1, int **nArrIniFs1, int *nArrFinial1, int **nArrIniPop2, int **nArrIniFs2, int *nArrFinial2, int nIniPopNum);
void fnNew2(int l, int m, int ***arr);
template <class T> T fnMax(T * nArr, T nL);
void fnDel2(int l, int m, int ***arr);
void fnIterativeLocalSearch(int nArrP[], int nArrT[], int &nGlobalBest);
void fnFindMin(int *arr, int nL, int& pos);

void fnFindMin(int *arr, int nL, int& pos)
{
	/*input: arr, nL,output: the position of the minimum element pos*/
	int nMin = INF;
	for (int i = 0; i < nL; i++)
	{
		if (arr[i] < nMin)
		{
			nMin = arr[i];
			pos = i;
		}
	}
}

void fnIterativeLocalSearch(int nArrP[], int nArrT[],  int &nGlobalBest)
{
	ofstream fileB;
	fileB.open("fileB.txt", ios::app);
	int nSearchJs = 0;
	int nMaxKeep = 10;
	int  **nArrIniPop1 = new int*[c_nPopNum];
	int  **nArrIniFs1 = new int *[c_nPopNum];
	int  **nArrIniPop2 = new int*[c_nPopNum];
	int  **nArrIniFs2 = new int *[c_nPopNum];
	int nArrFinial1[c_nPopNum];
	int nArrFinial2[c_nPopNum];
	for (int i = 0; i < c_nPopNum; i++)
	{
		nArrIniPop1[i] = new int[c_nGeneCam];
		nArrIniFs1[i] = new int[c_nGeneCam];
		nArrIniPop2[i] = new int[c_nGeneCam];
		nArrIniFs2[i] = new int[c_nGeneCam];
	}
	fnRandomInitialize(nArrIniPop1, c_nPopNum, c_nGeneCam);
	fnCompleteSequence(nArrIniPop1, nArrIniFs1, c_nPopNum, c_nGeneCam);
	fnLocalSearch(nArrIniPop1, nArrIniFs1, c_nMaxKeep, nArrFinial1, c_nPopNum);
	int nPos;
	fnFindMin(nArrFinial1, c_nPopNum, nPos);
	for (int i = 0; i < c_nGeneCam; i++)
	{
		nArrP[i] = nArrIniPop1[nPos][i];
		nArrT[i] = nArrIniFs1[nPos][i];
	}
	fileB << nArrFinial1[nPos] << '\t';
	nGlobalBest = nArrFinial1[nPos];
	int nFinal = 0;
	int generationCount = 0;
	while (generationCount < c_nMaxGeneration)
	{
		generationCount++;
		fnBigDisturbance(nArrIniPop1, c_nPopNum, c_nGeneCam, nArrIniPop2);
		fnCompleteSequence(nArrIniPop2, nArrIniFs2, c_nPopNum, c_nGeneCam);
		fnLocalSearch(nArrIniPop2, nArrIniFs2, c_nMaxKeep, nArrFinial2, c_nPopNum);
		int nPos;
		fnFindMin(nArrFinial2, c_nPopNum, nPos);
		if (nArrFinial2[nPos] < nGlobalBest)
		{
			for (int i = 0; i < c_nGeneCam; i++)
			{
				nArrP[i] = nArrIniPop2[nPos][i];
				nArrT[i] = nArrIniFs2[nPos][i];
			}
			nGlobalBest = nArrFinial2[nPos];
		}
		fileB << nGlobalBest << '\t';
		fnRoundBet(nArrIniPop1, nArrIniFs1, nArrFinial1, nArrIniPop2, nArrIniFs2, nArrFinial2, c_nPopNum);
	}
	fileB << endl;
	fileB.close();
	for (int i = 0; i < c_nPopNum; i++)
	{
		delete[] nArrIniPop1[i];
		delete[] nArrIniFs1[i];
		delete[] nArrIniPop2[i];
		delete[] nArrIniFs2[i];
	}
	delete[] nArrIniPop1;
	delete[] nArrIniFs1;
	delete[] nArrIniPop2;
	delete[] nArrIniFs2;
}

void fnNew2(int l, int m, int ***arr)
{
	*arr = new int*[l];
	for (int i = 0; i < l; i++)
	{
		(*arr)[i] = new int[m];
	}
}

template <class T> T fnMax(T * nArr, T nL)
{
	T tmax = -1;
	for (int i = 0; i < nL; i++)
	{
		if (nArr[i]>tmax)
		{
			tmax = nArr[i];
		}
	}
	return tmax;
}

void fnDel2(int l, int m, int ***arr)
{
	for (int i = 0; i < l; i++)
	{
		delete[](*arr)[i];// = new int[m];
	}
	delete[](*arr);
}

void fnRoundBet(int **nArrIniPop1, int **nArrIniFs1, int *nArrFinial1, int **nArrIniPop2, int **nArrIniFs2, int *nArrFinial2, int nIniPopNum)
{
	int **nArrIniPop, **nArrIniFs, *nArrFinial;
	fnNew2(nIniPopNum, c_nGeneCam, &nArrIniPop);
	fnNew2(nIniPopNum, c_nGeneCam, &nArrIniFs);
	nArrFinial = new int[nIniPopNum];
	int *dFitness = new int[2 * nIniPopNum];
	double *dSelecP = new double[2 * nIniPopNum];
	double *dProSel = new double[nIniPopNum];
	int maxF = max(fnMax(nArrFinial1, nIniPopNum), fnMax(nArrFinial2, nIniPopNum));
	int sumF = 0;
	for (int i = 0; i < nIniPopNum; i++)
	{
		dFitness[i] = maxF - nArrFinial1[i];
		sumF += dFitness[i];
	}
	for (int i = nIniPopNum; i < 2 * nIniPopNum; i++)
	{
		dFitness[i] = maxF - nArrFinial2[i - nIniPopNum];
		sumF += dFitness[i];
	}
	for (int i = 0; i < 2 * nIniPopNum; i++)
	{
		dSelecP[i] = dFitness[i] / (1.0*sumF);
	}
	for (int i = 0; i < nIniPopNum; i++)
	{
		dProSel[i] = (rand() % 10000) / 10000.0;
	}
	double *dS = new double[2 * nIniPopNum + 1];
	dS[0] = 0;
	for (int i = 1; i < 2 * nIniPopNum + 1; i++)
	{
		dS[i] = dS[i - 1] + dSelecP[i];
	}
	for (int i = 1; i < nIniPopNum + 1; i++)
	{
		for (int j = 1; j < 2 * nIniPopNum + 1; j++)
		{
			if (dProSel[i - 1] >= dS[j - 1] && dProSel[i - 1] <= dS[j])
			{
				if (j - 1 < nIniPopNum)
				{
					for (int k = 0; k < c_nGeneCam; k++)
					{
						nArrIniPop[i - 1][k] = nArrIniPop1[j - 1][k];
						nArrIniFs[i - 1][k] = nArrIniFs1[j - 1][k];
					}
					nArrFinial[i - 1] = nArrFinial1[j - 1];
				}
				else
				{
					for (int k = 0; k < c_nGeneCam; k++)
					{
						nArrIniPop[i - 1][k] = nArrIniPop2[j - 1 - nIniPopNum][k];
						nArrIniFs[i - 1][k] = nArrIniFs2[j - 1 - nIniPopNum][k];
					}
					nArrFinial[i - 1] = nArrFinial2[j - 1 - nIniPopNum];
				}
				break;
			}
			else if (dProSel[i - 1] > dS[2 * nIniPopNum])
			{
				for (int k = 0; k < c_nGeneCam; k++)
				{
					nArrIniPop[i - 1][k] = nArrIniPop2[nIniPopNum - 1][k];
					nArrIniFs[i - 1][k] = nArrIniFs2[nIniPopNum - 1][k];
				}
				nArrFinial[i - 1] = nArrFinial2[nIniPopNum - 1];
				break;
			}
		}
	}
	for (int i = 0; i < nIniPopNum; i++)
	{
		nArrFinial1[i] = nArrFinial[i];
		for (int j = 0; j < c_nGeneCam; j++)
		{
			nArrIniPop1[i][j] = nArrIniPop[i][j];
			nArrIniFs1[i][j] = nArrIniFs[i][j];
		}
	}
	delete[] dFitness;
	delete[] dSelecP;
	delete[] dProSel;
	fnDel2(nIniPopNum, c_nGeneCam, &nArrIniPop);
	fnDel2(nIniPopNum, c_nGeneCam, &nArrIniFs);
	delete[] nArrFinial;
}

void fnParetoDominanceAcceptanceCriterion(int **nArrIniPop1, int **nArrIniFs1, int *nArrFinial1, int **nArrIniPop2, int **nArrIniFs2, int *nArrFinial2, int nIniPopNum)
{
	// nArrIniPop1 is the set of chromosomes, nArrIniFs1 is the set of transition sequences according to nArrIniPop1, nArrFinial1 is the finish time of the chromosomes in nArrIniPop1
	// nArrIniPop2 is the set selected chromosomes by PD-AC on nArrIniPop1, nArrIniFs2 is the set of transition sequences according to nArrIniPop2, nArrFinial2 is the finish time of the chromosomes in nArrIniPop2
	// nIniPopNum is the number of chromosomes in nArrIniPop1
	// please refer to IV.B 3) to finish this part
}

void fnSearchLocal(int *nPmatrix, int nJss, int *nPmatrix1, int nC)
{
	int length = (int)(nC / 5);
	int rnum = (int)(nC*0.05);
	int *maxp = new int[rnum];
	int *maxo = new int[rnum];
	for (int i = 0; i < nC; i++)
	{
		nPmatrix1[i] = nPmatrix[i];
	}
	for (int i = 0; i < rnum; i++)
	{
		maxp[i] = 0;
	}
	for (int j = 0; j < rnum; j++)
	{
		int temp = length*(nJss - 1) + (int)(rand() % 100 / 100.0*(length - 1));
		if (j == 0)
		{
			maxp[j] = temp;
		}
		else if (j > 0)
		{
			int flag = 1;
			while (flag == 1)
			{
				flag = 0;
				for (int i = 0; i <= j - 1; i++)
				{
					if (temp == maxp[i])
					{
						flag = 1;
						break;
					}
				}
				if (flag == 1)
				{
					temp = length*(nJss - 1) + (int)(rand() % 100 / 100.0*(length - 1));
				}
			}
			maxp[j] = temp;
		}
	}
	for (int i = 0; i < rnum; i++)
	{
		maxo[i] = nPmatrix[maxp[i]];
	}
	int a, b, temp;
	for (int i = 0; i < rnum; i++)
	{
		a = rand() % rnum;
		b = rand() % rnum;
		temp = maxo[a];
		maxo[a] = maxo[b];
		maxo[b] = temp;
	}
	for (int i = 0; i < rnum; i++)
	{
		nPmatrix1[maxp[i]] = maxo[i];
	}
	delete[] maxo;
	delete[] maxp;
}

void fnEqual21L(int l, int m, int **arrL, int *arrR)
{
	for (int i = 0; i < m; i++)
	{
		arrR[i] = arrL[l][i];
	}
}

void fnEqual12(int l, int m, int **arrR, int *arrL)
{
	for (int i = 0; i < m; i++)
	{
		arrR[l][i] = arrL[i];
	}
}

void fnEDAPerturbation(int **nArrIniPop1, int nL, int nC, int **nArrIniPop2)
{
	// perturbation operations on chromosomes in nArrIniPop1, nL is the number of chromosomes in  nArrIniPop1, nC is the number of genes of a chromosome, nArrIniPop2 is the chromosomes after perturbation operation on nArrIniPop1
	// please refer to IV.B 1) to finish this part
}

void fnBigDisturbance(int **nArrIniPop1, int nL, int nC, int **nArrIniPop2)
{
	int rNum = (int)(c_nGeneCam*0.1);
	int *maxp = new int[rNum];
	int *maxo = new int[rNum];
	for (int j = 0; j < rNum; j++)
	{
		maxp[j] = 0;
		maxo[j] = 0;
	}
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nC; j++)
		{
			nArrIniPop2[i][j] = nArrIniPop1[i][j];
		}
		for (int j = 0; j < rNum; j++)
		{
			maxp[j] = 0;
		}
		for (int j = 0; j < rNum; j++)
		{
			int temp = (int)(rand() % 100 / 100.0*(nC - 1));
			if (j == 0)
			{
				maxp[j] = temp;
			}
			else if (j > 0)
			{
				int flag = 1;
				while (flag == 1)
				{
					flag = 0;
					for (int kk = 0; kk < j; kk++)
					{
						if (temp == maxp[kk])
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
					{
						temp = (int)(rand() % 100 / 100.0*(nC - 1));
					}
				}
				maxp[j] = temp;
			}
		}
		for (int k = 0; k < rNum; k++)
		{
			maxo[k] = nArrIniPop1[i][maxp[k]];
		}
		int a, b, temp;
		for (int k = 0; k < rNum; k++)
		{
			a = rand() % rNum;
			b = rand() % rNum;
			temp = maxo[a];
			maxo[a] = maxo[b];
			maxo[b] = temp;
		}
		for (int k = 0; k < rNum; k++)
		{
			nArrIniPop2[i][maxp[k]] = maxo[k];
		}
	}
	delete[] maxp;
	delete[] maxo;
}

void fnLocalSearch(int **nArrIniPop1, int **nArrIniFs1, const int nMaxKeep, int *nArrFinial1, int nL)
{
	int nPmatrix[c_nGeneCam];
	int nTmatrix[c_nGeneCam];
	int nPmatrix1[c_nGeneCam];
	int nTmatrix1[c_nGeneCam];
	for (int i = 0; i < nL; i++) nArrFinial1[i] = 0;
	for (int i = 0; i < nL; i++)
	{
		int nJsKeep = 0;
		int nJss = 0;
		fnEqual21L(i, c_nGeneCam, nArrIniPop1, nPmatrix);
		fnEqual21L(i, c_nGeneCam, nArrIniFs1, nTmatrix);
		fnComputeTime(c_nGeneCam, nPmatrix, nTmatrix, nArrFinial1[i]);
		while (nJsKeep < nMaxKeep)
		{
			if (nJss == 5)
			{
				nJss = 0;
			}
			nJss += 1;
			fnEqual21L(i, c_nGeneCam, nArrIniPop1, nPmatrix);
			fnSearchLocal(nPmatrix, nJss, nPmatrix1, c_nGeneCam);
			fnCompleteSequence1(nPmatrix1, nTmatrix1, c_nGeneCam);
			int nFinal1;
			fnComputeTime(c_nGeneCam, nPmatrix1, nTmatrix1, nFinal1);
			if (nFinal1 < nArrFinial1[i])
			{
				nArrFinial1[i] = nFinal1;
				nJsKeep = 0;
				fnEqual12(i, c_nGeneCam, nArrIniPop1, nPmatrix1);
				fnEqual12(i, c_nGeneCam, nArrIniFs1, nTmatrix1);
			}
			else
			{
				nJsKeep += 1;
			}
		}
	}
}

void fnRandomInitialize(int **nArrIniPop, int nL, int nC)
{
	int nTotal[c_nProNum];
	int nJs[c_nProNum];
	for (int i = 0; i < c_nProNum; i++)
	{
		if (i < c_nP1)
		{
			nTotal[i] = c_nPlaceNum1;
		}
		else if (c_nP1 <= i && i < c_nP1 + c_nP2)
		{
			nTotal[i] = c_nPlaceNum2;
		}
		else if (c_nP1 <= c_nP1 + c_nP2 && i < c_nP1 + c_nP2 + c_nP3)
		{
			nTotal[i] = c_nPlaceNum3;
		}
		else
		{
			nTotal[i] = c_nPlaceNum4;
		}
	}
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nC; j++)
		{
			nArrIniPop[i][j] = 1 + rand() % c_nProNum;
		}
		for (int j = 0; j < c_nProNum; j++)
		{
			nJs[j] = 0;
		}
		for (int j = 0; j < nC; j++)
		{
			nJs[nArrIniPop[i][j] - 1] += 1;
			if (nJs[nArrIniPop[i][j] - 1] > nTotal[nArrIniPop[i][j] - 1])
			{
				nArrIniPop[i][j] = 0;
			}
		}
		for (int j = 0; j < c_nProNum; j++)
		{
			if (nJs[j] < nTotal[j])
			{
				for (int t = 0; t < nC; t++)
				{
					if (nArrIniPop[i][t] == 0)
					{
						nArrIniPop[i][t] = j + 1;
						nJs[j] += 1;
						if (nJs[j] == nTotal[j])
						{
							break;
						}
					}
				}
			}
		}
	}
}

void fnComputeTime(int nC, int *pmatrix, int* tmatrix, int &finaltime)
{
	int nProFinish[c_nProNum][7] = { 0 };
	int tNow = 0;
	int nJs = 0;
	int preT = 0, start;
	for (int i = 0; i < nC; i++)
	{
		int hh = pmatrix[i] - 1;
		int lh = nProFinish[hh][6];
		if (lh == 0)
		{
			preT = 0;
		}
		else
		{
			preT = nProFinish[hh][lh - 1];
		}
		start = max(tNow, preT);
		tNow = start;
		if (pmatrix[i] <= c_nP1)
		{
			nProFinish[hh][lh] = start + c_nTime[lh];
		}
		else if (c_nP1 < pmatrix[i] && pmatrix[i] <= c_nP1 + c_nP2)
		{
			nProFinish[hh][lh] = start + c_nTime[c_nPlaceNum1 + lh];
		}
		else if (c_nP1 + c_nP2 < pmatrix[i] && pmatrix[i] <= c_nP1 + c_nP2 + c_nP3)
		{
			nProFinish[hh][lh] = start + c_nTime[c_nPlaceNum1 + c_nPlaceNum2 + lh];
		}
		else
		{
			nProFinish[hh][lh] = start + c_nTime[c_nPlaceNum1 + c_nPlaceNum2 + c_nPlaceNum3 + lh];
		}
		nProFinish[hh][6] += 1;
	}
	finaltime = 0;
	for (int i = 0; i < c_nProNum; i++)
	{
		finaltime = max(finaltime, nProFinish[i][nProFinish[i][6] - 1]);
	}
}

int fnComputeT(int **nArrIniPop1, int nl, int nc, int *nJs)
{
	int nT1 = 0;
	int nT = nJs[nArrIniPop1[nl][nc] - 1] + 1;
	if (nArrIniPop1[nl][nc] <= c_nP1)
	{
		nT1 = nT;
	}
	else if (c_nP1 < nArrIniPop1[nl][nc] && nArrIniPop1[nl][nc] <= c_nP1 + c_nP2)
	{
		nT1 = nT + c_nPlaceNum1;
	}
	else if (c_nP1 + c_nP2 < nArrIniPop1[nl][nc] && nArrIniPop1[nl][nc] <= c_nP1 + c_nP2 + c_nP3)
	{
		nT1 = nT + c_nPlaceNum1 + c_nPlaceNum2;
	}
	else
	{
		nT1 = nT + c_nPlaceNum1 + c_nPlaceNum2 + c_nPlaceNum3;
	}
	return nT1;
}

int fnComputeT1(int *nArrIniPop1, int nc, int *nJs)
{
	int nT1 = 0;
	int nT = nJs[nArrIniPop1[nc] - 1] + 1;
	if (nArrIniPop1[nc] <= c_nP1)
	{
		nT1 = nT;
	}
	else if (c_nP1 < nArrIniPop1[nc] && nArrIniPop1[nc] <= c_nP1 + c_nP2)
	{
		nT1 = nT + c_nPlaceNum1;
	}
	else if (c_nP1 + c_nP2 < nArrIniPop1[nc] && nArrIniPop1[nc] <= c_nP1 + c_nP2 + c_nP3)
	{
		nT1 = nT + c_nPlaceNum1 + c_nPlaceNum2;
	}
	else
	{
		nT1 = nT + c_nPlaceNum1 + c_nPlaceNum2 + c_nPlaceNum3;
	}
	return nT1;
}

void fnTranFir(int *M, int nT)
{
	for (int i = 0; i < c_nPlace; i++) //Ms=M+petrinet(:,T_num)';//%M[t'>Ms
	{
		M[i] = M[i] + c_nPetriNet[i][nT - 1];
	}
}

void fnCompleteSequence(int** nArrIniPop1, int **nArrIniFs1, int nL, int nC)
{
	// repair the chromosomes in nArrIniPop1, nArrIniFs1 is the transition sequences according to nArrIniPop1,nL is the number of chromosomes, nC is the number of genes of a chromosome
	// please refer to repair algorithm to finish this part
}

void fnCompleteSequence1(int* nArrIniPop1, int *nArrIniFs1, int nC)
{
	// repair a chromosome in nArrIniPop1, nArrIniFs1 is the transition sequence according to nArrIniPop1, nC is the number of genes of a chromosome
	// please refer to repair algorithm to finish this part
}

void  checkSequence(int** nArrIniPop1, int** nArrIniFs1)
{
	int nJS[c_nProNum];
	for (int i = 0; i < c_nPopNum; i++)
	{
		for (int j = 0; j < c_nProNum; j++)
		{
			nJS[j] = 0;
		}
		for (int j = 0; j < c_nGeneCam; j++)
		{
			nJS[nArrIniPop1[i][j] - 1]++;
		}
		bool checkResult = true;
		for (int j = 0; j < c_nProNum; j++)
		{
			if (j < c_nP1)
			{
				if (nJS[j] != c_nPlaceNum1)
				{
					checkResult = false;
					cout << "wrong generation" << endl;
				}
			}
			else if (c_nP1 <= j && j < c_nP1 + c_nP2)
			{
				if (nJS[j] != c_nPlaceNum2)
				{
					checkResult = false;
					cout << "wrong generation" << endl;
				}
			}
			else if (c_nP1 <= c_nP1 + c_nP2 && j < c_nP1 + c_nP2 + c_nP3)
			{
				if (nJS[j] != c_nPlaceNum3)
				{
					checkResult = false;
					cout << "wrong generation" << endl;
				}
			}
			else
			{
				if (nJS[j] != c_nPlaceNum4)
				{
					checkResult = false;
					cout << "wrong generation" << endl;
				}
			}
		}
		cout << "checkResult is: " << checkResult << endl;
	}
	for (int i = 0; i < c_nPopNum; i++)
	{
		int  M[c_nPlace];
		for (int j = 0; j < c_nPlace; j++)
		{
			M[j] = c_aInitialMark[j];
		}
		for (int j = 0; j < c_nGeneCam; j++)
		{
			fnTranFir(M, nArrIniFs1[i][j]);
		}
		for (int j = 0; j < c_nPlace; j++)
		{
			if (M[j] != c_aFinalMark[j])
			{
				cout << "transition sequence is wrong!" << endl;
			}
		}
		int finalTime;
		fnComputeTime(c_nGeneCam, nArrIniPop1[i], nArrIniFs1[i], finalTime);
		cout << finalTime << '\t';
	}
	cout << endl;
}

void changeCode2Transition()
{
	string changeTransition[c_nTransition] = {
		"t11", "t12", "t13", "t14", "t15", "t21", "t22", "t31", "t32", "t33", "t34", "t41", "t42", "t43" };
	int trans[c_nGeneCam];
	FILE *fp = fopen("fileT.txt", "r");
	if (!fp)
	{
		cout << "cant open file!" << endl;
	}
	int i = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%d", &trans[i]);
		i++;
	}
	fclose(fp);
	ofstream fileTStr;
	fileTStr.open("fileTStr.txt");
	for (i = 0; i < c_nGeneCam; i++)
	{
		fileTStr << changeTransition[trans[i] - 1] << ", ";
	}
}

void changePart2Transition()
{
	int  *nArrIniPop2 = new int[c_nGeneCam];
	int *nArrIniFs2 = new int[c_nGeneCam];
	FILE *fp = fopen("fileP.txt", "r");
	if (!fp)
	{
		cout << "cant open file!" << endl;
	}
	int i = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%d", &nArrIniPop2[i]);
		i++;
	}
	fclose(fp);
	fnCompleteSequence1(nArrIniPop2, nArrIniFs2, c_nGeneCam);
	ofstream fileT;
	fileT.open("fileT.txt");
	for (i = 0; i < c_nGeneCam; i++)
	{
		fileT << nArrIniFs2[i] << '\t';
	}
}

int main()
{
	srand((unsigned)time(NULL));
	ofstream fileP, fileT, fileRunTime;
	fileP.open("fileP.txt", ios::app);
	fileT.open("fileT.txt", ios::app);
	fileRunTime.open("fileRunTime.txt", ios::app);
	int nBest;
	int nArrP[c_nGeneCam], nArrT[c_nGeneCam];
	time_t tSta, tEnd;
	int retry = 10;
	while (retry > 0)
	{
		cout << retry-- << endl;
		tSta = time(NULL);
		fnIterativeLocalSearch(nArrP, nArrT, nBest);
		tEnd = time(NULL);
		cout << nBest << endl;
		for (int t = 0; t < c_nGeneCam; t++)
		{
			if (t == c_nGeneCam - 1)
			{
				fileT << nArrT[t] << endl;
				fileP << nArrP[t] << endl;
			}
			else
			{
				fileP << nArrP[t] << '\t';
				fileT << nArrT[t] << '\t';
			}
		}
		fileRunTime << (tEnd - tSta) << endl;
	}
	fileP.close();
	fileT.close();
	fileRunTime.close();
}
