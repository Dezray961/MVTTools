/*
g++ -fPIC -shared -std=c++17 -fopenmp findMVT.cpp MVTClass.hpp MVTClass.cpp -o findMVT.so
*/


# include "MVTClass.hpp"

// function to allocate memory for the MVTAnalysis class
extern "C" {
    MVTAnalysis* allocateMVTAnalysis(
        double *rate,
        double *time,
        double *rateErr,
        int lenghtOfData,
        bool generateEvenKSetFlag
    )
    {
        return new MVTAnalysis(rate, time, rateErr, lenghtOfData, generateEvenKSetFlag);
    }

    double* getVTSetArray(MVTAnalysis* analysis)
    {
        return analysis->VTSetArray;
    }

    int getVTSetSize(MVTAnalysis* analysis)
    {
        return analysis->kSetSize;
    }

    int* getKSetArray(MVTAnalysis* analysis)
    {
        return analysis->kSet.data();
    }

    double* getDeltaTArray(MVTAnalysis* analysis)
    {
        return analysis->deltaT.data();
    }

    double* getDeltaTErrorArray(MVTAnalysis* analysis)
    {
        return analysis->deltaTError.data();
    }
}

int main()
{
    return 0;
}