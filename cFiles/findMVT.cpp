/*
g++ -fPIC -shared findMVT.cpp MVTClass.hpp -o findMVT.so
*/


# include "MVTClass.hpp"

// function to allocate memory for the MVTAnalysis class
extern "C" {
    MVTAnalysis* allocateMVTAnalysis(
        double *rate,
        double *time,
        double *rateErr,
        int lenghtOfData
    )
    {
        return new MVTAnalysis(rate, time, rateErr, lenghtOfData);
    }

    double* getVTSetArray(MVTAnalysis* analysis)
    {
        return analysis->VTSetArray;
    }

    int getVTSetSize(MVTAnalysis* analysis)
    {
        return analysis->kSetSize;
    }
}

int main()
{
    return 0;
}