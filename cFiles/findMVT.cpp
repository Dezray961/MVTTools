/*
g++ -fPIC -shared -std=c++17 -fopenmp findMVT.cpp MVTClass.hpp MVTClass.cpp -o findMVT.so
*/


# include "permuteAnalysis.hpp"

// function to allocate memory for the MVTAnalysis class
extern "C" {
    PermuteAnalysis* allocatePermuteAnalysis(
        double *rate,
        double *time,
        double *rateErr,
        int lenghtOfData,
        int numberOfTimeBins
    )
    {
        return new PermuteAnalysis
        (
            rate,
            time,
            rateErr,
            lenghtOfData,
            numberOfTimeBins
        );
    }

    void runPermuteAnalysis(PermuteAnalysis* analysis)
    {
        analysis->runAnalysis();
    }

    int getPermutationsCompleted(PermuteAnalysis* analysis)
    {
        return analysis->getPermutationsCompleted();
    }

    int getTotalPermutations(PermuteAnalysis* analysis)
    {
        return analysis->getTotalPermutations();
    }

    int isAnalysisComplete(PermuteAnalysis* analysis)
    {
        return analysis->isAnalysisComplete() ? 1 : 0;
    }

    void freePermuteAnalysis(PermuteAnalysis* analysis)
    {
        delete analysis;
    }

    double* getLogBinEdges(PermuteAnalysis* analysis)
    {
        return analysis->logBinEdges.data();
    }

    double* getLogBinCenters(PermuteAnalysis* analysis)
    {
        return analysis->logBinCenters.data();
    }

    double* getPowerSetAverages(PermuteAnalysis* analysis, int index)
    {
        return &analysis->powerSetAverages[index];
    }

    double* getPowerSetStdDevs(PermuteAnalysis* analysis, int index)
    {
        return &analysis->powerSetStdDevs[index];
    }
}

int main()
{
    return 0;
}