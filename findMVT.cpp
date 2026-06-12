#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <numeric>

// class to perform MVT calculations
class MVTAnalysis
{
    public:
        // data members
        std::vector<double> rate;
        std::vector<double> time;
        std::vector<double> rateErr;
        const int lenghtOfData = rate.size();
        const int windowSize = pow(2, log2(lenghtOfData));
        const int kMax = windowSize / 2;
        std::vector<int> kSet;



        // constructor
        MVTAnalysis
        (
            std::vector<double>& rate,
            std::vector<double>& time,
            std::vector<double>& rateErr
        )
        {
            this->rate = rate;
            this->time = time;
            this->rateErr = rateErr;
        };
        

    private:
        // private member functions
        // function to find the k set
        void findKSet()
        {
            std::vector<int> kSet;
            int number = 1;
            while (number <= kMax)
            {
                kSet.push_back(number);
                number *= 2;
            }
        };


        // function to find the local average over k samples
        double localAverageOverKBins
        (
            int index,
            int k
        )
        {
            std::vector<double> summationList;
            for (int i = 0; i < k; i++)
            {
                int interiorIndex = index - i;
                summationList.push_back(rate[interiorIndex]);
            }
            return std::accumulate(summationList.begin(), summationList.end(), 0.0) / k;
        };


        // function to find the local average over k samples for the entire data set
        std::vector<double> localAverageOverKBinsForDataSet
        (
            int k
        )
        {
            std::vector<double> localAverages;
            for (int i = k; i < lenghtOfData; i++)
            {
                localAverages.push_back(localAverageOverKBins(i, k));
            }
            return localAverages;
        };

        // function to find the square of the difference between local averages
        std::vector<double> squaredDifference
        (
            int k
        )
        {
            std::vector<double> localAverages = localAverageOverKBinsForDataSet(k);
            std::vector<double> squaredDifferences;
            for (int i = k; i < lenghtOfData; i++)
            {
                double difference = rate[i] - localAverages[i - k];
                squaredDifferences.push_back(pow(difference, 2));
            }
            return squaredDifferences;
        };
    

        // function to find the VT from the square of the difference between local averages
        double findVT
        (
            int k
        )
        {
            std::vector<double> squaredDifferences = squaredDifference(k);
            double sumOfSquaredDifferences = std::accumulate(squaredDifferences.begin(), squaredDifferences.end(), 0.0);
            return sqrt(sumOfSquaredDifferences / (lenghtOfData - k));
        };


        // function to find the VT for all k in the k set
        std::vector<double> findVTForAllK()
        {
            std::vector<double> VTSet;
            for (int k : kSet)
            {
                VTSet.push_back(findVT(k));
            }
            return VTSet;
        };
};

// main function
int main()
{
    // this needs to beable to take a ctypes input from python. I will need to make a function to convert
    // the input from python to vectors of doubles. I will also need to return the output to python.
}
