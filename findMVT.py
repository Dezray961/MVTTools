#### This needs to be converted to c++ as it is too slow in python.
#### I am also not sure how the paper has so many data points for the VT vs k plot.






from importGRB import GRBData
import numpy as np
import pandas as pd


# class to perform the MVT analysis on a GRB
class MVTAnalysis(GRBData):
    # constructor for the MVTAnalysis class
    def __init__(self,
                 name: str,
                 CSVfilePath: str,
                 ) -> None:
        super().__init__(name, CSVfilePath)
        self.lengthOfData: int = len(self.data)
        # find the window size from the length of the data and the divisor (must be a power of 2)
        self.windowSize: int = 2 ** int(np.log2(self.lengthOfData))
        self.kMax: int = self.windowSize // 2
        self.__findKValues()
        

    # method to find the set of k values to use for the analysis
    def __findKValues(
            self
            ) -> None:
        kValues: list = []
        number: int = 1
        while number <= self.kMax:
            kValues.append(number)
            number *= 2
        # remove the first value of kValues (which is 1) as it is not useful for the analysis
        #kValues = kValues[1:]
        self.kSet: list = kValues


    # method to find the local average over k samples
    def __localAverageOverKBins(
            self,
            index: int,
            k: int
            ) -> np.ndarray:
        summationList: list = []
        for i in range(0, k):
            interiorIndex: int = index - i
            summationList.append(self.data['logRate'].to_numpy()[interiorIndex])
        return 1 / k * np.sum(summationList)
    

    # method to find the local average over k samples for the entire data set
    def __localAverageOverKBinsForDataSet(
            self,
            k: int
            ) -> np.ndarray:
        localAverage: np.ndarray = np.zeros(self.lengthOfData)
        for i in range(k, self.lengthOfData):
            localAverage[i] = self.__localAverageOverKBins(i, k)
        return localAverage

    
    # method to find the square of the difference between local averages
    def __squaredDifference(
            self,
            k: int
            ) -> np.ndarray:
        localAverage: np.ndarray = self.__localAverageOverKBinsForDataSet(k)
        squaredDifferences: np.ndarray = np.zeros(self.lengthOfData)
        for i in range(k, self.lengthOfData):
            squaredDifferences[i] = (localAverage[i] - localAverage[i - k]) ** 2
        return squaredDifferences


    # method to find the VT from the square of the difference between local averages
    def __findVT(
            self,
            k: int
            ) -> float:
        squaredDifferences: np.ndarray = self.__squaredDifference(k)
        return np.sqrt(np.mean(squaredDifferences))
    

    # method to find the VT for all possible values of k
    def findVTForAllK(
            self
            ) -> np.ndarray:
        vtValues: np.ndarray = np.zeros(len(self.kSet))
        for i in range(len(self.kSet)):
            vtValues[i] = self.__findVT(self.kSet[i])
            print(f"VT for k = {self.kSet[i]}: {vtValues[i]}")
        return vtValues


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    grbName: str = "GRB080319B"
    csvFilePath: str = f"data/processed/{grbName}LC.csv"
    analysis = MVTAnalysis(grbName, csvFilePath)
    #vt: np.ndarray = analysis.findVTForAllK()
    fig, ax = plt.subplots()
    ax.scatter(analysis.kSet, vt)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('k')
    ax.set_ylabel('VT')