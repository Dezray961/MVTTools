from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
import numpy as np

class MVTFinder:
    """
    Class to find the MVT of a given light curve data. It takes in a LightCurveData object and a time window, and provides methods to find the MVT of the light curve data within that time window.
    """
    def __init__(
            self,
            data: LightCurveData,
            timeWindow: tuple[int, int]
        ) -> None:
        """
        Initializes the MVTFinder class with the given light curve data and time window size.

        Args:
            data (LightCurveData): The light curve data to analyze. This should already be denoised and rebinned to a constant SNR.
            timeWindow (tuple[int, int]): The time window to consider for MVT calculation. Should be a tuple of two integers representing the start and end indices of the time window in the rebinned data.
        """
        # store the light curve data and time window as instance variables
        self.__data: LightCurveData = data
        self.__timeWindow: tuple[int, int] = timeWindow
        # window the data to the specified time window
        self.__windowData()
        # extract the log rate  from the rebinned data and store it as an instance variable
        self.__logRate: np.ndarray = np.log(data.rebinnedData['rate'].to_numpy())
        # calculate the propagated error of the log rate and store it as an instance variable
        self.__propagatedError: np.ndarray = data.rebinnedData['error'].to_numpy() / data.rebinnedData['rate'].to_numpy()
        # get the time in bins from the rebinned data and store it as an instance variable
        self.__timeInBins: np.ndarray = data.rebinnedData['timeInBin'].to_numpy()
        # find the undecimated Haar transform of the log rate and store the detail and scale coefficients as instance variables
        self.__detailCoefficients, self.__scaleCoefficients = MODWT(self.__logRate)
        # find the Allan variance of the detail coefficients and store it as an instance variable
        self.__allanVariance: np.ndarray = self.__getAllanVariance()
        self.haarWaveletSF: np.ndarray = np.sqrt(self.__allanVariance)
        print(self.haarWaveletSF)
        
    
    def __windowData(self) -> None:
        """
        Windows the light curve data to the specified time window. This method updates the data.rebinnedData attribute to only include the data within the specified time window. This is to save unnecessary calculations on data outside the time window of interest.
        """
        # extract the start and end indices of the time window
        startIndex, endIndex = self.__timeWindow
        # update the data.rebinnedData attribute to only include the data within the specified time window
        self.__data.rebinnedData = self.__data.rebinnedData.iloc[startIndex:endIndex + 1]


    @staticmethod
    def __haarSupport(
            level: int,
            timeIndex: int,
            N: int
            ) -> int:
        """
        Calculates the support of the Haar wavelet at a given level and time index. The support is the range of time indices over which the Haar wavelet is non-zero. This method is used to determine the range of time indices that contribute to the calculation of the MVT at a given time index.

        Args:
            level (int): The level of the Haar wavelet. Level 0 corresponds to the coarsest scale, and higher levels correspond to finer scales.
            timeIndex (int): The time index for which to calculate the support of the Haar wavelet.
            N (int): The total number of time indices in the data. This is used to wrap around the time indices when calculating the support of the Haar wavelet.

        Returns:
            int: The support of the Haar wavelet at the given level and time index.
        """
        L: int = 2 ** level
        return [
            (timeIndex + k) % N
            for k in range(L)
            ]
    


    def __getAllanVariance(self) -> np.ndarray:
        """
        Calculates the Allan variance of the detail coefficients obtained from the MODWT Haar transform of the np.log rate.

        Returns:
            np.ndarray: The Allan variance of the detail coefficients.
        """
        # get the number of levels in the detail coefficients
        numLevels: int = len(self.__detailCoefficients)
        allanVariance: np.ndarray = np.zeros(numLevels)
        for level in range(numLevels):
            # get the detail coefficients for the current level
            detailCoefficients: np.ndarray = self.__detailCoefficients[level]
            N: int = len(detailCoefficients) # length of the detail coefficients at the current level
            sigmaW2: np.ndarray = np.zeros(N) # array to store the Allan variance for the current level

            # loop over each time index in the detail coefficientsj
            for timeIndex in range(N):
                indices: list[int] = self.__haarSupport(level, timeIndex, N)

                sigmaW2[timeIndex] = (
                    np.sum(
                        self.__propagatedError[indices] ** 2
                    ) / (2 ** level)
                )
            allanVariance[level] = np.mean(detailCoefficients ** 2 - sigmaW2)

        return allanVariance







if __name__ == "__main__":
    from analysisTools.rebinLightCurve import rebinLightCurve
    from analysisTools.parametricMCUncertainty import MonteCarloUncertainty
    data = LightCurveData("GRB080319B")
    MonteCarloUncertainty(
        data = data,
        denoisedDataArgs = {
            "thresholdMethod": "hard",
            "thresholdScaleFactor": 0.5
        },
        numSimulations = 500,
        highRAMSystem = False
    )
    rebinLightCurve(data, "swiftBAT", snrThreshold=5.0)
    mvtFinder = MVTFinder(data, timeWindow=(0, len(data.rebinnedData)))