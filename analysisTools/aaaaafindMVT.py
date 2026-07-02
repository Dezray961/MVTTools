from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
import numpy as np
from tqdm import tqdm
from scipy.ndimage import uniform_filter1d

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
        print("Finding the undecimated Haar transform of the log rate...")
        self.__detailCoefficients, self.__scaleCoefficients = MODWT(self.__logRate)
        # find the Allan variance of the detail coefficients and store it as an instance variable
        print("Calculating the Allan variance of the detail coefficients...")
        self.__allanVariance: np.ndarray = self.__getAllanVariance()
        self.haarWaveletSF: np.ndarray = np.sqrt(self.__allanVariance)
        print(self.haarWaveletSF)
        print(len(self.haarWaveletSF))
        
    
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
        numLevels: int = len(self.__detailCoefficients)
        allanVariance: np.ndarray = np.zeros(numLevels)
        
        # 1. Pre-square the array in place (or create one copy, which is memory-safe)
        propagatedErrorSquared = self.__propagatedError ** 2

        # Simple loop tracker for the progress bar
        with tqdm(total=numLevels, desc="Processing Levels", unit="level") as pbar:
            for level in range(numLevels):
                detailCoefficients: np.ndarray = self.__detailCoefficients[level]
                
                # 2. Determine the window size for the Haar support at this level.
                # In MODWT Haar, the support width is typically 2^level.
                windowSize = 2 ** level 
                
                sigmaW2 = uniform_filter1d(
                    propagatedErrorSquared, 
                    size=windowSize, 
                    mode='wrap' # Matches MODWT periodic/circular boundary conditions
                ) * windowSize / 2
                
                # 4. Apply the remaining math
                sigmaW2 /= (2 ** level)
                allanVariance[level] = np.mean(detailCoefficients ** 2 - sigmaW2)
                
                pbar.update(1)

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
        numSimulations = 100,
        highRAMSystem = False
    )
    rebinLightCurve(data, "swiftBAT", snrThreshold=5.0)
    mvtFinder = MVTFinder(data, timeWindow=(0, len(data.rebinnedData)))