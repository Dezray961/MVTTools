from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
from analysisTools.haarDenoise import haarDenoise
from numpy import ceil, log, ndarray, arange, zeros, newaxis
from scipy.stats import chi2
import pandas as pd

import matplotlib.pyplot as plt


class HaarMTVFinder:
    """
    A class to find the Haar MTV (Maximum Time Variation) in a light curve.
    Attributes:
        data (LightCurveData): The light curve data for which to find the Haar MTV.


    """
    def __init__(
            self,
            data: LightCurveData,
            chi2CriticalTablePointer: ndarray,
            plot1: bool = False,
            SNRThreshold: float = 5.0,
            logFile: bool = False,              #TODO: implement log file functionality
            silent: bool = True                 #TODO: implement silent/verbose functionality
    ) -> None:
        """docstring here"""

        # initialise class attributes from constructor arguments
        self.data: LightCurveData = data # The light curve data for which to find the Haar MTV.
        self.__plot1: bool = plot1 # Whether to plot the Haar MTV finding process.
        self.__SNRThreshold: float = SNRThreshold # The threshold for signal-to-noise ratio.
        if logFile: # Whether to create a log file for the Haar MTV finding process.
            self.__logFile: str = f"{self.data.name}_haarMTVFinder_log.txt"
        self.__silent: bool = silent # Whether to print verbose output to the console.
        self.__chi2CriticalTable: ndarray = chi2CriticalTablePointer # The chi-squared critical values for the Haar MTV finding process. 

        # Set the geometric bins for the Haar MTV finding process.
        self.__setGeometricBins()

        # initialise arrays
        self.__initialiseArrays()

        # initialise scalars
        self.__initialiseScalars()

        



    def __setGeometricBins(
            self
            ) -> None:
        """Sets the geometric bins for the Haar MTV finding process."""
        # define the minimum and maximum delta time values
        minimumDeltaTime: float = 1.0e-4
        maximumDeltaTime: float = 1.0e3

        # define the number of bins and the bin factor
        self.__binningFactor: float = 2.0
        self.__numberOfBins: int = (
            self.__binningFactor * ceil(log(maximumDeltaTime / minimumDeltaTime) / log(2.0)))

        # log the minimum and maximum delta time values
        logMinimumDeltaTime: float = log(minimumDeltaTime) / log(2.0)
        logMaximumDeltaTime: float = log(maximumDeltaTime) / log(2.0)

        # geometrically spaced bins
        self.__timeBinStart: ndarray = 2.0 ** (
            logMinimumDeltaTime + (logMaximumDeltaTime - logMinimumDeltaTime)
            * arange(self.__numberOfBins) / (self.__numberOfBins - 1)
        )
        self.__timeBinEnd: ndarray = self.__shift(self.__timeBinStart, 1)

        # subtract 1 from the number of bins because we are using deltaTime rather than bin edges
        self.__numberOfBins -= 1

        # truncate the timeBinStart/End arrays to the number of bins
        self.__timeBinStart = self.__timeBinStart[:self.__numberOfBins]
        self.__timeBinEnd = self.__timeBinEnd[:self.__numberOfBins]


    def __initialiseArrays(
            self
    ) -> None:
        """Initialises the arrays for the Haar MTV finding process. Seperate method to keep the constructor clean. Using float32 to save memory and improve performance."""
        self.__binChi2Sum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of chi-squared values for each bin.
        self.__binAdjustedWeightSum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of adjusted weights for each bin.
        self.__binRawWeightSum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of raw weights for each bin.
        self.__noiseBaseline: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # noise baseline for each bin.
        self.__powerSpectrum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # power spectrum for each bin.
        self.__powerSpectrumError: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # error in the power spectrum for each bin
        self.binTermCounts: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # number of terms in each bin.
        
        # extract the data from the LightCurveData object and store it in arrays
        self.__time: ndarray = self.data.burstData['time'].to_numpy(dtype='float32') # time values for the light curve data.
        self.__timeBinDuration: ndarray = self.data.burstData['timeInBin'].to_numpy(dtype='float32') # duration of each time bin for the light curve data.
        self.__rate: ndarray = self.data.burstData['rate'].to_numpy(dtype='float32') # rate values for the light curve data.
        self.__deltaRate: ndarray = self.data.burstData['error'].to_numpy(dtype='float32') # error in the rate values for the light curve data.
        self.__logRate: ndarray = log(self.__rate) # logarithm of the rate values for the light curve data.
        self.__deltaLogRate: ndarray = self.__deltaRate / self.__rate # error in the logarithm of the rate values for the light curve data.


    def __initialiseScalars(
            self
    ) -> None:
        """Initialises the scalars for the Haar MTV finding process. Seperate method to keep the constructor clean."""
        self.__minimumDeltaTimeArray: float = self.__timeBinStart.max() # initalised to the maximum value of the timeBinStart array, will be updated during the Haar MTV finding process.
        self.__maximumDeltaTimeArray: float = 0.0 # initalised to 0, will be updated during the Haar MTV finding process.
        self.__degreesOfFreedom: int = len(self.__chi2CriticalTable) # initalised to the length of the chi2CriticalTable, will be updated during the Haar MTV finding process.
        self.__numberOfRepetitions: int = 1 # techincally this is only 1 for the Haar transform, if DB2, or some other WT, is implimented then it will need to be changed.











    @staticmethod
    def __shift(
        list: list,
        n: int
    ) -> list:
        """ Shift the elements of a list by n positions to the left.

        Args:
            list list[any]: The list to be shifted.
            n int: The number of positions to shift.

        Returns:
            list[any]: The shifted list.
        """
        return list[n:] + list[:n]



    def __padSignalTimeline(
            self
    ) -> None:
        
        # denoise the signal timeline using the MODWT
        self.__logRate = haarDenoise(self.__logRate) # NOTE Golkhu's code also hands in the error to this step

        # calculate the span of a single repetition of the signal timeline
        maximumTimeDifference: float = self.__time.max() - self.__time.min()
        totalCycles: int = self.__numberOfRepetitions + 1
        
        # tile the measurement arrays
        self.__timeBinDuration = self.__timeBinDuration.tile(totalCycles)
        self.__logRate = self.__logRate.tile(totalCycles)
        self.__deltaLogRate = self.__deltaLogRate.tile(totalCycles)

        # construct an advancing time array for the tiled signal timeline
        timeOffset: float = maximumTimeDifference * arange(totalCycles)
        self.__time = (self.__time[:, newaxis] + timeOffset).ravel(order = 'F')

        # perform the Haar wavelet transform on the tiled signal timeline
        # TODO: write a wrapper for the MODWT to calculate all the extra stuff needed for Golkhu's method




    
    def __processHaarWaveletTransform(
            self
    ) -> None:
        ...
    

    def __adaptiveRebinBySNR(
            self
    ) -> None:
        ...
    

    def __evaluteSignificanceAndFit(
            self
    ) -> None:
        ...


def getChi2LookupTable(
        data: LightCurveData,
    ) -> ndarray:
    """
    Generates a chi-squared lookup table for a given LightCurveData object. The lookup table is used to determine the critical chi-squared value for a given number of degrees of freedom.

    Args:
        data (LightCurveData): The light curve data for which to generate the lookup table.

    Returns:
        ndarray: The chi-squared lookup table.
    """
    expectedDegreesOfFreedom: int = len(data.burstData) + 1
    return chi2.ppf(0.95, arange(0, expectedDegreesOfFreedom))


if __name__ == "__main__":
    data: LightCurveData = LightCurveData("GRB080319B")
    chi2LookupTable: ndarray = getChi2LookupTable(data)