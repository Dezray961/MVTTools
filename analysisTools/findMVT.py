#### I am also not sure how the paper has so many data points for the VT vs k plot.


from analysisTools.importLightCurve import LightCurveData
import numpy as np
import pandas as pd
import ctypes
import os
from contextlib import chdir
from threading import Thread
import time
import matplotlib.pyplot as plt

from tqdm import tqdm


"""
THIS NEEDS EXTENSIVE WORK TO MAKE IT COMPATIBLE WITH THE PRE-PROCESSING DONE IN analysisPipe.py.
The C++ files:
    findMVT.cpp
    HaarCoefficient.hpp
    HaarCoefficient.cpp
    permuteAnalysis.cpp
    permuteAnalysis.hpp
all need to be commented and adjusted.
✅
"""




# time decorator
def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time} seconds to execute.")
        return result
    return wrapper


# class to hand data to the C++ code and get the results back
class analyseLightCurve:
    # constructor for the MVTAnalysis class
    def __init__(
            self,
            data: LightCurveData,
            numberOfTimeBins: int = 64,
            timeWindow: tuple[int, int] = None,
            ) -> None:
        super().__init__()
        self.__data = data
        if timeWindow is not None:
            self.__timeWindow = timeWindow
            self.__windowData()
        self.__rebinnedData = self.__data.rebinnedData.copy() # extract the rebinned data from the LightCurveData object
        self.__numberOfTimeBins = numberOfTimeBins
        # keep only valid bins before building the ctypes buffers
        validDataMask = np.isfinite(self.__rebinnedData['time'].to_numpy())
        validDataMask &= np.isfinite(self.__rebinnedData['rate'].to_numpy())
        validDataMask &= np.isfinite(self.__rebinnedData['error'].to_numpy())
        validDataMask &= self.__rebinnedData['rate'].to_numpy() > 0.0
        self.__rebinnedData = self.__rebinnedData.loc[validDataMask].copy()
        self.__lengthOfData = len(self.__rebinnedData)
        targetPermutationCount = max(256, self.__numberOfTimeBins * 8)
        self.__permutationStep = max(1, self.__lengthOfData // targetPermutationCount)
        if self.__lengthOfData < 2:
            raise ValueError("rebinned light curve must contain at least two positive-rate bins before MVT analysis")
        # get the log of the rate and the propagated log-rate variance for each time bin
        rateArray = self.__rebinnedData['rate'].to_numpy(dtype=np.double, copy=False)
        errorArray = self.__rebinnedData['error'].to_numpy(dtype=np.double, copy=False)
        self.__rebinnedData['logRate'] = np.log(rateArray)
        self.__rebinnedData['logRateErrorSquared'] = (errorArray / rateArray)**2
        # call the C++ wrapper to set up the analysis
        self.cWrapper()
        # run the analysis
        self.runAnalysis()


    def __windowData(self) -> None:
        """
        Windows the light curve data to the specified time window. This method updates the data.rebinnedData attribute to only include the data within the specified time window. This is to save unnecessary calculations on data outside the time window of interest.
        """
        # extract the start and end indices of the time window
        startIndex, endIndex = self.__timeWindow
        # update the data.rebinnedData attribute to only include the data within the specified time window
        self.__data.rebinnedData = self.__data.rebinnedData.iloc[startIndex:endIndex + 1]




    # method to hand the data to the C++ code and get the results back
    def cWrapper(
            self
    ) -> None:
        # tell the interpreter where to find the C++ shared library
        self.libPath: str = os.path.abspath(os.path.join(os.getcwd(), "cFiles", "findMVT.so"))
        self.lib: ctypes.CDLL = ctypes.CDLL(self.libPath)
        # declare the argument and return types for the C++ function
        self.lib.allocatePermuteAnalysis.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int
        ]
        self.lib.allocatePermuteAnalysis.restype = ctypes.c_void_p
        self.lib.getLogBinEdges.argtypes = [ctypes.c_void_p]
        self.lib.getLogBinEdges.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.getLogBinCenters.argtypes = [ctypes.c_void_p]
        self.lib.getLogBinCenters.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.getPowerSetAverages.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self.lib.getPowerSetAverages.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.getPowerSetStdDevs.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self.lib.getPowerSetStdDevs.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.runPermuteAnalysis.argtypes = [ctypes.c_void_p]
        self.lib.runPermuteAnalysis.restype = None
        self.lib.getPermutationsCompleted.argtypes = [ctypes.c_void_p]
        self.lib.getPermutationsCompleted.restype = ctypes.c_int
        self.lib.getTotalPermutations.argtypes = [ctypes.c_void_p]
        self.lib.getTotalPermutations.restype = ctypes.c_int
        self.lib.isAnalysisComplete.argtypes = [ctypes.c_void_p]
        self.lib.isAnalysisComplete.restype = ctypes.c_int
        self.lib.freePermuteAnalysis.argtypes = [ctypes.c_void_p]
        self.lib.freePermuteAnalysis.restype = None


    # method to run the C++ code to get the results back
    @timeit
    def runAnalysis(
            self
            ) -> None:
        if self.__lengthOfData < 2:
            raise ValueError("rebinned light curve must contain at least two positive-rate bins before MVT analysis")
        self.analysis = self.lib.allocatePermuteAnalysis(
            self.__rebinnedData['logRate'].to_numpy(dtype=np.double, copy=False).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.__rebinnedData['time'].to_numpy(dtype=np.double, copy=False).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.__rebinnedData['logRateErrorSquared'].to_numpy(dtype=np.double, copy=False).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.__lengthOfData,
            ctypes.c_int(self.__numberOfTimeBins),
            ctypes.c_int(self.__permutationStep)
        )

        worker = Thread(target=self.lib.runPermuteAnalysis, args=(self.analysis,))
        worker.start()

        totalPermutations = self.lib.getTotalPermutations(self.analysis)
        lastCompleted = 0
        progressBar = tqdm(total=totalPermutations, desc='Permutations', unit='perm')

        while worker.is_alive():
            completed = self.lib.getPermutationsCompleted(self.analysis)
            if completed > lastCompleted:
                if progressBar is not None:
                    progressBar.update(completed - lastCompleted)
                lastCompleted = completed
            time.sleep(0.05)

        worker.join()

        completed = self.lib.getPermutationsCompleted(self.analysis)
        if completed > lastCompleted and progressBar is not None:
            progressBar.update(completed - lastCompleted)

        if progressBar is not None:
            progressBar.close()

        self.logBinEdges: np.ndarray = self.getLogBinEdgesAsArray()
        self.logBinCentres: np.ndarray = self.getLogBinCentres()
        self.powerSetAverages: np.ndarray = self.getPowerSetAverages()
        self.powerSetStdDevs: np.ndarray = self.getPowerSetStdDevs()


    # method to conbert the log bin edges set to an array of doubles
    def getLogBinEdgesAsArray(
            self
            ) -> np.ndarray:
        logBinEdgesArrayPtr = self.lib.getLogBinEdges(self.analysis)
        logBinEdgesSize = self.__numberOfTimeBins + 1
        LogBinEdgesArray: np.ndarray = np.ctypeslib.as_array(
            logBinEdgesArrayPtr,
            shape=(logBinEdgesSize,))
        return np.array(LogBinEdgesArray, copy=True)


    # method to convert the log bin centres set to an array of doubles
    def getLogBinCentres(
            self
            ) -> np.ndarray:
        logBinCentresArrayPtr = self.lib.getLogBinCenters(self.analysis)
        logBinCentresSize = self.__numberOfTimeBins
        LogBinCentresArray: np.ndarray = np.ctypeslib.as_array(
            logBinCentresArrayPtr,
            shape=(logBinCentresSize,))
        return np.array(LogBinCentresArray, copy=True)
    

    # method to convert the power set averages to an array of doubles
    def getPowerSetAverages(
            self
            ) -> np.ndarray:
        PowerSetAveragesArray: np.ndarray = np.array(
            [self.lib.getPowerSetAverages(self.analysis, index)[0] for index in range(self.__numberOfTimeBins)]
        )
        return PowerSetAveragesArray
    

    # method to convert the power set standard deviations to an array of doubles
    def getPowerSetStdDevs(
            self
            ) -> np.ndarray:
        PowerSetStdDevsArray: np.ndarray = np.array(
            [self.lib.getPowerSetStdDevs(self.analysis, index)[0] for index in range(self.__numberOfTimeBins)]
        )
        return PowerSetStdDevsArray
    

    # method to find the uncertainty in the time bins
    def getTimeBinUncertainty(
            self
            ) -> np.ndarray:
        self.timeBinUncertaintyArray: list = []
        for i, edge in enumerate(self.logBinEdges):
            if i == 0:
                continue
            self.timeBinUncertaintyArray.append((edge - self.logBinEdges[i-1]) / 2)
        self.timeBinUncertaintyArray = np.array(self.timeBinUncertaintyArray)


    def __del__(self) -> None:
        if hasattr(self, "analysis") and self.analysis:
            if hasattr(self, "lib"):
                self.lib.freePermuteAnalysis(self.analysis)
            self.analysis = None


# class to analyse the light curves for pre-burst and burst data for a given GRB.
class analyseGRB:
    def __init__(
            self,
            GRBName: str,
            numberOfTimeBins: int = 64,
            timeWindow: tuple[int, int] = None
            ) -> None:
        self.GRBName: str = GRBName
        self.__numberOfTimeBins: int = numberOfTimeBins
        self.__timeWindow: tuple[int, int] = timeWindow
        self.data: LightCurveData = LightCurveData(GRBName)
        MonteCarloUncertainty(
            data = self.data,
            denoisedDataArgs = {
                "thresholdMethod": "hard",
                "thresholdScaleFactor": 0.5
            },
            numSimulations = 100,
            highRAMSystem = False
        )
        rebinLightCurve(self.data,
                        instrument = "swiftBAT",
                        snrThreshold = 5.0
                        )
        self.burstAnalysis: analyseLightCurve = analyseLightCurve(
            data = self.data,
            numberOfTimeBins = self.__numberOfTimeBins,
            timeWindow = self.__timeWindow
        )
        self.logBinEdges: np.ndarray = self.burstAnalysis.logBinEdges
        self.logBinCentres: np.ndarray = self.burstAnalysis.logBinCentres
        self.powerSetAverages: np.ndarray = self.burstAnalysis.powerSetAverages

    # method to filter the values that are less than 3σ from the background
    def filterValues(
            self
            ) -> None:
        pass

    def plotVTvsDeltaT(
            self
            ) -> None:
        fig, ax = plt.subplots()
        ax.errorbar(
            10**self.logBinCentres,
            self.powerSetAverages,
#            xerr=self.burstAnalysis.timeBinUncertaintyArray,
            ls='none',
            marker='x',
            color='blue')
        ax.set_ylabel('VT')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('$\\Delta t$')
        plt.savefig(f"data/processed/{self.GRBName}/{self.GRBName}VTvsDeltaT.png", dpi=300)
        plt.show()





if __name__ == "__main__":
    """
    this script won't generate the  plot if it is run in a terminal. In the interactive 
    window of VSCode, it generates the plot. However, the interpreter in the interactive
    window needs to be restarted after each run of the script. Otherwise, it will throw an
    error allocateMVTAnalysis() is an undefined symbol. I don't know why this happens. 
    """

    from analysisTools.parametricMCUncertainty import MonteCarloUncertainty
    from analysisTools.rebinLightCurve import rebinLightCurve


    # compile the C++ code to a shared library
    with chdir("cFiles"):
        os.system("g++ -O3 -fPIC -shared -std=c++17 -fopenmp findMVT.cpp HaarCoefficient.hpp HaarCoefficient.cpp permuteAnalysis.cpp permuteAnalysis.hpp -o findMVT.so")
    
    # run the analysis on a specific GRB
    grbName: str = "GRB080319B"
    analysis = analyseGRB(grbName, numberOfTimeBins=64)

