#### I am also not sure how the paper has so many data points for the VT vs k plot.


from importLightCurve import lightCurveData
import numpy as np
import pandas as pd
import ctypes
import os
from contextlib import chdir
from threading import Thread
import time

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


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
class analyseLightCurve(lightCurveData):
    # constructor for the MVTAnalysis class
    def __init__(
            self,
            name: str,
            CSVfilePath: str,
            numberOfTimeBins: int = 64
            ) -> None:
        super().__init__(name, CSVfilePath)
        self.lengthOfData: int = len(self.data)
        self.numberOfTimeBins: int = numberOfTimeBins
        # create a pointer to the MVTAnalysis class
        self.cWrapper()
        # run the C++ code to get the results back
        self.runAnalysis()
        self.getTimeBinUncertainty()


    # method to hand the data to the C++ code and get the results back
    def cWrapper(
            self
    ) -> None:
        # tell the interpreter where to find the C++ shared library
        self.libPath: str = os.path.join(os.path.dirname(__file__), "cFiles", "findMVT.so")
        self.lib: ctypes.CDLL = ctypes.CDLL(self.libPath)
        # declare the argument and return types for the C++ function
        self.lib.allocatePermuteAnalysis.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
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
        self.analysis = self.lib.allocatePermuteAnalysis(
            self.data['logRate'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.data['time'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.data['logRateErrorSquared'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.lengthOfData,
            ctypes.c_int(self.numberOfTimeBins)
            )

        worker = Thread(target=self.lib.runPermuteAnalysis, args=(self.analysis,))
        worker.start()

        totalPermutations = self.lib.getTotalPermutations(self.analysis)
        lastCompleted = 0
        progressBar = tqdm(total=totalPermutations, desc='Permutations', unit='perm') if tqdm is not None else None

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
        logBinEdgesSize = self.numberOfTimeBins + 1
        LogBinEdgesArray: np.ndarray = np.ctypeslib.as_array(
            logBinEdgesArrayPtr,
            shape=(logBinEdgesSize,))
        return LogBinEdgesArray


    # method to convert the log bin centres set to an array of doubles
    def getLogBinCentres(
            self
            ) -> np.ndarray:
        logBinCentresArrayPtr = self.lib.getLogBinCenters(self.analysis)
        logBinCentresSize = self.numberOfTimeBins
        LogBinCentresArray: np.ndarray = np.ctypeslib.as_array(
            logBinCentresArrayPtr,
            shape=(logBinCentresSize,))
        return LogBinCentresArray
    

    # method to convert the power set averages to an array of doubles
    def getPowerSetAverages(
            self
            ) -> np.ndarray:
        PowerSetAveragesArray: np.ndarray = np.array(
            [self.lib.getPowerSetAverages(self.analysis, index)[0] for index in range(self.numberOfTimeBins)]
        )
        return PowerSetAveragesArray
    

    # method to convert the power set standard deviations to an array of doubles
    def getPowerSetStdDevs(
            self
            ) -> np.ndarray:
        PowerSetStdDevsArray: np.ndarray = np.array(
            [self.lib.getPowerSetStdDevs(self.analysis, index)[0] for index in range(self.numberOfTimeBins)]
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








if __name__ == "__main__":
    """
    this script won't generate the  plot if it is run in a terminal. In the interactive 
    window of VSCode, it generates the plot. However, the interpreter in the interactive
    window needs to be restarted after each run of the script. Otherwise, it will throw an
    error allocateMVTAnalysis() is an undefined symbol. I don't know why this happens. 
    """
    import matplotlib.pyplot as plt
    import sys


    # plotting function
    def plotVTvsDeltaT()-> None:
        fig, ax = plt.subplots()
        ax.errorbar(
            10**analysis.logBinCentres,
            analysis.powerSetAverages,
            ls='none',
            marker='x',
            color='blue')
        ax.set_ylabel('VT')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('$\\Delta t$')
        plt.show()


    # compile the C++ code to a shared library
    with chdir(os.path.join(os.path.dirname(__file__), "cFiles")):
        os.system("g++ -O3 -fPIC -shared -std=c++17 -fopenmp findMVT.cpp HaarCoefficient.hpp HaarCoefficient.cpp permuteAnalysis.cpp permuteAnalysis.hpp -o findMVT.so")
    
    # run the analysis on a specific GRB
    grbName: str = "GRB080319B"
    csvFilePath: str = f"data/processed/{grbName}LC.csv"
    analysis = analyseLightCurve(grbName, csvFilePath)

    #plot the results
    plotVTvsDeltaT()