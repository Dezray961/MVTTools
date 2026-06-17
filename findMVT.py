#### I am also not sure how the paper has so many data points for the VT vs k plot.


from importGRB import GRBData
import numpy as np
import pandas as pd
import ctypes
import os
from contextlib import chdir
import time


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
class analyseGRB(GRBData):
    # constructor for the MVTAnalysis class
    def __init__(self,
                 name: str,
                 CSVfilePath: str,
                 ) -> None:
        super().__init__(name, CSVfilePath)
        self.lengthOfData: int = len(self.data)
        # create a pointer to the MVTAnalysis class
        self.cWrapper()
        rate: np.ndarray = np.ascontiguousarray(self.data['rate'].to_numpy(), dtype=np.float64)
        time: np.ndarray = np.ascontiguousarray(self.data['time'].to_numpy(), dtype=np.float64)
        rateErr: np.ndarray = np.ascontiguousarray(self.data['error'].to_numpy(), dtype=np.float64)
        # run the C++ code to get the results back
        self.runAnalysis()


    # method to hand the data to the C++ code and get the results back
    def cWrapper(
            self
    ) -> None:
        # tell the interpreter where to find the C++ shared library
        self.libPath: str = os.path.join(os.path.dirname(__file__), "cFiles", "findMVT.so")
        self.lib: ctypes.CDLL = ctypes.CDLL(self.libPath)
        # declare the argument and return types for the C++ function
        self.lib.allocateMVTAnalysis.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
        ]
        self.lib.allocateMVTAnalysis.restype = ctypes.c_void_p
        self.lib.getVTSetArray.argtypes = [ctypes.c_void_p]
        self.lib.getVTSetArray.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.getVTSetSize.argtypes = [ctypes.c_void_p]
        self.lib.getVTSetSize.restype = ctypes.c_int
        self.lib.getKSetArray.argtypes = [ctypes.c_void_p]
        self.lib.getKSetArray.restype = ctypes.POINTER(ctypes.c_int)
        self.lib.getDeltaTArray.argtypes = [ctypes.c_void_p]
        self.lib.getDeltaTArray.restype = ctypes.POINTER(ctypes.c_double)
        self.lib.getDeltaTErrorArray.argtypes = [ctypes.c_void_p]
        self.lib.getDeltaTErrorArray.restype = ctypes.POINTER(ctypes.c_double)


    # method to run the C++ code to get the results back
    @timeit
    def runAnalysis(
            self
            ) -> None:
        self.analysis = self.lib.allocateMVTAnalysis(
            self.data['rate'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.data['time'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.data['error'].to_numpy().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.lengthOfData
            )
        self.VTSet: np.ndarray = self.getVTSetAsArray()
        self.kSet: np.ndarray = self.getKSetAsArray()
        self.deltaT: np.ndarray = self.getDeltaTAsArray()
        self.deltaTError: np.ndarray = self.getDeltaTErrorAsArray()


    # method to conbert the VT set to an array of doubles
    def getVTSetAsArray(
            self
            ) -> np.ndarray:
        vtSetArrayPtr = self.lib.getVTSetArray(self.analysis)
        vtSetSize = self.lib.getVTSetSize(self.analysis)
        VTSetArray: np.ndarray = np.ctypeslib.as_array(
            vtSetArrayPtr,
            shape=(vtSetSize,))
        return VTSetArray


    # method to convert the k set to an array of integers
    def getKSetAsArray(
            self
            ) -> np.ndarray:
        kSetArrayPtr = self.lib.getKSetArray(self.analysis)
        kSetSize = self.lib.getVTSetSize(self.analysis)
        KSetArray: np.ndarray = np.ctypeslib.as_array(
            kSetArrayPtr,
            shape=(kSetSize,))
        return KSetArray
    

    # method to convert the deltaT set to an array of doubles
    def getDeltaTAsArray(
            self
            ) -> np.ndarray:
        deltaTArrayPtr = self.lib.getDeltaTArray(self.analysis)
        deltaTSize = self.lib.getVTSetSize(self.analysis)
        DeltaTArray: np.ndarray = np.ctypeslib.as_array(
            deltaTArrayPtr,
            shape=(deltaTSize,))
        return DeltaTArray


    # method to convert the deltaTError set to an array of doubles
    def getDeltaTErrorAsArray(
            self
            ) -> np.ndarray:
        deltaTErrorArrayPtr = self.lib.getDeltaTErrorArray(self.analysis)
        deltaTErrorSize = self.lib.getVTSetSize(self.analysis)
        DeltaTErrorArray: np.ndarray = np.ctypeslib.as_array(
            deltaTErrorArrayPtr,
            shape=(deltaTErrorSize,))
        return DeltaTErrorArray


if __name__ == "__main__":
    """
    this script won't generate the  plot if it is run in a terminal. In the interactive 
    window of VSCode, it generates the plot. However, the interpreter in the interactive
    window needs to be restarted after each run of the script. Otherwise, it will throw an
    error allocateMVTAnalysis() is an undefined symbol. I don't know why this happens. 
    """
    import matplotlib.pyplot as plt
    import sys
    # compile the C++ code to a shared library
    with chdir(os.path.join(os.path.dirname(__file__), "cFiles")):
        os.system("g++ -O3 -fPIC -shared -std=c++17 -fopenmp findMVT.cpp MVTClass.hpp MVTClass.cpp -o findMVT.so")
    grbName: str = "GRB080319B"
    csvFilePath: str = f"data/processed/{grbName}LC.csv"
    analysis = analyseGRB(grbName, csvFilePath)
    fig, ax = plt.subplots()
    ax.scatter(
        analysis.deltaT,
        analysis.VTSet,
        color='blue')
    ax.set_ylabel('VT')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$\\Delta t$')
    plt.show()