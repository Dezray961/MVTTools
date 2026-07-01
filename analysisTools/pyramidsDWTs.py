"""
ctypes wrapper for the pyramidDWTs C++ library.

This module is meant to be imported by other Python scripts. It wraps the `pyramidDWTs` export
from `cFiles/pyramidDWTs.cpp` and converts NumPy arrays to the flattened output format used by
the C++ API.
"""

from ctypes import CDLL, c_int
from pathlib import Path
from os import system
from numpy import ndarray, ctypeslib, float64, ascontiguousarray, asarray, zeros, floor, log2

class PyramidDWTs:
    """
    Class to wrap the pyramidDWTs C++ library for MODWT and inverse MODWT operations.
    This class provides methods to perform the MODWT decomposition and inverse reconstruction
    using the Haar and Daubechies 2 wavelet filters. It handles the conversion between NumPy arrays and the C++ API, and manages the shared library loading and function signatures.
    """
    def __init__(
            self,
            libraryPath: str | Path | None = None
            ) -> None:
        """Initializes the PyramidDWTs class.

        Args:
            libraryPath (str | Path | None, optional): Path to shared object library. Defaults to None.
        """
        self.libraryPath = self.__resolveLibraryPath(libraryPath)
        self.lib = CDLL(str(self.libraryPath))
        self.__configureSignatures()


    @staticmethod
    def __resolveLibraryPath(
        libraryPath: str | Path | None
    ) -> Path:
        """Resolves the path to the shared library. If a path is provided, it is used. Otherwise, the function looks for the library in the cFiles directory. If the library is not found, it compiles the C++ code into a shared library.

        Args:
            libraryPath (str | Path | None): Path to shared object library. If None, the function looks for the library in the cFiles directory. 
        """
        # If a library path is provided, use it. 
        if libraryPath is not None:
            return Path(libraryPath).expanduser().resolve()

        # Otherwise, look for the library in the cFiles directory
        module_dir = Path(__file__).resolve().parent
        candidates = [
            module_dir.parent / "cFiles" / "pyramidDWTs.so",
            module_dir.parent / "cFiles" / "libpyramidDWTs.so",
            module_dir.parent / "cFiles" / "pyramidDWTs.dylib",
        ]

        for candidate in candidates:
            if candidate.exists():
                return candidate
        # g++ -std=c++20 -openmp -shared -fPIC pyramidDWTs.cpp -o pyramidDWTs.so
        # If we reach this point, no valid library was found so compile the C++ code into a shared library
        system(f"g++ -std=c++20 -openmp -shared -fPIC {module_dir.parent / 'cFiles' / 'pyramidDWTs.cpp'} -o {candidates[0]}")
        return candidates[0]

    
    def __configureSignatures(
            self
            ) -> None:
        """Configures the argument and return types for the C++ functions."""
        # configure the argument and return types for the pyramidDWTs function MODWTTransform
        self.lib.MODWTTransform.argtypes = [
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # signal
            c_int, # signalLength
            c_int, # transformType
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # output
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # scalingCoeffs
        ]
        # configure the return type for the pyramidDWTs function MODWTTransform
        self.lib.MODWTTransform.restype = c_int # exit status

        # configure the argument and return types for the pyramidDWTs function MODWTInverseTransform
        self.lib.MODWTInverseTransform.argtypes = [
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # waveletCoeffs
            c_int, # waveletCoeffsLength
            c_int, # waveletLevels
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # scalingCoeffs
            c_int, # transformType
            ctypeslib.ndpointer(
                dtype = float64,
                ndim = 1,
                flags = "C_CONTIGUOUS"
            ), # output
        ]
        # configure the return type for the pyramidDWTs function MODWTInverseTransform
        self.lib.MODWTInverseTransform.restype = c_int # exit status


    @staticmethod
    def __asDoubleArray(
        values: ndarray | list[float]
        ) -> ndarray:
        """Converts an iterable of floats to a contiguous NumPy array of type float64."""
        return ascontiguousarray(asarray(values, dtype=float64))


    def MODWT(
        self,
        signal: ndarray | list[float],
        transformType: str
        ) -> tuple[ndarray, ndarray]:
        """Performs the MODWT decomposition using the provided wavelet filters and returns the wavelet and scaling coefficients.

        Args:
            signal (ndarray | list[float]): The input signal to be decomposed.
            transformType (str): The type of transform to perform. "haar" for Haar, "db2" for Daubechies 2.
        
        Raises:
            ValueError: If the transformType is not "haar" or "db2".
            RuntimeError: If the MODWT decomposition fails.
        
        Returns:
            tuple[ndarray, ndarray]: The wavelet coefficient pyramid with shape
            (jMax + 1, signalLength) and the final scaling coefficients.
        """
        signalArray: ndarray = self.__asDoubleArray(signal)
        signalLength: int = signalArray.size

        # get the integer code for the transform type
        match transformType: # using match-case so that other transform types can be added in the future
            case "haar":
                transformTypeCode: int = 0
            case "db2":
                transformTypeCode: int = 1
            case _:
                raise ValueError(f"Invalid transform type: {transformType}. Must be 'haar' or 'db2'.")

        waveletLevels: int = max(int(floor(log2(signalLength))) - 1, 0)

        # prepare the output arrays
        waveletCoeffs: ndarray = zeros(waveletLevels * signalLength, dtype=float64)
        scalingCoeffs: ndarray = zeros(signalLength, dtype=float64)

        # call the C++ function
        status: int = self.lib.MODWTTransform(
            signalArray,
            signalLength,
            transformTypeCode,
            waveletCoeffs,
            scalingCoeffs
        )

        # check the status
        if status != 0:
            raise RuntimeError(f"MODWT decomposition failed with status code {status}.")
        
        return waveletCoeffs.reshape(waveletLevels, signalLength), scalingCoeffs


    def inverseMODWT(
        self,
        transformCoefficients: ndarray | list[float],
        scalingCoefficients: ndarray | list[float],
        transformType: str
        ) -> ndarray:
        """Performs the inverse MODWT reconstruction using the provided wavelet and scaling coefficients.
        
        Args:
            transformCoefficients (ndarray | list[float]): The wavelet coefficients from the MODWT decomposition.
            scalingCoefficients (ndarray | list[float]): The scaling coefficients from the MODWT decomposition.
            transformType (str): The type of transform to perform. "haar" for Haar, "db2" for Daubechies 2.
        
        Raises:
            ValueError: If the transformType is not "haar" or "db2".
            RuntimeError: If the inverse MODWT reconstruction fails.

        Returns:
            ndarray: The reconstructed signal.
        """
        transformCoeffsArray: ndarray = asarray(transformCoefficients, dtype=float64)
        if transformCoeffsArray.ndim == 2:
            waveletLevels: int = transformCoeffsArray.shape[0]
            transformCoeffsArray = ascontiguousarray(transformCoeffsArray.reshape(-1))
        elif transformCoeffsArray.ndim == 1:
            waveletLevels = 1
            transformCoeffsArray = ascontiguousarray(transformCoeffsArray)
        else:
            raise ValueError("transformCoefficients must be a 1D or 2D array-like object.")

        scalingCoeffsArray: ndarray = self.__asDoubleArray(scalingCoefficients)
        transformCoeffsLength: int = transformCoeffsArray.size
        signalLength: int = transformCoeffsLength // waveletLevels

        # get the integer code for the transform type
        match transformType: # using match-case so that other transform types can be added in the future
            case "haar":
                transformTypeCode: int = 0
            case "db2":
                transformTypeCode: int = 1
            case _:
                raise ValueError(f"Invalid transform type: {transformType}. Must be 'haar' or 'db2'.")
        
        # prepare the output array
        output: ndarray = zeros(signalLength, dtype=float64)

        # call the C++ function
        status: int = self.lib.MODWTInverseTransform(
            transformCoeffsArray,
            transformCoeffsLength,
            waveletLevels,
            scalingCoeffsArray,
            transformTypeCode,
            output
        )

        # check the status
        if status != 0:
            raise RuntimeError(f"Inverse MODWT reconstruction failed with status code {status}.")

        return output


def MODWT(
    signal: ndarray | list[float],
    transformType: str = "haar"
    ) -> tuple[ndarray, ndarray]:
    """Performs the MODWT decomposition using the provided wavelet filters and returns the wavelet and scaling coefficients.

    Args:
        signal (ndarray | list[float]): The input signal to be decomposed.
        transformType (str, optional): The type of transform to perform. "haar" for Haar, "db2" for Daubechies 2. Defaults to "haar".
    
    Raises:
        ValueError: If the transformType is not "haar" or "db2".
        RuntimeError: If the MODWT decomposition fails.

    Returns:
        tuple[ndarray, ndarray]: The wavelet and scaling coefficients.
    """
    return PyramidDWTs().MODWT(signal, transformType)


def inverseMODWT(
    transformCoefficients: ndarray | list[float],
    scalingCoefficients: ndarray | list[float],
    transformType: str = "haar"
    ) -> ndarray:
    """Performs the inverse MODWT reconstruction using the provided wavelet and scaling coefficients.
    
    Args:
        transformCoefficients (ndarray | list[float]): The wavelet coefficients from the MODWT decomposition.
        scalingCoefficients (ndarray | list[float]): The scaling coefficients from the MODWT decomposition.
        transformType (str, optional): The type of transform to perform. "haar" for Haar, "db2" for Daubechies 2. Defaults to "haar".
    
    Raises:
        ValueError: If the transformType is not "haar" or "db2".
        RuntimeError: If the inverse MODWT reconstruction fails.

    Returns:
        ndarray: The reconstructed signal.
    """
    return PyramidDWTs().inverseMODWT(transformCoefficients, scalingCoefficients, transformType)