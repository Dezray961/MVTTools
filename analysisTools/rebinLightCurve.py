"""
Rebins the light curve data to a constant SNR. The rebinned data is added to the LightCurveData object
as a new dataframe.

This will require a method for Swift BAT data and a method for Fermi GBM data. It will likely require a
method for SVOM data as well. 

Swift BAT data:


"""

from analysisTools.importLightCurve import LightCurveData
import numpy as np
import pandas as pd

from time import time


def swiftBATRebin(
        data: LightCurveData,
        snrThreshold: float
        ) -> None:
    """Rebins the Swift BAT light curve data to a constant SNR. The rebinned data is added to the LightCurveData object as a new dataframe.

    Args:
        data (LightCurveData): light curve data to rebin
        snrThreshold (float): the SNR threshold for rebinning

    Returns:
        None
    """

    # Extract the relevant data from the LightCurveData object
    burstRate: np.ndarray = data.burstData['denoisedRate'].to_numpy()
    burstError: np.ndarray = data.burstData['error'].to_numpy()


    # convert the burst rate and error to counts per bin
    deltaTime: float = 100e-6  # 100 microseconds in seconds
    burstCounts: np.ndarray = burstRate * deltaTime
    burstCountsError: np.ndarray = burstError * deltaTime

    # square the burst count error to get the variance
    burstCountsVariance: np.ndarray = burstCountsError ** 2


    # rebin the data to a constant SNR
    rebinnedCounts: list[float] = []
    rebinnedVariance: list[float] = []
    rebinnedTime: list[float] = []
    rebinnedTimeInBin: list[float] = []
    accumulatedCounts: float = 0.0
    accumulatedVariance: float = 0.0
    accumulatedTime: float = 0.0
    accumulatedTimeInBin: float = 0.0
    for i in range(len(burstCounts)):
        accumulatedCounts += burstCounts[i]
        accumulatedVariance += burstCountsVariance[i]
        accumulatedTime += deltaTime
        accumulatedTimeInBin += deltaTime

        # calculate the SNR for the current bin
        if accumulatedVariance > 0:
            currentSNR: float = accumulatedCounts / np.sqrt(accumulatedVariance)
        else:
            currentSNR: float = 0.0

        # if the SNR exceeds the threshold, save the rebinned data and reset the counters (not time)
        if currentSNR >= snrThreshold:
            rebinnedCounts.append(accumulatedCounts)
            rebinnedVariance.append(accumulatedVariance)
            rebinnedTime.append(accumulatedTime)
            rebinnedTimeInBin.append(accumulatedTimeInBin)
            accumulatedCounts = 0.0
            accumulatedVariance = 0.0
            accumulatedTimeInBin = 0.0
    
    # convert the rebinned counts and variance back to rates and errors
    rebinnedRate: np.ndarray = np.array(rebinnedCounts) / np.array(rebinnedTimeInBin)
    rebinnedError: np.ndarray = np.sqrt(np.array(rebinnedVariance)) / np.array(rebinnedTimeInBin)

    # create a new dataframe for the rebinned data
    rebinnedData: pd.DataFrame = pd.DataFrame(
        {
        'time': rebinnedTime,
        'rate': rebinnedRate,
        'error': rebinnedError,
        'timeInBin': rebinnedTimeInBin
        }
        )
    data.rebinnedData = rebinnedData




def fermiGBMRebin(
        data: LightCurveData,
        snrThreshold: float
        ) -> None:
    """Rebins the Fermi GBM light curve data to a constant SNR. The rebinned data is added to the LightCurveData object as a new dataframe.

    Args:
        data (LightCurveData): light curve data to rebin
        snrThreshold (float): the SNR threshold for rebinning
    
    Returns:
        None
    """
    pass


def svomRebin(
        data: LightCurveData,
        snrThreshold: float
        ) -> None:
    """Rebins the SVOM light curve data to a constant SNR. The rebinned data is added to the LightCurveData object as a new dataframe.

    Args:
        data (LightCurveData): light curve data to rebin
        snrThreshold (float): the SNR threshold for rebinning

    Returns:
        None
    """
    pass


def rebinLightCurve(
        data: LightCurveData,
        instrument: str,
        snrThreshold: float = 5.0
        ) -> None:
    """Rebins the light curve data to a constant SNR. The rebinned data is added to the LightCurveData object as a new dataframe.

    Args:
        data (LightCurveData): light curve data to rebin
        instrument (str): the instrument for which to rebin data
        snrThreshold (float): the SNR threshold for rebinning
    
    Raises:
        ValueError: if the instrument is not supported for rebinning

    Returns:
        None
    """
    print(f"Rebinning light curve data for instrument: {instrument} with SNR threshold: {snrThreshold}")
    match instrument:
        case "swiftBAT":
            swiftBATRebin(data, snrThreshold)
        case "fermiGBM":
            fermiGBMRebin(data, snrThreshold)
        case "SVOM":
            svomRebin(data, snrThreshold)
        case _:
            raise ValueError(f"Instrument {instrument} not supported for rebinning.")
    print("Rebinning complete.")


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from analysisTools.haarDenoise import haarDenoise
    def quickAnalysis(
            thresholdMethod: str,
            thresholdScaleFactor: float
            ) -> None:
        def plotlightCurve(
            data: LightCurveData
            ) -> None:


            plt.figure(figsize=(10, 8))
            plt.plot(
                data.rebinnedData['time'],
                data.rebinnedData['rate'],
                label='Denoised Light Curve'
                )
            plt.xlabel('Time (s)')
            plt.ylabel('Rate (counts/s)')
            plt.ylim(bottom=0, top=16)
            plt.legend()
            plt.show()
        

        grbName: str = "GRB080319B"
        data: LightCurveData = LightCurveData(grbName)
        haarDenoise(data, thresholdMethod = thresholdMethod, thresholdScaleFactor = thresholdScaleFactor)
        # swiftBAT test
        rebinLightCurve(data, "swiftBAT", 5.0)
        plotlightCurve(data)
    
    
    quickAnalysis("soft", 1.0)
    print("Soft thresholding with scale factor 1.0")
    
    quickAnalysis("soft", 0.5)
    print("Soft thresholding with scale factor 0.5")
    
    quickAnalysis("hard", 1.0)
    print("Hard thresholding with scale factor 1.0")
    
    quickAnalysis("hard", 0.5)
    print("Hard thresholding with scale factor 0.5")
