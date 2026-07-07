from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def findDeltaT(
        data: LightCurveData
    ) -> tuple[np.ndarray, np.ndarray]:

    # unpack the rebinned data into numpy arrays for processing
    rate: np.ndarray = data.rebinnedData['rate'].to_numpy()
    rateErr: np.ndarray = data.rebinnedData['error'].to_numpy()
    timeBinEdges: np.ndarray = data.rebinnedData['time'].to_numpy()
    timeInBins: np.ndarray = data.rebinnedData['timeInBin'].to_numpy()

    # perform the undecimated Haar transform on the rate and time bin edges
    rateTransformCoefficients, _ = MODWT(rate)
    timeInBinsTransformCoefficients, timeInBinsScalingCoefficients = MODWT(timeInBins)

    # find the set of decomposition levels
    jMax: int = len(rateTransformCoefficients)
    levels: list[int] = list(range(jMax))
    normalisationCoefficients: list[float] = [2 ** (j) for j in levels]

    # find the average bin width
    averageBinWidth: float = timeBinEdges[-1] / (len(timeBinEdges) - 1)

    # find DeltaT for each level of the decomposition
    deltaT: np.ndarray = np.zeros_like(rateTransformCoefficients, dtype=float)
    for j in levels:
        # get the normalisation coefficient for this level
        normalisationCoefficient: float = normalisationCoefficients[j]
        # normalise the time bin edges transform coefficients to find DeltaT
        deltaT[j] = np.abs(timeInBinsTransformCoefficients[j] * (averageBinWidth + normalisationCoefficient))

    # find the power of the rate transform coefficients for each level
    power: np.ndarray = np.square(rateTransformCoefficients)

    # find the standard deviation of the pre-burst data to use as a threshold for the power
    preBurstStdDev: float = np.std(data.preBurstData['denoisedRate'].to_numpy())
    power *= (power > 3 * preBurstStdDev)  # apply a 3 sigma threshold to the power

    # flatten the deltaT and power arrays to 1D arrays
    deltaTFlat: np.ndarray = deltaT.flatten()
    powerFlat: np.ndarray = power.flatten()


    # remove any entries where deltaT is zero or power is zero
    nonZeroIndices: np.ndarray = np.where((deltaTFlat > 1e-2) & (powerFlat > 0))
    deltaTFlat: np.ndarray = deltaTFlat[nonZeroIndices]
    powerFlat: np.ndarray = powerFlat[nonZeroIndices]

    # return the flattened deltaT and power arrays
    return deltaTFlat, powerFlat


def findAveragePower(
    
    deltaT: np.ndarray,
    power: np.ndarray,
    numBins: int = 50
) -> tuple[np.ndarray, np.ndarray]:

    deltaTMin: float = np.min(deltaT)
    deltaTMax: float = np.max(deltaT)
    
    bins = np.logspace(np.log10(deltaTMin), np.log10(deltaTMax), num=numBins + 1)
    counts, _ = np.histogram(deltaT, bins=bins)
    sums, _ = np.histogram(deltaT, bins=bins, weights=power)

    # Prevent division by zero for empty bins
    averagePower = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)


    # Use geometric mean for proper log-spaced bin centers
    binCenters: np.ndarray = np.sqrt(bins[:-1] * bins[1:])
    
    return binCenters, averagePower






if __name__ == "__main__":
    from analysisTools.rebinLightCurve import rebinLightCurve
    from analysisTools.parametricMCUncertainty import MonteCarloUncertainty

    data: LightCurveData = LightCurveData("GRB080319B")
    MonteCarloUncertainty(
        data = data,
        denoisedDataArgs = {
            "thresholdMethod": "hard",
            "thresholdScaleFactor": 0.5
        },
        numSimulations = 100,
        highRAMSystem = False
    )
    rebinLightCurve(
        data,
        instrument = "swiftBAT",
        snrThreshold = 5.0
        )
    deltaT, power = findDeltaT(data)
    averageDeltaT, averagePower = findAveragePower(deltaT, power, numBins=50)

    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.scatter(averageDeltaT, averagePower, alpha=0.5)
    plt.axvline(x=0.04, color='r', linestyle='-', label='Threshold')
    plt.axvline(x=0.03, color='g', linestyle='--', label='Threshold - 0.01')
    plt.axvline(x=0.05, color='g', linestyle='--', label='Threshold + 0.01')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\Delta t$ (s)')
    plt.ylabel('Power')
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(data.rebinnedData['time'], data.rebinnedData['rate'], alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Rate')
    plt.show()
