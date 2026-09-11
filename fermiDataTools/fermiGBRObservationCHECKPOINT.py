"""
This will NOT run on a laptop due to the large memory requirements of the Fermi GBM data!
"""


from analysisTools.GRBData import GRBData
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.binning.unbinned import bin_by_edges as binByEdges
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.background.primitives import BackgroundRates
import numpy as np
import gc
from more_itertools import distinct_combinations as combinations
from analysisTools.haarMethods import haarPowerMod


# standard logging/configuration setup
from loadConfig import config
import logging
logger = logging.getLogger(__name__)
if __name__ == "__main__":
    from loggerSetup import initialiseLogging
    initialiseLogging()



class ProcessFermiData:
    """
    Class to process Fermi GBM data for a given GRB. This class will process a given GRBData object into a CSV file with the binned photon counts. 
    
    Follows the method outlined in Bala et al. (2026)
    """
    def __init__(
            self,
            grb: GRBData,
            timeRange: tuple[float, float] = (-10, 60)
        ) -> None:
        """Constructor 

        Args:
            grb (GRBData): Instance of the GRBData class for the GRB to be processed.
        """
        self._grb: GRBData = grb
        self._detectorList: list[str] = ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb"]
        self._detectorData: dict[str, Phaii] = {}
        self._binSize: float = float(config.preProcessingConfig.fermiGMBConfig.processing.initialBinSize)
        self._timeRange: tuple[float, float] = timeRange
        self._energyRange: tuple[float, float] = (
            float(config.preProcessingConfig.fermiGMBConfig.download.energyRange.split("-")[0]),
            float(config.preProcessingConfig.fermiGMBConfig.download.energyRange.split("-")[1])
        )

        # initial read in
        self._readInTTEData()


        # find the optimal detector combination for the GRB
        self._generateDetectorCombinations()
        optimalDetectorCombination: tuple[str] = self._findOptimalDetectorCombination()
        logger.debug(f"Optimal detector combination: {optimalDetectorCombination}")


#        for detector, data in self._detectorData.items():
#            logger.debug(f"Detector: {detector}, Counts: {data['counts']}, Background CPS: {data['backgroundCPS']}, Bin Width: {data['binWidth']}")


        # bin to 100 μs


        # create a CSV file with the binned photon counts


    def _readInTTEData(
            self
        ) -> None:
        for detector in self._detectorList:
            try:
                self._detectorData[detector] = self._getSingleDetectorData(detector)
            except Exception as e:
                logger.error(f"Error processing detector {detector}: {e}")
                continue


    def _getSingleDetectorData(
            self,
            detector: str
        ) -> dict:
        """
        Retrieves, bins, fits background, and collapses data for a given detector.

        Returns:
            dict: Contains 'counts' (1D array), 'countUncertainties' (1D array), 'backgroundCPS' (float), and 'binWidth' (float).
        """
        # generate path for the TTE file
        filePath: str = f"{self._grb.folderPath}/current/glg_tte_{detector}_bn{self._grb.observationID}_v00.fit"

        # generate a strict time grid to ensure all detectors match shape
        startTime: float = self._timeRange[0]
        endTime: float = self._timeRange[1]
        timeGrid: np.ndarray = np.arange(startTime, endTime + self._binSize, self._binSize)

        # read and process the TTE event list within local scope
        logger.debug(f"Processing data for GRB {self._grb.name} from detector {detector} with file path {filePath}")
        tteData: GbmTte = GbmTte.open(filePath)

        logger.debug(f"Binning data for GRB {self._grb.name} from detector {detector}")
        tempPhaii: Phaii = tteData.to_phaii(
            binByEdges,
            timeGrid,
            time_range=self._timeRange,
            energy_range=self._energyRange
        )

        # convert the Phaii object into a lightcurve and extract the counts array
        lightCurve: Lightcurve = tempPhaii.to_lightcurve(energy_range=self._energyRange)
        extractedCounts: np.ndarray = np.array(lightCurve.counts)
        extractedUncertainties: np.ndarray = np.array(lightCurve.count_uncertainty)
        
        # get the scalar width of a single time bin
        scalarBinWidth: float = float(tempPhaii.data.time_widths[0])

        # Fit the background directly to the active Phaii object
        logger.debug(f"Fitting background to data for detector {detector}")
        fitter = BackgroundFitter.from_phaii(tempPhaii, Polynomial) 
        fitter.fit(order=1) 
        
#        logger.debug(f"Fitting statistic: {fitter.statistic_name} = {fitter.statistic}")
#        logger.debug(f"Fitting parameters: {fitter.parameters}")
#        logger.debug(f"Fitting degrees of freedom: {fitter.dof}")
        
        # interpolate the background to the exact time boundaries of the file
        backgroundRatesObject = fitter.interpolate_bins(
            tempPhaii.data.tstart,
            tempPhaii.data.tstop
        )
        
        # Extract the average background count rate (CPS) across your active window
        backgroundRates: np.ndarray = backgroundRatesObject.rates
        averageBackgroundCPS: float = float(np.mean(backgroundRates))

        # clear all references and wipe the heavy objects from RAM
        logger.debug(f"Purging GDT tracking structures from RAM for detector {detector}")
        tteData.close()
        
        del tteData, tempPhaii, fitter, backgroundRatesObject, lightCurve
        gc.collect()

        return {
            "counts": extractedCounts,
            "countUncertainties": extractedUncertainties,
            "backgroundCPS": averageBackgroundCPS,
            "binWidth": scalarBinWidth
        }


    def _generateDetectorCombinations(
            self
        ) -> None:
        """
        Generates detector combinations for the fermi mission. This will generate combinations up to a maximum specified in the config file.
        """
        maxCombinations: int = config.preProcessingConfig.fermiGMBConfig.processing.maxDetectorCombinations
        self._combinationsList: list[str] = []
        for i in range(1, maxCombinations + 1):
            self._combinationsList += list(combinations(self._detectorList, i))


    def _combineDetectorData(
            self,
            detectorCombination: tuple[str, ...]
        ) -> dict:
        """
        Combines datasets from multiple detectors into a single dataset.

        Args:
            detectorCombination (tuple[str]): The detector combination to combine.
        
        Returns:
            dict: The combined data containing total counts, combined count uncertainties, summed background CPS, and the baseline bin width.
        """
        # get the first detector's data as the baseline
        firstDetector: str = detectorCombination[0]
        combinedCounts: np.ndarray = np.array(self._detectorData[firstDetector]["counts"], dtype=float)
        combinedBackgroundCPS: float = float(self._detectorData[firstDetector]["backgroundCPS"])
        baselineBinWidth: float = float(self._detectorData[firstDetector]["binWidth"])
        combinedVariance: np.ndarray = np.array(self._detectorData[firstDetector]["countUncertainties"], dtype=float) ** 2

        # loop through and add the remaining detectors
        for detector in detectorCombination[1:]:
            combinedCounts += self._detectorData[detector]["counts"]
            combinedBackgroundCPS += self._detectorData[detector]["backgroundCPS"]
            combinedVariance += self._detectorData[detector]["countUncertainties"] ** 2

        # convert the combined variance back to standard deviation
        combinedUncertainties: np.ndarray = np.sqrt(combinedVariance)

        return {
            "counts": combinedCounts,
            "countUncertainties": combinedUncertainties,
            "backgroundCPS": combinedBackgroundCPS,
            "binWidth": baselineBinWidth
        }


    def _findPeakSNR(
            self,
            detectorCombination: tuple[str, ...],
            data: dict
        ) -> tuple[float, float, float, float]:
        """
        Calculates the peak signal-to-noise ratio for a given combined detector dataset.

        Args:
            detectorCombination (tuple[str]): The detector combination being analyzed.
            data (dict): A dictionary containing 'counts', 'countUncertainties', 'backgroundCPS', and 'binWidth'.
        
        Returns:
            tuple: (calculatedSNR, peakCounts, expectedBackgroundInBin, backgroundCPS)
        """
        try:
            backgroundCPS: float = data["backgroundCPS"]
            binWidthSec: float = data["binWidth"]
            countsArray: np.ndarray = data["counts"]
            countUncertainties: np.ndarray = data["countUncertainties"]

            # calculate the expected background counts inside a single bin
            expectedBackgroundInBin: float = backgroundCPS * binWidthSec
            
            # safety validation check
            if expectedBackgroundInBin <= 0 or not countsArray.any():
                return 0.0, 0.0, 0.0, 0.0, np.array([]), np.array([])
            
            # calculate the peak SNR
            peakCounts: float = np.max(countsArray)
            calculatedSNR: float = peakCounts / np.sqrt(expectedBackgroundInBin)

            # clean up values
            peakCounts = round(peakCounts, 2)
            expectedBackgroundInBin = round(expectedBackgroundInBin, 2)
            backgroundCPS = round(backgroundCPS, 2)

            return calculatedSNR, peakCounts, expectedBackgroundInBin, backgroundCPS, countsArray, countUncertainties

        except Exception as errorTrace:
            logger.error(f"An error occurred during final SNR calculation for detector combination {detectorCombination}: {errorTrace}")
            return 0.0, 0.0, 0.0, 0.0, np.array([]), np.array([])


    def _findHigestSNRCombination(
            self
        ) -> tuple[tuple[str, ...], float, float, float, float]:
        """
        Finds the detector combination with the highest peak SNR.

        Returns:
            tuple: (bestDetectorCombination, bestSNR, peakCounts, expectedBackgroundInBin, backgroundCPS)
        """
        bestSNR: float = 0.0
        bestDetectorCombination: tuple[str, ...] = ()
        bestPeakCounts: float = 0.0
        bestExpectedBackgroundInBin: float = 0.0
        bestBackgroundCPS: float = 0.0

        for detectorCombination in self._combinationsList:
            combinedData: dict = self._combineDetectorData(detectorCombination)
            calculatedSNR, peakCounts, expectedBackgroundInBin, backgroundCPS, countsArray, countUncertainties = self._findPeakSNR(detectorCombination, combinedData)

            if calculatedSNR > bestSNR:
                bestSNR = calculatedSNR
                bestDetectorCombination = detectorCombination
                bestPeakCounts = peakCounts
                bestExpectedBackgroundInBin = expectedBackgroundInBin
                bestBackgroundCPS = backgroundCPS
                bestCountsArray = countsArray
                bestCountUncertainties = countUncertainties

        return bestDetectorCombination, bestSNR, bestPeakCounts, bestExpectedBackgroundInBin, bestBackgroundCPS, bestCountsArray, bestCountUncertainties


    def _findOptimalDetectorCombination(
            self
        ) -> tuple[tuple[str, ...], float, float, float, float]:
        """
        Finds the optimal detector combination for the GRB based on the highest peak SNR.

        Returns:
            tuple[str]: The optimal detector combination
        """
        currentSNRValue: float = 0.0
        currentDetectorCombination: tuple[str, ...] = None
        while True:
            bestDetectorCombination, bestSNR, _, _, _, countsArray, countUncertainties = self._findHigestSNRCombination()
            SNRDifference: float = abs(bestSNR - currentSNRValue)
            logger.debug(f"Best SNR: {bestSNR} for detector combination: {bestDetectorCombination} with SNR difference: {SNRDifference}")
            if round(SNRDifference, 2) < 0.001:
                logger.info(f"Converged SNR: {bestSNR} for detector combination: {bestDetectorCombination}")
                return bestDetectorCombination
            currentSNRValue = bestSNR
            currentDetectorCombination = bestDetectorCombination
            logger.debug(f"Current SNR: {currentSNRValue} for detector combination: {currentDetectorCombination}")

            # find the MVT for the GRB using the optimal detector combination
            minimumVariabilityTimescale = haarPowerMod(
                countsArray,
                countUncertainties
            )[3]
            logger.debug(f"Calculated MVT: {minimumVariabilityTimescale} for detector combination: {currentDetectorCombination}")
            # set the bin size to the MVT
            self._binSize = minimumVariabilityTimescale
            logger.debug(f"Setting bin size to MVT: {minimumVariabilityTimescale}")

            # read in the TTE data again with the new bin size
            self._readInTTEData()
            # repeat until the SNR converges to a stable value

            # TODO This is really slow on the second pass, probably due to slow I/O on the network drive
            # Look into how big the tteData object in _getSingleDetectorData is and if it can be held in
            # memory. It also might not be. Have a look at np.histogram, might be able to do the binning
            # using a vectorized approach rather than using the GDT library.



from fermiDataTools.fermiGBMCatalogueGRB import getFermiGRBData
testGRBData = getFermiGRBData("GRB230307A")
testProcessFermiData = ProcessFermiData(testGRBData)