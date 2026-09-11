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
        Retrieves raw event timestamps within energy bounds and fits a continuous background model. Discards heavy objects to allow fast custom rebinning later.

        Args:
            detector (str): The name of the detector to process e.g. 'n0', 'n1', 'n2', etc.

        Returns:
            dict: Contains 'rawTimes' (1D array), 'backgroundModel' (GDT object) and 'initialBinWidth' (float).
        """
        # generate path for the TTE file
        filePath: str = f"{self._grb.folderPath}/current/glg_tte_{detector}_bn{self._grb.observationID}_v00.fit"

        # read and process the TTE event list within local scope
        logger.debug(f"Processing data for GRB {self._grb.name} from detector {detector} with file path {filePath}")
        tteData: GbmTte = GbmTte.open(filePath)

        # filter the TTE data by energy
        tteDataFiltered = tteData.slice_energy(self._energyRange)

        # extract the raw event arrival times (unbinned)
        extractedRawTimes: np.ndarray = np.array(tteDataFiltered.data.time)

        # create a coarse time grid for background fitting
        coarseStartTime: float = self._timeRange[0]
        coarseEndTime: float = self._timeRange[1]
        coarseGrid: np.ndarray = np.arange(coarseStartTime, coarseEndTime + self._binSize, self._binSize)

        logger.debug(f"Binning coarse temporary data for background fit on detector {detector}")
        tempPhaii: Phaii = tteDataFiltered.to_phaii(
            binByEdges,
            coarseGrid,
            time_range=self._timeRange
        )

        # fit the background model to the temporary coarse structure
        logger.debug(f"Fitting continuous background model for detector {detector}")
        fitter = BackgroundFitter.from_phaii(tempPhaii, Polynomial) 
        fitter.fit(order=1) 

        # interpolate the fitted background model to create a continuous background model
        continuousBackgroundModel = fitter.interpolate()

        # clear all references and wipe the heavy objects from RAM
        logger.debug(f"Purging GDT tracking structures from RAM for detector {detector}")
        tteData.close()
        tteDataFiltered.close()

        del tteData, tteDataFiltered, tempPhaii, fitter
        gc.collect()

        return {
            "rawTimes": extractedRawTimes,
            "backgroundModel": continuousBackgroundModel,
            "initialBinWidth": float(self._binSize)
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
        Dynamically bins and combines raw datasets from multiple detectors based on the current class iteration bin width.

        Args:
            detectorCombination (tuple[str, ...]): The detector combination to combine e.g. ('n0', 'n1', 'n2').
        
        Returns:
            dict: The combined data containing total counts, combined count uncertainties, summed expected background counts, and the iteration bin width.
        """
        # generate the strict, updated time boundaries based on the current bin size
        startTime: float = self._timeRange[0]
        endTime: float = self._timeRange[1]
        binEdges: np.ndarray = np.arange(startTime, endTime + self._binSize, self._binSize)
        
        # calculate bin centers and widths
        binCenters: np.ndarray = (binEdges[:-1] + binEdges[1:]) / 2.0
        binWidths: np.ndarray = binEdges[1:] - binEdges[:-1]

        # initialise blank arrays matching the active timeline resolution
        totalBins: int = len(binEdges) - 1
        combinedCounts: np.ndarray = np.zeros(totalBins, dtype=float)
        combinedBackgroundCounts: np.ndarray = np.zeros(totalBins, dtype=float)
        combinedVariance: np.ndarray = np.zeros(totalBins, dtype=float)

        # loop across all targeted detectors
        for detector in detectorCombination:
            detectorData:  = self._detectorData[detector]
            detectorCounts, _ = np.histogram(detectorData["rawTimes"], bins=binEdges)
            combinedCounts += detectorCounts
            # re-calculate individual variance contribution (Poisson variance = counts)
            combinedVariance += detectorCounts
            # continuous background model rates (CPS)
            detectorRates: np.ndarray = detectorData["backgroundModel"].rate(binCenters)
            # convert rates to expected counts
            combinedBackgroundCounts += (detectorRates * binWidths)

        # convert the total variance array back to standard deviation
        combinedUncertainties: np.ndarray = np.sqrt(combinedVariance)

        return {
            "counts": combinedCounts,
            "countUncertainties": combinedUncertainties,
            "background": combinedBackgroundCounts,
            "binWidth": float(self._binSize)
        }


    def _findPeakSNR(
            self,
            detectorCombination: tuple[str, ...],
            data: dict
        ) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
        """
        Calculates the peak signal-to-noise ratio for a given combined detector dataset using time-resolved background count arrays.

        Args:
            detectorCombination (tuple[str, ...]): The detector combination being analyzed e.g. ('n0', 'n1', 'n2').
            data (dict): A dictionary containing 'counts', 'countUncertainties', 'background', and 'binWidth'.
        
        Returns:
            tuple: (calculatedSNR, peakCounts, expectedBackgroundInPeakBin, averageBackgroundCPS, countsArray, countUncertainties)
        """
        try:
            countsArray: np.ndarray = data["counts"]
            countUncertainties: np.ndarray = data["countUncertainties"]
            backgroundArray: np.ndarray = data["background"]
            binWidthSec: float = data["binWidth"]

            # safety validation check to avoid dividing by zero or empty datasets
            if not countsArray.any() or backgroundArray.min() <= 0:
                return 0.0, 0.0, 0.0, 0.0, np.array([]), np.array([])
            
            # calculate the SNR timeline
            snrTimeline: np.ndarray = countsArray / np.sqrt(backgroundArray)
            
            # find the peak SNR index
            peakIndex: int = int(np.argmax(snrTimeline))
            calculatedSNR: float = float(snrTimeline[peakIndex])
            peakCounts: float = float(countsArray[peakIndex])
            expectedBackgroundInPeakBin: float = float(backgroundArray[peakIndex])
            
            # reconstruct the average background counts per second (CPS) across the entire time range
            totalTimeDuration: float = len(backgroundArray) * binWidthSec
            averageBackgroundCPS: float = float(np.sum(backgroundArray) / totalTimeDuration)

            # clean up values
            peakCounts = round(peakCounts, 2)
            expectedBackgroundInPeakBin = round(expectedBackgroundInPeakBin, 2)
            averageBackgroundCPS = round(averageBackgroundCPS, 2)

            return calculatedSNR, peakCounts, expectedBackgroundInPeakBin, averageBackgroundCPS, countsArray, countUncertainties

        except Exception as errorTrace:
            logger.error(f"An error occurred during final SNR calculation for detector combination {detectorCombination}: {errorTrace}")
            return 0.0, 0.0, 0.0, 0.0, np.array([]), np.array([])


    def _findHigestSNRCombination(
            self
        ) -> tuple[tuple[str, ...], float, float, float, float, np.ndarray, np.ndarray]:
        """
        Finds the detector combination with the highest peak SNR.

        Returns:
            tuple: (bestDetectorCombination, bestSNR, bestPeakCounts, 
                    bestExpectedBackgroundInBin, bestBackgroundCPS, 
                    bestCountsArray, bestCountUncertainties)
        """
        bestSNR: float = 0.0
        bestDetectorCombination: tuple[str, ...] = ()
        bestPeakCounts: float = 0.0
        bestExpectedBackgroundInBin: float = 0.0
        bestBackgroundCPS: float = 0.0
        
        # initialise empty arrays
        bestCountsArray: np.ndarray = np.array([])
        bestCountUncertainties: np.ndarray = np.array([])

        for detectorCombination in self._combinationsList:
            combinedData: dict = self._combineDetectorData(detectorCombination)
            
            # compute the time-resolved SNR
            (calculatedSNR, 
            peakCounts, 
            expectedBackgroundInBin, 
            backgroundCPS, 
            countsArray, 
            countUncertainties) = self._findPeakSNR(detectorCombination, combinedData)

            # tracking step
            if calculatedSNR > bestSNR:
                bestSNR = calculatedSNR
                bestDetectorCombination = detectorCombination
                bestPeakCounts = peakCounts
                bestExpectedBackgroundInBin = expectedBackgroundInBin
                bestBackgroundCPS = backgroundCPS
                bestCountsArray = countsArray
                bestCountUncertainties = countUncertainties

        return (
            bestDetectorCombination, 
            bestSNR, 
            bestPeakCounts, 
            bestExpectedBackgroundInBin, 
            bestBackgroundCPS, 
            bestCountsArray, 
            bestCountUncertainties
        )


    def _findOptimalDetectorCombination(
            self
        ) -> tuple[str, ...]:
        """
        Finds the optimal detector combination for the GRB based on the highest peak SNR.
        Converges the bin size using the Minimum Variability Timescale (MVT) completely in memory.

        Returns:
            tuple[str, ...]: The optimal detector combination.
        """
        currentSNRValue: float = 0.0
        currentDetectorCombination: tuple[str, ...] = None
        loopCounter: int = 0
        maxIterations: int = 15  # Safety threshold to prevent infinite loops

        # --- ONE-TIME FILE INITIALISATION ---
        # Read the raw files EXACTLY ONCE at the start. This builds your raw timestamps 
        # and continuous background model cache inside self._detectorData.
        logger.info(f"Initialising raw data cache and background models for GRB {self._grb.name}...")
        self._readInTTEData() 

        while loopCounter < maxIterations:
            loopCounter += 1
            
            # This now executes in milliseconds because everything is processed via in-memory math
            (bestDetectorCombination, 
            bestSNR, _, _, _, 
            countsArray, 
            countUncertainties) = self._findHigestSNRCombination()
            
            snrDifference: float = abs(bestSNR - currentSNRValue)
            logger.debug(f"Iteration {loopCounter} - Best SNR: {bestSNR} for combo: {bestDetectorCombination} with difference: {snrDifference}")
            
            # Check for mathematical convergence
            if round(snrDifference, 3) < 0.001:
                logger.info(f"🎉 Converged successfully after {loopCounter} passes! SNR: {bestSNR} for combination: {bestDetectorCombination}")
                return bestDetectorCombination
            
            # Update state values for the next convergence evaluation step
            currentSNRValue = bestSNR
            currentDetectorCombination = bestDetectorCombination
            logger.debug(f"Current state set - SNR: {currentSNRValue} for combo: {currentDetectorCombination}")

            # Find the MVT for the GRB using the winner of this pass
            # Note: Confirm that HaarPowerMod accepts your flat NumPy variance arrays directly.
            minimumVariabilityTimescale = haarPowerMod(
                countsArray,
                countUncertainties
            )[3]
            logger.debug(f"Calculated MVT: {minimumVariabilityTimescale} for detector combination: {currentDetectorCombination}")
            
            # Set the class variable bin size to the new MVT width
            self._binSize = minimumVariabilityTimescale
            logger.debug(f"Setting class active bin size to MVT: {minimumVariabilityTimescale}")

            # --- THE SPEED FIX ---
            # DO NOT call self._readInTTEData() here! 
            # Because _combineDetectorData reads self._binSize dynamically on every loop execution,
            # np.histogram will automatically re-slice the cached data at this new resolution next pass.

        logger.warning(f"Reached maximum iteration limit ({maxIterations}) without perfect convergence. Returning best guess.")
        return currentDetectorCombination




from fermiDataTools.fermiGBMCatalogueGRB import getFermiGRBData
testGRBData = getFermiGRBData("GRB230307A")
testProcessFermiData = ProcessFermiData(testGRBData)