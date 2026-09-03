"""
This will NOT run on a laptop due to the large memory requirements of the Fermi GBM data!
"""


from analysisTools.GRBData import GRBData
from gdt.core import data_path
from gdt.core.binning.unbinned import bin_by_time as binByTime
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.background.primitives import BackgroundRates
import numpy as np
import gc
from more_itertools import distinct_combinations as combinations


# standard logging/configuration setup
from loadConfig import config
import logging
logger = logging.getLogger(__name__)
if __name__ == "__main__":
    from loggerSetup import initialiseLogging
    initialiseLogging()


def getDetectorData(
        grb: GRBData,
        detector: str,
        binSize: float,
        timeRange: tuple[float, float],
        energyRange: tuple[float, float]
    ) -> Phaii:
    """
    Retrieves the data for a given GRB and detector and bins it into a Phaii object.

    Args:
        grb (GRBData): Instance of the GRBData class for the GRB to be processed.
        detector (str): The detector for which to process the data.
        binSize (float): The bin size to use for the data.
        timeRange (tuple[float, float]): The time range to use for the data.
        energyRange (tuple[float, float]): The energy range to use for the data.

    Returns:
        Phaii: The binned data for the given GRB and detector.
    """
    # generate the file path for the data file
    filePath: str = f"{grb.folderPath}/current/glg_tte_{detector}_bn{grb.observationID}_v00.fit"

    # read the data from the data file
    logger.debug(f"Processing data for GRB {grb.name} from detector {detector} with file path {filePath}")
    tteData: GbmTte = GbmTte.open(filePath)

    # bin the data into the given bin size for the given time range
    logger.debug(f"Binning data for GRB {grb.name} from detector {detector} with bin size {binSize}, time range {timeRange}, and energy range {energyRange}")
    phaii: Phaii = tteData.to_phaii(
        binByTime,
        binSize,
        time_range=timeRange,
        energy_range=energyRange
    )

    # clean the tte instance and trigger garbage collection to free up memory
    logger.debug(f"Cleaning up tteData instance for GRB {grb.name} from detector {detector}")
    tteData.close()
    del tteData
    gc.collect()

    return phaii




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
        self._getDetectorData()

        # find the optimal detector combination for the GRB
        self._generateDetectorCombinations()
        print(self._findPeakSNR(self._detectorData[self._detectorList[0]]))


        # bin to 100 μs


        # create a CSV file with the binned photon counts


    def _getDetectorData(
            self
    ) -> None:
        """
        Retrieves the data for each detector and bins it into a Phaii object.
        """
        for detector in self._detectorList:
            try:
                self._detectorData[detector] = getDetectorData(
                    self._grb,
                    detector,
                    self._binSize,
                    self._timeRange,
                    self._energyRange
                )
            except FileNotFoundError:
                logger.warning(f"Data file for detector {detector} not found for GRB {self._grb.name}. Skipping this detector.")
                continue


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
            detectorCombination: tuple[str]
        ) -> Phaii:
        """
        Combines the data from the given detector combination into a single Phaii object.

        Args:
            detectorCombination (tuple[str]): The detector combination to combine
        
        Returns:
            Phaii: The combined data for the given detector combination.
        """
        return Phaii.merge([self._detectorData[detector] for detector in detectorCombination])


    def _findPeakSNR(
            self,
            phaii: Phaii
        ) -> float:
        """
        Finds the peak signal-to-noise ratio for the given Phaii object.

        Args:
            phaii (Phaii): The Phaii object to find the SNR for.
        """
        # fit the background to the data
        fitter: BackgroundFitter = BackgroundFitter.from_phaii(phaii, Polynomial) # initialise the fitter
        fitter.fit(order = 1) # fit the background
        logger.debug("Fitted background to data for SNR calculation")
        logger.debug(f"Fitting statistic: {fitter.statistic_name} = {fitter.statistic}")
        logger.debug(f"Fitting parameters: {fitter.parameters}")
        logger.debug(f"Fitting degrees of freedom: {fitter.dof}")
        # interpolate the background to the data
        background: BackgroundRates = fitter.interpolate_bins(
            phaii.data.tstart,
            phaii.data.tstop
        )
        # calculate background level (CPS) from the rates model
        backgroundCPS = np.mean(background.rates)
        # extract the lightcurve
        lightCurveTotal = phaii.to_lightcurve(energy_range = self._energyRange)
        # get the bin width
        binWidthSec = phaii.data.time_widths[0]
        # calculate the SNR
        expectedBackgroundInBin = backgroundCPS * binWidthSec
        
        if expectedBackgroundInBin <= 0 or not lightCurveTotal.counts.any():
            return 0.0
        
        peakCounts = np.max(lightCurveTotal.counts)
        calculatedSnr = peakCounts / np.sqrt(expectedBackgroundInBin)

        peakCounts = round(peakCounts, 2)
        expectedBackgroundInBin = round(expectedBackgroundInBin, 2)
        backgroundCPS = round(backgroundCPS, 2)

        return calculatedSnr, peakCounts, expectedBackgroundInBin, backgroundCPS





from fermiDataTools.fermiGBMCatalogueGRB import getFermiGRBData
testGRBData = getFermiGRBData("GRB230307A")
testProcessFermiData = ProcessFermiData(testGRBData)