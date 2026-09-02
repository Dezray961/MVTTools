from analysisTools.GRBData import GRBData
from gdt.core import data_path
from gdt.core.binning.unbinned import bin_by_time as binByTime
from gdt.missions.fermi.gbm.tte import GbmTte
import numpy as np
from gdt.core.phaii import Phaii
import gc


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
        timeRange: tuple[float, float]
    ) -> Phaii:
    """
    Retrieves the data for a given GRB and detector and bins it into a Phaii object.

    Args:
        grb (GRBData): Instance of the GRBData class for the GRB to be processed.
        detector (str): The detector for which to process the data.
        timeRange (tuple[float, float]): The time range to use for the data.
    
    Returns:
        Phaii: The binned data for the given GRB and detector.
    """
    binSize: float = float(config.preProcessingConfig.fermiGMBConfig.processing.initialBinSize)
    energyRangeStrings: tuple[str, str] = config.preProcessingConfig.fermiGMBConfig.download.energyRange.split("-")
    energyRange: tuple[float, float] = (float(energyRangeStrings[0]), float(energyRangeStrings[1]))
    timeRange: tuple[float, float] = timeRange

    # generate the file path for the data file
    fileName: str = f"glg_tte_{detector}_bn{grb.observationID}_v00.fit"
    filePath: str = f"{grb.folderPath}/current/{fileName}"

    # read the data from the data file
    tteData: GbmTte = GbmTte.open(filePath)

    # bin the data into the given bin size for the given time range
    phaii: Phaii = tteData.to_phaii(
        binByTime,
        binSize,
        time_range=timeRange
    )

    # clean the tte instance to free up memory
    tteData.close()
    del tteData

    # trigger garbage collection to free up memory
    gc.collect()

    return phaii




class ProcessFermiData:
    """
    Class to process Fermi GBM data for a given GRB. This class will process a given GRBData object into a CSV file with the binned photon counts. 
    
    Follows the method outlined in Bala et al. (2026)
    """
    def __init__(
            self,
            grb: GRBData
        ) -> None:
        """Constructor 

        Args:
            grb (GRBData): Instance of the GRBData class for the GRB to be processed.
        """
        self._grb: GRBData = grb
        self._detectorList: list[str] = ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb"]
        self._detectorData: dict[str, Phaii] = {}

        # read in the data from the files
        for detector in self._detectorList:
            try:
                self._detectorData[detector] = getDetectorData(self._grb, detector, (-10, 60))
            except FileNotFoundError:
                logger.warning(f"Data file for detector {detector} not found for GRB {self._grb.name}. Skipping this detector.")
                continue


        # find the optimal detector combination for the GRB


        # bin to 100 μs


        # create a CSV file with the binned photon counts



from fermiDataTools.fermiGBMCatalogueGRB import getFermiGRBData
testGRBData = getFermiGRBData("GRB230307A")
testProcessFermiData = ProcessFermiData(testGRBData)