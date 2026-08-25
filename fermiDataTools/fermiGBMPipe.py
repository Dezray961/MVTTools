from analysisTools.GRBData import GRBData
from astroquery.heasarc import Heasarc
from fermiDataTools.fermiGBMCatalogueGRB import getFermiGRBData
from astropy.table import Table
from loadConfig import config
from shutil import rmtree
from pathlib import Path

# standard logging/configuration setup
from loadConfig import config
import logging
logger = logging.getLogger(__name__)
if __name__ == "__main__":
    from loggerSetup import initialiseLogging
    initialiseLogging()

def downloadFermiGBMData(
        grb: GRBData
    )-> None:
    """
    Downloads data from the Fermi GBM archive for a given GRB. Assumes the existance of a `data/fermiGBM` directory to download the data into. Initially downloads a tar file and then extracts the contents of the tar file into `data/fermiGBM/{grb.name}/current` and `data/fermiGBM/{grb.name}/quicklook` folders.

    Args:
        grb (GRBData): The GRB data for which to download the data.
    """
    # get the links to the data files for the GRB from the HEASARC database
    linksTable: Table = Heasarc.locate_data(
        grb.table,
        'fermigbrst'
    )
    # generate the folder path to download the data into
    folderPath: str = f"{config.generalSettings.directories.dataPath}/fermiGBM/{grb.name}"
    # if the folder path does not exist, create it
    if not Path(folderPath).exists():
        Path(folderPath).mkdir(parents=True, exist_ok=True)
    # download the data files into the folder path
    Heasarc.download_data(
        linksTable,
        host = 'heasarc',
        location = folderPath
    )


def processFermiGBMData(
        grb: GRBData
    ) -> None:
    """
    Processes the data for a given GRB name into a usable format for analysis. This includes downloading the data from the Fermi GBM archive, unzipping the data if it is compressed, and using the HEASoft tools to process the data into light curves for different energy ranges.

    Args:
        grb (GRBData): The GRB data for which to process the data.
    """
    # download the data using wget
    if config.preProcessingConfig.general.downloadData:
        logger.info(f"Downloading data for {grb.name}...")
        # download the data using the wget statement
        downloadFermiGBMData(grb)
        logger.info("Download complete.")

    # process the data
    logger.info(f"Processing data for {grb.name}...")
    # ... (processing logic here)
    logger.info("Processing complete.")

    # delete the original data file if specified
    if config.preProcessingConfig.general.deleteWhenDone:
        logger.info(f"Deleting original data for {grb.name}...")
        rmtree(f"{config.generalSettings.directories.dataPath}/fermiGBM/{grb.name}")
        logger.info("Deletion complete.")


if __name__ == "__main__":
    testGRBData = getFermiGRBData("GRB230307A")
    processFermiGBMData(testGRBData)
