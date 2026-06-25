"""
This needs to be split into a terminal running script with a class that can be imported into other
scripts. It also needs to use 'with chdir()' to change the working directory rather than using 
os as it will automatically change back to the original working directory when the block is exited.
Further, the SwiftBAT tools should be in their own directory.
"""

import swiftDataTools.swiftBATCatalogueGRB as catalogue
import swiftDataTools.swiftBATDataFetcher as fetcher
from swiftDataTools.swiftProcessor import ProcessSwiftData
from contextlib import chdir
from pathlib import Path
import pandas as pd
import shutil
import os

# function to process the data into a usable format for analysis
def processSwiftBATData(
        GRBName: str,
        energyRange: str,
        download: bool = False,
        deleteOriginal: bool = False,
        ) -> None:
    """Processes the data for a given GRB name into a usable format for analysis. This includes downloading the data from the swift.ac.uk archive, unzipping the data if it is compressed, and using the HEASoft tools to process the data into light curves for different energy ranges.

    Args:
        GRBName (str): The name of the GRB to be processed. Fomatted as "GRBYYMMDDX", where X is the letter assigned to the GRB. 
        SNRThreshold (float): The signal-to-noise ratio threshold to be used in the light curves bins. This will determine the size of the time bins in the light curves, with higher SNR thresholds resulting in larger time bins. If timeBinSize is provided, this will be ignored and the time bins will be of the specified size instead.
        energyRange (str): The energy range to be used in the light curves. This should be in the format "min-max", where min and max are the minimum and maximum energies in keV. The BAT sensor has a energy range of 15-350 keV.
        timeBinSize (float): The size of the time bins, in seconds, to be used in the light curves.
        download (bool): Whether to download the data from the swift.ac.uk archive. If False, it is assumed that the data has already been downloaded and is available in the "data/reproc" directory.
        deleteOriginal (bool): Whether to delete the original data file after processing. If True, the original data file will be deleted after the processed data has been saved.

    Returns:
        None: This function does not return anything, but it will save the processed data to a file that can be used for analysis. The processed data will include light curves for the specified energy range. The file will be saved in the same directory as the original data, with the same name but with "_processed" appended to the end of the file name.
    """

    # generate the wget statement to download the data for the given GRB name
    wgetStatement: str = fetcher.generateWgetStatement(GRBName)

    # download the data using the generated wget statement
    if download:
        print(f"Downloading data for {GRBName}...")
        fetcher.downloadSwiftBATData(wgetStatement)
        print("Download complete.")

    # generate the file path for the downloaded data
    data, _ = catalogue.importData('summary_general.csv')
    observationID: str = catalogue.getObservationID(GRBName, data)

    # process the data using the HEASoft tools to generate light curves for the specified energy range and time bin size
    print(f"Processing data for {GRBName}...")
    ProcessSwiftData(
        GRBName = GRBName,
        energyBins = energyRange
    )

    print("Processing complete.")

    # delete the original data file if specified
    if deleteOriginal:
        print(f"Deleting original data file for {GRBName}...")
        shutil.rmtree(f'data/reproc/{observationID}')
        print("Deletion complete.")

if __name__ == "__main__":
    grbName: str = "GRB080319B"
    processSwiftBATData(
        GRBName = grbName,
        energyRange = "15-350",
        download = True,
        deleteOriginal = False,
        )
