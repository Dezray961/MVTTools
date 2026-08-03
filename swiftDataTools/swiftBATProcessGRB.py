"""
Complete pipe for processing Swift BAT data for a given GRB name. Steps:
1. Download the data from the swift.ac.uk archive using wget.
2. Process the data using SwiftProcessor to generate photon counts.
3. Calculate the time-resolved MVT using the photon counts and errors.
4. Calculate the time-resolved Epeak.
"""

from swiftDataTools.swiftBATCatalogueGRB import importData, getObservationID
from swiftDataTools.swiftBATDataFetcher import generateWgetStatement, downloadSwiftBATData
from swiftDataTools.swiftProcessor import ProcessSwiftData
from shutil import rmtree

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
        energyRange (str): The energy range to be used in the light curves. This should be in the format "min-max", where min and max are the minimum and maximum energies in keV. The BAT sensor has a energy range of 15-350 keV.
        download (bool): Whether to download the data from the swift.ac.uk archive. If False, it is assumed that the data has already been downloaded and is available in the "data/reproc" directory.
        deleteOriginal (bool): Whether to delete the original data file after processing. If True, the original data file will be deleted after the processed data has been saved.

    Returns:
        None: This function does not return anything, but it will save the processed data to a file that can be used for analysis. The processed data will include light curves for the specified energy range. The file will be saved in the same directory as the original data, with the same name but with "_processed" appended to the end of the file name.
    """

    # generate the wget statement to download the data for the given GRB name
    wgetStatement: str = generateWgetStatement(GRBName)

    # download the data using the generated wget statement
    if download:
        print(f"Downloading data for {GRBName}...")
        downloadSwiftBATData(wgetStatement)
        print("Download complete.")

    # generate the file path for the downloaded data
    data, _ = importData('summary_general.csv')
    observationID: str = getObservationID(GRBName, data)

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
        rmtree(f'data/reproc/{observationID}')
        print("Deletion complete.")

if __name__ == "__main__":
    grbName: str = "GRB080319B"
    processSwiftBATData(
        GRBName = grbName,
        energyRange = "15-350",
        download = False,
        deleteOriginal = False,
        )
