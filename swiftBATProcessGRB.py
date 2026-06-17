"""
This needs to be split into a terminal running script with a class that can be imported into other
scripts. It also needs to use 'with chdir()' to change the working directory rather than using 
os as it will automatically change back to the original working directory when the block is exited.
Further, the SwiftBAT tools should be in their own directory.
"""

import swiftBATCatalogueGRB as catalogue
import swiftBATDataFetcher as fetcher
import swiftBATInitalProcessing as processing
import pexpect
from contextlib import chdir
from pathlib import Path
import pandas as pd
import shutil

# function to process the data into a usable format for analysis
def processSwiftBATData(
        GRBName: str,
        SNRThreshold: float,
        energyRange: str,
        timeBinSize: float = 0.0,
        download: bool = False,
        deleteOriginal: bool = False,
        fullDataSet: bool = False
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

    # function to convert the light curve data into a CSV file
    def convertLightCurveToCSV(
            filename: str,
            GRBName: str
            ) -> None:
        # split the filename to get the directory
        directory: str = filename.rsplit("/", 1)[0]
        # open a shell
        with chdir(directory):
            shell = pexpect.spawn("sh", encoding = 'utf-8')
            logFile = open('shellOutput.txt', 'w')
            shell.logfile_read = logFile
            # initilaise the HEASoft tools
            shell.sendline('source $CALDB/software/tools/caldbinit.sh')
            shell.expect('CALDB/software/tools/caldbinit.sh') 
            shell.sendline('source $HEADAS/headas-init.sh')
            shell.expect('headas-init.sh') 
            # run the command to convert the light curve data into a CSV file
            shell.sendline('fdump output.lc outfile=output.txt prhead=no clobber=yes && echo Done')
            shell.expect('fdump output', timeout = None)
            shell.expect('(?i)Names')
            shell.send('\n')
            shell.expect('(?i)Lists')
            shell.send('\n')
            shell.expect('Done', timeout = None)
            shell.close()
        
        # create a new file path for the CSV file
        Path('data/processed').mkdir(parents=True, exist_ok=True)
        csvFilePath: str = f"data/processed/{GRBName}LC.csv"
        # read the output.txt file and write the data to the new CSV file
        with open(f'{directory}/output.txt', 'r') as outputFile, open(csvFilePath, 'w') as csvFile:
            lines: list[str] = outputFile.readlines()

        # find the gap lines
        gapLines: list[int] = [i for i, line in enumerate(lines) if not line.strip()]
        # The first two are blank, the data rows are between the 2nd and 3rd, and the 3rd and 4th
        block1Lines: list[str] = lines[gapLines[1]:gapLines[2]]
        block2Lines: list[str] = lines[gapLines[2]:gapLines[3]]
        # clean the data
        def linesToDataFrame(lines: list[str]) -> list[str]:
            data = [line.strip().split() for line in lines if line.strip()]
            # change the column names to lower case and remove the units
            data[0] = [col.split('(')[0].lower() for col in data[0]]
            data.pop(1)
            # seperate the index column into a separate list
            indexList = [row[0] for row in data[1:]]
            for i, list in enumerate(data):
                if i == 0:
                    continue
                list.pop(0)
                    
            return pd.DataFrame(data[1:], columns=data[0], index=indexList)
        

        cleanedBlock1Lines: list[str] = linesToDataFrame(block1Lines)
        cleanedBlock2Lines: list[str] = linesToDataFrame(block2Lines)

        # combine the two blocks of data side by side, using the index to align the rows
        mergedData: pd.DataFrame = pd.concat([cleanedBlock1Lines, cleanedBlock2Lines], axis=1)
        mergedData.to_csv(csvFilePath)


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
    dataFilePath: str = f"data/reproc/{observationID}/bat/event/sw{observationID}bevshsp_uf.evt.gz"

    # process the data using the HEASoft tools to generate light curves for the specified energy range and time bin size
    print(f"Processing data for {GRBName}...")
    processing.processSwiftBATData(
        filename = dataFilePath,
        SNRThreshold = SNRThreshold,
        energyBins = energyRange,
        GRBName = GRBName,
        timeDeliniation = timeBinSize,
        fullDataSet = fullDataSet
    )
    print("Processing complete.")

    # convert the processed data into a CSV file for analysis
    print(f"Converting processed data for {GRBName} into CSV format...")
    convertLightCurveToCSV(dataFilePath, GRBName)
    print("Conversion complete.")

    # delete the original data file if specified
    if deleteOriginal:
        print(f"Deleting original data file for {GRBName}...")
        shutil.rmtree(f'data/reproc/{observationID}')
        print("Deletion complete.")

if __name__ == "__main__":
    grbName: str = "GRB080319B"
    processSwiftBATData(
        GRBName = grbName,
        SNRThreshold = 5.0,
        energyRange = "15-350",
        download = False,
        deleteOriginal = False,
        #timeBinSize = 100e-6,
        fullDataSet = True
        )
