"""
download currently commented out.
"""

from swiftBATCatalogueGRB import getObservationID, importData
from swiftBATDataFetcher import downloadSwiftBATData, generateWgetStatement
from swiftBATInitalProcessing import processSwiftBATData as processSwiftBATEventData
import pexpect
from contextlib import chdir
from pathlib import Path

# function to convert the light curve data into a CSV file
def convertLightCurveToCSV(
        filename: str,
        GRBName: str
        ) -> None:
    # split the filename to get the directory
    directory: str = filename.rsplit("/", 1)[0]
    # open a shell
    with chdir(directory):
        shell = pexpect.spawn("sh", encoding='utf-8')
        logFile = open('shellOutput.txt', 'w')
        shell.logfile_read = logFile
        # initilaise the HEASoft tools
        shell.sendline('source $CALDB/software/tools/caldbinit.sh')
        shell.expect('CALDB/software/tools/caldbinit.sh') 
        shell.sendline('source $HEADAS/headas-init.sh')
        shell.expect('headas-init.sh') 
        # run the command to convert the light curve data into a CSV file
        shell.sendline('fdump output.lc outfile=output.txt prhead=no clobber=yes && echo Done')
        shell.expect('fdump output')
        shell.expect('(?i)Names')
        shell.send('\n')
        shell.expect('(?i)Lists')
        shell.send('\n')
        shell.expect('Done')
        shell.close()
    
    # create a new file path for the CSV file
    Path('data/processed').mkdir(parents=True, exist_ok=True)
    csvFilePath: str = f"data/processed/{GRBName}LC.csv"
    # read the output.txt file and write the data to the new CSV file
    with open(f'{directory}/output.txt', 'r') as outputFile, open(csvFilePath, 'w') as csvFile:
        lines: list[str] = outputFile.readlines()
        for line in lines:
            if line.strip() == "":
                continue
            line = line.strip()
            lineSplit: list[str] = line.split(" ")
            lineSplit = [field for field in lineSplit if field != ""]
            lineSplit.pop(0) # remove the first column, which is just the row number
            line = ",".join(lineSplit) + "\n"
            

            csvFile.write(line)



# function to process the data into a usable format for analysis
def processSwiftBATData(
        GRBName: str,
        timeBinSize: float,
        energyRange: str
        ) -> None:
    """Processes the data for a given GRB name into a usable format for analysis. This includes downloading the data from the swift.ac.uk archive, unzipping the data if it is compressed, and using the HEASoft tools to process the data into light curves for different energy ranges.

    Args:
        GRBName (str): The name of the GRB to be processed. Fomatted as "GRBYYMMDDX", where X is the letter assigned to the GRB. 
        timeBinSize (float): The size of the time bins, in seconds, to be used in the light curves.
        energyRange (str): The energy range to be used in the light curves. This should be in the format "min-max", where min and max are the minimum and maximum energies in keV. The BAT sensor has a energy range of 15-350 keV.

    Returns:
        None: This function does not return anything, but it will save the processed data to a file that can be used for analysis. The processed data will include light curves for the specified energy range. The file will be saved in the same directory as the original data, with the same name but with "_processed" appended to the end of the file name.
    """
    # generate the wget statement to download the data for the given GRB name
    wgetStatement: str = generateWgetStatement(GRBName)

    # download the data using the generated wget statement
    print(f"Downloading data for {GRBName}...")
#    downloadSwiftBATData(wgetStatement)
    print("Download complete.")

    # generate the file path for the downloaded data
    data, _ = importData('summary_general.csv')
    observationID: str = getObservationID(GRBName, data)
    dataFilePath: str = f"data/reproc/{observationID}/bat/event/sw{observationID}bevshsp_uf.evt.gz"

    # process the data using the HEASoft tools to generate light curves for the specified energy range and time bin size
    print(f"Processing data for {GRBName}...")
    processSwiftBATEventData(dataFilePath, timeBinSize, energyRange, GRBName)
    print("Processing complete.")

    # convert the processed data into a CSV file for analysis
    print(f"Converting processed data for {GRBName} into CSV format...")
    convertLightCurveToCSV(dataFilePath, GRBName)
    print("Conversion complete.")

if __name__ == "__main__":
    grbName: str = "GRB080319B"
    processSwiftBATData(grbName, 100e-6, "15-350")
