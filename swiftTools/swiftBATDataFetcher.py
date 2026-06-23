"""
To be implimented later, dumping some things I have found.
Data needs to be downloaded from the www.swift.ac.uk archive.
They have the ability to generate a wget statement to download the 
data from a shell.
The statements are generated for each trigger ID that is included in the
address bar of the webpage. E.g.
https://www.swift.ac.uk/archive/download.sh?reproc=1&tid=00145675&source=obs&subdir=bat
I can work supply an address for the data I want to download, and use
their wget statement to download the data. 

The data is downloaded in a directory structure and the wget will need
to be run in the data folder to be compatible with initialProcessing.py

Only the initial 000 file has the event data

"""
import pexpect
from swiftTools.swiftBATCatalogueGRB import importData, getObservationID
from contextlib import chdir


# function to generate the wget statement for a given trigger ID or GRB name
def generateWgetStatement(ID: str) -> str:
    """Generates a wget statement to download data from the swife.ac.uk archive

    Args:
        ID (str): The trigger ID, observation ID or GRB name for the data to be downloaded

    Returns:
        str: A wget statement to download the data for the given trigger ID or GRB name
    """
    # read in the catalogue data to get the observation ID for the given GRB name
    data, _ = importData('summary_general.csv')
    
    observationID: str = ""
    # Check if the ID is a GRB name or a trigger ID
    if ID.startswith("GRB"):
        # Get the observation ID for the given GRB name
        for row in data:
            if row[0] == ID:
                observationID = row[18]
    else:
        observationID = getObservationID(ID, data, isTrigID=True)

    # Generate the wget statement
    wgetStatement = f'wget -e \'robots=off\' -nv -w 2 -nH --cut-dirs=1 -r --no-parent --reject "index.html*" https://www.swift.ac.uk/archive/reproc/{observationID}/bat/'
    return wgetStatement


# function to execute the wget statement to download the data
def downloadSwiftBATData(wgetStatement: str) -> None:
    """Executes the given wget statement to download data from the swift.ac.uk archive. Assumes the existance
        of a data directory to download the data into.

    Args:
        wgetStatement (str): The wget statement to be executed to download the data.
    """
    # change the working directory to the data folder. Returns to the original working directory after.
    with chdir("data"): 
        shell = pexpect.run(wgetStatement, timeout=None)
    

if __name__ == "__main__":
    grbName = "GRB080319B"
    wgetStatement = generateWgetStatement(grbName)
    print(wgetStatement)
    downloadSwiftBATData(wgetStatement)