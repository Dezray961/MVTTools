import pandas as pd
import requests
import os
from contextlib import chdir

# decorator to change the current working directory
def chdirDecorator(function):
    def wrapper(*args, **kwargs):
        currentDir: str = os.getcwd()
        if not currentDir.endswith("data"):
            # slit the current directory into a list
            currentDirSplit: list[str] = currentDir.split("/")
            # find the index of the "MVTTools" directory
            dataIndex: int = currentDirSplit.index("MVTTools")
            # join the list back into a string up to the "MVTTools" directory
            chdirString: str = "/".join(currentDirSplit[:dataIndex + 1])
        else:
            chdirString: str = currentDir
        with chdir(chdirString):
            return function(*args, **kwargs)
    return wrapper
        



# Function to return a observation id for a given GRB name
def getObservationID(grbName: str, data: list[list[str]], isTrigID: bool = False) -> str:
    if isTrigID:
        trigID: str = grbName
    else:
        # find the Trig_ID for the given GRB name
        trigID: str = ""
        for row in data:
            if row[0] == grbName:
                trigID = row[1]
    # change the trigID to a string
    match trigID:
        case int():
            trigID = str(trigID)
        case float():
            trigID = str(int(trigID))
        case str():
            pass
        case _:
            raise ValueError(f"Trig_ID is of type {type(trigID)}, expected int, float, or str.")
    # return the observation id for the given Trig_ID
    if len(trigID) < 11:
        trigID = str(trigID) + "000"
        while len(trigID) < 11:
            trigID = "0" + trigID
    return trigID


def convertSummaryGeneralToList(filename: str) -> tuple[list[list[str]], list[list[str]]]:
    # function to read in text file
    def readFile(filename: str) -> tuple[list[str], list[str], list[str]]:
        with open(filename, 'r') as f:
            lines: list[str] = f.readlines()

        discriptionLines: list[str] = []
        dataLines: list[str] = []
        # seperate the file discrption lines from the data lines
        for i, line in enumerate(lines):
            if line.startswith('#'):
                discriptionLines.append(line)
            else:
                dataLines.append(line)

        # the last line of the discription lines is the column names
        columnNames: list[str] = [discriptionLines.pop(-1).strip("## ")]
        return discriptionLines, dataLines, columnNames


    # function to clean the data lines
    def cleanData(dataLines: list[str]) -> list[list[str]]:
        cleanedData: list[list[str]] = []
        for line in dataLines:
            line = line.split("|")
            line = [field.strip() for field in line]
            cleanedData.append(line)
        return cleanedData


    # function to append the observation id to the data
    def appendObservationID(data: list[list[str]], columnNames: list[list[str]]) -> tuple[list[list[str]], list[list[str]]]:
        columnNames[0].append("Observation_ID")
        for row in data:
            observationID = getObservationID(row[0], data)
            row.append(observationID)
        return data, columnNames


    discriptionLines, data, columnNames = readFile(filename)
    data = cleanData(data)
    columnNames = cleanData(columnNames)
    data, columnNames = appendObservationID(data, columnNames)
    return  data, columnNames


# function to updatet the data from a new summary_general.txt file
@chdirDecorator
def exportDataToCSV(filename: str, csvFilename: str) -> None:
    # function to write the data to a csv file
    def writeToCSV(data: list[list[str]], columnNames: list[list[str]], filename: str) -> None:
        df = pd.DataFrame(data, columns=columnNames[0])
        df.to_csv(filename, index=False)

    
    data, columnNames = convertSummaryGeneralToList(filename)
    writeToCSV(data, columnNames, csvFilename)


# function to import the data from a csv file
def importData(filename: str) -> tuple[list[list[str]], list[str]]:
    currentDir: str = os.getcwd()
    if currentDir.endswith("event"):
        with chdir("../../../../../"):
            df = pd.read_csv(filename).fillna('').convert_dtypes()
            data = df.values.tolist()
            columnNames = df.columns.tolist()
    else:
        df = pd.read_csv(filename).fillna('').convert_dtypes()
        data = df.values.tolist()
        columnNames = df.columns.tolist()
    return data, columnNames


# function to get the coordinates of a GRB given its Trig_ID
def getCoordinates(trigID: str, data: list[list[str]]) -> tuple[float, float]:
    for row in data:
        if row[1] == int(trigID):
            return float(row[4]), float(row[5])
        

# function to get the start/stop time of a GRB given its GRB name
def getStartStopTime(grbName: str, data: list[list[str]]) -> tuple[float, float, float]:
    triggerTime: float = 0.0
    t90Time: float = 0.0
    t90Error: float = 0.0
    for row in data:
        if row[0] == grbName:
            triggerTime = float(row[2])
            t90Time = float(row[8])
            t90Error = float(row[9])
    stopTime: float = triggerTime + t90Time
    return triggerTime, stopTime, t90Error


# function to download the summary_general.txt file from the Swift BAT website
def downloadSummaryGeneralFile() -> None:
    url: str = "https://swift.gsfc.nasa.gov/results/batgrbcat/summary_cflux/summary_general_info/summary_general.txt"
    response = requests.get(url)
    with open("summary_general.txt", "wb") as f:
        f.write(response.content)




if __name__ == "__main__":
    filename = 'summary_general.txt'
    data, columnNames = convertSummaryGeneralToList(filename)
    exportDataToCSV(filename, 'summary_general.csv')
    for i, columnName in enumerate(columnNames[0]):
        print(f"{i}: {columnName}")
    print(data)

    data, columnNames = importData('summary_general.csv')
    print(getCoordinates('145675', data))
    print(getStartStopTime('GRB080319B', data))
