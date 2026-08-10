"""
Tools for working with the Swift BAT GRB catalogue. Currently set up to use a "data" list of lists as
a catalogue, but should really be modified to use a class instead. The downloadSummaryGeneralFile
function should probably be moved to initial setup. As should exportDataToCSV.
"""
from pathlib import Path
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
GRB_CATALOGUE_DIR = REPO_ROOT / "data" / "GRBCatalogue"


def _resolveCataloguePath(filename: str) -> Path:
    """Resolve a catalogue file without depending on the current working directory."""
    path = Path(filename)
    candidates = [path]

    if not path.is_absolute():
        candidates.extend([
            REPO_ROOT / path,
            GRB_CATALOGUE_DIR / path.name,
            Path(__file__).resolve().parent / path,
        ])

    for candidate in candidates:
        if candidate.exists():
            return candidate

    return candidates[-1]

# function to import the data from a csv file
def importData(
        filename: str
        ) -> tuple[list[list[str]], list[str]]:
    """
    Imports the data from a CSV file and returns it as a list of lists along with the column names.

    Args:
        filename (str): The path to the CSV file.

    Returns:
        tuple[list[list[str]], list[str]]: A tuple containing the data and column names.
    """
    cataloguePath = _resolveCataloguePath(filename)
    df = pd.read_csv(cataloguePath).fillna('').convert_dtypes()
    data = df.values.tolist()
    columnNames = df.columns.tolist()
    return data, columnNames


# function to get the coordinates of a GRB given its Trig_ID
def getCoordinates(
        GRBName: str,
        data: list[list[str]]
        ) -> tuple[float, float]:
    """
    Gets the coordinates (RA, Dec) for a given GRB name from the data.

    Args:
        GRBName (str): The name of the GRB to get coordinates for.
        data (list[list[str]]): The data containing the GRB information.

    Returns:
        tuple[float, float]: The RA and Dec for the given GRB name.
    """
    for row in data:
        if row[0] == GRBName:
            return float(row[4]), float(row[5])
        

# function to get the start/stop time of a GRB given its GRB name
def getStartStopTime(
        grbName: str,
        data: list[list[str]]
        ) -> tuple[float, float, float]:
    """Gets the start, stop and t90 times for a given GRB name

    Args:
        grbName (str): The name of the GRB to get the start and stop times for
        data (list[list[str]]): The data containing the GRB information

    Returns:
        tuple[float, float, float]: The start, stop, and t90 times for the GRB
    """
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


if __name__ == "__main__":
    data, columnNames = importData('summary_general.csv')
    print(getCoordinates('145675', data))
    print(getStartStopTime('GRB080319B', data))
