"""
Tools for working with the Swift BAT GRB catalogue. Currently set up to use a "data" list of lists as
a catalogue, but should really be modified to use a class instead.
"""
from pathlib import Path
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
GRB_CATALOGUE_DIR = REPO_ROOT / "data" / "GRBCatalogue"


class GRBData:
    """
    A class to represent a GRB data entry. Use SwiftGRBCatalogue.getGRBData(grbName) to get an instance of this class for a given GRB name.

    Attributes:
        name (str): The name of the GRB.
        triggerTime (float): The trigger time of the GRB.
        ra (float): The right ascension of the GRB.
        dec (float): The declination of the GRB.
        t90 (float): The T90 duration of the GRB.
        t90Error (float): The error in the T90 duration of the GRB.
    """

    def __init__(
            self,
            name: str,
            triggerTime: float,
            ra: float,
            dec: float,
            t90: float,
            t90Error: float
            ) -> None:
        self.name = name
        self.triggerTime = triggerTime
        self.ra = ra
        self.dec = dec
        self.t90 = t90
        self.t90Error = t90Error


    def __repr__(self) -> str:
        return f"GRBData(name={self.name}, triggerTime={self.triggerTime}, ra={self.ra}, dec={self.dec}, t90={self.t90}, t90Error={self.t90Error})"

class SwiftGRBCatalogue:
    """
    A class to represent a GRB catalogue.

    Attributes:
        filename (str): The path to the GRB catalogue file.
        data (list[list[str]]): The data from the GRB catalogue.
        columnNames (list[str]): The column names from the GRB catalogue.
    """
    def __init__(
            self,
            filename: str
            ) -> None:
        self.filename = filename
        self.data, self.columnNames = self._importData()


    def _resolveCataloguePath(
            self
            ) -> Path:
        """Resolve a catalogue file without depending on the current working directory."""
        path = Path(self.filename)
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


    def _importData(
            self
            ) -> tuple[list[list[str]], list[str]]:
        """
        Imports the data from a CSV file and returns it as a list of lists along with the column names.

        Args:
            filename (str): The path to the CSV file.

        Returns:
            tuple[list[list[str]], list[str]]: A tuple containing the data and column names.
        """
        cataloguePath = self._resolveCataloguePath()
        df = pd.read_csv(cataloguePath).fillna('').convert_dtypes()
        data = df.values.tolist()
        columnNames = df.columns.tolist()
        return data, columnNames


    def _getCoordinates(
            self,
            GRBName: str
            ) -> tuple[float, float]:
        """
        Gets the coordinates (RA, Dec) for a given GRB name from the data.

        Args:
            GRBName (str): The name of the GRB to get coordinates for.
            data (list[list[str]]): The data containing the GRB information.

        Returns:
            tuple[float, float]: The RA and Dec for the given GRB name.
        """
        for row in self.data:
            if row[0] == GRBName:
                return float(row[4]), float(row[5])
            

    def _getStartStopTime(
            self,
            grbName: str
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
        for row in self.data:
            if row[0] == grbName:
                triggerTime = float(row[2])
                t90Time = float(row[8])
                t90Error = float(row[9])
        stopTime: float = triggerTime + t90Time
        return triggerTime, stopTime, t90Error


    def getGRBData(
            self,
            grbName: str
            ) -> GRBData:
        """
        Gets the GRB data for a given GRB name.

        Args:
            grbName (str): The name of the GRB to get data for.
        
        Returns:
            GRBData: An instance of the GRBData class containing the data for the given GRB name.
        """
        triggerTime, stopTime, t90Error = self._getStartStopTime(grbName)
        ra, dec = self._getCoordinates(grbName)
        t90 = stopTime - triggerTime
        return GRBData(grbName, triggerTime, ra, dec, t90, t90Error)



if __name__ == "__main__":
    catalogue = SwiftGRBCatalogue("summary_general.csv")
    GRB080319A = catalogue.getGRBData("GRB080319A")
    print(GRB080319A)