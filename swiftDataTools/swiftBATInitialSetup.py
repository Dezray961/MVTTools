from swiftDataTools import swiftBATCatalogueGRB
from pathlib import Path



# download summary_general.txt file from the Swift BAT website and update the data in the catalogue
def updateSummaryGeneralData() -> None:
    # download the summary_general.txt file from the Swift BAT website
    swiftBATCatalogueGRB.downloadSummaryGeneralFile()
    # update the data in the catalogue
    swiftBATCatalogueGRB.exportDataToCSV("summary_general.txt", 'summary_general.csv')


# create the directory structure for the data
def createDirectoryStructure() -> None:
    # create the directory structure for the data
    Path('data').mkdir(parents=True, exist_ok=True)
    Path('data/processed').mkdir(parents=True, exist_ok=True)




# main function
def main() -> None:
    createDirectoryStructure()
    updateSummaryGeneralData()


if __name__ == "__main__":
    main()