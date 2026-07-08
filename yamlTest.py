# to install yaml: conda install -c anaconda pyyaml


import yaml

configFilePath = "config.yaml"
directories = yaml.safe_load(open(configFilePath, 'r'))['generalSettings']['directories']

print("Data Path:", directories['dataPath'])
print("Processed Data Path:", directories['processedDataPath'])
print("Shared Objects Path:", directories['sharedObjectsPath'])
