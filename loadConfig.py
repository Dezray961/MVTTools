# to install yaml: conda install -c anaconda pyyaml

def importConfiguration(configFilePath: str = "config.yaml"):
    """Load the configuration from a YAML file and return it as a DictToClass object.
    
    Args:
        configFilePath (str): Path to the YAML configuration file. Defaults to "config.yaml".

    Returns:
        DictToClass: An object that allows attribute-style access to the configuration settings.

    Note:
        Example usage:
        ```
            config = importConfiguration()
            dataPath: str = config.generalSettings.directories.dataPath
        ```
    """
    import yaml

    
    class DictToClass(dict):
        def __getattr__(self, name):
            try:
                value = self[name]
                if isinstance(value, dict):
                    value = DictToClass(value)
                return value
            except KeyError:
                raise AttributeError(f"'DictToClass' object has no attribute '{name}'")


    configYAML = yaml.safe_load(open(configFilePath, 'r'))
    return DictToClass(configYAML)





if __name__ == "__main__":
    config = importConfiguration()
    print("Configuration loaded successfully.")
    print(config.generalSettings.directories.dataPath)