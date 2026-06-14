"""
Need to find which shell to use for the current OS.
    def configure_shell(self) -> None:
        """Configures non-login/login shell to initialize with heainit command"""

        installdir = os.getcwd()
        config = 0

        # Loops through each specified shell initialization file
        for conf in self.shell_config:
            path = os.path.expandvars(os.path.join(*conf.split()))

            # Checks if shell configuration file is present and writes to file if present
            if os.path.exists(path):
                self.__write_script(path, installdir)
                config += 1
                break

        # Creates login configuration file if not present and writes there
        if not config:
            path = os.path.expandvars(os.path.join(*self.shell_config[-1].split()))
            os.makedirs(os.path.dirname(path), exist_ok=True)
            self.__write_script(path, installdir)
            print(
                f"\n{path} was not found. "
                f"{os.path.basename(path)} has been created and can be "
                f"found at {os.path.dirname(path)}. Make sure that {os.path.basename(path)} "
                f"is read by {self.def_shell} before execution"
            )

        return

This is from the HEAInstaller (see readme.md for reference). It looks like it finds the shell by checking the shell 
configuration files for the current OS. If it finds one, it writes a script to that file. 

Use config.json for a list of shells and conme commands for them. The "source" key is really the one that I need to use. 
Unsure if the rest is useful. 

Config.json is a heavily modified version of the original config.json file from HEAInstaller. I will need to assign 
attribution to the original author of HEAInstaller in the readme.md file. Do I need to ask permission to do this?
What is the license for HEAInstaller? I will need to check this.

"""