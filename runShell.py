import subprocess
import sys

# class to open and run shell commands
class ShellRunner:
    # constructor to initialize the shell runner
    def __init__(self) -> None:
        self.__caldbinitString: str = ". $CALDB/software/tools/caldbinit.sh"
        self.__headasinitString: str = ". $HEADAS/software/tools/headasinit.sh"
        self.__openShell()  # open the shell process
        tempShellCommandOutput: list[str] = self.tempShellCommand()  # run the temp shell command

        
    def shellOutputWrapper(self, function: callable) -> list[str]:
        # wrapper function to run a shell command and get the output
        def wrapper(*args, **kwargs) -> list[str]:
            # run the shell command
            function(*args, **kwargs)
            # get the output from the shell process
            output: list[str] = []
            for line in iter(self.__shellProcess.stdout.readline, ''):
                output.append(line.strip())
            sys.stdout.flush()  # flush the output buffer
            return output
        return wrapper


    @shellOutputWrapper
    def tempShellCommand(self) -> None:
        # example shell command to run
        self.__shellProcess.stdin.write("echo 'Running shell command...'\n")
        self.__shellProcess.stdin.write("ls -l\n")
        self.__shellProcess.stdin.write("echo 'Shell command finished.'\n")
        self.__shellProcess.stdin.flush()  # flush the input buffer


    def __openShell(self) -> None:
        # open a new shell process as an atrubute of the class
        self.__shellProcess: subprocess.Popen = subprocess.Popen(
            ["sh"], # run a shell
            shell=True, # use the shell
            stdin=subprocess.PIPE, # pipe the input
            stdout=subprocess.PIPE, # pipe the output
            stderr=subprocess.PIPE, # pipe the error
            text=True, # use text mode for input/output
            bufsize=1, # line-buffered
        )

        






if __name__ == "__main__":
    shell: ShellRunner = ShellRunner()

