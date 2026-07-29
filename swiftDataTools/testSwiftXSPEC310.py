import os, time, subprocess

# =====================================================================
# GLOBAL CONFIGURATION - Edit these paths once for your system
# =====================================================================
PYTHON_310_PATH = "/home/derekpinkett/.conda/envs/XspecEnv/bin/python"
SCRIPT_310_PATH = "/home/derekpinkett/coding/MVTTools/swiftDataTools/swiftXSPEC310.py"
TARGET_FOLDER   = "data/reproc/00306757000/bat/event"

# NEW: Explicitly specify your preferred output directory and filename
OUTPUT_DIR      = "/home/derekpinkett/coding/MVTTools/results" 
OUTPUT_FILENAME = "timeResolvedEpeakResults.csv"

TIMEOUT_SECONDS = 30
# =====================================================================

def testSingleFolderWithActiveTimeout():
    print("=== Python 3.13 Subprocess Test: Absolute Path Routing ===")
    
    # ensure the target output directory exists before spawning anything
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    # combine the directory and filename into an absolute path
    finalCsvResult = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)
    
    if os.path.exists(finalCsvResult):
        os.remove(finalCsvResult)

    # clone environment variables and inject the path attributes
    testEnv = os.environ.copy()
    testEnv["SLICE_ID"] = "1"
    testEnv["TARGET_DIR"] = TARGET_FOLDER
    testEnv["OUTPUT_CSV"] = finalCsvResult
    testEnv["PYTHONUNBUFFERED"] = "1"

    commandArray = [PYTHON_310_PATH, SCRIPT_310_PATH]
    print(f"\n--- Spawning Process (Hard Timeout Safety: {TIMEOUT_SECONDS}s) ---")
    
    startTime = time.time()
    try:
        process = subprocess.Popen(
            commandArray,
            env=testEnv,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        # make the read loop non-blocking to prevent kernel stalls
        os.set_blocking(process.stdout.fileno(), False)
        
        while True:
            pollStatus = process.poll()
            
            try:
                outputLine = process.stdout.readline()
                if outputLine:
                    cleanLine = outputLine.strip()
                    print(f"[3.10 Output]: {cleanLine}")
                    
                    # instant kill rules
                    if "Error:" in cleanLine or "FITTING EXCEPTION" in cleanLine:
                        print("\n[ABORT]: Caught direct failure flag. Terminating...")
                        process.kill()
                        return
            except Exception:
                pass

            if pollStatus is not None:
                break 
                
            # runtime verification loop
            if (time.time() - startTime) > TIMEOUT_SECONDS:
                print(f"\n[TIMEOUT EXPIRED]: Process exceeded {TIMEOUT_SECONDS}s. Force killing...")
                process.kill()
                print("-> Process successfully killed. Kernel protected.")
                return
                
            time.sleep(0.1) 
            
        print("\n" + "="*50)
        print(f"SUBPROCESS CLOSED (Exit Code: {process.returncode})")
        print("="*50)
        print(f"Verified: Output saved to {finalCsvResult}")

    except Exception as systemError:
        print(f"\nFailed to launch process pool: {str(systemError)}")

if __name__ == "__main__":
    testSingleFolderWithActiveTimeout()
