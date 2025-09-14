import time
import os
import psutil
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import getpass
import argparse # For handling command-line arguments

CURRENT_USER = getpass.getuser()
TARGET_SUBDIRS = ["GENIE3", "GRNBOOST2"]
TARGET_FILENAME = "outFile.txt" 

class GRNFileHandler(FileSystemEventHandler):
    def __init__(self, monitor_root):
        super().__init__()
        self.monitor_root = monitor_root

    def on_created(self, event):
        if event.is_directory:
            return

        rel_path = os.path.relpath(event.src_path, self.monitor_root)

        for method in TARGET_SUBDIRS:
            expected_path = os.path.join(method, TARGET_FILENAME)
            if rel_path.endswith(expected_path):
                dataset_path = os.path.dirname(os.path.dirname(event.src_path))
                dataset_rel_path = os.path.relpath(dataset_path, self.monitor_root)
                print(f"[Watcher] Detected {method} outFile at {event.src_path}")
                time.sleep(10)
                kill_grn_process_for_dataset(dataset_rel_path, method)
                break


def kill_grn_process_for_dataset(dataset_path, method):
    print(f"[Watcher] Attempting to kill {method} processes for dataset: {dataset_path}")
    killed_any = False

    for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'cwd', 'username']):
        try:
            if proc.info['username'] == CURRENT_USER:
                cmdline_str = ' '.join(proc.info.get('cmdline', [])) or ''
                cwd = proc.info.get('cwd') or ''
                if method.lower() in cmdline_str.lower():
                    if dataset_path in cmdline_str or dataset_path in cwd:
                        if "time" not in cmdline_str.lower():
                            print(f"[Watcher] Killing {method} process {proc.pid} for dataset {dataset_path}")
                            print(f"    └─ CMD: {cmdline_str}")
                            print(f"    └─ CWD: {cwd}")
                            proc.kill()
                            time.sleep(1)
                            killed_any = True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass

    if not killed_any:
        print(f"[Watcher] No {method} process found for dataset {dataset_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Monitor a specified output directory and kill GRN processes upon completion."
    )
    parser.add_argument(
        "monitor_directory",
        help="The root output directory to monitor recursively (e.g., '/home/user/Beeline/outputs')."
    )
    args = parser.parse_args()

    event_handler = GRNFileHandler(monitor_root=args.monitor_directory)
    observer = Observer()
    
    observer.schedule(event_handler, path=args.monitor_directory, recursive=True)
    observer.start()

    print(f"[Watcher] Monitoring GRN method outputs under {args.monitor_directory}...")

    try:
        while True:
            time.sleep(2)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()