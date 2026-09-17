"""THROWAWAY INSTRUMENTATION — NEVER MERGE. RECORDS WHICH EXCEPTION HANDLERS FIRE."""
import json
import os
import threading

_LOCK = threading.Lock()
_DEFAULT_PATH = "/tmp/soxspipe_exc_probe.jsonl"


def record(
        filePath,
        lineNumber,
        exception):
    """RECORD ONE HANDLER FIRING. MUST NEVER RAISE AND MUST NEVER CHANGE CONTROL FLOW."""
    try:
        outputPath = os.environ.get("SOXSPIPE_EXC_PROBE", _DEFAULT_PATH)
        row = {
            "file": filePath,
            "line": lineNumber,
            "pid": os.getpid(),
            "exc_type": type(exception).__name__ if exception is not None else None,
            "exc_module": type(exception).__module__ if exception is not None else None,
            "exc_msg": str(exception)[:300] if exception is not None else None,
        }
        payload = json.dumps(row) + "\n"
        with _LOCK:
            with open(outputPath, "a") as probeFile:
                probeFile.write(payload)
    except Exception:
        pass
