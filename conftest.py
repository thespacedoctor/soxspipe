from types import SimpleNamespace
import sys


sys.modules.setdefault("line_profiler", SimpleNamespace(profile=lambda func: func))
