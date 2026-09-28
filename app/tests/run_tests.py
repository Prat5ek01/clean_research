"""Unit and integration tests for cleanx (plain Python; pytest is not needed).

Run from app\\:   .\\cx.bat tests\\run_tests.py            (unit + integration, ~2-4 min)
                 .\\cx.bat tests\\run_tests.py --unit-only
Prints [PASS]/[FAIL] per test and ALL <n> TESTS PASSED at the end; exits 1 on any failure. The integration test
builds a synthetic app in a temporary folder and never writes into data/ or results/ of this repository."""
import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

MODULES = ["test_metrics_ext", "test_evalcard", "test_certify", "test_integration"]


def main():
    unit_only = "--unit-only" in sys.argv
    n_pass, failed = 0, []
    t_all = time.time()
    for name in MODULES:
        if unit_only and name == "test_integration":
            continue
        mod = importlib.import_module(name)
        for fn in sorted(k for k in dir(mod) if k.startswith("test_") and callable(getattr(mod, k))):
            t0 = time.time()
            try:
                getattr(mod, fn)()
                print(f"[PASS] {name}.{fn} ({time.time() - t0:.1f}s)", flush=True)
                n_pass += 1
            except Exception:
                print(f"[FAIL] {name}.{fn}", flush=True)
                traceback.print_exc()
                failed.append(f"{name}.{fn}")
    print()
    if failed:
        print(f"{len(failed)} of {n_pass + len(failed)} TESTS FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"ALL {n_pass} TESTS PASSED in {time.time() - t_all:.0f}s")


if __name__ == "__main__":
    main()
