"""Require committed analysis code for confirmatory runs and record its digests."""
import hashlib
import subprocess
from pathlib import Path


def _git(args, cwd):
    r = subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip() or f"git {' '.join(args)} failed")
    return r.stdout


def _app():
    return Path(__file__).resolve().parents[1]


def code_state(paths, root=None):
    """{'commit': HEAD hash or None, 'dirty': [uncommitted, untracked or outside-repository files], 'git': bool}.
    root: the folder the paths belong to (default: app/ of this repository; tests pass a temporary repository)."""
    app = Path(root).resolve() if root else _app()
    inside, outside = [], []
    for p in paths:
        p = Path(p).resolve()
        try:
            inside.append(p.relative_to(app).as_posix())
        except ValueError:
            outside.append(str(p))
    try:
        commit = _git(["rev-parse", "HEAD"], app).strip()
        status = _git(["status", "--porcelain", "--", *inside], app) if inside else ""
    except (OSError, RuntimeError):
        return {"commit": None, "dirty": inside + outside, "git": False}
    dirty = [ln[3:].strip() for ln in status.splitlines() if ln.strip()]
    return {"commit": commit, "dirty": dirty + [f"{p} (outside the repository)" for p in outside], "git": True}


def require_frozen(paths, allow_dirty, what, root=None):
    """Returns a label for the outputs: 'confirmatory @<commit>' or 'NOT CONFIRMATORY (...)'.

    Stops the run when the code is not committed and allow_dirty is False."""
    st = code_state(paths, root)
    if st["git"] and not st["dirty"]:
        return f"confirmatory @{st['commit'][:12]}"
    why = ("git unavailable" if not st["git"] else "uncommitted changes in " + ", ".join(st["dirty"]))
    if not allow_dirty:
        raise SystemExit(f"{what}: the confirmatory run needs the pre-registered analysis code committed ({why}). "
                         f"Commit it first; --allow-dirty is for development and tests only.")
    label = f"NOT CONFIRMATORY ({why})"
    print(f"  WARNING: {label}", flush=True)
    return label


def digests(paths, root=None):
    """{path relative to app/ (or absolute): SHA-256} of the analysis files."""
    app, out = (Path(root).resolve() if root else _app()), {}
    for p in paths:
        p = Path(p).resolve()
        try:
            key = p.relative_to(app).as_posix()
        except ValueError:
            key = str(p)
        out[key] = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "missing"
    return out


def app_paths(*names):
    """Absolute paths of analysis files inside app/ (the real repository, not CLEANX_APP_DIR)."""
    return [_app() / n for n in names]
