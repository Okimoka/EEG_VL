"""Small file helpers for dataset views that must never write through aliases."""
import os
from pathlib import Path
import tempfile


def real_directory(path):
    """Create directories while refusing every symlinked path component."""
    path = Path(os.path.abspath(path))
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"Destination directory must not be a symlink: {current}")
        current.mkdir(exist_ok=True)
    return path


def write_independent(path, text):
    """Atomically replace sidecars; unchanged independent files retain mtime.

    Replacing the directory entry breaks source hardlinks/file symlinks without
    changing the source. No existing destination file is opened for writing.
    """
    path = Path(path)
    real_directory(path.parent)
    content = text.encode("utf-8")
    if (not path.is_symlink() and path.is_file() and path.stat().st_nlink == 1
            and path.read_bytes() == content):
        return
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("wb", dir=path.parent, prefix=".view-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def link_source(source, target):
    """Link unchanged source files, preserving correct existing links."""
    source, target = Path(source), Path(target)
    real_directory(target.parent)
    if target.is_symlink():
        if target.resolve() == source.resolve():
            return
        target.unlink()
    elif target.exists():
        raise FileExistsError(f"Expected a source symlink, found an existing file: {target}")
    target.symlink_to(os.path.relpath(source, target.parent))


def remove_linked_locks(root):
    """Discard only inherited lock symlinks, never a process's real lock file.

    Runtime lock files must be local to a view. The reader opens them with
    O_NOFOLLOW, so a symlink cannot be used as its active file lock.
    """
    for path in Path(root).rglob("*.lock"):
        if path.is_symlink():
            path.unlink()
