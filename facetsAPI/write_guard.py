"""Stops a Python process from changing anything outside the folders it was given.

Hard rule for generated jobs: they read FACETS repositories and never write, modify or delete anything in
them. File permissions can't enforce that, because the account running jobs can also write to the
repositories for other work. So a job installs this guard before doing anything else.

It uses an audit hook, which Python calls before each file operation and which can't be removed once added.
Writes are refused unless the real path, with symlinks resolved, is inside an allowed folder and outside
every protected folder. Starting other programs is refused too, as they would not be covered by the hook,
except the ASCETS R script when asked for: that code has been read and writes only to its output folder.

Standard library only.
"""
import os
import sys
import tempfile

_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND

# Audit event -> positions of the arguments that are paths being changed.
_PATH_EVENTS = {
    "os.mkdir": (0,), "os.remove": (0,), "os.rmdir": (0,), "os.rename": (0, 1), "os.link": (1,), "os.symlink": (1,),
    "os.truncate": (0,), "os.chmod": (0,), "os.chown": (0,), "os.utime": (0,), "os.chflags": (0,),
    "os.setxattr": (0,), "os.removexattr": (0,), "shutil.rmtree": (0,), "shutil.move": (0, 1),
    "shutil.copyfile": (1,), "shutil.copymode": (1,), "shutil.copystat": (1,), "shutil.copytree": (1,),
    "shutil.chown": (0,), "shutil.make_archive": (0,),
}
# Position of the dir_fd argument: a path relative to an open folder can't be checked, so it is refused.
_DIR_FD_ARG = {"os.mkdir": 2, "os.remove": 1, "os.rmdir": 1, "os.rename": 2, "os.chmod": 2, "os.chown": 3, "os.symlink": 2}

_PROCESS_EVENTS = ("os.system", "os.exec", "os.posix_spawn", "os.spawn", "os.fork", "os.forkpty", "os.startfile", "pty.spawn")

ASCETS_SCRIPT = "run_ascets.R"


class WriteBlocked(PermissionError):
    """A job tried to change something outside its own folders."""


def _real(path):
    return os.path.realpath(os.fsdecode(path))


def _inside(path, folder):
    return path == folder or path.startswith(folder.rstrip("/") + "/")


def enable_write_guard(allowed, protected=(), allow_ascets=False):
    """Allow writes only under the ``allowed`` folders and the temp folder, and never under ``protected``."""
    allowed = [_real(path) for path in allowed]
    protected = [_real(path) for path in protected]
    for folder in allowed:
        if any(_inside(folder, root) or _inside(root, folder) for root in protected):
            raise WriteBlocked("Refusing to run: %s overlaps a protected repository folder." % folder)
    writable = allowed + [_real(tempfile.gettempdir())]

    def check(path, event):
        if isinstance(path, int):
            return
        real = _real(path)
        if real == os.devnull:
            return
        if any(_inside(real, root) for root in protected) or not any(_inside(real, folder) for folder in writable):
            raise WriteBlocked("Blocked %s on %s: jobs may only write to their own output and job folders." % (event, real))

    def hook(event, args):
        if event == "open":
            path, mode, flags = args
            wants_write = (isinstance(mode, str) and any(c in mode for c in "wax+")) or (flags or 0) & _WRITE_FLAGS
            if wants_write and path is not None:
                check(path, event)
        elif event in _PATH_EVENTS:
            fd_position = _DIR_FD_ARG.get(event)
            if fd_position is not None and args[fd_position] not in (None, -1):
                raise WriteBlocked("Blocked %s relative to an open folder: the target can't be checked." % event)
            for position in _PATH_EVENTS[event]:
                check(args[position], event)
        elif event == "subprocess.Popen":
            command = args[1]
            is_ascets = (allow_ascets and isinstance(command, (list, tuple)) and len(command) == 3
                         and list(command[:2]) == ["bash", "-lc"] and ASCETS_SCRIPT in str(command[2]))
            if not is_ascets:
                raise WriteBlocked("Blocked starting a program: %s" % (command,))
        elif event in _PROCESS_EVENTS:
            raise WriteBlocked("Blocked %s: jobs may not start other programs." % event)

    sys.addaudithook(hook)
