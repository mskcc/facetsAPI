"""Exception hierarchy for facetsAPI.

Library code raises these instead of calling ``sys.exit()`` so that a long-running
caller (a dashboard listener, a notebook, another package) can catch and report
failures. Per-sample soft failures (a sample missing a file, a malformed fit) are
NOT exceptions; they are counted and skipped by the builders.
"""


class FacetsAPIError(Exception):
    """Base class for every error raised by facetsAPI."""


class FacetsConfigError(FacetsAPIError):
    """Bad constructor arguments or configuration (invalid fit_class, bad mode, missing config file)."""


class FacetsRepositoryError(FacetsAPIError):
    """The FACETS repository or the clinical sample file cannot be used (missing path, unreadable file)."""


class ManifestError(FacetsRepositoryError):
    """A facets_review.manifest is unreadable or lacks required columns."""


class ParseError(FacetsAPIError):
    """A FACETS output file we chose to load has malformed content."""


class ExternalToolError(FacetsAPIError):
    """An external command (Rscript via bash) exited nonzero."""

    def __init__(self, cmd, returncode, stderr=""):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr or ""
        tail = self.stderr.strip().splitlines()[-20:]
        msg = "command exited with status %s: %s" % (returncode, cmd)
        if tail:
            msg += "\n--- stderr (tail) ---\n" + "\n".join(tail)
        super().__init__(msg)


__all__ = [
    "FacetsAPIError",
    "FacetsConfigError",
    "FacetsRepositoryError",
    "ManifestError",
    "ParseError",
    "ExternalToolError",
]
