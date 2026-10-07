"""Runtime configuration for the external (R) tools facetsAPI can drive.

facetsAPI itself is pure Python. ``ExtTools.runAscets`` and ``FPTools.runGenerateGenomicAnnotations``
call ``Rscript``; where R lives differs per system, so the invocation is configurable:

    r_prefix          shell snippet run before Rscript, in the SAME shell. On the IRIS cluster:
                      "module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3"
                      (Lmod's ``module`` is an exported shell function, available to ``bash -lc``.)
    rscript           the Rscript executable (default "Rscript")
    fp_wrapper_script path to R/fp_generate_genomic_anno_wrapper.R
    fp_config_file    facets-preview config json (required by FPTools; no default)
    r_lib_loc         R library directory holding facetsPreview, passed to the wrapper
    ascets_dir        directory holding run_ascets.R and the genomic_arm_coordinates_*.txt files

Precedence: explicit constructor argument > environment variable > default.

    FACETSAPI_R_PREFIX, FACETSAPI_RSCRIPT, FACETSAPI_FP_WRAPPER, FACETSAPI_FP_CONFIG,
    FACETSAPI_R_LIB_LOC, FACETSAPI_ASCETS_DIR

Every R call is composed into ONE command string and run with ``bash -lc`` so the prefix and
Rscript share an environment. (The old ``os.system("module load ...")`` ran in its own throwaway
shell and never affected the following Rscript call.)
"""

import os
import shlex
import subprocess
from dataclasses import dataclass, fields, replace

from .exceptions import ExternalToolError, FacetsConfigError

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_DIR = os.path.normpath(os.path.join(_PKG_DIR, ".."))
DEFAULT_FP_WRAPPER = os.path.join(_REPO_DIR, "R", "fp_generate_genomic_anno_wrapper.R")
DEFAULT_ASCETS_DIR = os.path.join(_REPO_DIR, "tools", "ascets")

ENV_VARS = {
    "r_prefix":          "FACETSAPI_R_PREFIX",
    "rscript":           "FACETSAPI_RSCRIPT",
    "fp_wrapper_script": "FACETSAPI_FP_WRAPPER",
    "fp_config_file":    "FACETSAPI_FP_CONFIG",
    "r_lib_loc":         "FACETSAPI_R_LIB_LOC",
    "ascets_dir":        "FACETSAPI_ASCETS_DIR",
}


@dataclass
class FacetsConfig:
    r_prefix: str = ""
    rscript: str = "Rscript"
    fp_wrapper_script: str = DEFAULT_FP_WRAPPER
    fp_config_file: str = None
    r_lib_loc: str = None
    ascets_dir: str = DEFAULT_ASCETS_DIR

    @classmethod
    def from_env(cls, environ=None, **overrides):
        """Build a config from defaults, then environment variables, then explicit overrides."""
        environ = os.environ if environ is None else environ
        values = {}
        for f in fields(cls):
            env_name = ENV_VARS.get(f.name)
            if env_name and environ.get(env_name) not in (None, ""):
                values[f.name] = environ[env_name]
        for k, v in overrides.items():
            if k not in ENV_VARS:
                raise FacetsConfigError("FacetsConfig has no setting %r" % k)
            if v is not None:
                values[k] = v
        return cls(**values)

    def with_(self, **changes):
        return replace(self, **changes)

    def compose_r_command(self, args):
        """The exact shell string run for an Rscript call (arguments are shell-quoted)."""
        if not self.rscript:
            raise FacetsConfigError("FacetsConfig.rscript is empty.")
        cmd = " ".join(shlex.quote(str(a)) for a in [self.rscript] + list(args))
        if self.r_prefix and self.r_prefix.strip():
            cmd = self.r_prefix.strip() + " && " + cmd
        return cmd

    def run_r(self, args, cwd=None, check=True):
        """Run Rscript with ``args`` through ``bash -lc``. Raises ExternalToolError on nonzero exit when check."""
        cmd = self.compose_r_command(args)
        try:
            result = subprocess.run(["bash", "-lc", cmd], cwd=cwd, capture_output=True, text=True, check=False)
        except OSError as e:
            raise ExternalToolError(cmd, -1, str(e)) from e
        if check and result.returncode != 0:
            raise ExternalToolError(cmd, result.returncode, result.stderr)
        return result


__all__ = ["FacetsConfig", "ENV_VARS", "DEFAULT_FP_WRAPPER", "DEFAULT_ASCETS_DIR"]
