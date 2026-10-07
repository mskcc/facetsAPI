## Dependencies

FacetsAPI requires Python 3.9 or newer with the pandas, numpy and matplotlib libraries (installed automatically by pip).

R is only needed for the optional external tools (`ExtTools.runAscets`, `FPTools.runGenerateGenomicAnnotations`):
R 3.6 with tidyverse/data.table for ascets, and facetsPreview for genomic annotations. See [Configuration](configuration.md)
for telling facetsAPI how to reach `Rscript` on your system (for example through Lmod modules on IRIS).

## Installation

Facets API can be installed using pip:

```
pip install git+https://github.com/mskcc/facetsAPI#facetsAPI
```

or, from a checkout, as an editable install (recommended on the cluster so the R wrapper and ascets files are found automatically):

```
git clone https://github.com/mskcc/facetsAPI.git
pip install -e facetsAPI
```

## Loading the Package

The easiest way to load the package once it is installed is the following:

```
from facetsAPI import *
```

It is also possible to import the package by:

```
import facetsAPI
```

but this approach requires prefacing facetsAPI functions.  For instance `facetsAPI.FacetsMeta(...)` instead of simply `FacetsMeta(...)`.
Either approach is valid, but the method of import will determine the syntax necessary for calling facetsAPI functions.

## Running the tests

The unit tests build a synthetic FACETS repository (standard and 2n samples) in a temporary directory and need no cluster access:

```
python -m unittest discover -s tests -p 'test_*.py'
```

`tests/smoke_real_mount.py` is a manual, read-only check against a real repository:

```
python tests/smoke_real_mount.py /path/to/impact/facets/all/ /path/to/impact_2n/facets/all/
```
