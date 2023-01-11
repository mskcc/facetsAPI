## Dependencies

FacetsAPI requires the Pandas python library.  
If you do not have it installed, see instructions [here](https://pandas.pydata.org/docs/getting_started/install.html).

FacetsAPI also requires matplotlib. 

If you do not have it installed, see instructions [here](https://matplotlib.org/stable/users/installing/index.html).


## Installation

Facets API can be easilly installed using pip.

```
pip install git+https://github.com/mskcc/facetsAPI#facetsAPI
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

but this approach requires prefacing facetsAPI functions.  For instance `facetsAPI.buildFacetsSample()` instead of simply `buildFacetsSample()`. 
Either approach is valid, but the method of import will determine the syntax necessary for calling facetsAPI functions.
