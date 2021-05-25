# Global Variables

In `Xtals.jl`, global variables are stored in a dictionary with which the user may interact by way of a simple API.
This is intended to improve the user experience by reducing the number of independent environment variables to track, as well as assist in writing and using 
libraries dependent upon `Xtals`.

# Docstrings

```@docs
    get_global
    set_global
    get_global_dict
```
