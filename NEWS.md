# wdnet (development version)

+ Updated default values of the arguments in `rpanet`.
  + `seednetwork = NULL` to 
  `seednetwork = list(edgelist = matrix(c(1, 2), nrow = 1))`.
  + `control = NULL` to `control = list()`.
+ Changed `rpactl` to `rpacontrol`.

# wdnet 0.0.5

## Minor changes

+ Updated returned items from `dprewire` and `dprewire.range`.
  + Removed solved `eta` and corresponding assortativity levels.
  + Added solver results from `CVXR`.
+ Allowed user-defined preference functions in `rpactl.preference`.

# wdnet 0.0.4

## Minor changes

+ Renamed internal C++ functions.
+ Edited `rpanet` with `binary` approach: renamed node structures (fix LTO issues).


# wdnet 0.0.3

+ First version of wdnet submitted to CRAN.
