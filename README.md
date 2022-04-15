# OUenv: Climatic/Environment dependent Ornstein-Uhlenbeck model of trait evolution

version 1.1  J. Clavel & A. Brinkworth
 
This model fits a generalized Ornstein-Uhlenbeck process (or Hull-White model) of trait evolution on phylogenetic trees where the optimum <img src="https://render.githubusercontent.com/render/math?math=\theta(t) "> is a function of an environmental or climatic variable through time (Eq. 1).

<img src="https://render.githubusercontent.com/render/math?math=dX(t) = \alpha (\theta(t) - X(t))dt + \sigma dB(t)"> Eq. (1)


Use of 'fit_t_env' function:
==

"tree" - is a phylo object (phylogenetic tree) from "ape" package.

"data" - a named vector of traits values for each species in "tree".

"fun" - is a three parameters function specifying the relationship between the optimum and the climatic curve. For instance, fun <- function(x, par, theta0) theta0 + par*Temperatures(x) ; where "x" is a time component (with x=0 is present days), "par" is a parameter controling the relationship to some envitonmental/climatic variable (Temperatures() here) to be estimated by the model, and "theta0" is an intercept like parameter that will be estimated by the model.

"startvalues" - is a vector containing some starting values for the optimization of the model - i.e. some guesses provided by the user for the parameter(s) "par" used in "fun".

"echo" - logical, whether messages should be prompted.

"method" - Optimization method used by the function (see ?optim). If method="spg",  large scale Barzilai-Borwein algorithm (from BB package) is used.

"nuisance" - logical, if TRUE a nuisance parameter is estimated for accounting for intraspecific/measurement error variance during model fit.

"mserr" - (optional) a vector with standard error of the mean for the traits values.

"control" - a list with control parameters for optim' methods.

Notes
==

This model requires (is identifiable on) a non-ultrametric time-tree. For instance, a phylogenetic tree with fossil species.

Example (tutorial)
==

TODO


References
==

Clavel & Morlon - Proc. Nat. Acad. Sci., 114(16): 4183-4188
Troyer et al. - in review
Clavel et al. - in prep.

