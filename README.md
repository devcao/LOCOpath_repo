# LOCOpath project code repo
This github repo contains all project code for the LOCO path high dimensional inference.
Apart from the code listed here, you may also need to install the R package [```LOCOpath```](https://github.com/devcao/LOCOpath).

## Install the R package

To install, please first install R package ```devtools``` and then 
```R
devtools::install_github("devcao/LOCOpath")
```

## Road map for each file.

[NetTS.R](./NetTS.R):  LOCO path statistic calculations for linear/logistic/Poisson regression, backend R package ```glmnet```.

[NetResampleTS.R](./NetResampleTS.R):  Bootstrapping and power simulation code based on LOCO path statistic for **linear** regression, backend R package ```glmnet```.

[NetResampleLogisticTS.R](./NetResampleLogisticTS.R):  Bootstrapping and power simulation code based on LOCO path statistic for **logistic/Poisson** regression, backend R package ```glmnet```, ```gglasso``` and [Logistic_Enet.R](./Logistic_Enet.R) .

[Logistic_Enet.R](./Logistic_Enet.R): Modified coordinate descent for logistic regression, enabling constraint like beta_1=1 while looping the coordinate descent, backend C++ code [lognet.cpp](./lognet.cpp) .

[graphLASSO.R](./graphLASSO.R): All the graphical LASSO code, including our wrapper of ```glasso``` package, our own constraint graphical LASSO code and variable screening code, backend R package ```glasso```.

[compare_power.R](./compare_power.R): All the power simulation codes for other method we need to compare, including desparsified LASSO, T/F/Wald test.

## Power simulation functions road map and some examples
### high-dimensional linear regression
```R
Path.Resample.Power
```

### high-dimensional logistic/Poisson regression
```R
Path.Resample.Power
```

### sparse graphical models
```R
Path.Resample.Power
```

## Real data analysis

## Some code to run on Slurm cluster




