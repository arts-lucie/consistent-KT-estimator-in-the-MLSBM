# Code for consistent model selection in the multi-layer SBM with the penalized Krichevsky-Trofimov estimator

## Description  
 
The codes presented here were used for the experiments in my article *Consistent Model Selection in a Collection of Stochastic Block Models* <!--([preprint link]())-->. They are used to estimate the number of communities in a multi-layer SBM based on the Krichevsky-Trofimov penalized estimator.

These codes are inspired by the [**mixer**](https://cran.r-project.org/web/packages/mixer/index.html) package, developed by Christophe Ambroise, Gilles Grasseau, Mark Hoebeke, Pierre Latouche, Vincent Miele, Franck Picard, Alexander Smith, the LAPACK authors (copyrights apply to src/*.f), the Laboratoire Statistique & Génome, and Carter T. Butts. In this project, these codes have been adapted specifically for a Multi-Layer SBM.
 

## Installation & Dependencies  

To run the code, make sure you have the following R packages installed:  

```r
install.packages(c("latex2exp", "ggplot2", "randnet", "igraph"))
```

## File structure  

- **code_MLSBM.R**: Contains all the functions used for the simulations.  
- **accuracy.R**: Compares the accuracy of different methods (KT-BHMC-NCV-PML).  
- **computation_time.R**: Measures the computational time required for the different methods (KT-BHMC-NCV-PML).  
- **rate_CV_edges.R**: Experiment attempting to observe the rate of convergence in the case where the maximum number of interactions is constant.  
- **rate_CV_nodes.R**: Experiment attempting to observe the rate of convergence in the case where the number of nodes is constant.  
- **sparse.R**: Experiments in a sparse regime context. 
  


## License

The entire source code of this product is available under the [CeCILL](http://www.cecill.info/) license, which is compatible with the [GNU GPL](https://www.gnu.org/licenses/gpl-3.0.html).
