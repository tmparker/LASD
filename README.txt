Note about computation

All the empirical and simulation computation was done using a combination of R
and C++ code, relying on the Rcpp package.  More information and the package are
available at

https://cran.r-project.org/web/packages/Rcpp/index.html  

Using Rcpp requires (in addition to R) a C++ compiler, for example using the
Rtools executable (MS Windows), Xcode command line tools (Mac OS) or gcc (Linux).

For those who prefer an all-R implementation, the commands contained in the file
lasd.cpp have been copied into R parlance and are contained in the file pdomR.R.
The pure R code is, however, much slower than when using lasd.cpp and Rcpp, and
is probably not suited for very intense computation.

