PathwayAnalyser (PA) is a tool for systems biologists for analysis of metabolic pathways, particularly by Flux Balance Analysis and ODE simulation. PA interfaces with the [GNU Linear Programming Toolkit (GLPK)](http://www.gnu.org/software/glpk/glpk.html) for linear programming, as well as with [Taylor](http://www.maia.ub.es/~angel/taylor/) for performing high precision simulations (e.g. using [GMP](http://www.swox.com/gmp/)). [OOQP](http://www.cs.wisc.edu/~swright/ooqp/) is used for quadratic programming/MoMA. It can give a comprehensive report on gene deletions from the [Systems Biology Markup Language (SBML)](http://sbml.org/) Model input. It is currently a command-line tool but a GUI may be added soon.

It is written in C++ and uses [libSBML](http://www.sbml.org/software/libsbml/) and GLPK. PA is meant for linux (excl. [cygwin](http://www.cygwin.com/) at the moment) but it is planned that it will be usable in all platforms in the future.

# Documents #
  * PathwayAnalyser poster: http://precedings.nature.com/documents/1868/version/1

# External Links #
  * PathwayAnalyser, http://sourceforge.net/projects/pathwayanalyser
  * GNU Linear Programming Toolkit (GLPK), http://www.gnu.org/software/glpk/glpk.html
  * Taylor, http://www.maia.ub.es/~angel/taylor/
  * GMP, http://www.swox.com/gmp/
  * OOQP, http://www.cs.wisc.edu/~swright/ooqp/
  * Systems Biology Markup Language (SBML), http://sbml.org/
  * libSBML, http://www.sbml.org/software/libsbml/