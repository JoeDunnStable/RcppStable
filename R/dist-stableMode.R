## Part of R package 'stabledist' (part of the Rmetrics project).

## The stabledist R package is free software; you can redistribute it and/or
## modify it under the terms of the GNU Library General Public
## License as published by the Free Software Foundation; either
## version 2 of the License, or (at your option) any later version.
##
## This R package is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Library General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/


################################################################################
# FUNCTIONS:		DESCRIPTION:
#  stableMode		 Computes the mode of the stable DF
################################################################################

##' Computes the mode of the alpha stable distribution
##' @title Mode of the stable distribution
##' @param alpha
##' @param beta
##' @param tol numerical tolerance used to find mode
##' @return a number, the stable mode
##' @author Diethelm Wuertz and Martin Maechler
stableMode <- function(alpha, beta, tol = 64*.Machine$double.eps)
{
  verbose=getOption("stableMode.debug", default=F)
  stopifnot(0 <= alpha, alpha <= 2, length(alpha) == 1,
	          -1 <= beta, beta <= 1, length(beta) == 1)
  sdstableMode(alpha,beta,tol,tol,1000,verbose)
}
