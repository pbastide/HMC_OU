# Copyright (C) {2020} {PB, LSTH, GB, PL, MAS}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##
#' @title Add iid error to each tip
#'
#' @param dat a data frame with ntaxa row and nTrait columns
#' @param sigma_e a nTrait x nTrait variance matrix
#' 
#' @return Data frame with added errors
##
add_errors <- function(dat, sigma_e) {
  dat <- dat + MASS::mvrnorm(nrow(dat), mu = rep(0, ncol(dat)), Sigma = sigma_e)
}