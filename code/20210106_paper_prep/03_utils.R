## 'utils.R' contains useful functions I copied from somewhere or written by my collaborators (see each function for sources).
##
## Copyright (C) 2014 Ester Pantaleo (for the function get.counts.h5)
##
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.




library(rhdf5)

#**************************************************************************
#  function to Get counts from hdf5
#  written by Ester Pantaleo
#**************************************************************************
get.counts.h5 <- function(list_path, chr, locus.start, locus.end, list_name = NULL, print_message=FALSE){
  M <-  NULL
  for (h5file in list_path){
    if(print_message){
      print(paste0("Loading ", h5file))
      print(paste0("h5read(", h5file, ",", chr, ", index=list(", locus.start, ":", locus.end, ")"))
    }
    M        <- rbind(M, h5read(h5file, chr, index=list(locus.start:locus.end)))
  }
  row.names(M) = list_name
  return(M)
}

