#' Cumulative sum over nearest neighbors
#'
#' \code{nn.cumsum} computes the cumulative sum of \code{y}
#' for the sequences of indices in each element of the list
#' contained in \code{nn}.
#'
#' @param nn A list of nearest neighbors in the format
#'   produced by \code{\link{nnpop}}.
#' @param y A numeric vector of values to be summed over.
#' @param simplify A logical value indicating whether the
#'   results should be simplified to a numeric vector.  The
#'   default is \code{TRUE}.
#'
#' @return A vector or list, depending on the value of
#'   \code{simplify}.
#' @export
#'
#' @examples
#' # show nn.cumsum example for a circular scan setting
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' d = sp::spDists(coords, longlat = TRUE)
#' # compute circular nearest neigbhors
#' nn = nnpop(d, pop = nydf$pop, ubpop = 0.1)
#' # compute cumulative sums over all nn
#' cnn = nn.cumsum(nn, cases)
#' # compute cumulative sums over just the first set of nn
#' cnn1 = cumsum(cases[nn[[1]]])
#' # check equality
#' all.equal(cnn1, cnn[seq_along(cnn1)])
nn.cumsum = function(nn, y, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  if (simplify) {
    unlist(lapply(nn, function(x) cumsum(y[x])), use.names = FALSE)
  } else {
    lapply(nn, function(x) cumsum(y[x]))
  }
}

#' @export
#' @rdname nn.cumsum
nn.cumlen = function(nn, y, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  if (simplify) {
    unlist(lapply(nn, function(x) seq_len(length(x))))
  } else {
    lapply(nn, function(x) seq_len(length(x)))
  }
}

#' @export
#' @rdname nn.cumsum
nn.cumvar = function(nn, y, out = FALSE, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  dseq = nn.cumlen(nn, y, simplify = FALSE)
  dseqv = unlist(dseq)
    if (out) {
      cumvars = unlist(mapply(function(nni, indx) nn.zonepvar(y, nni, indx, TRUE), nn, dseq, SIMPLIFY = simplify, USE.NAMES = FALSE))
      # mapply(function(nni, indx) cumstats::cumvar(y[-nni[seq(1, indx, 1)]]) * (indx-1) / indx, nn, dseq, SIMPLIFY = TRUE)
    } else {
      cumvars = unlist(sapply(nn, function(nni) as.vector(roll_var(y[nni]))))
      #unlist(mapply(function(nni, indx) nn.zonepvar(y, nni, indx, FALSE), nn, dseq, SIMPLIFY = simplify, USE.NAMES = FALSE))
      # mapply(function(nni, indx) cumstats::cumvar(y[nni[seq(1, indx, 1)]]) * (indx-1) / indx, nn, dseq, SIMPLIFY = TRUE)
    }
    cumvars * (dseqv - 1) / dseqv
}



#' @export
#' @rdname nn.cumsum
nn.cumsumsq = function(nn, y, out = FALSE, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  dseq = nn.cumlen(nn, y, simplify = FALSE)
    if (out) {
      res = unlist(mapply(function(nni, indx) sum(y^2) - as.vector(roll_psum(y[nni], 2)), nn, dseq, SIMPLIFY = simplify, USE.NAMES = FALSE))
      # mapply(function(nni, indx) cumsum(y[-nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = TRUE)
    } else {
      res = unlist(mapply(function(nni, indx) as.vector(roll_psum(y[nni], 2)), nn, dseq, SIMPLIFY = simplify, USE.NAMES = FALSE))
      # mapply(function(nni, indx) cumsum(y[nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = TRUE)
    }
    res
  
}



#' @export
#' @rdname nn.cumsum
nn.cumsumsqold = function(nn, y, out = FALSE, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  dseq = nn.cumlen(nn, y, simplify = FALSE)
  if (simplify) {
    if (out) {
      res = mapply(function(nni, indx) nn.zonepsum(y, nni, indx, 2, TRUE), nn, dseq, SIMPLIFY = TRUE)
      # mapply(function(nni, indx) cumsum(y[-nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = TRUE)
    } else {
      res = mapply(function(nni, indx) nn.zonepsum(y, nni, indx, 2, FALSE), nn, dseq, SIMPLIFY = TRUE)
      # mapply(function(nni, indx) cumsum(y[nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = TRUE)
    }
    unlist(res)
  } else {
    if (out) {
      res = mapply(function(nni, indx) nn.zonepsum(y, nni, indx, 2, TRUE), nn, dseq, SIMPLIFY = FALSE)
      # mapply(function(nni, indx) cumsum(y[-nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = FALSE)
    } else {
      res = mapply(function(nni, indx) nn.zonepsum(y, nni, indx, 2, FALSE), nn, dseq, SIMPLIFY = FALSE)
      # mapply(function(nni, indx) cumsum(y[nni[seq(1, indx, 1)]]^2), nn, dseq, SIMPLIFY = FALSE)
    }
    res
  }
}

#'
#' @param y A numeric vector of values to be summed over.
#' @param nni A numeric vector of nearest neighbors in the format
#'   produced by \code{\link{nnpop}}.
#' @param indices A numeric vector of sequential values that represent the indices
#' @export
#' @rdname nn.cumsum
nn.zonepsum = function(y, nni, indices, power = 1, out = FALSE) {
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(nni)) stop("nni must be a numeric vector")
  if (!is.numeric(indices)) stop("indices must be a numeric vector")
  if (!is.logical(out)) stop("out must be a logical value")
  if (out) {
    sapply(indices, function(index) sum(y[-nni[seq(1, index, 1)]]^power), USE.NAMES = FALSE)
  } else {
    sapply(indices, function(index) sum(y[nni[seq(1, index, 1)]]^power), USE.NAMES = FALSE)
  }
}

#' @export
#' @rdname nn.cumsum
nn.zonepvar = function(y, nni, indices, out = FALSE) {
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(nni)) stop("nni must be a numeric vector")
  if (!is.numeric(indices)) stop("indices must be a numeric vector")
  if (!is.logical(out)) stop("out must be a logical value")
  if (out) {
    sapply(indices, function(index) var(y[-nni[seq(1, index, 1)]]), USE.NAMES = FALSE)
  } else {
    sapply(indices, function(index) var(y[nni[seq(1, index, 1)]]), USE.NAMES = FALSE)
  }
}

