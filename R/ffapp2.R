ffapp2 = function (x, y, adjustvmode = TRUE, ...) 
{
    if (is.null(x)) {
        if (is.ff(y)) {
            return(ff:::clone.ff(y))  # trouble with clone.default dispatching
        }
        else {
            return(if (length(y)) as.ff(y))
        }
    }
    len <- length(x)
    to <- length(y)
    if (!to) 
        return(x)
    length(x) <- len + to
    if (is.factor.ff(x)) {
        levels(x) <- appendLevels(levels(x), levels(y))
    }
    if (adjustvmode == TRUE) {
        x <- ffbase:::coerce_to_highest_vmode(x = x, y = y, onlytest = FALSE)
    }
    for (i in bit::chunk(x, from = 1, to = to, ...)) {
        if (is.atomic(y)) {
            i <- as.which(i)
        }
        x[(i + len)] <- y[i]
    }
    x
}
