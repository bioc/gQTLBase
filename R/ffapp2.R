# > ffbase:::coerce_to_highest_vmode
coerce_to_highest_vmode = function (x, y, onlytest = TRUE) 
{
    test <- data.frame(x.vmode = vmode(x), y.vmode = vmode(y), 
        stringsAsFactors = FALSE)
    test$maxffmode <- apply(test[, , drop = FALSE], MARGIN = 1, 
        FUN = function(x) names(maxffmode(x)))
    needtocoerce <- list(coerce = test$x.vmode != test$maxffmode, 
        coerceto = test$maxffmode)
    if (onlytest) {
        return(needtocoerce)
    }
    if (sum(needtocoerce$coerce) > 0) {
        if (inherits(x, "ffdf")) {
            for (i in which(needtocoerce$coerce == TRUE)) {
                column <- names(x)[i]
                x[[column]] <- clone(x[[column]], vmode = needtocoerce$coerceto[i])
            }
            x <- x[names(x)]
        }
        else {
            x <- clone(x, vmode = needtocoerce$coerceto)
        }
    }
    x
}

ffapp2 = function (x, y, adjustvmode = TRUE, ...) 
{
    if (is.null(x)) {
        if (is.ff(y)) {
            return(clone.ff(y))  # trouble with clone.default dispatching
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
        x <- coerce_to_highest_vmode(x = x, y = y, onlytest = FALSE)
    }
    for (i in bit::chunk(x, from = 1, to = to, ...)) {
        if (is.atomic(y)) {
            i <- as.which(i)
        }
        x[(i + len)] <- y[i]
    }
    x
}
