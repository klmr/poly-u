modules::import('klmr/functional/lambda')

# FIXME: This MUST exist in R. But I can only find `runmed`.
runmean = function (data, k) {
    shift = indices ~ range -> {
        if (head(indices, 1) < range[1])
            indices = indices + range[1] - head(indices, 1)
        if (tail(indices, 1) > range[2])
            indices = indices + range[2] - tail(indices, 1)
        indices
    }

    window = data ~ i ~ k -> {
        w = k %/% 2
        shift(seq((i - w), (i + w)), c(1, length(data)))
    }

    if (k %% 2 == 0) {
        k = k + 1
        warning(sprintf('%s must be odd. Changing to %s', sQuote('k'), k))
    }
    unlist(lapply(seq_along(data),
                  i -> mean(data[window(data, i, k)])))
}
