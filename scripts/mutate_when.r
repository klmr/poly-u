lazyeval = import_package('lazyeval')

mutate_when = function (.data, .filter, ...) {
    mutate_when_(.data, lazyeval$lazy(.filter), .dots = lazyeval$lazy_dots(...))
}

mutate_when_ = function (.data, .filter, ..., .dots) {
    dots = lazyeval$all_dots(.dots, ..., all_named = TRUE)
    selected_rows = lazyeval$lazy_eval(.filter, .data)

    .data[selected_rows, names(dots)] =
        lapply(dots, lazyeval$lazy_eval, data = .data[selected_rows, ])
    .data
}
