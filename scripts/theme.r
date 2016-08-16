gg = modules::import_package('ggplot2', attach = TRUE)

theme = theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          rect = element_blank(),
          axis.line.x = element_line('grey50'),
          axis.line.y = element_line('grey50'),
          legend.title = element_blank())

theme_set(theme)

.color_values = c('gray50', 'tomato', 'cornflowerblue', 'tomato', 'cornflowerblue')
.color_limits = c('All', 'Germline', 'Infection response', 'WT', 'KO')

scale_color_discrete = function ()
    gg$scale_color_manual(values = setNames(.color_values, .color_limits))

scale_fill_discrete = function ()
    gg$scale_fill_manual(values = setNames(.color_values, .color_limits))

ggplot = function (data = NULL, mapping = aes(), ..., environment = parent.frame())
    gg$ggplot(data = data, mapping = mapping, ..., environment = environment) +
        scale_color_discrete() +
        scale_fill_discrete()
