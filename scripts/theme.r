gg = modules::import_package('ggplot2', attach = TRUE)

theme = theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          rect = element_blank(),
          axis.line.x = element_line('grey50'),
          axis.line.y = element_line('grey50'),
          legend.title = element_blank())

theme_set(theme)

scale_color_discrete = function ()
    gg$scale_color_grey()
