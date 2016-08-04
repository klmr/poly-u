taginfo_data = taginfo %>%
    group_by(Type, Treatment, Sample, Mod) %>%
    summarize(ModCount = n()) %>%
    mutate(Frequency = ModCount / sum(ModCount)) %>%
    ungroup()

taginfo_test_treatment = taginfo_data %>%
    group_by(Mod) %>%
    do(p = t.test(.[.$Treatment == 'Infected', ]$Frequency,
                  .[.$Treatment == 'Uninfected', ]$Frequency,
                  var.equal = TRUE)$p.value) %>%
    mutate(p = unlist(p)) %>%
    ungroup()

taginfo_test_type = taginfo_data %>%
    group_by(Mod) %>%
    do(p = t.test(.[.$Type == 'WT', ]$Frequency,
                  .[.$Type == 'KO', ]$Frequency,
                  var.equal = TRUE)$p.value) %>%
    mutate(p = unlist(p)) %>%
    ungroup()

none = '#00000000'

ggplot(taginfo_data) +
    aes(paste(Type, Treatment), Frequency, color = as.factor(Sample)) +
    geom_point(shape = 21, fill = none, size = 2) +
    facet_wrap(~ Mod, scales = 'free_y') +
    labs(x = 'Condition', color = 'Sample') +
    theme_bw()

ggsave('data/plots/mod-frequency-change.pdf', width = 10, height = 7)
