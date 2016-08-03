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
