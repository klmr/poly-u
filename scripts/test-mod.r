taginfo_stat = taginfo %>%
    group_by(Type, Treatment, Sample, Mod) %>%
    summarize(ModCount = n()) %>%
    mutate(Frequencies = ModCount / sum(ModCount)) %>%
    group_by(Type, Treatment, Mod) %>%
    summarize(Mean = mean(Frequencies), SD = sd(Frequencies), Frequencies = list(Frequencies)) %>%
    ungroup()

taginfo_test = taginfo_stat %>%
    select(Type, Treatment, Mod, Frequencies) %>%
    tidyr$spread(Treatment, Frequencies)

taginfo_test = taginfo_test %>%
    rowwise() %>%
    do(p = t.test(unlist(.$Infected), unlist(.$Uninfected), var.equal = TRUE)$p.value) %>%
    mutate(p = unlist(p)) %>%
    ungroup() %>%
    bind_cols(taginfo_test, .)
