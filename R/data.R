## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for common internal data structures generation
## Keywords: SimuRg, internal data

MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d")
funSum_sim <- list(mean   = ~mean(., na.rm = T),
                   median = ~median(., na.rm = T),
                   min    = ~min(., na.rm = T),
                   max    = ~max(., na.rm = T),
                   sd     = ~sd(., na.rm = T),
                   P025   = ~quantile(., 0.025, na.rm = T),
                   P05    = ~quantile(., 0.05, na.rm = T),
                   P10    = ~quantile(., 0.10, na.rm = T),
                   P15   = ~quantile(., 0.15, na.rm = T),
                   P25    = ~quantile(., 0.25, na.rm = T),
                   P75    = ~quantile(., 0.75, na.rm = T),
                   P85   = ~quantile(., 0.85, na.rm = T),
                   P90    = ~quantile(., 0.90, na.rm = T),
                   P95    = ~quantile(., 0.95, na.rm = T),
                   P975   = ~quantile(., 0.975, na.rm = T),
                   geom_mean = ~exp(mean(log(.), na.rm = T)),
                   CV     = ~sd(., na.rm = T)/mean(., na.rm = T)*100)
