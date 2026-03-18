#' Plot random effects or individual parameters vs covariates
#'
#' Generates ggplot2 visualizations for either random effects (RE) or individual parameters (IndPar)
#' versus continuous or categorical covariates from a Simurg object.
#'
#' @inheritParams sg_dummy
#' @param ptype Character. Type of plot: `"REvsCov"` for random effects or `"IndParvsCov"` for individual parameters.
#' @param cat_cov Optional tibble with categorical covariates. Must have columns `COV` and optionally `COVNAME` for labels.
#' @param cont_cov Optional tibble with continuous covariates. Must have columns `COV` and optionally `COVNAME` for labels.
#' @param color_palette Character vector. Colors used in plots.
#'
#' @return A list of ggplot objects. Returns versus continious covariates
#' `vs_contcov` and versus categorical covariates `vs_catcov` plots.
#'
#' @examples
#' library(tibble)
#' cont_cov <- tibble(
#'   COV = c("AGE", "WEIGHT"),
#'   COVNAME = c("Age, years", "Body weight, kg")
#' )
#' cat_cov <- tibble(
#'   COV = c("SEX", "VKORC1_gentyp"),
#'   COVNAME = c("Sex, M/F", "VKORC1 genotype")
#' )
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' p <- sg_gof_par_cov(
#'   fpath_i = fpath_i,
#'   ptype = "IndParvsCov",
#'   cont_cov = cont_cov,
#'   cat_cov = cat_cov
#' )
#' p$vs_contcov
#' p$vs_catcov
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom grid unit
#' @export
sg_gof_par_cov <- function(fpath_i,
                           ptype = "REvsCov",
                           cat_cov = NULL,
                           cont_cov = NULL,
                           color_palette = MSDcol) {
  smrg_obj <- read_smrg_obj(fpath_i)
  if (is.null(smrg_obj$PATAB) | is.null(smrg_obj$COTAB) | is.null(smrg_obj$CATAB) |
      is.null(smrg_obj$SUMTAB)) {
    stop("sg_fit object must contain PATAB, COTAB, CATAB and SUMTAB components")
  }
  patab <- smrg_obj$PATAB
  cotab <- smrg_obj$COTAB
  catab <- smrg_obj$CATAB
  sumtab <- smrg_obj$SUMTAB

  if (is.data.frame(patab) && nrow(patab) == 0) {
    stop("PATAB is empty (no rows)")
  } else if (is.data.frame(cotab) && nrow(cotab) == 0) {
    stop("COTAB is empty (no rows)")
  } else if (is.data.frame(catab) && nrow(catab) == 0) {
    stop("CATAB is empty (no rows)")
  } else if (is.data.frame(sumtab) && nrow(sumtab) == 0) {
    stop("SUMTAB is empty (no rows)")
  }

  if (is.list(patab) && length(patab) == 0) {
    stop("PATAB is empty (no elements)")
  } else if (is.list(cotab) && length(cotab) == 0) {
    stop("COTAB is empty (no elements)")
  } else if (is.list(catab) && length(catab) == 0) {
    stop("CATAB is empty (no elements)")
  } else if (is.list(sumtab) && length(sumtab) == 0) {
    stop("SUMTAB is empty (no elements)")
  }

  if (is.data.frame(patab)) { patab <- patab
  } else if (is.list(patab)) { patab <- as.data.frame(do.call(rbind, patab))
  } else stop("PATAB must be a data frame or a list of data frames")

  if (is.data.frame(cotab)) { cotab <- cotab
  } else if (is.list(cotab)) {cotab <- as.data.frame(do.call(rbind, cotab))
  } else stop("COTAB must be a data frame or a list of data frames")

  if (is.data.frame(catab)) { catab <- catab
  } else if (is.list(catab)) {catab <- as.data.frame(do.call(rbind, catab))
  } else stop("CATAB must be a data frame or a list of data frames")

  if (is.data.frame(sumtab)) { sumtab <- sumtab
  } else if (is.list(sumtab)) {sumtab <- as.data.frame(do.call(rbind, sumtab))
  } else stop("SUMTAB must be a data frame or a list of data frames")



  # IndPar columns
  typical_pars <- sumtab$PAR[sumtab$TYPE == "Typical values"]
  indpar_cols <- intersect(names(patab), gsub("_pop$", "", typical_pars, ignore.case = TRUE))
  # RE columns
  re_pars <- sumtab$PAR[grepl("Random effects", sumtab$TYPE, ignore.case = TRUE)]
  eta_cols <- intersect(names(patab), paste0("eta_", gsub("omega_|gamma_", "", re_pars)))

  # Covariate columns
  contcov_cols <- if (!is.null(cont_cov)) intersect(cont_cov$COV, setdiff(names(cotab), "ID")) else setdiff(names(cotab), "ID")
  catcov_cols <- if (!is.null(cat_cov)) intersect(cat_cov$COV, setdiff(names(catab), "ID")) else setdiff(names(catab), "ID")

  # Merge covariates
  cov_df <- catab %>%
    left_join(cotab, by = "ID")

  # Merge with parameters
  df <- patab %>%
    left_join(cov_df, by = "ID") %>%
    mutate(COHORTC = "All subjects")

  # Handle covariate labels
  if (!is.null(suppressWarnings(cont_cov$COVNAME))) {
    contcov_labels <- cont_cov
  } else {
    contcov_labels <- tibble(COV = contcov_cols, COVNAME = contcov_cols)
  }

  if (!is.null(suppressWarnings(cat_cov$COVNAME))) {
    catcov_labels <- cat_cov
  } else {
    catcov_labels <- tibble(COV = catcov_cols, COVNAME = catcov_cols)
  }

  # Helper functions: correlation + labeling
  add_corr_labels <- function(data) {
    data %>%
      group_by(PNAME, COV) %>%
      summarise(
        r = suppressWarnings(cor(COVVAL, VALUE, use = "complete.obs")),
        p = suppressWarnings(cor.test(COVVAL, VALUE)$p.value),
        .groups = "drop"
      ) %>%
      mutate(
        PEARS = sprintf("r = %.2f\np = %.3f", r, p),
        SIGNPEARS = if_else(p < 0.05, "firebrick", "black")
      )
  }

  add_anova_labels <- function(data) {
    data %>%
      group_by(PNAME, COV) %>%
      mutate(
        ANOVA = suppressWarnings(summary(aov(VALUE ~ as.factor(COVVAL))) [[1]][["Pr(>F)"]][1])
      ) %>%
      ungroup() %>%
      mutate(
        PVAL = ifelse(ANOVA < 0.001, "p < 0.001", paste0("p = ", round(ANOVA, 2))),
        LARGPVAL = ifelse(ANOVA < 0.05, "2", "1") # "2" significant, "1" not
      )
  }

  # select RE or IndPar columns
  if (ptype == "REvsCov") {
    # Random effects vs Continuous Covariates
    if (length(contcov_cols) > 0) {
    re_long <- df %>%
      select(ID, COHORTC, all_of(eta_cols), all_of(contcov_cols)) %>%
      pivot_longer(cols = all_of(eta_cols), names_to = "PNAME", values_to = "VALUE") %>%
      pivot_longer(cols = all_of(contcov_cols), names_to = "COV", values_to = "COVVAL")

    corr_stats <- add_corr_labels(re_long) %>%
      left_join(contcov_labels, by = "COV")

    plot_data <- re_long %>%
      left_join(corr_stats %>% select(PNAME, COV, COVNAME, PEARS, SIGNPEARS), by = c("PNAME", "COV"))

    re_vs_cont_p <- unique(select(corr_stats, PNAME, COV, COVNAME, PEARS, SIGNPEARS))
    # Plot REvsContCov
    p_cont <- ggplot(data = plot_data, aes(x = COVVAL, y = VALUE)) +
      geom_label(
        data = re_vs_cont_p,
        aes(x = -Inf, y = Inf, label = PEARS), colour = re_vs_cont_p$SIGNPEARS,
        hjust = -0.1, vjust = 1.1, show.legend = FALSE, size = 2.5
      ) +
      geom_point(aes(color = COHORTC), size = 1.5, alpha = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0, col = "grey25", lty = "dotted", linewidth = 0.5) +
      geom_smooth(method = "lm", formula = y ~ x, color = color_palette[3], se = FALSE) +
      facet_grid(PNAME ~ COVNAME, scales = "free", switch = "y") +
      scale_x_continuous(breaks = pretty_breaks(7), name = "Covariate value") +
      scale_y_continuous(breaks = pretty_breaks(7), name = "Random effect", position = "right") +
      scale_colour_manual(values = color_palette) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0.4, "lines"),
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(title = "Cohort")) +
      labs(title = "Random Effects vs Continious Covariates")
  } else p_cont <- NULL
    # Random effects vs Categorical Covariates
    if (length(catcov_cols) > 0) {
      re_cat <- df %>%
        select(ID, COHORTC, all_of(eta_cols), all_of(catcov_cols)) %>%
        mutate(across(all_of(catcov_cols), as.character))  %>%
        pivot_longer(cols = all_of(eta_cols), names_to = "PNAME", values_to = "VALUE") %>%
        pivot_longer(cols = all_of(catcov_cols), names_to = "COV", values_to = "COVVAL")

      re_cat <- add_anova_labels(re_cat) %>%
        left_join(catcov_labels, by = "COV")

      # Plot REvsCatCov
      p_cat <- ggplot(re_cat, aes(x = as.factor(COVVAL), y = VALUE)) +
        geom_boxplot(fill = color_palette[1], alpha = 0.5, outlier.colour = color_palette[3], outlier.shape = 3) +
        geom_hline(yintercept = 0, col = "grey25", lty = "dotted", linewidth = 0.5) +
        geom_label(data = unique(select(re_cat, PNAME, COV, COVNAME, PVAL, LARGPVAL)),
                   aes(x = -Inf, y = Inf, label = PVAL, col = LARGPVAL),
                   hjust = -0.1, vjust = 1.1, show.legend = FALSE, size = 2.5) +
        facet_grid(PNAME ~ COVNAME, scales = "free", switch = "y") +
        scale_y_continuous(breaks = pretty_breaks(7), name = "Random effect", position = "right") +
        scale_color_manual(values = c("1" = "black", "2" = "firebrick")) +
        labs(x = "Category", title = "Random Effects vs Categorical Covariates") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
              panel.grid.minor = element_blank())
    } else p_cat <- NULL

    return(list(vs_contcov = p_cont, vs_catcov = p_cat))

  } else if (ptype == "IndParvsCov") {

    indpar_cols <- intersect(indpar_cols, gsub("^eta_", "", eta_cols))

    # Individual parameter vs Continuous Covariates
    if ((length(contcov_cols) > 0) && (length(indpar_cols) > 0)) {
    ind_long <- df %>%
      select(ID, COHORTC, all_of(indpar_cols), all_of(contcov_cols)) %>%
      pivot_longer(cols = all_of(indpar_cols), names_to = "PNAME", values_to = "VALUE") %>%
      pivot_longer(cols = all_of(contcov_cols), names_to = "COV", values_to = "COVVAL")

    corr_stats <- add_corr_labels(ind_long) %>%
      left_join(contcov_labels, by = "COV")

    plot_data <- ind_long %>%
      left_join(corr_stats %>% select(PNAME, COV, COVNAME, PEARS, SIGNPEARS), by = c("PNAME", "COV"))

    ip_vs_cont_p <- unique(select(corr_stats, PNAME, COV, COVNAME, PEARS, SIGNPEARS))
    # Plot IndParvsContCov
    p_cont <- ggplot(plot_data, aes(x = COVVAL, y = VALUE)) +
      geom_label(
        data = ip_vs_cont_p,
        aes(x = -Inf, y = Inf, label = PEARS), colour = ip_vs_cont_p$SIGNPEARS,
        hjust = -0.1, vjust = 1.1, show.legend = FALSE, size = 2.5
      ) +
      geom_point(aes(color = COHORTC), size = 1.5, alpha = 0.7, show.legend = FALSE) +
      geom_smooth(method = "lm", formula = y ~ x, color = color_palette[3], se = FALSE) +
      facet_grid(PNAME ~ COVNAME, scales = "free", switch = "y") +
      scale_x_continuous(breaks = pretty_breaks(7), name = "Covariate value") +
      scale_y_continuous(breaks = pretty_breaks(7), name = "Individual parameter", position = "right") +
      scale_colour_manual(values = color_palette) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "top",
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0.4, "lines")
      ) +
      guides(color = guide_legend(title = "Cohort")) +
      labs(title = "Individual Parameters vs Continious Covariates")
  } else p_cont <- NULL
  # IndPar vs Categorical Covariates
  if ((length(catcov_cols) > 0) && (length(indpar_cols) > 0)) {
    ind_cat <- df %>%
      select(ID, COHORTC, all_of(indpar_cols), all_of(catcov_cols)) %>%
      mutate(across(all_of(catcov_cols), as.character))  %>%
      pivot_longer(cols = all_of(indpar_cols), names_to = "PNAME", values_to = "VALUE") %>%
      pivot_longer(cols = all_of(catcov_cols), names_to = "COV", values_to = "COVVAL")

    ind_cat <- add_anova_labels(ind_cat) %>%
      left_join(catcov_labels, by = "COV")

    # Plot IndParvsCatCov
    p_cat <- ggplot(ind_cat, aes(x = as.factor(COVVAL), y = VALUE)) +
      geom_boxplot(fill = color_palette[1], alpha = 0.5, outlier.colour = color_palette[3], outlier.shape = 3) +
      geom_label(data = unique(select(ind_cat, PNAME, COV, COVNAME, PVAL, LARGPVAL)),
                 aes(x = -Inf, y = Inf, label = PVAL, col = LARGPVAL),
                 hjust = -0.1, vjust = 1.1, show.legend = FALSE, size = 2.5) +
      facet_grid(PNAME ~ COVNAME, scales = "free", switch = "y") +
      scale_y_continuous(breaks = pretty_breaks(7), name = "Individual parameter", position = "right") +
      scale_color_manual(values = c("1" = "black", "2" = "firebrick")) +
      labs(x = "Category", title = "Individual Parameters vs Categorical Covariates") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
  } else p_cat <- NULL

  return(list(vs_contcov = p_cont, vs_catcov = p_cat))

  }
}
