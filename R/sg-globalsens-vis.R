#' Visualize global sensitivity analysis results
#'
#' @description
#' Generates visualization plots for results produced by
#' `sg_globalsens_sim()`. Supported visualizations:
#'
#' - PRCC: heatmap or barplot
#' - eFAST: barplot of first and total order indices
#'
#' @param x Result object returned by `sg_globalsens_sim()`
#' @param type Plot type for PRCC: `"heatmap"` or `"bar"`
#' @param params Optional vector of parameters to display
#' @param vars Optional vector of output variables
#' @param stats Optional vector of statistics to display
#'
#' @return ggplot object
#'
#' @examples
#' \donttest{
#' # Example model
#' library(rxode2)
#' library(tibble)
#' source(system.file("extdata", "RxODE_model", "example_rxode_model.R", package = "SimuRg")) # mod_ex
#'
#'
#' # Set up event table
#' et_base <- tribble(
#'   ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
#'   1,   0,     1,     1,    10,   2,     24,  112,   1
#' )
#'
#' # Define parameter bounds
#' inits <- rxInits(mod_ex)
#' par_bounds <- tibble::tibble(
#'   PAR = c("POPCL", "POPVC"),
#'   LB  = inits[c("POPCL", "POPVC")] * (1 - 0.9),
#'   UB  = inits[c("POPCL", "POPVC")] * (1 + 0.9)
#' )
#'
#' # Run eFAST sensitivity analysis
#' res_efast <- sg_globalsens_sim(
#'   method = c("eFAST"),
#'   model = mod_ex,
#'   params = c("POPCL"),
#'   par_bounds = par_bounds,
#'   n_sim = 100,
#'   stimes = seq(0, 168, 10),
#'   output = "Cc",
#'   cov = c("IGFR", "POPN"),
#'   et = et_base,
#'   stat_comp = c("mean")
#' )
#'
#'   efast_p <- sg_globalsens_vis(
#'     res_efast,
#'     params = c("POPCL"),
#'     vars = "Cc",
#'     stats = c("mean")
#'   )
#'
#'   efast_p
#' # Run PRCC sensitivity analysis
#' res_prcc <- sg_globalsens_sim(
#'   method = c("PRCC"),
#'   model = mod_ex,
#'   params = c("POPCL", "POPVC"),
#'   par_bounds = par_bounds,
#'   n_sim = 100,
#'   stimes = seq(0, 168, 10),
#'   output = "Cc",
#'   cov = c("IGFR", "POPN"),
#'   et = et_base,
#'   stat_comp = c("mean")
#' )
#'
#'   prcc_p <- sg_globalsens_vis(
#'     res_prcc,
#'     type = c("heatmap"),
#'     params = c("POPCL", "POPVC"),
#'     vars = "Cc",
#'     stats = c("mean")
#'   )
#'
#'   prcc_p
#' }
#'
#' @import ggplot2 dplyr tidyr
#' @export
sg_globalsens_vis <- function(x,
                              type = c("heatmap", "bar"),
                              params = NULL,
                              vars = NULL,
                              stats = NULL) {

  type <- match.arg(type)

  if (!"result" %in% names(x))
    stop("Input must be an object returned by sg_globalsens_sim()")

  df <- x$result

  # Detect method
  method <- if ("TYPE" %in% names(df) && "VALUE" %in% names(df)) {
    "eFAST"
  } else {
    "PRCC"
  }

  if (!is.null(params)) df <- df %>% dplyr::filter(PAR %in% params)
  if (!is.null(vars))   df <- df %>% dplyr::filter(VAR %in% vars)
  if (!is.null(stats))  df <- df %>% dplyr::filter(STAT %in% stats)

  # Create combined label
  if ("STAT" %in% names(df))
    df <- df %>%
    mutate(VAR_STAT = paste(VAR, STAT, sep = "_"))
  else
    df$VAR_STAT <- df$VAR


  # PRCC VISUALISATIONS
  if (method == "PRCC") {

    if (type == "heatmap") {

      p <- ggplot(df,
                  aes(x = PAR, y = VAR_STAT, fill = estimate)) +
        geom_tile(color = "grey25") +
        geom_text(aes(label = round(estimate, 2)), size = 4) +
        scale_fill_gradient2(
          low = MSDcol[3],
          mid = MSDcol[2],
          high = MSDcol[1],
          midpoint = 0,
          limits = c(-1, 1),
          name = "PRCC"
        ) +
        labs(
          x = "Parameters",
          y = "Output / statistic"
        ) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_minimal(base_size = 14)

    } else {

      p <- ggplot(df,
                  aes(x = PAR, y = estimate, fill = VAR_STAT)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_bar(stat = "identity",
                 position = position_dodge()) +
        geom_text(aes(label = round(estimate,2),
                      vjust = ifelse(estimate > 0, -0.6, 1.3)),
                  position = position_dodge(0.9),
                  size = 4) +
        scale_fill_manual(values = MSDcol) +
        ylim(-1.1, 1.1) +
        labs(
          x = "Parameters",
          y = "Partial rank correlation coefficient",
          fill = "Output"
        ) +
        theme_minimal(base_size = 14)

    }

  }


  # EFAST VISUALISATION
  if (method == "eFAST") {

    p <- ggplot(df,
                aes(x = PAR, y = VALUE, fill = TYPE)) +
      geom_bar(stat = "identity",
               position = "dodge",
               color = "black") +
      geom_hline(yintercept = 0) +
      facet_grid(STAT ~ VAR) +
      scale_fill_manual(values = c(MSDcol[2],MSDcol[3])) +
      scale_y_continuous(
        limits = c(-0.1, 1),
        breaks = seq(0, 1, 0.1),
        expand = c(0, 0)
      ) +
      labs(
        x = "Parameters",
        y = "eFAST sensitivity index",
        fill = "Index type"
      ) +
      theme_minimal(base_size = 14)

  }

  return(p)
}
