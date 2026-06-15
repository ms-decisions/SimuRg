#' @description
#' \figure{logo.png}{options: style="float: right" alt="logo" width="120"}
#' SimuRg is an R toolkit for PK/PD modeling and simulation workflows, with a focus on importing, organizing, and exploring nonlinear mixed-effects model results.
#' @docType package
#' @keywords internal

"_PACKAGE"

rlang::on_load(local_use_cli())

# Singletons
the <- rlang::new_environment()
