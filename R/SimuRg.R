#' @description
#' SimuRg is an R toolkit for PK/PD modeling and simulation workflows, with a focus on importing, organizing, and exploring nonlinear mixed-effects model results.
#' @keywords internal

"_PACKAGE"

rlang::on_load(local_use_cli())

# Singletons
the <- rlang::new_environment()
