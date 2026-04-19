#' Install CLI launchers for project12pkg
#'
#' Installs command-line launchers for the package CLI app.
#'
#' @return Invisibly returns the result of \code{Rapp::install_pkg_cli_apps()}.
#' @export
install_project12pkg_cli <- function() {
  Rapp::install_pkg_cli_apps("project12pkg")
}
