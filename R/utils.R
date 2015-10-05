## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_STABLEDIST_CHECK_EXTRA")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
