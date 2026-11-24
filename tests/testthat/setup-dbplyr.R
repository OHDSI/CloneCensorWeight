# Test-only workaround for a slow spot in dbplyr <= 2.5.2: every variable of a summarise() on a lazy
# table (so every count() and tally()) formats an error message with cli, about 5 ms each, in case a
# later variable refers to it; the message is used only then. That took about 15% of the test time.
# The development version of dbplyr formats the message only when it is used, and the patch below does
# the same (the message is formatted when the error is raised). It changes no results. It is applied
# only while the installed dbplyr formats the message up front, so it switches itself off once dbplyr
# is updated.
local({
  namespace <- asNamespace("dbplyr")
  if (!exists("summarise_bind_error1", envir = namespace, inherits = FALSE)) {
    return()
  }
  bindError <- get("summarise_bind_error1", envir = namespace)
  code <- body(bindError)
  # The body is { msg <- cli::format_message(c(...)); env_bind_lazy(error_env, !!dot_name := { abort(msg, call = NULL) }) }:
  if (length(code) != 3 || !identical(code[[2]][[1]], as.name("<-")) || !identical(deparse(code[[2]][[3]][[1]]), "cli::format_message")) {
    return()
  }
  formatMessage <- code[[2]][[3]]
  # msg <- c(...), and the message formatted inside the lazy binding:
  code[[2]][[3]] <- formatMessage[[2]]
  formatMessage[[2]] <- as.name("msg")
  code[[3]] <- do.call(substitute, list(code[[3]], list(msg = formatMessage)))
  body(bindError) <- code
  utils::assignInNamespace("summarise_bind_error1", bindError, ns = "dbplyr")
})
