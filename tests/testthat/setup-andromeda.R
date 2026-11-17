# Test-only workaround for a bug in Andromeda <= 1.2.1: flushAndromeda() starts with
# `if (isAndromeda(andromeda()))`, which calls the function Andromeda::andromeda() instead of using
# its argument, so every flush creates (and leaks) a new temporary duckdb database. Andromeda flushes
# on every `x$name <- <tbl>`, so the tests spent about half their time opening databases. The patch
# below uses the argument, as intended; it changes no results. It is applied only while the
# installed Andromeda still has the bug, so it switches itself off once Andromeda is fixed.
local({
  flush <- Andromeda::flushAndromeda
  if (any(grepl("isAndromeda(andromeda())", deparse(body(flush)), fixed = TRUE))) {
    # body(flush)[[2]] is the `if`, [[2]] its condition and [[2]] the call andromeda():
    body(flush)[[2]][[2]][[2]] <- quote(andromeda)
    utils::assignInNamespace("flushAndromeda", flush, ns = "Andromeda")
  }
})
