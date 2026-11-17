# Test-only workaround for a bug in Andromeda 1.2.0 and 1.2.1: flushAndromeda() starts with
# `if (isAndromeda(andromeda()))`, which calls the function Andromeda::andromeda() instead of using
# its argument, so every flush creates a new temporary duckdb database that lives until garbage
# collection. Andromeda flushes on every `x$name <- <tbl>`, so the tests spent about half their time
# opening databases. The patch below tests the argument, as intended. Only the condition changes: the
# CHECKPOINT and the cache eviction already ran on the argument, so for an Andromeda argument (which
# every caller in this package's dependencies passes) the patched and unpatched versions flush the same
# database and give the same results. (The fix is on the local Andromeda
# branch fix-flush-andromeda, as version 1.2.2.) It is applied only while the installed Andromeda still
# has the bug, so it switches itself off once Andromeda is fixed.
local({
  flush <- Andromeda::flushAndromeda
  if (any(grepl("isAndromeda(andromeda())", deparse(body(flush)), fixed = TRUE))) {
    # body(flush)[[2]] is the `if`, [[2]] its condition and [[2]] the call andromeda():
    body(flush)[[2]][[2]][[2]] <- quote(andromeda)
    utils::assignInNamespace("flushAndromeda", flush, ns = "Andromeda")
  }
})
