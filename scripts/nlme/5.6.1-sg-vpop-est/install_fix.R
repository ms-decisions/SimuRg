# FIX FOR SYNTHPOP / RSOLNP INSTALLATION
# ----------------------------------------

# 1. Define a binary-friendly CRAN mirror
repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos)

# 2. Clean renv state for these specific packages to avoid cache errors
if (requireNamespace("renv", quietly = TRUE)) {
  message("Cleaning renv cache for Rsolnp and synthpop...")
  try(renv::purge("Rsolnp"), silent = TRUE)
  try(renv::purge("synthpop"), silent = TRUE)
}

# 3. Install Rsolnp specifically as a BINARY
# This bypasses the "compilation failed" error on Windows
message("\nAttempting to install Rsolnp (binary)...")
tryCatch({
  # Try via renv first with explicit binary type
  renv::install("Rsolnp", type = "binary", prompt = FALSE)
  message("SUCCESS: Rsolnp installed via renv.")
}, error = function(e) {
  message("renv install failed, trying standard install.packages...")
  # Fallback to standard install
  install.packages("Rsolnp", type = "binary")
  # Record it in renv
  renv::snapshot(packages = "Rsolnp", prompt = FALSE)
  message("SUCCESS: Rsolnp installed via standard method.")
})

# 4. Install synthpop (now that dependency is fixed)
message("\nAttempting to install synthpop...")
tryCatch({
  renv::install("synthpop", prompt = FALSE)
  message("SUCCESS: synthpop installed!")
}, error = function(e) {
  message("FAILED to install synthpop. Please check the error message above.")
})

