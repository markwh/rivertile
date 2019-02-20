# reticulate.R
# Manual fix(es?) of function(s?) in Reticulate package


#' Fix for reticulate::use_condaenv 
#' 
#' Previously ran into a warning about condition having length > 1
use_condaenv <- function(condaenv, conda = "auto", required = FALSE) {

# list all conda environments
conda_envs <- conda_list(conda)

# look for one with that name
conda_env_python <- subset(conda_envs, conda_envs$name == condaenv)$python
if (length(conda_env_python) == 0 && required)
  stop("Unable to locate conda environment '", condaenv, "'.")

if (!is.null(condaenv))
  use_python(unique(conda_env_python), required = required)

invisible(NULL)
}
