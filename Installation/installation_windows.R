rm(list = ls())
cat("\014")

# INSTRUCTIONS - GO THROUGH THIS SCRIPT AFTER OPENING R IN A FRESH R SESSION 
# AND NOT WITHIN ANY R PROJECT (.Rproj)

# Make sure necessary packages are installed
if(!require("tidyverse", character.only = TRUE)){
  install.packages("tidyverse", dependencies = TRUE)
  library("tidyverse",  character.only = TRUE)
}

if(any(grepl("tidyverse", search()))){
  cat("'tidyverse' is attached. Continue." )
}else{
  warning(strwrap("'tidyverse' is not attached. Go back and make sure it is 
                  installed and attached."))
}

## TODO: We can add more packages if necessary
packages <- c("bayesplot", "brms", "collapsibleTree", "patchwork", "posterior",
              "rstanarm", "tidybayes", "ggforce", "gganimate", "gifski", 
              "ggpubr", "latex2exp")

walk(packages, .f = function(.x){
  if(!require(.x, character.only = TRUE)){
    install.packages(.x, dependencies = TRUE)
    library(.x, character.only = TRUE)
  }
})

# Install cmdstanr
## we recommend running this in a fresh R session or restarting your current 
## session
if(!require("cmdstanr", character.only = TRUE)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", 
                                         getOption("repos")))
  library("cmdstanr",  character.only = TRUE)
}

walk(c(packages, "cmdstanr"), 
     .f = function(.x) {
       if(any(str_detect(search(), .x))){
         cat(str_c("'", .x, "'", " is attached. Continue.\n"))
       }else{
         warning(str_wrap(str_c(.x, " is not attached. Go back and make sure it 
         is installed and attached.")))
       }
     })

# Make sure C++ toolchain is installed and working properly. For Windows this
# means RTools. It should be fine in Mac and Linux
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

# Install CmdStan  (optional)
## Change parallel::detectCores() to whatever number you want to use. 
## More cores = faster
install_cmdstan(cores = parallel::detectCores()) 

# Check if CmdStan is working properly
model_cmdstan <- cmdstan_model(file.path(cmdstan_path(), "examples", 
                                         "bernoulli", "bernoulli.stan"))
model_cmdstan$print()
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit_cmdstan <- model_cmdstan$sample(data = data_list, 
                                    seed = 345, 
                                    chains = 4, 
                                    parallel_chains = 4,
                                    refresh = 500)

fit_cmdstan$summary()

# Install Torsten (mandatory)
system("git clone https://github.com/metrumresearchgroup/Torsten.git")
shell("cd Torsten/cmdstan && mingw32-make build")

set_cmdstan_path("~/Torsten/cmdstan/")
cmdstan_version()


# Check to see if Torsten is working properly (mandatory)
## First with a simple model in pure Stan code
model_torsten <- cmdstan_model(file.path(cmdstan_path(), "examples", 
                                         "bernoulli", "bernoulli.stan"))
model_torsten$print()
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit_torsten <- model_torsten$sample(data = data_list, 
                                    seed = 345, 
                                    chains = 4, 
                                    parallel_chains = 4,
                                    refresh = 500)

fit_torsten$summary()

## Now with a simple two-compartment model that actually uses Torsten functions
file_path_base <- file.path("~", "Torsten", "example-models", "pk2cpt")
model_pk2cpt <- cmdstan_model(file.path(file_path_base, "pk2cpt.stan"))

fit_pk2cpt <- model_pk2cpt$sample(data = file.path(file_path_base, 
                                                   "pk2cpt.data.R"),
                                  seed = 345, 
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  chains = 4,
                                  parallel_chains = 4,
                                  refresh = 500)

fit_pk2cpt$summary()

