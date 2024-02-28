# This file configures the virtualenv and Python paths differently depending on
# the environment the app is running in (local vs remote server).

VIRTUALENV_NAME = 'r-reticulate'
# Edit this name if desired when starting a new 


# ------------------------- Settings (Do not edit) -------------------------- #

if (Sys.info()[['user']] == 'shiny'){
  
  # Running on shinyapps.io
  # The PYTHON_PATH environment variable gets fed into the python argument of
  # the virtualenv_create function.  The python argument seems to be a little 
  # strange.  It seems to want an executable (or a symlink to an executalble)
  # such as python3.  But it will take a version.  On my linux box, an 
  # executable worked.  After a lot of pain and agony, it seems like to get it 
  # to work from my Mac, that I need to use a version.  This might be a little
  # problematic when shinyapps updates their version of python, but whatever.
  # Another solution is to point it to the actual executable - which on 
  # shinyapps is /usr/bin/python3.10 but that seems to have the same problem
#  #  Sys.setenv(PYTHON_PATH = 'python3') - semes ideal.  Worked from linux
#  #  Sys.setenv(PYTHON_PATH = '/usr/bin/python3.10') - another solution-for now
  Sys.setenv(PYTHON_PATH = '3.10')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # Installs into default shiny virtualenvs dir
  # It seems like reticulate and virtual environments really like having the 
  # RETICULATE_PYTHON environment variable set, but despite setting it, shinyapps
  # doesn't seem to use it. Just complains that it isn't set, despite the fact 
  # that I set it.  So I unset it.
#  # Sys.setenv(RETICULATE_PYTHON = paste0('/home/shiny/.virtualenvs/', VIRTUALENV_NAME, '/bin/python'))
  
} else if (Sys.info()[['user']] == 'rstudio-connect'){
  
  # Running on remote server
  Sys.setenv(PYTHON_PATH = '/opt/python/3.7.6/bin/python')
  Sys.setenv(VIRTUALENV_NAME = paste0(VIRTUALENV_NAME, '/')) # include '/' => installs into rstudio-connect/apps/
  Sys.setenv(RETICULATE_PYTHON = paste0(VIRTUALENV_NAME, '/bin/python'))
  
} else {
  
  # Running locally
  options(shiny.port = 7450)
  Sys.setenv(PYTHON_PATH = 'python3')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # exclude '/' => installs into ~/.virtualenvs/
  # RETICULATE_PYTHON is not required locally, RStudio infers it based on the ~/.virtualenvs path
}
