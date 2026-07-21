# Initialize renv in the current project
renv::init()


#take sapshot everytime that a package is installed or updated
renv::snapshot()


#restore renv in a new computer for reproducibility
renv::restore()
