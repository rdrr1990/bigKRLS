.pkgenv <- new.env(parent = emptyenv())
# this code is adapted from Dirk Eddelbuettel's https://github.com/eddelbuettel/anytime/blob/ef8b1e52b80a99e96f46232dfe29180686327887/R/init.R#L49-L52
# in light of this discussion: http://stackoverflow.com/questions/43247649/rcpparmadillo-bigmemory-crashes-windows-rstudio-but-no-other-gui-os-type
.onLoad <- function(libname, pkgname){
  
  packageStartupMessage("\n\nCheck out vignette(\"bigKRLS_basics\") for a brief explanation of the statistics, references, and syntax.")
  
  boostable <- if(eval(parse(text=".Platform$OS.type == \"unix\" | .Platform$GUI != \"RStudio\""))){
    if (Sys.getenv("RSTUDIO", unset="0") == "1" &&
        exists("RStudio.Version") &&
        ## the following is evil but keeps R CMD check off our back
        eval(parse(text=paste("RStudio.Version()$version", ">=", "\"1.0.136\"")))) TRUE else FALSE
  }else{
    if (Sys.getenv("RSTUDIO", unset="0") == "1" &&
        exists("RStudio.Version") &&
        ## the following is evil but keeps R CMD check off our back
        eval(parse(text=paste("RStudio.Version()$version", ">=", "\"1.1.129\"")))) TRUE else FALSE
    
  }
    
  .pkgenv[["boostable"]] <- boostable
  
}