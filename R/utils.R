RStudio_outdated <- function(){
  .pkgenv$RStudio_outofdate
}

check_platform <- function(){
  
  problem <- RStudio_outdated()
  if(is.null(problem)){
    problem <- if(.Platform$GUI == "RStudio"){
      threshold <- if(.Platform$OS.type == "unix") "1.0.136" else "1.1.129"
      !(Sys.getenv("RSTUDIO", unset="0") == "1" &&
          exists("RStudio.Version") &&
          eval(parse(text=paste0("RStudio.Version()$version ",  " >= ", "\"", threshold, "\""))))
    } else FALSE
  }
  
  if(problem){
    if(.Platform$OS.type == "unix"){
      stop("bigKRLS requires RStudio 1.0.136 or higher.\n       To use bigKRLS, switch to RGui or check the following webpages:\n       https://www.rstudio.com/products/rstudio/download/ \n") 
    }else{
      stop("bigKRLS requires RStudio 1.1.129 or higher on Windows.\n       To use bigKRLS with Windows, switch to RGui or check the following webpages:\n       https://www.rstudio.com/products/rstudio/download/\n       https://dailies.rstudio.com/ \n") 
    }
  }
}
