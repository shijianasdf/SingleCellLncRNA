#' @description   导入程序库和数据集
#' @author shi jian

#' @export
QuickStart <- function(script.path,
                       data.path,
                       result.path){
  ## create result path
  if(!file.exists(result.path)){
    cat(paste0("Create ",result.path))
    dir.create(result.path,recursive = T)
  }
  result.path <<- result.path #设置全局变量
  data.path <<- data.path
  ## run scripts
  s <- list.files(path = script.path,pattern = ".R",full.names = T)
  for(i in s) source(i,encoding = "UTF-8")
}

#' @export
## load1 help load local .rda files by a pattern
load1 <- function(pattern,path=".",envir=parent.frame()){
  for(i in pattern){
    file = list.files(path = path,pattern = i,full.names = T)
    if(length(file) > 1){
      print("Attention!Multiple files.")
    }
    for(j in file){
      load(j,envir = envir)
    }
    
  }
}
QuickStart("D:/Rsources/Rcode/BioinforRCode-master/R","D:/Rsources/Rcode/BioinforRCode-master/data","D:Rsources/Rcode/results")
load1(".rda","D:/Rsources/Rcode/BioinforRCode-master/data")
