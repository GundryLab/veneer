processZip <- function(zipfile) {
  td = tempdir()
  ud = (paste0(td,'/unzipped'))
  if( dir.exists(ud) ){
    unlink( paste0(ud,"/*") )
  } else {
    dir.create(ud)
  }
  unzip(zipfile, exdir=ud)
  files <- list.files(ud, full.names=TRUE)
  
  return(files)
}  

processSingle <- function(datapath, filename) {
  td = tempdir()
  ud = (paste0(td,'/unzipped'))
  if( dir.exists(ud) ){
    unlink( paste0(ud,"/*") )
  } else {
    dir.create(ud)
  }
  file.copy(datapath, paste0(ud, "/", filename))
  files <- list.files(ud, full.names=TRUE)
  
  return(files)
  
}


