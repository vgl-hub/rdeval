library("bit64")

input_file <- "m54306Ue_230222_164948.bc1015--bc1015.hifi_reads.fastq.gz.rd"

# rdFile class definition to store object data
rdFileClass <- setRefClass("rdFileClass", fields = list(input_file = "character", 
                                                        md5s = "list", 
                                                        lengths = "vector", 
                                                        qualities = "vector", 
                                                        A_count = "integer64", 
                                                        C_count = "integer64", 
                                                        G_count = "integer64", 
                                                        T_count = "integer64", 
                                                        N_count = "integer64"))
generateRdFile <- function(input_file) {
  # open connection to file to read the header
  rd_path <- gzfile(input_file)
  open(rd_path, "r+b")
  
  # read number of md5s stored in file
  md5sN <- readBin(rd_path, integer(), n = 1, size = 4, endian = "little")
  md5s <- list()
  # header size
  header_size <- 4
  
  # read each filename and md5 in rdFile object
  for (x in 1:md5sN) {
    n_char <- readBin(rd_path, integer(), n = 1, size = 2, endian = "little", signed = FALSE)
    filename <- rawToChar(readBin(rd_path, raw(), n = n_char, endian = "little"))
    header_size <- header_size + 2 + n_char
    n_char <- readBin(rd_path, integer(), n = 1, size = 2, endian = "little", signed = FALSE)
    md5 <- rawToChar(readBin(rd_path, raw(), n = n_char, endian = "little"))
    header_size <- header_size + 2 + n_char
    md5s <- c(md5s, list(c(filename, md5)))
  }
  
  # read uncompressed size
  uncompressed_size <- readBin(rd_path, numeric(), n = 1, size = 8, endian = "little")
  class(uncompressed_size)<-"integer64";
  header_size <- header_size + 8
  
  # read gzip compressed information and decompress it
  data <- memDecompress(readBin(rd_path, "raw", n = file.info(input_file)$size-header_size, size = 1, endian = "little"), "gzip")
  
  # we are done with the file
  close(rd_path)
  
  # read ACGTN
  A <- readBin(data[1:8], numeric(), n = 1, size = 8, endian = "little")
  C <- readBin(data[9:16], numeric(), n = 1, size = 8, endian = "little")
  G <- readBin(data[17:24], numeric(), n = 1, size = 8, endian = "little")
  T <- readBin(data[25:32], numeric(), n = 1, size = 8, endian = "little")
  N <- readBin(data[33:40], numeric(), n = 1, size = 8, endian = "little")
  class(A)<-"integer64"; class(C)<-"integer64"; class(G)<-"integer64"; class(T)<-"integer64"; class(N)<-"integer64"
  
  # read counts for types
  len8 <- readBin(data[41:48], numeric(), n = 1, size = 8, endian = "little")
  len16 <- readBin(data[49:56], numeric(), n = 1, size = 8, endian = "little")
  len64 <- readBin(data[57:64], numeric(), n = 1, size = 8, endian = "little")
  class(len8)<-"integer64"; class(len16)<-"integer64"; class(len64)<-"integer64";
  
  if (len8 != as.integer(len8) | len16 != as.integer(len16) | len64 != as.integer(len64)) {
    print("Error: too many reads to display (>2^32)")
  }else{
    len8<-as.integer(len8); len16<-as.integer(len16); len64<-as.integer(len64)
  }
  
  # read length and quality for types
  len8_matrix  <- matrix(data[65:(64+(8*len8))], ncol = 8, byrow = TRUE)
  len16_matrix <- matrix(data[(65+(8*len8)):(64+(8*len8)+(8*len16))], ncol = 8, byrow = TRUE)
  len64_matrix <- matrix(data[(65+(8*len8)+(8*len16)):(64+(8*len8)+(8*len16)+(16*len64))], ncol = 16, byrow = TRUE)
  
  lengths <- vector()
  qualities <- vector()
  # append to rdFile object
  lengths <- c(lengths, readBin(as.raw(t(cbind(len64_matrix[,1:8]))), integer(), n = len64, size = 8, endian = "little")) # note that reads >2^32 will be trimmed
  qualities <- c(qualities, readBin(as.raw(t(cbind(len64_matrix[,9:12]))), numeric(), n = len64, size = 4, endian = "little"))
  lengths <- c(lengths, readBin(as.raw(t(cbind(len16_matrix[,1:2]))), integer(), n = len16, size = 2, endian = "little", signed = FALSE))
  qualities <- c(qualities, readBin(as.raw(t(cbind(len16_matrix[,5:8]))), numeric(), n = len16, size = 4, endian = "little"))
  lengths <- c(lengths, readBin(len8_matrix[,1], integer(), n = len8, size = 1, endian = "little", signed = FALSE))
  qualities <- c(qualities, readBin(as.raw(t(cbind(len8_matrix[,5:8]))), numeric(), n = len8, size = 4, endian = "little"))
  
  # new rdFile instance
  rdFile <- rdFileClass$new(input_file = input_file,
                            md5s = md5s,
                            A_count = A, 
                            C_count = C, 
                            G_count = G, 
                            T_count = T, 
                            N_count = N,
                            lengths = lengths,
                            qualities = qualities)
  return(rdFile)
}

printRdSummary <- function(rdFile) {
  if (class(rdFile)[1] == 'rdFileClass') {
    gc_content = round(
      (rdFile$C_count + rdFile$G_count)/(rdFile$A_count + rdFile$C_count + rdFile$T_count + rdFile$G_count) * 100, 4)
    cat(paste0("### ", rdFile$md5s[[1]][1], "\n"))
    cat(paste0("Number of reads: ", length(rdFile$lengths), "\n\n"))
    cat(paste0("Total read length: ", sum(rdFile$lengths), "\n\n"))
    cat(paste0("Average read length: ", round(mean(rdFile$lengths), 1), "\n\n"))
    cat(paste0("Read N50: ", quantile(rdFile$lengths, .5)[['50%']], "\n\n"))
    cat(paste0("Smallest read length: ", tail(rdFile$lengths, n=1), "\n\n"))
    cat(paste0("Largest read length: ", rdFile$lengths[1], "\n\n"))
    cat(paste0("GC content %: ", gc_content, "\n\n"))
    cat(paste0("Base composition (A:C:T:G): ", 
                 rdFile$A_count, ":", 
                 rdFile$C_count, ":", 
                 rdFile$T_count, ":",
                 rdFile$G_count, "\n\n"))
    cat(paste0("Average read quality: ", round(mean(rdFile$qualities)), "\n\n"))
    cat()
  } else {
    print('Input is not of class rdFileClass.')
  }
}
