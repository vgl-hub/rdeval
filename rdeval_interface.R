library("bit64")

input_file <- "bTaeGut7-dna3_1.fastq.gz.rd"

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

# usage examples:
tail(rdFile$lengths, n=1) # shortest read (they are sorted)
rdFile$lengths[1] # longest read
mean(rdFile$lengths) # average read length
mean(rdFile$qualities) # average read quality (not weighted by length)
rdFile$A_count # As count etc
input_file # name of the input
rdFile$md5s # md5s of files that were used to generate the .rd file
