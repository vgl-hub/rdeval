input_file <- "../Documents/GitHub/rdeval/bTaeGut7_PAG68436_20230103.cat.fastq.gz.rd"

# rdFile class definition to store object data
rdFileClass <- setClass("rdFileClass", slots = c(md5s = "list", lengths = "list", qualities = "list"))
setGeneric("md5s", function(x) standardGeneric("md5s"))
setGeneric("md5s<-", function(x, value) standardGeneric("md5s<-"))
setMethod("md5s", "rdFileClass", function(x) x@md5s)
setMethod("md5s<-", "rdFileClass", function(x, value) { x@md5s <- value; x })

setGeneric("lengths<-", function(x, value) standardGeneric("lengths<-"))
setMethod("lengths<-", "rdFileClass", function(x, value) { x@lengths <- value; x })

setGeneric("qualities<-", function(x, value) standardGeneric("qualities<-"))
setMethod("qualities<-", "rdFileClass", function(x, value) { x@qualities <- value; x })

# new rdFile instance
rdFile <- rdFileClass()

# open connection to file to read the header
rd_path <- gzfile(input_file)
open(rd_path, "r+b")

# read number of md5s stored in file
md5sN <- readBin(rd_path, integer(), n = 1, size = 4, endian = "little")
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
  md5s(rdFile) <- c(rdFile@md5s, list(c(filename, md5)))
}

# read uncompressed size
uncompressed_size <- readBin(rd_path, integer(), n = 1, size = 8, endian = "little")
header_size <- header_size + 8

# read gzip compressed information and decompress it
data <- memDecompress(readBin(rd_path, "raw", n = file.info(input_file)$size-header_size, size = 1, endian = "little"), "gzip")

# we are done with the file
close(rd_path)

# read ACGTN
A <- readBin(data[1:8], integer(), n = 1, size = 8, endian = "little")
C <- readBin(data[9:16], integer(), n = 1, size = 8, endian = "little")
G <- readBin(data[17:24], integer(), n = 1, size = 8, endian = "little")
T <- readBin(data[25:32], integer(), n = 1, size = 8, endian = "little")
N <- readBin(data[33:40], integer(), n = 1, size = 8, endian = "little")

# read counts for types
len8 <- readBin(data[41:48], integer(), n = 1, size = 8, endian = "little")
len16 <- readBin(data[49:56], integer(), n = 1, size = 8, endian = "little")
len64 <- readBin(data[57:64], integer(), n = 1, size = 8, endian = "little")

# read length and quality for types
len8_matrix  <- matrix(data[65:(64+(8*len8))], ncol = 8, byrow = TRUE)
len16_matrix <- matrix(data[(65+(8*len8)):(64+(8*len8)+(8*len16))], ncol = 8, byrow = TRUE)
len64_matrix <- matrix(data[(65+(8*len8)+(8*len16)):(64+(8*len8)+(8*len16)+(16*len64))], ncol = 16, byrow = TRUE)

# append to rdFile object
lengths(rdFile) <- c(rdFile@lengths, readBin(as.raw(t(cbind(len64_matrix[,1:2]))), integer(), n = len64, size = 8, endian = "little"))
qualities(rdFile) <- c(rdFile@qualities, readBin(as.raw(t(cbind(len64_matrix[,5:8]))), numeric(), n = len64, size = 4, endian = "little"))
lengths(rdFile) <- c(rdFile@lengths, readBin(as.raw(t(cbind(len16_matrix[,1:2]))), integer(), n = len16, size = 2, endian = "little", signed = FALSE))
qualities(rdFile) <- c(rdFile@qualities, readBin(as.raw(t(cbind(len16_matrix[,5:8]))), numeric(), n = len16, size = 4, endian = "little"))
lengths(rdFile) <- c(rdFile@lengths, readBin(len8_matrix[,1], integer(), n = len8, size = 1, endian = "little", signed = FALSE))
qualities(rdFile) <- c(rdFile@qualities, readBin(as.raw(t(cbind(len8_matrix[,5:8]))), numeric(), n = len8, size = 4, endian = "little"))
