PRJEB36534 <- read.csv("input_ena_lists/filereport_read_run_PRJEB36534.tsv", sep = "\t")
download_list = c(unlist(strsplit(PRJEB36534$fastq_ftp[23:27], ";")))

PRJNA597671 <- read.csv("input_ena_lists/filereport_read_run_PRJNA597671.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJNA597671$fastq_ftp, ";")))

PRJNA823845 <- read.csv("input_ena_lists/filereport_read_run_PRJNA823845.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJNA823845$fastq_ftp[order(PRJNA823845$sra_bytes, decreasing = T)[1:5]], ";")))

PRJNA630300 <- read.csv("input_ena_lists/filereport_read_run_PRJNA630300.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJNA630300$fastq_ftp[order(PRJNA630300$sra_bytes, decreasing = T)[1:5]], ";")))

PRJEB38557 <- read.csv("input_ena_lists/filereport_read_run_PRJEB38557.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJEB38557$fastq_ftp, ";")))

PRJEB41174 <- read.csv("input_ena_lists/filereport_read_run_PRJEB41174.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJEB41174$fastq_ftp[order(PRJEB41174$sra_bytes, decreasing = T)[1:10]], ";")))

PRJEB31111 <- read.csv("input_ena_lists/filereport_read_run_PRJEB31111.tsv", sep = "\t")
download_list = c(download_list, unlist(strsplit(PRJEB31111$fastq_ftp, ";")))

download_list <- unique(gsub("_.*", "", gsub(".*/", "", download_list)))

write.table(data.frame(download_list), col.names = F, "soil_files_to_download.tsv", sep = "\t", row.names = F, quote = F)

