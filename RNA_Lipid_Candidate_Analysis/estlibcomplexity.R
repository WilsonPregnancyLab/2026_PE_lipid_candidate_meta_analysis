list_est_comp <- list.files(path = ".", pattern = ".txt", full.names = FALSE)     

for (file in list_est_comp) {
    lines <- readLines(file)[7:8]
    headers <- strsplit(lines[1], "\t") [[1]]
    values <- strsplit(lines[2], "\t")[[1]]
    named_values <- setNames(as.numeric(values), headers)
    data_list[[basename(file)]] <- named_values
}

combined_est_lib <- do.call(cbind, data_list)

