library(tidyverse)

### Import
# Path
path <- "C:/.../"
# Merged classes
data <- read_tsv(str_glue("{path}merged_classes__autoDescribed.tsv")) |>
					separate_longer_delim(basic_class, delim = ", ")
nrows <- nrow(data)

### Check
results <- tibble(merged_class = rep(NA, nrows), basic_class = rep(NA, nrows), identity = rep(NA, nrows))
for (i in seq(1:nrows)) {
	query_vec <- data[i, 2] |> pull() |> str_split(", ") |> unlist()
	target_class <- data[i, 3] |> pull()
	results[i, 1] <- data[i, 1] |> pull()
	results[i, 2] <- target_class
	if(all(query_vec == get(target_class))) {
		results[i, 3] <- TRUE
	} else {
		results[i, 3] <- FALSE
	}
}
if (all(results |> pull(identity))) {
	print("OK")
}