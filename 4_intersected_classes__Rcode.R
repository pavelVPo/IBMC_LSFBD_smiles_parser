library(tidyverse)

## Files to work with
in_f <- ".../chars_&_symbs.tsv"
out_f <- ".../chars_co-oc.tsv"

## Import the data
data_raw <- read_tsv(in_f)

## Process the data
# Select characters
chars <- data_raw |> select(symbType, symbClass, charClass, chars)
# Create tsv file with the intersections of characters to be further decorated in R, Excel, Google Sheets or somewhere else
string <- chars |> pull(symbType)  |> str_c(collapse = "\t")
write(str_c("\t", "\t", "\t", string, collapse = ""), file=out_f, append=TRUE)
string <- chars |> pull(symbClass) |> str_c(collapse = "\t")
write(str_c("\t", "\t", "\t", string, collapse = ""), file=out_f, append=TRUE)
string <- chars |> pull(charClass) |> str_c(collapse = "\t")
write(str_c("\t", "\t", "\t", string, collapse = ""), file=out_f, append=TRUE)
for (i in seq(1:nrow(chars))) {
	# Reset the string
	string <- paste(chars[i,1] |> pull(), "\t", chars[i,2] |> pull(), "\t", chars[i,3] |> pull(), "\t")
	for (k in seq(1:nrow(chars))) {
		# Get the intersection for the character set corresponding to this row and all other character sets.
		# Represent it as a string
		# This string's length could be used latter to characterize the intersection
		intersection <- intersect(chars[i,4] |> pull() |> str_split(", ") |> unlist(), chars[k,4] |> pull() |> str_split(", ") |> unlist()) |> unlist() |> sort() |> str_c(collapse = ", ")
		if (intersection |> str_length() == chars[i,4] |> pull() |> str_length() & intersection |> str_length() == chars[k,4] |> pull() |> str_length()) {
			intersection <- 10
		} else {
			if (intersection |> str_length() > 0) {
				intersection <- 5
			} else {
				intersection <- 0
			}
		}
		string <- paste0(string, intersection, "\t")
	}
	write(string, file=out_f, append=TRUE)
}