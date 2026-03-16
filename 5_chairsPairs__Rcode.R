library(tidyverse)
# Input chars, classes and types
data <- read_tsv("C:/.../chars_&_symbs.tsv") |>
				select(chars, charClass, symbClass, symbType) |>
				separate_longer_delim(chars, delim = ", ")
# Get unique characters
unique_chars <- data |> pull(chars) |> unique()
# Generate all the pairs possible in theory
pairs <- expand.grid(left_char = unique_chars, right_char = unique_chars)

# Assess the number of theoretically possible pairs
nrow(pairs)