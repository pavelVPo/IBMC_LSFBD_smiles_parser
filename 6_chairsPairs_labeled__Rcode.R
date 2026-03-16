library(tidyverse)
# Input chars, classes and types
data <- read_tsv("C:/.../chars_&_symbs.tsv") |>
                select(chars, charClass, symbClass, symbType) |>
                separate_longer_delim(chars, delim = ", ")


###
# Get unique characters
unique_chars <- data |> pull(chars) |> unique()
# Generate all the pairs possible in theory
pairs <- expand.grid(left_char = unique_chars, right_char = unique_chars)
# Get all the theoretically possible pairs considering classes and types
pairs_labeled <- pairs |> inner_join(data, by = c("left_char" = "chars"), relationship = "many-to-many") |>
                          rename(left_charClass = charClass, left_symbClass = symbClass, left_symbType = symbType) |>
                          inner_join(data, by = c("right_char" = "chars"), relationship = "many-to-many") |>
                          rename(right_charClass = charClass, right_symbClass = symbClass, right_symbType = symbType)
# Assess the number of theoretically possible pairs
nrow(pairs_labeled)


###
# Get the distinct theoretically possible pairs of symbol types
pairs_symbType <- pairs_labeled |> select(left_symbType, right_symbType) |> distinct()
# Assess the number of theoretically possible symbol pairs
nrow(pairs_symbType)
# Export these results
write_tsv(pairs_symbType, "C:/.../theory_pairs_symbType.tsv")