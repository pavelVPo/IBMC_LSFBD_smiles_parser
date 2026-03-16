library(tidyverse)
# Input chars, classes and types
data <- read_tsv(".../chars_&_symbs.tsv") |>
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

###
# Exclude the symbol type pairs, which are not allowed
pairs_symbTypes_not <- read_tsv(".../theory_pairs_symbType.tsv") |>
                            filter(allowed == "no")
pairs_labeled <- pairs_labeled |> anti_join(pairs_symbTypes_not)
# Assess the number of theoretically possible pairs
nrow(pairs_labeled)

###
# Get the distinct theoretically possible pairs of symbol types
pairs_symbClass <- pairs_labeled |> select(left_symbClass, right_symbClass) |> distinct()
# Assess the number of theoretically possible symbol pairs
nrow(pairs_symbClass)
# Export these results
write_tsv(pairs_symbClass, ".../theory_pairs_symbClass.tsv")

###
# 1, 2. Isotope symbols could only be paired with the bracket on the left and bracket atom on the right, also isotope_m could be paired with itself on both sides.
isotope_allowed <- pairs_labeled |> filter(
                                                (left_symbClass == "isotope" & str_detect(right_symbClass, "atom_b")) |
                                                (right_symbClass == "isotope" & left_symbClass == "bracket") |
                                                (left_symbClass == "isotope_m" & (str_detect(right_symbClass, "atom_b") | right_symbClass == "isotope_m")) |
                                                (right_symbClass == "isotope_m" & (left_symbClass == "bracket" | left_symbClass == "isotope_m"))
                                                ) |>
                                 select(left_symbClass, right_symbClass) |> distinct()
# 3. Features symbols, besides isotope, could only be paired with the other feature symbols or bracket symbols on the right and with the other feature symbols or bracket atom symbols on the left.
# Or they can go in pairs with themselves
fs_allowed <- pairs_labeled |> filter(
                                      (str_detect(left_symbClass, "chiral.*")   & (str_detect(right_symbClass, "hydro.*") | str_detect(right_symbClass, "charge.*") | str_detect(right_symbClass, "class.*") | right_symbClass == "bracket")) |
                                      (str_detect(left_symbClass, "hydro.*")    & (str_detect(right_symbClass, "charge.*") | str_detect(right_symbClass, "class.*") | right_symbClass == "bracket")) |
                                      (str_detect(left_symbClass, "charge.*")   & (str_detect(right_symbClass, "class.*") | right_symbClass == "bracket")) |
                                      (str_detect(left_symbClass, "class.*")    & right_symbClass == "bracket") |
                                      (str_detect(right_symbClass, "chiral.*")  & str_detect(left_symbClass, "atom_b.*")) |
                                      (str_detect(right_symbClass, "hydro.*")   & (str_detect(left_symbClass, "chiral.*") | str_detect(left_symbClass, "atom_b.*"))) |
                                      (str_detect(right_symbClass, "charge.*")  & (str_detect(left_symbClass, "chiral.*") | str_detect(left_symbClass, "hydro.*") | str_detect(left_symbClass, "atom_b.*"))) |
                                      (str_detect(right_symbClass, "class.*")   & (str_detect(left_symbClass, "chiral.*") | str_detect(left_symbClass, "hydro.*") | str_detect(left_symbClass, "charge.*") | str_detect(left_symbClass, "atom_b.*"))) |
                                      (str_detect(left_symbClass,  "chiral_.*|hydro_.*|charge_.*|class_.*") & left_symbClass == right_symbClass) |
                                      (str_detect(right_symbClass, "chiral_.*|hydro_.*|charge_.*|class_.*") & left_symbClass == right_symbClass)
                               ) |>
                               select(left_symbClass, right_symbClass) |> distinct()
# 5. Bracket atoms cannot be paired with the symbols contained outside the brackets and they can only be paired with the isotope symbols if those symbols are on the left and they can only be paired with the other bracket atom symbols on condition that those symbols has the same length, which is greater than 1.
ba_allowed <- pairs_labeled |> filter(
                                        (str_detect(left_symbClass, "atom_b.*")     & str_detect(right_symbClass, "chiral.*|hydro.*|charge.*|class.*|chiral.*|bracket")) |
                                        (str_detect(right_symbClass, "atom_b.*")    & str_detect(left_symbClass,  "bracket|isotope.*")) |
                                        (str_detect(left_symbClass, "atom_b.*_.*")  & right_symbClass == left_symbClass) |
                                        (str_detect(right_symbClass, "atom_b.*_.*") & right_symbClass == left_symbClass)
                                ) |>
                                select(left_symbClass, right_symbClass) |> distinct()
# 6. Organic atom symbols can not be paired with the symbols contained inside the brackets.
oa_allowed <- pairs_labeled |> filter(
                                        (str_detect(left_symbClass, "atom_o.*")    & !str_detect(right_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*")) |
                                        (str_detect(right_symbClass, "atom_o.*")   & !str_detect(left_symbClass,  "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*"))
                                ) |>
                                select(left_symbClass, right_symbClass) |> distinct()
# 7. Bond symbols and bond modifying symbols can not be paired with the symbols contained inside the brackets.
bbm_allowed <- pairs_labeled |> filter(
                                        (str_detect(left_symbClass, ".*_bond") &  !str_detect(right_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*")) |
                                        (str_detect(right_symbClass, ".*_bond") &  !str_detect(left_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*")) |
                                        (str_detect(left_symbClass, "bm.*") &  !str_detect(right_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*")) |
                                        (str_detect(right_symbClass, "bm.*") &  !str_detect(left_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*"))
                                ) |>
                                select(left_symbClass, right_symbClass) |> distinct()
# 8. Cis/trans symbols are only allowed on the left side from the organic atoms and brackets.
ct_allowed <- pairs_labeled |> filter(
                                        (left_symbClass == "ct" & (str_detect(right_symbClass, "atom_o.*_.*") | right_symbClass == "bracket")) |
                                        (right_symbClass == "ct" & (!str_detect(left_symbClass, "atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*"))) 
                               ) |>
                               select(left_symbClass, right_symbClass) |> distinct()
# Delete records, to which the rules were applied from the main set of pairs
pairs_unknown <- pairs_labeled |> filter(
                                        !str_detect(left_symbClass,  "ct|bm.*|.*_bond|atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*|chiral.*") &
                                        !str_detect(right_symbClass, "ct|bm.*|.*_bond|atom_b.*|isotope.*|chiral.*|hydro.*|charge.*|class.*|chiral.*")
                                ) |>
                                select(left_symbClass, right_symbClass) |> distinct()
# Gather labeled pairs
pairs_rulesApplied <- bind_rows(isotope_allowed, fs_allowed, ba_allowed, oa_allowed, bbm_allowed, pairs_unknown) |> distinct()
# Export these prefiltered results
write_tsv(pairs_rulesApplied, ".../theory_pairs_symbClass_rules.tsv")