library(tidyverse)
### Functions
# SEE: https://stackoverflow.com/questions/71777528/how-to-filter-not-in
`%nin%` <- Negate(`%in%`)

### Input classes of symbols
data <- read_tsv("C:/.../chars_&_symbs.tsv")
data_class   <- data |> select(symbClass, symbType) |> distinct()
allowed_type <- read_tsv("C:/.../theory_pairs_symbType.tsv") |>
                    filter(allowed == "yes") |>
                    select(left_symbType, right_symbType) |>
                    distinct()


### Process
##  Generate all the possible pairs of symbol classes
data_combs <- data_class |> select(symbClass) |>
                      summarize(symbClass = str_c(symbClass, collapse = ", ")) |>
                      ungroup() |>
                      mutate(left = symbClass, right = symbClass) |>
                      select(left, right)

##  Filter according to the list of allowed pairs of symbol types
data_combs <- data_combs |> separate_longer_delim(left, delim = ", ") |>
                    separate_longer_delim(right, delim = ", ") |>
                    distinct() |>
                    inner_join(data_class, by =c("left" = "symbClass")) |>
                    rename(left_type = symbType) |>
                    inner_join(data_class, by =c("right" = "symbClass")) |>
                    rename(right_type = symbType) |>
                    inner_join(allowed_type, by = c("left_type" = "left_symbType", "right_type" = "right_symbType"))

##  Filter according to the rules of symbol classes pairings
data_rules <- data_combs |> filter(
                                # aromatic bond
                                !(left == "aromatic_bond"  & right %in% c("atom_bal", "atom_bal_2", "atom_oal", "atom_oal_2")) &
                                !(right == "aromatic_bond" & left  %in% c("atom_bal", "atom_bal_2", "atom_oal", "atom_oal_2")) &
                                # isotopes
                                !(left  == "isotope"    & right %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(left  == "isotope_m"  & right %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "isotope"    & left %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "isotope_m"  & left %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                # feature symbols, chiral
                                !(left  == "chiral"     & right %nin% c("hydro", "hydro_2", "charge", "charge_2", "charge_m", "class_m", "bracket")) &
                                !(left  == "chiral_2"   & right %nin% c("hydro", "hydro_2", "charge", "charge_2", "charge_m", "class_m", "bracket")) &
                                !(left  == "chiral_m"   & right %nin% c("hydro", "hydro_2", "charge", "charge_2", "charge_m", "class_m", "bracket")) &
                                !(right == "chiral"     & left  %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "chiral_2"   & left  %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "chiral_m"   & left  %nin% c("atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                # feature symbols, hydrogen count
                                !(left  == "hydro"      & right %nin% c("charge", "charge_2", "charge_m", "class_m", "bracket")) &
                                !(left  == "hydro_2"    & right %nin% c("charge", "charge_2", "charge_m", "class_m", "bracket")) &
                                !(right == "hydro"      & left  %nin% c("chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "hydro_2"    & left  %nin% c("chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                # feature symbols, charge
                                !(left  == "charge"     & right %nin% c("class_m", "bracket")) &
                                !(left  == "charge_2"   & right %nin% c("class_m", "bracket")) &
                                !(left  == "charge_m"   & right %nin% c("class_m", "bracket")) &
                                !(right == "charge"     & left  %nin% c("hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "charge_2"   & left  %nin% c("hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                !(right == "charge_m"   & left  %nin% c("hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                # feature symbols, class
                                !(left  == "class_m"    & right %nin% c("bracket")) &
                                !(right == "class_m"    & left  %nin% c("charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2")) &
                                # bracket atoms
                                !(left   %in% c("atom_bal", "atom_bar") & right %nin% c("class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "bracket")) &
                                !(left   == "atom_bal_2" & right %nin% c("atom_bal_2", "class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "bracket")) &
                                !(left   == "atom_bar_2" & right %nin% c("atom_bar_2", "class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "bracket")) &
                                !(right  %in% c("atom_bal", "atom_bar") & left  %nin% c("isotope_m", "isotope", "bracket")) &
                                !(right  == "atom_bal_2" & left  %nin% c("atom_bal_2", "isotope_m", "isotope", "bracket")) &
                                !(right  == "atom_bar_2" & left  %nin% c("atom_bar_2", "isotope_m", "isotope", "bracket")) &
                                # organic atoms lacking additional properties
                                !(left  %in% c("atom_oal", "atom_oal_2", "atom_oar", "atom_oar_2") & right %in% c("class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2", "isotope", "isotope_m")) &
                                !(right %in% c("atom_oal", "atom_oal_2", "atom_oar", "atom_oar_2") & left  %in% c("class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2", "isotope", "isotope_m")) &
                                # bonds and bond modifying symbols should not be paired with the atoms and features inside the brackets
                                !(left  %in% c("aromatic_bond", "double_bond", "no_bond", "quadruple_bond", "single_bond", "triple_bond", "bm_ibi", "bm_iri", "bm_tbi", "bm_tri", "bm_ibe_2", "bm_ire_2", "bm_iri_3", "bm_ire_4", "bm_tbe_2", "bm_tre_2", "bm_tre_4", "bm_tri_3") &
                                  right %in% c("class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2", "isotope", "isotope_m")) &
                                !(right  %in% c("aromatic_bond", "double_bond", "no_bond", "quadruple_bond", "single_bond", "triple_bond", "bm_ibi", "bm_iri", "bm_tbi", "bm_tri", "bm_ibe_2", "bm_ire_2", "bm_iri_3", "bm_ire_4", "bm_tbe_2", "bm_tre_2", "bm_tre_4", "bm_tri_3") &
                                  left %in% c("class_m", "charge", "charge_2", "charge_m", "hydro", "hydro_2", "chiral", "chiral_2", "chiral_m", "atom_bal", "atom_bal_2", "atom_bar", "atom_bar_2", "isotope", "isotope_m")) &
                                # bond modifying symbols, initiators and terminators of branching
                                !(left %in% c("bm_ibi", "bm_ibe_2") & right %in% c("bm_tbi", "bm_tbe_2")) &
                                # two character (explicit) bond modifying symbols on the left, initiators (terminator) should not have new branching iniator (terminator) on the right
                                !(left %in% c("bm_ibe_2") & right %in% c("bm_ibi")) &
                                !(left %in% c("bm_ibe_2") & right %in% c("bm_ibe_2")) &
                                !(left %in% c("bm_tbe_2") & right %in% c("bm_tbi")) &
                                !(left %in% c("bm_tbe_2") & right %in% c("bm_tbe_2"))
                            )


write_tsv(data_rules, "C:/.../theory_pairs_symbClass_rules.tsv")