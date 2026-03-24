library(tidyverse)

# NOT IN function, SEE: https://stackoverflow.com/questions/67140696/what-does-in-functionx-yinx-y-mean
'%!in%' <- function(x,y)!('%in%'(x,y))

## Input additional data on character classes
char_classes <- read_tsv("C:/.../chars_&_symbs.tsv")

## Prefiltering
# Input chars, classes and types
data_raw <- read_tsv("C:/.../chars_&_symbs.tsv") |>
                select(charClass, symbClass) |>
                distinct()
# Input the list of the allowed symbol classes
sc_allowed <- read_tsv("C:/.../theory_pairs_symbClass_rules.tsv")
# Prepare the data for furhter filtering on the level of character classes
data_all <- sc_allowed |> inner_join(data_raw, by = c("left_symbClass" = "symbClass"), relationship = "many-to-many") |>
					  rename(left_charClass = charClass) |>
					  inner_join(data_raw, by = c("right_symbClass" = "symbClass"), relationship = "many-to-many") |>
					  rename(right_charClass = charClass) |>
					  distinct() |>
					  select(left_charClass, right_charClass, left_symbClass, right_symbClass)
# Export pairs
write_tsv(data_all, "C:/.../theory_pairs_charClasses_rules.tsv")

## Filtering
# Get the correct subset of single character symbols
whole_allowed <- data_all |> filter( (left_charClass %in% c("w_anything", "w_aromatic_bond", "w_atom_bal", "w_atom_bar", "w_atom_oal", "w_atom_oar",
															"w_bm_ibi", "w_bm_iri", "w_bm_tbi", "w_bm_tri", "w_charge", "w_chiral", "w_double_bond",
															"w_hydro", "w_isotope", "w_no_bond", "w_quadruple_bond", "w_single_bond", "w_triple_bond") &
									  right_charClass %!in% c("e_atom_bal", "e_atom_bar", "e_atom_oal", "e_bm_ibe", "e_bm_ire", "e_bm_tbe", "e_bm_tre",
									  						  "e_charge", "e_chiral", "e_hydro", "m_chiral", "n_bm_ire", "n_bm_tre", "r_bm_ire",
									  						  "r_bm_iri", "r_bm_tre", "r_bm_tri", "r_charge", "r_chiral", "r_class", "r_isotope")) |
									  (right_charClass %in% c("w_anything", "w_aromatic_bond", "w_atom_bal", "w_atom_bar", "w_atom_oal", "w_atom_oar",
															"w_bm_ibi", "w_bm_iri", "w_bm_tbi", "w_bm_tri", "w_charge", "w_chiral", "w_double_bond",
															"w_hydro", "w_isotope", "w_no_bond", "w_quadruple_bond", "w_single_bond", "w_triple_bond") &
									   left_charClass %!in% c("m_chiral", "n_bm_ire", "n_bm_tre", "s_atom_bal", "s_atom_bar", "s_atom_oal", "s_bm_ibe", "s_bm_ire",
									   							"s_bm_ire_4", "s_bm_iri", "s_bm_tbe", "s_bm_tre", "s_bm_tre_4", "s_bm_tri", "s_charge", "s_charge_m",
									   							"s_chiral", "s_chiral_m", "s_class", "s_hydro", "s_isotope")) )
data_inProgress <- data_all |> anti_join(whole_allowed, by = c("left_charClass", "right_charClass"))
# Get the correct subset including first character of the multicharacter class
start_allowed <- data_inProgress |> filter(
								# left start, 2 character
								(left_charClass %in% c("s_atom_bal", "s_atom_bar", "s_hydro", "s_chiral", "s_charge", "s_atom_oal", "s_bm_ibe", "s_bm_ire", "s_bm_tbe", "s_bm_tre") & right_charClass == str_replace(left_charClass, "^s_", "e_")) |
								# left start, multiple
								(left_charClass == "s_chiral_m" & right_charClass == "m_chiral_m") |
								# left start, 4-character
								(left_charClass %in% c("s_bm_ire_4", "s_bm_tre_4") & right_charClass == str_replace(left_charClass, "^s_", "n_")) |
								# left start, rest
								(left_charClass %in% c("s_isotope", "s_charge_m", "s_class", "s_bm_iri", "s_bm_tri") & right_charClass == str_replace(left_charClass, "^s_", "r_"))
							)
data_inProgress <- data_all |> anti_join(start_allowed, by = c("left_charClass", "right_charClass")) |> filter(left_charClass %!in% c("s_isotope", "s_charge_m", "s_class", "s_bm_iri", "s_bm_tri", "s_bm_ire_4", "s_bm_tre_4", "s_chiral_m", "s_atom_bal", "s_atom_bar", "s_hydro", "s_chiral", "s_charge", "s_atom_oal", "s_bm_ibe", "s_bm_ire", "s_bm_tbe", "s_bm_tre"))
# Get the correct subset including intermediate character of the multicharacter class
inter_allowed <- data_inProgress |> filter (
										# left, m_chiral case
										(left_charClass == "m_chiral" & (right_charClass == "m_chiral" | right_charClass == "r_chiral")) |
										# left, n_bm cases
										(left_charClass == "n_bm_ire" & (right_charClass == "r_bm_ire")) |
										(left_charClass == "n_bm_tre" & (right_charClass == "r_bm_tre")) |
										# right, m_chiral case
										(right_charClass == "m_chiral" & (left_charClass == "s_chiral_m" | left_charClass == "m_chiral"))
									)
# Get the correct subset including last character of the multicharacter class
end_allowed <- data_inProgress |> filter(
										# left end
										(left_charClass %in% c("e_atom_bal", "e_atom_bar", "e_atom_oal", "e_bm_ibe", "e_bm_ire", "e_bm_tbe", "e_bm_tre", "e_bracket", "e_charge", "e_chiral", "e_hydro",
											"r_bm_ire", "r_bm_iri", "r_bm_tre", "r_bm_tri", "r_charge", "r_chiral", "r_class", "r_ct", "l_ct", "r_isotope") &
										(right_charClass %in% c("l_ct", "r_ct", "s_atom_bal", "s_atom_bar", "s_atom_oal", "s_bm_ibe", "s_bm_ire", "s_bm_ire_4", "s_bm_iri", "s_bm_tbe", "s_bm_tre", "s_bm_tre_4", "s_bm_tri", "s_bracket", "s_charge", "s_charge_m", "s_chiral", "s_chiral_m", "s_class", "s_hydro", "s_isotope", "w_anything", "w_aromatic_bond", "w_atom_bal", "w_atom_bar", "w_atom_oal", "w_atom_oar", "w_bm_ibi", "w_bm_iri", "w_bm_tbi", "w_bm_tri", "w_charge", "w_chiral", "w_double_bond", "w_hydro", "w_isotope", "w_no_bond", "w_quadruple_bond",  "w_single_bond", "w_triple_bond")))
								  ) |> filter( !(left_charClass == "r_isotope" & right_charClass == "s_isotope") & 
									  		   !(left_charClass == "r_chiral" & right_charClass == "s_chiral_m") &
									  		   !(left_charClass == "e_hydro" & right_charClass == "s_hydro") &
									  		   !(left_charClass == "e_charge" & right_charClass == "s_charge") &
								   		  	   !(left_charClass == "e_charge" & right_charClass == "s_charge_m") &
								   		  	   !(left_charClass == "e_chiral" & right_charClass == "s_chiral") &
								   		  	   !(left_charClass == "e_atom_bal" & right_charClass == "s_atom_bal") &
								   		  	   !(left_charClass == "e_atom_bar" & right_charClass == "s_atom_bar")
								   	 )
# Process further only pairs, which have the correct order
data_inProgress <- bind_rows(whole_allowed, start_allowed, inter_allowed, end_allowed) |> distinct()
# Get the correct subset for l_ct and r_ct
ct_allowed <- data_inProgress |> filter(
									# lct is on the left
									(left_charClass == "l_ct"  & right_charClass %in% c("s_bracket", "w_anything", "w_atom_oal", "w_atom_oar", "s_atom_oal", "s_bm_ire", "s_bm_ire_4", "s_bm_iri", "s_bm_tre", "s_bm_tre_4", "s_bm_tri", "w_bm_iri", "s_bm_tri")) |
									# rct is on the left
									(left_charClass == "r_ct"  & right_charClass %in% c("s_bracket", "w_anything", "w_atom_oal", "w_atom_oar", "s_atom_oal", "s_bm_ire", "s_bm_ire_4", "s_bm_iri", "s_bm_tre", "s_bm_tre_4", "s_bm_tri", "w_bm_iri", "s_bm_tri")) |
									# lct is on the right
									(right_charClass == "l_ct" & left_charClass  %in% c("w_atom_oal", "w_atom_oar", "w_bm_ibi", "w_bm_iri", "w_bm_tri", "w_bm_tbi")) |
									# rct is on the right
									(right_charClass == "r_ct" & left_charClass  %in% c("w_atom_oal", "w_atom_oar", "w_bm_ibi", "w_bm_iri", "w_bm_tri", "w_bm_tbi"))
								)
data_inProgress <- data_inProgress |> anti_join(ct_allowed, by = c("left_charClass", "right_charClass")) |> filter(left_charClass != "l_ct" & right_charClass != "l_ct") |> filter(left_charClass != "r_ct" & right_charClass != "r_ct")
# Get the correct subset including start and end of brackets
bracket_allowed   <- data_inProgress |> filter(
											# bracket start is on the left
											(left_charClass == "s_bracket"   & right_charClass %in% c("w_isotope", "s_isotope", "w_anything", "w_atom_bar", "s_atom_bar", "w_atom_bal", "s_atom_bal")) |
											# bracket start is on the right
											( right_charClass == "s_bracket" & left_charClass %in% c("e_bracket", "l_ct", "r_ct", "w_anything", "w_aromatic_bond", "w_atom_oal", "w_atom_oar", "w_bm_ibi", "w_bm_iri", "w_bm_tbi", "w_bm_tri", "w_no_bond", "w_quadruple_bond", "w_single_bond", "w_triple_bond")) |
											# bracket end is on the right
											(right_charClass == "e_bracket"  & left_charClass %in% c("w_anything", "w_atom_bar", "w_atom_bal", "e_atom_bar", "e_atom_bal", "e_charge", "e_chiral", "e_hydro", "r_charge", "r_chiral", "r_class", "w_hydro")) |
											# bracket end is on the left
											(left_charClass == "e_bracket"   & right_charClass %in% c("s_bracket", "l_ct", "r_ct", "w_anything", "w_aromatic_bond", "w_atom_oal", "w_atom_oar", "w_bm_ibi", "w_bm_iri", "w_bm_tbi", "w_bm_tri", "w_no_bond", "w_quadruple_bond", "w_single_bond", "w_triple_bond"))
										)
data_inProgress <- data_inProgress |> anti_join(bracket_allowed, by = c("left_charClass", "right_charClass")) |> filter(left_charClass != "s_bracket") |> filter(right_charClass != "s_bracket") |> filter(left_charClass != "e_bracket") |> filter(right_charClass != "e_bracket")
# Filter to re-check the rest
data_inProgress <- data_inProgress |> anti_join(start_allowed, by = c("left_charClass", "right_charClass"))
data <- bind_rows(data_inProgress, ct_allowed, bracket_allowed)
# Data examples to check current results
data_examples <- data |> inner_join(char_classes, by = c("left_charClass" = "charClass")) |>
						 select(left_charClass, right_charClass, left_symbClass, right_symbClass, chars) |>
						 rename(left_chars = chars) |>
						 inner_join(char_classes, by = c("right_charClass" = "charClass")) |>
						 select(left_charClass, right_charClass, left_symbClass, right_symbClass, left_chars, chars) |>
						 rename(right_chars = chars) |>
						 rowwise() |>
						 mutate(example_pair = str_c(str_sub(left_chars, 1, 1), str_sub(right_chars, 1, 1))) |>
						 select(left_charClass, right_charClass, example_pair)
# Export pairs
write_tsv(data, "C:/.../theory_pairs_charClasses_examples.tsv")