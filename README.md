# SMILES Parser (work in progress)

For the SMILES (Simplified Molecular Input Line Entry System) reference, please, SEE:

-   <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>

-   <http://opensmiles.org/opensmiles.html>

-   <https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System>

In short:

-   SMILES is a language, which is used to describe a molecule in the terms of and according to the rules of chemical valence and mathematical graph theories.

-   SMILES string is the linear notation of a spanning tree of a graph representing molecule.

-   SMILES string includes atom symbols, bonds between them and some characteristics of these entities.

-   Also, according to the Wikipedia, "From the view point of a formal language theory, SMILES is a word. A SMILES is parsable with a context-free parser".

-   To build such a parser it will be useful to understand, what are the symbols and corresponding characters allowed in SMILES and what are their allowed combinations.

-   Also, according to the Wikipedia, "In terms of a graph-based computational procedure, SMILES is a string obtained by printing the symbol nodes encountered in a depth-first tree traversal of a chemical graph. The chemical graph is first trimmed to remove hydrogen atoms and cycles are broken to turn it into a spanning tree."

-   Thus, via direct left to right parsing of the SMILES, it is possible to obtain readily useful representation of the chemical structure.

## Enumeration and classification of the symbols and characters allowed in SMILES

| *This part may require some further adjustments and corrections, but should be OK in general.*

| *Update from 19-02-2026*: using the older version of the materials in this repo it was possible to correctly classify each character in the SMILES from ChEMBL v36 in the terms of its possible belonging to the described classes. The following major flaw was identified: the classification of characters itself is not strict enough allowing some inconsistencies.
| Thus, more consistent approach is needed, terminology should be described.

The following categories (upper level of classification) of SMILES symbols (sequence of characters having specific meaning) could be enumerated:

1.  Atom symbols
2.  Symbol of any atom or basically anything
3.  Square bracket symbols
4.  Bond symbols
5.  Bond modifying (multiplying) symbols
6.  Cis/Trans symbols
7.  All the symbols inside the square brackets besides the main atom symbol (isotope symbols, chirality symbols, hydrogen symbols, charge symbols, atom class symbols)

Using faceted classification scheme (<https://en.wikipedia.org/wiki/Faceted_classification>) the following aspects meaningful for the parsing task could be used to describe symbols in SMILES further:

-   number of characters constituting the symbol

-   specific grammatical requirements (symbols of some atoms could only be valid when they are enclosed within the square brackets, symbols of other atoms do not require such enclosing)

-   aromaticity (for atoms only)

-   whether symbol marks the start or the end of something in case of symbols, which go in pairs (initiator (left), terminator (right))

-   whether symbol means the single additional bond (branch) or initiation / termination of the cycle (ring) - only for the bond multiplying symbols

-   whether symbol includes bonds explicitly or implicitly - only for the bond multiplying symbols

Using the information above, it is possible to

-   construct character classes describing symbols in SMILES or their parts providing some convenience for parsing

-   assess their intersections and frequency in the available data, which will be useful while selecting particular parsing approach

And then, select particular parsing approach and set of rules within it to hopefully finally come up with the pretty normal SMILES parser.

### Atom symbols category and corresponding character classes

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two facets by their length:

-   Symbols consisting of the single character

-   Symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two facets by their grammatical requirements:

-   Symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   Symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on the nature of their bonding:

-   Symbols of the aromatic atoms

-   Symbols of the aliphatic atoms

Thus, the following subcategories of atom symbols allowed in SMILES could be enumerated and labeled:

1.  Single character atom symbols of organic (from so called *organic* subset) aromatic atoms lacking the additional grammatical requirements and features:

| b, c, n, o, s, p

Corresponding characters could be designated as distinct character class, **w_atom_oar**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

2.  Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features:

| B, C, N, O, S, P, F, I

Corresponding characters could be designated as distinct character class, **w_atom_oal**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

3.  Single character atom symbols of aromatic atoms enclosed within brackets:

| b, c, n, o, s, p

Corresponding characters could be designated as **w_atom_bar**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic. As it can be seen, this class contains the same symbols as w_atom_oar, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

4.  Single character atom symbols of bracket aliphatic atoms:

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **w_atom_bal**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic. As it can be seen, this class contains the same symbols as w_atom_oal, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

5.  Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features:

| Cl, Br

Corresponding character classes could be designated as **s_atom_oal** & **e_atom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic. Further division of the characters describing symbol into two classes could be useful if the resulting parser will operate one character at time.

6.  Two character atom symbols of aromatic atoms, which should be enclosed within the square brackets:

| se, as, te

Corresponding characters could be designated as **s_atom_bar** & **e_atom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

7.  Two character atom symbols of aliphatic atoms, which should be enclosed within the square brackets:

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **s_atom_bal** & **e_atom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

Also, **\*** is an allowed symbol in SMILES, it corresponds to the any atom symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms:

8.  Single character symbol of any atom or basically anything:

| \*

Corresponding character could be designated as **anything**, since it could mean basically anything.

### Square bracket symbols and corresponding character classes

#### What are they?

9.  Square bracket symbols **[, ]** is the SMILES way to mark the start of the atom record including its various properties is the way to designate the end of the atom record.

| [, ]

Corresponding characters could be designated as **s_bracket & e_bracket** classes.

### Bond symbols and corresponding character classes

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, all of them are single character and do not have other peculiar aspects, five of them correspond to the conventional type of chemical bond:**

10. Single character bond symbol corresponding to the single bond:

| -

This single character symbol could be and typically is omitted, since by default all the atoms, which symbols are written side by side in SMILES string, are presumed to be connected by this type of bond. Corresponding character class will be designated as **single_bond**.

11. Single character bond symbol corresponding to the double bond:

| =

Corresponding character class will be designated as **double_bond**.

12. Single character bond symbol corresponding to the triple bond:

| \#

Corresponding character class will be designated as **triple_bond**.

13. Single character symbol corresponding to the quadruple bond:

| \$

Corresponding character class will be designated as **quadruple_bond**.

14. Single character symbol corresponding to the aromatic bond:

| :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications. Corresponding character class will be designated as **aromatic_bond_obsolete**.

15. Single character symbol corresponding to the absence of the bond between the two specific atoms:

| .

As it was said earlier, by default each atom in SMILES string is considered to be connected with its immediate neighbors via the single bond (**-**). Thus, symbol corresponding to the negation of the bond is needed sometimes, and here it is: **.** Corresponding character class will be designated as **no_bond**.

### Bond modifying (multiplying) symbols and corresponding characters

#### What are they?

Bond multiplying symbols is used in SMILES to extend the number of atoms, for which connections to the current atom could be written using linear notation (SMILES).

Bond multiplying symbols allowed in SMILES could be divided into four facets by their length:

-   Single character symbols

-   Two-character symbols

-   Three-character symbols

-   Four-character symbols

Bond multiplying symbols allowed in SMILES go in pairs and thus could be divided into two facets according to their role in completing the task of the symbols pair:

-   Symbols initiators

-   Symbols terminators

Bond multiplying symbols allowed in SMILES could be divided into two facets by their task:

-   Symbols used to indicate simple additional bond for the current atom (branch)

-   Symbols used to indicate additional bond, which allows for the cycle (ring) to be formed

Bond multiplying symbols allowed in SMILES could be divided into two categories according to whether additional bond is explicitly written or assumed:

-   Symbols explicitly including additional bond (any bond could be added)

-   Symbols implicitly including additional bond (only single bond could be added)

##### Thus, the following subcategories of bond multiplying symbols could be found in SMILES:

16. Single character bond multiplying symbols initiators of branching with implicit bond:

| (

Corresponding character class could be designated as **w_bm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

17. Single character bond multiplying symbols initiators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_bm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

18. Single character bond multiplying symbols terminators of branching with implicit bond:

| )

Corresponding characters could be designated as **w_bm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

19. Single character bond multiplying symbols terminators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **w_bm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

20. Two-character bond multiplying symbols initiators of branching with explicit bond:

| (-, (=, (#, (\$, (:, (.

Corresponding character classes could be designated as **s_bm_ibe** & **e_bm_ibe** where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, second suffix **e** - for explicit.

21. Two-character bond multiplying symbols initiators of rings with explicit bond:

| -0, =0, #0, \$0, :0, .0,
| -1, =1, #1, \$1, :1, .1,
| -2, =2, #2, \$2, :2, .2,
| -3, =3, #3, \$3, :3, .3,
| -4, =4, #4, \$4, :4, .4,
| -5, =5, #5, \$5, :5, .5,
| -6, =6, #6, \$6, :6, .6,
| -7, =7, #7, \$7, :7, .7,
| -8, =8, #8, \$8, :8, .8,
| -9, =9, #9, \$9, :9, .9

Corresponding characters could be designated as **s_bm_ire** & **e_bm_ire**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit.

22. Three-character bond multiplying symbols initiators of rings with implicit bond:

| %01, %02, %03, %04, %05, %06, %07, %08, %09, %10,
| %11, %12, %13, %14, %15, %16, %17, %18, %19, %20,
| %21, %22, %23, %24, %25, %26, %27, %28, %29, %30,
| %31, %32, %33, %34, %35, %36, %37, %38, %39, %40,
| %41, %42, %43, %44, %45, %46, %47, %48, %49, %50,
| %51, %52, %53, %54, %55, %56, %57, %58, %59, %60,
| %61, %62, %63, %64, %65, %66, %67, %68, %69, %70,
| %71, %72, %73, %74, %75, %76, %77, %78, %79, %80,
| %81, %82, %83, %84, %85, %86, %87, %88, %89, %90,
| %91, %92, %93, %94, %95, %96, %97, %98, %99

Corresponding character classes could be designated as **s_bm_iri** & **r_bm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, next suffix **i** - for implicit.

23. Four-character bond multiplying symbols initiators of rings with explicit bond:

| -%0[1-9], =%0[1-9], #%0[1-9], \$%0[1-9], :%0[1-9], .%0[1-9],
| -%[1-9][0-9], =%[1-9][0-9], #%[1-9][0-9], \$%[1-9][0-9], :%[1-9][0-9], .%[1-9][0-9]

Corresponding character classes could be designated as **s_bm_ire, n_bm_ire** & **r_bm_ire**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit.

24. Two-character bond multiplying symbols terminators of branching with explicit bond:

| )-, )=, )#, )\$, ):, ).

Corresponding character classes could be designated as **s_bm_tbe** & **e_bm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

25. Two-character bond multiplying symbols terminators of rings with explicit bond:

| -0, =0, #0, \$0, :0, .0,
| -1, =1, #1, \$1, :1, .1,
| -2, =2, #2, \$2, :2, .2,
| -3, =3, #3, \$3, :3, .3,
| -4, =4, #4, \$4, :4, .4,
| -5, =5, #5, \$5, :5, .5,
| -6, =6, #6, \$6, :6, .6,
| -7, =7, #7, \$7, :7, .7,
| -8, =8, #8, \$8, :8, .8,
| -9, =9, #9, \$9, :9, .9

Corresponding character classes could be designated as **s_bm_tre** & **e_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit.

26. Four-character bond multiplying symbols terminators of rings with explicit bond:

| -%0[1-9], =%0[1-9], #%0[1-9], \$%0[1-9], :%0[1-9], .%0[1-9],
| -%[1-9][0-9], =%[1-9][0-9], #%[1-9][0-9], \$%[1-9][0-9], :%[1-9][0-9], .%[1-9][0-9]

Corresponding characters could be designated as **s_bm_tre, n_bm_tre** & **r_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit.

27. Three-character bond multiplying symbols terminators of rings with implicit bond:

| %01, %02, %03, %04, %05, %06, %07, %08, %09, %10,
| %11, %12, %13, %14, %15, %16, %17, %18, %19, %20,
| %21, %22, %23, %24, %25, %26, %27, %28, %29, %30,
| %31, %32, %33, %34, %35, %36, %37, %38, %39, %40,
| %41, %42, %43, %44, %45, %46, %47, %48, %49, %50,
| %51, %52, %53, %54, %55, %56, %57, %58, %59, %60,
| %61, %62, %63, %64, %65, %66, %67, %68, %69, %70,
| %71, %72, %73, %74, %75, %76, %77, %78, %79, %80,
| %81, %82, %83, %84, %85, %86, %87, %88, %89, %90,
| %91, %92, %93, %94, %95, %96, %97, %98, %99

Corresponding characters could be designated as **s_bm_tri** & **r_bm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, next suffix **i** - for implicit.

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two facets of cis/trans symbols are allowed in SMILES:

28. Cis/trans single character symbols on the left side of the rotary non-permissive bond:

| /, \\

Corresponding single character characters could be designated as **l_ct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

29. Cis/trans symbols on the right side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **r_ct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\\F)=C/F are equivalent, since

| The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.
| *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****
| **\*** <http://opensmiles.org/opensmiles.html>

------------------------------------------------------------------------

#### Question

Just of curiosity, how many SMILES strings do contain cis/trans symbols, which could be misinterpreted by the *several SMILES interpreters*?

It is quite easy to approximate the answer to this question using

-   ChEMBL data (Zdrazil, Barbara. "Fifteen years of ChEMBL and its role in cheminformatics and drug discovery." *Journal of Cheminformatics* 17.1 (2025): 1-9.)

-   R ([https://cran.rstudio.org/bin/windows/](#0){style="font-size: 11pt;"})

-   Tidyverse ([https://tidyverse.org/](#0){style="font-size: 11pt;"})

-   DBI (<https://cran.r-project.org/web/packages/DBI>)

-   RMariaDB (<https://cran.r-project.org/web/packages/RMariaDB/index.html>)

Here is the code:

``` r
library(RMariaDB)
library(DBI)
library(tidyverse)

# 30, Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__lct <- c("/", "\\\\")
lct <- symb__lct
# 31, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__rct <- c("/", "\\\\")
rct <- symb__rct

## Patterns to search for
# Parsing SMILES is a task, which is harder than it may appears on the first glance
# Thus, at this stage the regexps will be used, which allows to extact substring containing only the first pair of cis/trans symbols
# ^ - stands for the start of the string
# [^\\\\/]* - means thath from the start of the string and up to the next meaningful part of regexp there should not be matches with the characters of cis/trans symbols
# The last meaningful (and variable) part of the regexps stands for the one of the variants of cis/trans symbols usage from http://opensmiles.org/opensmiles.html
# for example, in pattern_baseOne, [^\\(]/[^\\\\/]*=[^\\\\/]*[aA-zZ]/. matches F/C=C/F and does not match C(/F)=C/F
# for example, in pattern_hardOne, *C\\(/.\\)=C/. matches C(/F)=C/F and does not match F/C=C/F
pattern_baseOne <- "^[^\\\\/]*[^\\(]/[Cc]=[Cc]/."
pattern_baseTwo <- "^[^\\\\/]*[^\\(]\\\\.[Cc]=[Cc]/."
pattern_hardOne <- "^[^\\\\/]*[Cc]\\(/.\\)=[Cc]/."
pattern_hardTwo <- "^[^\\\\/]*[Cc]\\(\\\\.\\)=[Cc]/."
str_extract('C(/F)=C/F', pattern_baseOne)
str_extract('F/C=C/F', pattern_baseOne)
str_extract('C(/F)=C/F', pattern_hardOne)
str_extract('F/C=C/F', pattern_hardOne)
## Connect to DB
mysql_password = '***'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_36',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)
## Extract SMILES
cs__query <- dbSendQuery(con, 'SELECT canonical_smiles FROM compound_structures')
cs_smiles <- dbFetch(cs__query) |> distinct()
dbClearResult(cs__query)
# Close the connection
dbDisconnect(con)
## Check patterns against SMILES
cs_smiles_checked <- cs_smiles |> rowwise() |>
                    mutate(
                        pattern_baseOne = str_extract(canonical_smiles, pattern_baseOne),
                        pattern_baseTwo = str_extract(canonical_smiles, pattern_baseTwo),
                        pattern_hardOne = str_extract(canonical_smiles, pattern_hardOne),
                        pattern_hardTwo = str_extract(canonical_smiles, pattern_hardTwo)
                    ) |>
                    ungroup()
## Filter and Count
cs_smiles_matched <- cs_smiles_checked |> filter(if_any(starts_with("pattern"), ~ !is.na(.)))                       # 49 673  records matched one of the patterns
cs_smiles_matched_hard <- cs_smiles_matched |> filter(if_any(starts_with("pattern_hard"), ~ !is.na(.)))             # 0       records matched one of the hard patterns
cs_smiles_matched_base <- cs_smiles_matched |> filter(if_any(starts_with("pattern_base"), ~ !is.na(.)))             # 49 673  records matched one of the base patterns
```

The first pairs of cis/trans symbols in ChEMBL data (canonical smiles) do not include the hard cases similar to the ones described in <http://opensmiles.org/opensmiles.html>

Probably, it is safe to say that the percentage of such cases should be quite low at least in the curated databases.

From this, it is possible to assume that reliability of the SMILES as a form of representation of chemical structures comes not only from the basic rules of this language, but also from the standards of its usage adopted in the community.

Still, ability to parse SMILES using basic rules are essential to maintain this status.

------------------------------------------------------------------------

### All the symbols and corresponding character classes inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating speciffic isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

30. Single character isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_isotope**, where prefix **w** stands for the whole symbol.

31. Multicharacter (from 2 to 3 characters) isotope symbols:

| [0-9][0-9], [0-9][0-9][0-9]

Corresponding characters could be designated as **s_isotope_m & r_isotope_m**, where prefix **s** stands for the start and prefix **r** stands for the rest of the symbol and suffix **m** stands for the multi.

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

32. Single character chirality symbols:

| \@

Corresponding character class could be designated as **w_chiral**, where prefix **w** stands for the whole symbol.

33. Two-character chirality symbols:

| \@\@

Corresponding character classes could be designated as **s_chiral & e_chiral**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

34. Multi-character (four or five character) chirality symbols:

| \@, T, H, A, L, S, P, B, O, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character classes could be designated as **s_chiral_m, m_chiral_m & r_chiral_m**, where prefix **s** stands for the start, *prefix* **m** stands for the next from the medium (two characters), prefix **r** stands for the rest of the symbol and *suffix* **m** stands for the multi.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 facets by their length:

35. Single character hydrogen symbols:

| H

Corresponding character classes could be designated as **w_hydro**, where prefix **w** stands for the whole symbol.

36. Two-character hydrogen symbols:

| H0, H1, H2, H3, H4, H5, H6, H7, H8, H9

Corresponding character classes could be designated as **s_hydro & ehydro**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom.

Charge symbols allowed in SMILES could be divided into 2 facets by their length:

37. Single character charge symbols:

| +, -

Corresponding character classes could be designated as **w_charge**, where prefix **w** stands for the whole symbol.

38. Two-character charge obsolete symbols:

| ++, - -

Corresponding characters could be designated as **s_charge & e_charge**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

39. Multicharacter (two or three characters) charge symbols:

| +1, +2, +3, +4, +5, +6, +7, +8, +9, +10, +11, +12, +13, +14, +15,
|  -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15

Corresponding characters could be designated as **s_charge_m & r_charge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi.

##### Class symbols

Class symbols designate user-defined class of the atom.

Class symbols allowed in SMILES may have variable length, but there is no point to divide them into facets based on this aspect, so there is only:

40. Multicharacter (from 2 to 4 characters) class symbols:

| :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]

Corresponding character classes could be designated as **s_aclass & r_aclass**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

## Classification shortening

Following character classes describing chemically meaningful symbols occurring in SMILES could be listed:

1.  Character classes describing atom symbols:

    | e\_**atom**\_bal, e\_**atom**\_bar, e\_**atom**\_oal, s\_**atom**\_bal, s\_**atom**\_bar, s\_**atom**\_oal, w\_**atom**\_bal, w\_**atom**\_bar, w\_**atom**\_oal, w\_**atom**\_oar
    | 
    | Some of these original classes could be merged on condition that they are identical on the character level and having further convenience in mind:
    | 
    | w\_**atom**\_ar

2.  Character class describing any atom or basically anything:

    | **anything**

3.  Character classes describing bonds:

    | single\_**bond**, double\_**bond**, triple\_**bond**, quadruple\_**bond**, aromatic\_**bond**\_obsolete, no\_**bond**

4.  Character classes describing bond modifying (multiplying) symbols:

    | e\_**bm**\_ibe, e\_**bm**\_ire, e\_**bm**\_tbe, e\_**bm**\_tre, n\_**bm**\_ire, n\_**bm**\_tre, r\_**bm**\_ire, r\_**bm**\_iri, r\_**bm**\_tre, r\_**bm**\_tri, s\_**bm**\_ibe, s\_**bm**\_ire, s\_**bm**\_iri, s\_**bm**\_tbe, s\_**bm**\_tre, s\_**bm**\_tri, w\_**bm**\_ibi, w\_**bm**\_iri, w\_**bm**\_tbi, w\_**bm**\_tri
    | 
    | Some of these original classes could be merged on condition that they are identical on the character level and having further convenience in mind:
    | 
    | es\_**bm**\_e, s\_**bm**\_ri, n\_**bm**\_r, sw\_**bm**\_ib, sw\_**bm**\_tb, w\_**bm**\_ri, er\_**bm**\_r

5.  Character classes describing symbols of cis / trans features:

    | l\_**ct**, r\_**ct**
    | 
    | These original classes could be merged on condition that they are identical on the character level and having further convenience in mind:
    | 
    | **ct**

6.  Character classes describing symbols of square brackets:

    | **bracket**\_start, **bracket**\_end

7.  Character classes describing symbols of atomic features written inside the square brackets:

    | e\_**charge**, e\_**chiral**, e\_**hydro**, n\_**chiral**, r\_**aclass**, r\_**charge**, r\_**chiral**, r\_**isotope**, s\_**aclass**, s\_**charge**, s\_**chiral**, s\_**hydro**, s\_**isotope**, w\_**charge**, w\_**chiral**, w\_**hydro**
    | 
    | Some of these original classes could be merged on condition that they are identical on the character level and having further convenience in mind: 
    | 
    | esw\_**charge**

## Shortening's check up

Using the data provided in the *class_list.tsv* and *merged_classes\_\_autoDescribed.tsv*, it is possible to verify that characters from the merged character classes are identical to those from basic character classes (also it is possible to do the same thing using the code from *1_symbols&classes\_\_Rcode.R*, *2_classesMerged\_\_Rcode.R* and *3_checkMerged\_\_Rcode.R*).

## Intersection of character classes

As it can be seen without the further analysis, the aforementioned classes of symbols are highly interconnected via characters constituting them, even besides those, which are identical; the degree could be assessed numerically using previously introduced tools and visualized using UpSet plot from ggupset library (<https://cran.r-project.org/web/packages/ggupset/>).
