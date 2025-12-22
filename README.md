---
---
---

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

## Analysis of symbols and characters allowed in SMILES

### Atom symbols and corresponding characters

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two categories by their length:

-   Symbols consisting of the single character

-   Symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two categories by their grammatical requirements:

-   Symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   Symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on their nature and surroundings:

-   Symbols of the aromatic atoms

-   Symbols of the aliphatic atoms

##### Thus, the following type of atom symbols probably could be found in SMILES:

-   Single character atom symbols of organic aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could be designated as **watom_oar**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

-   Single character atom symbols of organic aliphatic atoms:

| B, C, N, O, S, P, F, I

Corresponding characters could be designated as **watom_oal**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

-   Single character atom symbols of bracket aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could be designated as **watom_bar**,where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic.

-   Single character atom symbols of bracket aliphatic atoms:

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **watom_bal**,where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic.

-   Two character atom symbols of organic aromatic atoms:

| Such symbols do not exist

-   Two character atom symbols of organic aliphatic atoms:

| Cl, Br

Corresponding characters could be designated as **satom_oal** & **eatom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic.

-   Two character atom symbols of bracket aromatic atoms:

| se, as, te

Corresponding characters could be designated as **satom_bar** & **eatom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

-   Two character atom symbols of bracket aliphatic atoms:

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **satom_bal** & **eatom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

Also, **\*** is an allowed symbol in SMILES, which corresponds to the any atom symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms:

-   Single character atom symbol of any atom:

| \*

Corresponding character could be designated as **anything**, since it could mean basically anything.

#### Question

So, "SMILES is parsable with a context-free parser", which in practice means that given the set of SMILES strings, each of them could be correctly parsed without considering any other string.

Maybe it is also possible to guess the atom and its properties expressed as the type of atom symbol in the single SMILES string by a single character, without considering any other character in this SMILES string?

To be honest, the answer is definitely not, but the details matter.

Trying to answer this question by checking the intersection between the sets of characters describing atom symbols of different types, it is possible to get some insights into what also should be considered, besides the one character at time, while parsing SMILES string.

This check up is quite easy to do using R (<https://cran.rstudio.org/bin/windows/>), Tidyverse (<https://tidyverse.org/>) and ggupset (<https://github.com/const-ae/ggupset>):

``` {#char_intersect__atomSymbols .R .illustration}
library(tidyverse)
library(ggupset)
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
# each charcter corresponds to the whole symbol
watom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
# each charcter corresponds to the whole symbol
watom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
# each charcter corresponds to the whole symbol
watom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
# each charcter corresponds to the whole symbol
watom_bal <- symb__atom_bal
#05, two character atom symbol of organic aromatic atoms
# symb_atom_oar IS EMPTY
# 06, two character atom symbol of organic aliphatic atoms
symb__atom_oal <- c('Cl', 'Br')
# Get the first character for each atom symbol
satom_oal <- map(symb__atom_oal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_oal <- map(symb__atom_oal, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 07, two character atom symbol of bracket aromatic atom
symb__atom_bar <- c('se', 'as', 'te')
# Get the first character for each atom symbol
satom_bar <- map(symb__atom_bar, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bar <- map(symb__atom_bar, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 08, two character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar',
                    'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb',
                    'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
                    'Te', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt',
                    'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
                    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl',
                    'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np',
                    'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
satom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 09, single character symbol of anything
anything <- c('*')

## Analyze characters constituting atom symbols from different classes
atom_symbol_characters <- tibble(
    character = c(watom_oar,
                  watom_oal,
                  watom_bar,
                  watom_bal,
                  satom_oal,
                  eatom_oal,
                  satom_bar,
                  eatom_bar,
                  satom_bal,
                  eatom_bal,
                  anything),
    class = c(rep("watom_oar", length(watom_oar)),
              rep("watom_oal", length(watom_oal)),
              rep("watom_bar", length(watom_bar)),
              rep("watom_bal", length(watom_bal)),
              rep("satom_oal", length(satom_oal)),
              rep("eatom_oal", length(eatom_oal)),
              rep("satom_bar", length(satom_bar)),
              rep("eatom_bar", length(eatom_bar)),
              rep("satom_bal", length(satom_bal)),
              rep("eatom_bal", length(eatom_bal)),
              rep("anything", length(anything)))
    ) |>
    group_by(character) |>
    summarise(class = list(class))
                    
# UpSet plot for classes
plot_class <- ggplot(atom_symbol_characters, aes(x=class)) +
                geom_bar() +
                scale_x_upset() +
                theme_minimal()
plot_class
```

#### Result

![**Figure 1.** Intersections of the character sets describing different atom symbol types.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_atomchars.png)

As it can be clearly seen from the **Figure 1:** yes, it will not be possible to guess the atom and its type operating single atom symbol character at time while parsing SMILES string, since there are intersections between the different sets of characters. However, on the positive side: intersections are not uniform:

| e.g.,
| satom_bal, watom_bal & watom_oal sets have the **largest** intersection, while intersection between eatom_bal, watom_bar, watom_oarm, satom_bar andeatom_bar sets is the **smallest**

Probably, this observations could be used latter to provide the future SMILES parser with some efficiency.

### Bond symbols and corresponding characters

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, each of which corresponds to the certain type of chemical bond:**

-   Single character bond symbol corresponding to the single bond:

| -

This single character symbol could be and often is omitted, since by default all the atoms, which symbols are written side by side in SMILES string are presumed to be connected by this type of bond, thus, this symbol is often omitted.

-   Single character bond symbol corresponding to the double bond:

| =

-   Single character bond symbol corresponding to the triple bond:

| \#

-   Single character symbol corresponding to the quadruple bond:

| \$

-   Single character symbol corresponding to the aromatic bond:

| :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications.

-   Single character symbol corresponding to the absence of the bond between the two specific atoms:

| .

As it was said earlier, by default each atom in SMILES string is considered to be connected with its immediate neighbors via the single bond (**-**). Thus, symbol corresponding to the negation of the bond is needed sometimes, and here it is: **.**

### Bond modifying (multiplying) symbols and corresponding characters

#### What are they?

Bond multiplying symbols is used in SMILES to extend the number of atoms, for which connections to the current atom could be written using linear notation (SMILES).

Bond multiplying symbols allowed in SMILES could be divided into two categories by their length:

-   Symbols consisting of the single character

-   Symbols consisting of more than one character

Bond multiplying symbols allowed in SMILES could be divided into two categories according to their role in completing the task:

-   Symbols initiators

-   Symbols terminators

Bond multiplying symbols allowed in SMILES could be divided into two categories by their task:

-   Symbols used to indicate simple additional bond for the current atom (branch)

-   Symbols used to indicate additional bond, which allows for the cycle (ring) to be formed

Bond multiplying symbols allowed in SMILES could be divided into two categories according to the quality of the additional bond:

-   Symbols explicitly including additional bond

-   Symbols implicitly including additional bond

##### Thus, the following types of bond multiplying symbols probably could be found in SMILES:

-   Single character bond multiplying symbols initiators of branching with explicit bond:

| Such symbols do not exist

-   Single character bond multiplying symbols initiators of branching with implicit bond:

| (

Corresponding characters could be designated as **wbm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols initiators of rings with explicit bond:

| Such symbols do not exist

-   Single character bond multiplying symbols initiators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wbm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols terminators of branching with explicit bond:

| Such symbols do not exist

-   Single character bond multiplying symbols terminators of branching with implicit bond:

| )

Corresponding characters could be designated as **wbm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols terminators of rings with explicit bond:

| Such symbols do not exist

-   Single character bond multiplying symbols terminators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wbm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

-   Multicharacter bond multiplying symbols initiators of branching with explicit bond:

| (-, (=, (#, (\$, (:, (.

Corresponding characters could be designated as **sbm_ibe** & **ebm_ibe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Multicharacter bond multiplying symbols initiators of branching with implicit bond:

| Such symbols do not exist

-   Multicharacter bond multiplying symbols initiators of rings with explicit bond:

| -[0-9], =[0-9], #[0-9], \$[0-9], :[0-9], .[0-9], -%[0-9][0-9], =%[0-9][0-9], #%[0-9][0-9], \$%[0-9][0-9], :%[0-9][0-9], .%[0-9][0-9]

Corresponding characters could be designated as **sbm_ire** & **ebm_ire**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and last suffix **e** - for explicit.

-   Multicharacter bond multiplying symbols initiators of rings with implicit bond:

| %[0-9][0-9]

Corresponding characters could be designated as **sbm_iri** & **ebm_iri**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and last suffix **i** - for implicit.

-   Multicharacter bond multiplying symbols terminators of branching with explicit bond:

| )-, )=, )#, )\$, ):, ).

Corresponding characters could be designated as **sbm_tbe** & **ebm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

-   Multicharacter bond multiplying symbols terminators of branching with implicit bond:

| Such symbols do not exist

-   Multicharacter bond multiplying symbols terminators of rings with explicit bond:

| -[0-9], =[0-9], #[0-9], \$[0-9], :[0-9], .[0-9], -%[0-9][0-9], =%[0-9][0-9], #%[0-9][0-9], \$%[0-9][0-9], :%[0-9][0-9], .%[0-9][0-9]

Corresponding characters could be designated as **sbm_tre** & **ebm_tre**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, and last suffix **e** - for explicit.

-   Multicharacter bond multiplying symbols terminators of rings with implicit bond:

| %[0-9][0-9]

Corresponding characters could be designated as **sbm_tri** & **ebm_tri**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, and last suffix i - for implicit.

#### Sort of the same question

From simple enumeration its quite clear that it will not possible to discriminate between the different classes of bond modifying (multiplying) symbols operating one character at time, also it will be not possible to discriminate between the some classes of bond multipliers and bonds while checking the first character of the symbol.

However, what are the intersections between the character sets describing different classes of bond multipliers?

Here is the code:

``` {#char_intersect__bmSymbols .R .illustration}
library(tidyverse)
library(ggupset)
# 02, single character bond multiplying symbols initiators of branch with implicit bond
wbm__ibi <- c('(')
# 03, single character bond multiplying symbols terminators of branch with implicit bond
wbm__tbi <- c(')')
# 04, single character bond multiplying symbols initiators of rings with implicit bond
wbm__iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
# 05, single character bond multiplying symbols terminators of rings with implicit bond
wbm__tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
# 06, multicharacter bond multiplying symbols initiators of branch with explicit bond 
sbm__ibe <- c('(')
ebm__ibe <- c('-', '=', '#', '$', ':', '.')
# 07, multicharacter bond multiplying symbols initiators of rings with explicit bond
sbm__ire  <- c('-', '=', '#', '$', ':', '.', '%')
ebm__ire  <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
# 08, multicharacter bond multiplying symbols initiators of rings with implicit bond
sbm__iri <- c('%')
ebm__iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
# 09, multicharacter bond multiplying symbols terminators of branch with explicit bond
sbm__tbe <- c(')')
ebm__tbe <- c('-', '=', '#', '$', ':', '.')
# 10, multicharacter bond multiplying symbols terminators of branch with explicit bond
sbm__tre <- c('-', '=', '#', '$', ':', '.')
ebm__tre <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
# 11, multicharacter bond multiplying symbols terminators of rings with explicit bond
sbm__tri <- c('%')
ebm__tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
## Analyze characters from different symbol classes
bm_symbol_characters <- tibble(character = c(wbm__ibi, wbm__tbi, wbm__iri, wbm__tri, sbm__ibe, ebm__ibe, sbm__ire, ebm__ire,
                                                sbm__iri, ebm__iri, sbm__tbe, ebm__tbe, sbm__tre, ebm__tre, sbm__tri, ebm__tri),
                                 class = c(rep("wbm__ibi", length(wbm__ibi)),
                                            rep("wbm__tbi", length(wbm__tbi)),
                                            rep("wbm__iri", length(wbm__iri)),
                                            rep("wbm__tri", length(wbm__tri)),
                                            rep("sbm__ibe", length(sbm__ibe)),
                                            rep("ebm__ibe", length(ebm__ibe)),
                                            rep("sbm__ire", length(sbm__ire)),
                                            rep("ebm__ire", length(ebm__ire)),
                                            rep("sbm__iri", length(sbm__iri)),
                                            rep("ebm__iri", length(ebm__iri)),
                                            rep("sbm__tbe", length(sbm__tbe)),
                                            rep("ebm__tbe", length(ebm__tbe)),
                                            rep("sbm__tre", length(sbm__tre)),
                                            rep("ebm__tre", length(ebm__tre)),
                                            rep("sbm__tri", length(sbm__tri)),
                                            rep("ebm__tri", length(ebm__tri)) )
                                 ) |>
                            group_by(character) |>
                            summarise(class = list(class))
# UpSet plot for classes
plot_class <- ggplot(bm_symbol_characters, aes(x=class)) +
                geom_bar() +
                scale_x_upset() +
                theme_minimal()
plot_class
```

![**Figure 2.** Intersections of the character sets describing different bond multiplying symbols types.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_bmchars.png)

As it can be clearly seen from the **Figure 2:** the sets of characters have large intersections in this case, it will be useful later while creating rules for parsing SMILES. At the moment these intersection suggest that the characters describing symbols from different classes were enumerated well: the largest intersection is between the set of characters describing end of the ring symbols or whole rings. As it should be ([0-9]).

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two categories of cis/trans symbols are allowed in SMILES:

-   Cis/trans symbols on the left side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **lct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

-   Cis/trans symbols on the right side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **rct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\F)=C/F are equivalent, since

> The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.
>
> *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****
>
> **\*** <http://opensmiles.org/opensmiles.html>

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
library(tidyverse)
library(RMariaDB)
library(DBI)
## Characters of cis/trans symbols
# 01.     char           Single character left side cis/trans              lct
# 02.     char           Single character right side cis/trans             rct
lct <- c('/', '\\')
rct <- c('/', '\\')
## Patterns to search for
# Parsing SMILES is a task, which may be harder than it seems on the first glance
# Thus, at this stage the regexps will be used, which allows to extact substring containing only the first pair of cis/trans symbols
# ^ - stands for the start of the string
# followed by
# [^\\\\/]* - means thath from the start of the string and up to the next meaningful part of regexp there should not be matches with the characters of cis/trans symbols
# The last meaningful (and variable) part of the regexps stands for the one of the variants of cis/trans symbols usage from http://opensmiles.org/opensmiles.html
# for example, in pattern_baseOne, [^\\(]/[Cc]=[Cc]/. matches F/C=C/F and does not match C(/F)=C/F
# for example, in pattern_hardOne, [Cc]\\(/.\\)=[Cc]/. matches C(/F)=C/F and does not match F/C=C/F
pattern_baseOne <- "^[^\\\\/]*[^\\(]/[Cc]=[Cc]/."
pattern_baseTwo <- "^[^\\\\/]*[^\\(]\\\\.[Cc]=[Cc]/."
pattern_hardOne <- "^[^\\\\/]*[Cc]\\(/.\\)=[Cc]/."
pattern_hardTwo <- "^[^\\\\/]*[Cc]\\(\\\\.\\)=[Cc]/."

## Connect to DB
mysql_password = '*****'
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
                ungroup() # 2 854 654 SMILES strings extracted

## Filter and Count
cs_smiles_matched <- cs_smiles_checked |>
                        filter(if_any(starts_with("pattern"), ~ !is.na(.)))
# 49 673  records matched one of the patterns
cs_smiles_matched_hard <- cs_smiles_matched |>
                        filter(if_any(starts_with("pattern_hard"), ~ !is.na(.)))
# 0       records matched one of the hard patterns
cs_smiles_matched_base <- cs_smiles_matched |>
                        filter(if_any(starts_with("pattern_base"), ~ !is.na(.)))
# 49 673  records matched one of the base patterns
```

The first pairs of cis/trans symbols in ChEMBL data (canonical smiles) do not include the hard cases similar to the ones described in <http://opensmiles.org/opensmiles.html>

Probably, it is safe to say that the percentage of such cases should be quite low at least in the curated databases.

From this, it is possible to assume that reliability of the SMILES as a form of representation of chemical structures comes not only from the basic rules of this language, but also from the standards of its usage adopted in the community.

Still, ability to parse SMILES using basic rules are essential to maintain this status.

### All the symbols and corresponding characters inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating speciffic isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

-   Single character isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wisotope**, where prefix **w** stands for the whole symbol.

-   Two-character isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sisotope & eisotope**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

-   Three-character isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sisotope, misotope, eisotope**, where prefix **s** stands for the start and prefix **m** stands for the middle and prefix **e** stands for the end of the symbol.

##### Chiral symbols

Chiral symbols are used to show that an atom is a stereocenter.

Chiral symbols allowed in SMILES could be divided into 5 categories by their length:

-   Single character chiral symbols:

| \@

Corresponding characters could be designated as **wchiral**, where prefix **w** stands for the whole symbol.

-   Two-character chiral symbols:

| \@

Corresponding characters could be designated as **schiral & echiral**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

-   Four-character chiral symbols:

| \@, T, H, A, L, S, P, B, O, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **schiral & mchiral & echiral**, where prefix **s** stands for the start, prefix **m** stands for the middle (2 characters) and prefix **e** stands for the end of the symbol.

-   Five-character chiral symbols:

| \@, T, B, O, H, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **schiral & mchiral & echiral**, where prefix **s** stands for the start, prefix **m** stands for the middle (3 characters) and prefix **e** stands for the end of the symbol.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 categories by their length:

-   Single character hydrogen symbols:

| H

Corresponding characters could be designated as **whydrogen**, where prefix **w** stands for the whole symbol.

-   Two-character hydrogen symbols:

| H, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **shydrogen & ehydrogen**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 categories by their length:

-   Single character charge symbols:

| +, -

Corresponding characters could be designated as **wcharge**, where prefix **w** stands for the whole symbol.

-   Two-character hydrogen symbols:

| +, -, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **scharge & echarge**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Class symbols

Class symbols designate the class of the atom, which is thing defined by the author of SMILES string.

Class symbols allowed in SMILES could be divided into 3 categories by their length:

-   Two-character class symbols:

| :, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sclass & eclass**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

-   Three-character class symbols:

| :, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sclass & mclass & eclass**, where prefix **s** stands for the start, **m** (1 character) stands for the middle and **e** stands for the end of the symbol.

-   Four-character class symbols:

| :, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sclass & mclass & eclass**, where prefix **s** stands for the start, **m** (2 characters) stands for the middle and **e** stands for the end of the symbol.

-   Five-character class symbols:

| :, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sclass & mclass & eclass**, where prefix **s** stands for the start, **m** (3 characters) stands for the middle and **e** stands for the end of the symbol.

As it can be seen without the further analysis, the aforementioned in-brackets categories of symbols are highly interconnected, but still it is quite easy to discriminate between the different symbols inside the square brackets, since they appear in fixed order.
