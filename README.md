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

-   Two-character bond multiplying symbols initiators of rings with explicit bond:

| -[0-9], =[0-9], #[0-9], \$[0-9], :[0-9], .[0-9]

Corresponding characters could be designated as **sbm_ire_2** & **ebm_ire_2**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **2** - for two-character.

-   Four-character bond multiplying symbols terminators of rings with explicit bond:

| -%[0-9][1-9], =%[0-9][1-9], #%[0-9][1-9], \$%[0-9][1-9], :%[0-9][1-9], .%[0-9][1-9]

Corresponding characters could be designated as **sbm_ire_4, mbm_ire_4** & **rbm_ire_4**,where prefix **s** stands for the start of the symbol, prefix **m** stands for the middle of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **4** - for four-character.

-   Multicharacter bond multiplying symbols initiators of rings with implicit bond:

| %[0-9][0-9]

Corresponding characters could be designated as **sbm_iri** & **rbm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and last suffix **i** - for implicit.

-   Multicharacter bond multiplying symbols terminators of branching with explicit bond:

| )-, )=, )#, )\$, ):, ).

Corresponding characters could be designated as **sbm_tbe** & **ebm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

-   Multicharacter bond multiplying symbols terminators of branching with implicit bond:

| Such symbols do not exist

-   Two-character bond multiplying symbols terminators of rings with explicit bond:

| -[0-9], =[0-9], #[0-9], \$[0-9], :[0-9], .[0-9]

Corresponding characters could be designated as **sbm_tre_2** & **ebm_tre_2**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **2** - for two-character.

-   Four-character bond multiplying symbols terminators of rings with explicit bond:

| -%[0-9][1-9], =%[0-9][1-9], #%[0-9][1-9], \$%[0-9][1-9], :%[0-9][1-9], .%[0-9][1-9]

Corresponding characters could be designated as **sbm_tre_4, mbm_tre_4** & **rbm_tre_4**,where prefix **s** stands for the start of the symbol, prefix **m** stands for the middle of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **4** - for four-character.

-   Multicharacter bond multiplying symbols terminators of rings with implicit bond:

| %[0-9][0-9]

Corresponding characters could be designated as **sbm_tri** & **rbm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, and last suffix i - for implicit.

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

-   Cis/trans single character symbols on the left side of the rotary non-permissive bond:

| /, \\

Corresponding single character characters could be designated as **lct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

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
lct <- c('/', '\\\\')
rct <- c('/', '\\\\')
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

-   Multicharacter (from 2 to 3 characters) isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sisotope_m & risotope_m**, where prefix **s** stands for the start and prefix r stands for the rest of the symbol and suffix **m** stands for the multi.

##### Chiral symbols

Chiral symbols are used to show that an atom is a stereocenter.

Chiral symbols allowed in SMILES could be divided into 5 categories by their length:

-   Single character chirality symbols:

| \@

Corresponding characters could be designated as **wchiral**, where prefix **w** stands for the whole symbol.

-   Two-character chirality symbols:

| \@

Corresponding characters could be designated as **schiral & echiral**, where prefix **e** stands for the start and prefix **r** stands for the end of the symbol.

-   Multicharacter chirality symbols:

| \@, T, H, A, L, S, P, B, O, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **schiral & rchiral**, where prefix **s** stands for the start, prefix **r** stands for the rest of the symbol.

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

-   Two-character charge obsolete symbols:

| +, -

Corresponding characters could be designated as **scharge_obsolete & echarge_obsolete**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

-   Multicharacter charge obsolete symbols:

| +, -, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **scharge_m & rcharge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi.

##### Class symbols

Class symbols designate the class of the atom, which is thing defined by the author of SMILES string.

Class symbols allowed in SMILES could be divided into 3 categories by their length:

-   Multicharacter (from 2 to 4 characters) class symbols:

| :, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **sclass & rclass**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

As it can be seen without the further analysis, the aforementioned in-brackets categories of symbols are highly interconnected, but still it is quite easy to discriminate between the different symbols inside the square brackets, since they appear in fixed order.

**This, *Analysis of symbols and characters allowed in SMILES*, will be updated to be in line with the further developments.**

## Pairs of characters in SMILES and in SMILES from ChEMBL (for example)

In short:

It is possible to read SMILES string character by character. In this case, each new character is the point where some questions should be answered and some decisions should be taken:

-   Is this character allowed in SMILES?

-   Is this character allowed after the previous character?

-   What does this character could mean in the current context (with the given previous characters)?

-   What are the questions, which should be addressed by the next character given this one and previous ones (in the updated context)?

Thus, it will be useful to prepare the list of all possible pairs of symbols, character classes and characters and check viability of their pairs in theory and in the real-world data to direct the further thinking.

Here is the code to enumerate the classes of characters described earlier:

``` r
library(RMariaDB)
library(DBI)
library(tidyverse)
library(ggupset)

### Enumeration of the character classes allowed in SMILES
## ATOMS
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
watom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
watom_bal <- symb__atom_bal
#06, two character atom symbol of organic aliphatic atoms
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
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
                          'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                          'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                          'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
                          'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
satom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()
## ANYTHING
# 09, single character symbol of anything
symb__anything <- c('*')
anything <- symb__anything
## BONDS
symb__single_bond <- c("-")
single_bond <- symb__single_bond
symb__double_bond <- c("=")
double_bond <- symb__double_bond
symb__triple_bond <- c("#")
triple_bond <- symb__triple_bond
symb__quadruple_bond <- c("$")
quadruple_bond <- symb__quadruple_bond
symb__aromatic_bond_obsolete <- c(":")
aromatic_bond_obsolete <- symb__aromatic_bond_obsolete
symb__no_bond <- c(".")
no_bond <- symb__no_bond
## BOND modifiers (multipliers)
# 10, single character bond multiplying symbols initiators of branch with implicit bond
symb__wbm_ibi <- c('(')
wbm_ibi <- symb__wbm_ibi
# 11, single character bond multiplying symbols terminators of branch with implicit bond
symb__wbm_tbi <- c(')')
wbm_tbi <- symb__wbm_tbi
# 12, single character bond multiplying symbols initiators of rings with implicit bond
symb__wbm_iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_iri <- symb__wbm_iri
# 13, single character bond multiplying symbols terminators of rings with implicit bond
symb__wbm_tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_tri <- symb__wbm_tri
# 14, multicharacter bond multiplying symbols initiators of branch with explicit bond
symb__bm_ibe <- c('(-', '(=', '(#', '($', '(:', '(.')
sbm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 15, two-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_2  <- c('-0', '=0', '#0', '$0', ':0', '.0',
					'-1', '=1', '#1', '$1', ':1', '.1',
					'-2', '=2', '#2', '$2', ':2', '.2',
					'-3', '=3', '#3', '$3', ':3', '.3',
					'-4', '=4', '#4', '$4', ':4', '.4',
					'-5', '=5', '#5', '$5', ':5', '.5',
					'-6', '=6', '#6', '$6', ':6', '.6',
					'-7', '=7', '#7', '$7', ':7', '.7',
					'-8', '=8', '#8', '$8', ':8', '.8',
					'-9', '=9', '#9', '$9', ':9', '.9' )
sbm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, 2)) |> unlist() |> unique()
# 16, four-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_4  <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
					'-%02', '=%02', '#%02', '$%02', ':%02', '.%02',
					'-%03', '=%03', '#%03', '$%03', ':%03', '.%03',
					'-%04', '=%04', '#%04', '$%04', ':%04', '.%04',
					'-%05', '=%05', '#%05', '$%05', ':%05', '.%05',
					'-%06', '=%06', '#%06', '$%06', ':%06', '.%06',
					'-%07', '=%07', '#%07', '$%07', ':%07', '.%07',
					'-%08', '=%08', '#%08', '$%08', ':%08', '.%08',
					'-%09', '=%09', '#%09', '$%09', ':%09', '.%09', 
					'-%10', '=%10', '#%10', '$%10', ':%10', '.%10',
					'-%11', '=%11', '#%11', '$%11', ':%11', '.%11',
					'-%12', '=%12', '#%12', '$%12', ':%12', '.%12',
					'-%13', '=%13', '#%13', '$%13', ':%13', '.%13',
					'-%14', '=%14', '#%14', '$%14', ':%14', '.%14',
					'-%15', '=%15', '#%15', '$%15', ':%15', '.%15',
					'-%16', '=%16', '#%16', '$%16', ':%16', '.%16',
					'-%17', '=%17', '#%17', '$%17', ':%17', '.%17',
					'-%18', '=%18', '#%18', '$%18', ':%18', '.%18',
					'-%19', '=%19', '#%19', '$%19', ':%19', '.%19',
					'-%20', '=%20', '#%20', '$%20', ':%20', '.%20',
					'-%21', '=%21', '#%21', '$%21', ':%21', '.%21',
					'-%22', '=%22', '#%22', '$%22', ':%22', '.%22',
					'-%23', '=%23', '#%23', '$%23', ':%23', '.%23',
					'-%24', '=%24', '#%24', '$%24', ':%24', '.%24',
					'-%25', '=%25', '#%25', '$%25', ':%25', '.%25',
					'-%26', '=%26', '#%26', '$%26', ':%26', '.%26',
					'-%27', '=%27', '#%27', '$%27', ':%27', '.%27',
					'-%28', '=%28', '#%28', '$%28', ':%28', '.%28',
					'-%29', '=%29', '#%29', '$%29', ':%29', '.%29',
					'-%30', '=%30', '#%30', '$%30', ':%30', '.%30',
					'-%31', '=%31', '#%31', '$%31', ':%31', '.%31',
					'-%32', '=%32', '#%32', '$%32', ':%32', '.%32',
					'-%33', '=%33', '#%33', '$%33', ':%33', '.%33',
					'-%34', '=%34', '#%34', '$%34', ':%34', '.%34',
					'-%35', '=%35', '#%35', '$%35', ':%35', '.%35',
					'-%36', '=%36', '#%36', '$%36', ':%36', '.%36',
					'-%37', '=%37', '#%37', '$%37', ':%37', '.%37',
					'-%38', '=%38', '#%38', '$%38', ':%38', '.%38',
					'-%39', '=%39', '#%39', '$%39', ':%39', '.%39',
					'-%40', '=%40', '#%40', '$%40', ':%40', '.%40',
					'-%41', '=%41', '#%41', '$%41', ':%41', '.%41',
					'-%42', '=%42', '#%42', '$%42', ':%42', '.%42',
					'-%43', '=%43', '#%43', '$%43', ':%43', '.%43',
					'-%44', '=%44', '#%44', '$%44', ':%44', '.%44',
					'-%45', '=%45', '#%45', '$%45', ':%45', '.%45',
					'-%46', '=%46', '#%46', '$%46', ':%46', '.%46',
					'-%47', '=%47', '#%47', '$%47', ':%47', '.%47',
					'-%48', '=%48', '#%48', '$%48', ':%48', '.%48',
					'-%49', '=%49', '#%49', '$%49', ':%49', '.%49',
					'-%50', '=%50', '#%50', '$%50', ':%50', '.%50',
					'-%51', '=%51', '#%51', '$%51', ':%51', '.%51',
					'-%52', '=%52', '#%52', '$%52', ':%52', '.%52',
					'-%53', '=%53', '#%53', '$%53', ':%53', '.%53',
					'-%54', '=%54', '#%54', '$%54', ':%54', '.%54',
					'-%55', '=%55', '#%55', '$%55', ':%55', '.%55',
					'-%56', '=%56', '#%56', '$%56', ':%56', '.%56',
					'-%57', '=%57', '#%57', '$%57', ':%57', '.%57',
					'-%58', '=%58', '#%58', '$%58', ':%58', '.%58',
					'-%59', '=%59', '#%59', '$%59', ':%59', '.%59',
					'-%60', '=%60', '#%60', '$%60', ':%60', '.%60',
					'-%61', '=%61', '#%61', '$%61', ':%61', '.%61',
					'-%62', '=%62', '#%62', '$%62', ':%62', '.%62',
					'-%63', '=%63', '#%63', '$%63', ':%63', '.%63',
					'-%64', '=%64', '#%64', '$%64', ':%64', '.%64',
					'-%65', '=%65', '#%65', '$%65', ':%65', '.%65',
					'-%66', '=%66', '#%66', '$%66', ':%66', '.%66',
					'-%67', '=%67', '#%67', '$%67', ':%67', '.%67',
					'-%68', '=%68', '#%68', '$%68', ':%68', '.%68',
					'-%69', '=%69', '#%69', '$%69', ':%69', '.%69',
					'-%70', '=%70', '#%70', '$%70', ':%70', '.%70',
					'-%71', '=%71', '#%71', '$%71', ':%71', '.%71',
					'-%72', '=%72', '#%72', '$%72', ':%72', '.%72',
					'-%73', '=%73', '#%73', '$%73', ':%73', '.%73',
					'-%74', '=%74', '#%74', '$%74', ':%74', '.%74',
					'-%75', '=%75', '#%75', '$%75', ':%75', '.%75',
					'-%76', '=%76', '#%76', '$%76', ':%76', '.%76',
					'-%77', '=%77', '#%77', '$%77', ':%77', '.%77',
					'-%78', '=%78', '#%78', '$%78', ':%78', '.%78',
					'-%79', '=%79', '#%79', '$%79', ':%79', '.%79',
					'-%80', '=%80', '#%80', '$%80', ':%80', '.%80',
					'-%81', '=%81', '#%81', '$%81', ':%81', '.%81',
					'-%82', '=%82', '#%82', '$%82', ':%82', '.%82',
					'-%83', '=%83', '#%83', '$%83', ':%83', '.%83',
					'-%84', '=%84', '#%84', '$%84', ':%84', '.%84',
					'-%85', '=%85', '#%85', '$%85', ':%85', '.%85',
					'-%86', '=%86', '#%86', '$%86', ':%86', '.%86',
					'-%87', '=%87', '#%87', '$%87', ':%87', '.%87',
					'-%88', '=%88', '#%88', '$%88', ':%88', '.%88',
					'-%89', '=%89', '#%89', '$%89', ':%89', '.%89',
					'-%90', '=%90', '#%90', '$%90', ':%90', '.%90',
					'-%91', '=%91', '#%91', '$%91', ':%91', '.%91',
					'-%92', '=%92', '#%92', '$%92', ':%92', '.%92',
					'-%93', '=%93', '#%93', '$%93', ':%93', '.%93',
					'-%94', '=%94', '#%94', '$%94', ':%94', '.%94',
					'-%95', '=%95', '#%95', '$%95', ':%95', '.%95',
					'-%96', '=%96', '#%96', '$%96', ':%96', '.%96',
					'-%97', '=%97', '#%97', '$%97', ':%97', '.%97',
					'-%98', '=%98', '#%98', '$%98', ':%98', '.%98',
					'-%99', '=%99', '#%99', '$%99', ':%99', '.%99')
sbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
mbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 3)) |> unlist() |> unique()
# 17, multicharacter bond multiplying symbols initiators of rings with implicit bond
symb__bm_iri <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
					'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
					'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
					'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
					'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
					'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
					'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
					'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
					'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' )
sbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique()
# 18, multicharacter bond multiplying symbols terminators of branch with explicit bond
symb__tbe <- c(')-', ')=', ')#', ')$', '):', ').')
sbm_tbe   <- map(symb__tbe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tbe   <- map(symb__tbe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 19, Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_2  <- symb__bm_ire_2
sbm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 20, Four-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_4  <- symb__bm_ire_4
sbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
mbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 3)) |> unlist() |> unique()
# 21, multicharacter bond multiplying symbols terminators of rings with implicit bond
symb__bm_tri <- symb__bm_iri
sbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique()
# 22, Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__lct <- c("/", "\\\\")
lct <- symb__lct
# 23, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__rct <- c("/", "\\\\")
rct <- symb__rct
## All the things inside the square brackets besides the main atom symbol
# 24, single character isotope symbols:
symb__isotope_1 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
wisotope <- symb__isotope_1
# 25, two-character isotope symbols:
symb__isotope_2 <- map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__isotope_3 <- map(symb__isotope_2, \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
sisotope_m <- map(symb__isotope_3, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
risotope_m <- map(symb__isotope_3, \(x) str_sub(x, 2)) |> unlist() |> unique()
# 26, single character chirality symbol:
symb__chiral_1 <- c('@')
wchiral <- symb__chiral_1
# 27, two-character chirality symbol:
symb__chiral_2 <- c('@@')
schiral_2 <- map(symb__chiral_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
echiral_2 <- map(symb__chiral_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 28, multicharacter chirality symbol:
symb__chir_m <- c("@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3", "@TB1", "@TB2", "@TB3", "@TB4",
					"@TB5", "@TB6", "@TB7", "@TB8", "@TB9", "@TB10", "@TB11", "@TB12", "@TB13", "@TB14", "@TB15", "@TB16", "@TB17", "@TB18", "@TB19", "@TB20",
					"@OH1", "@OH2", "@OH3", "@OH4", "@OH5", "@OH6", "@OH7", "@OH8", "@OH9", "@OH10", "@OH11", "@OH12", "@OH13", "@OH14", "@OH15", "@OH16", "@OH17", "@OH18", "@OH19", "@OH20",
					"@OH21", "@OH22", "@OH23", "@OH24", "@OH25", "@OH26", "@OH27", "@OH28", "@OH29", "@OH30")
schiral_m <- map(symb__chir_m, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rchiral_m <- map(symb__chir_m, \(x) str_sub(x, 2)) |> unlist() |> unique()
# 29, single character hydrogen symbol:
symb__hydrogen_1 <- c("H")
whydrogen <- symb__hydrogen_1
# 30, two-character hydrogen symbol:
symb__hydrogen_2 <- c("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9")
shydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ehydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 31, single character charge symbol:
symb__charge_1 <- c("+", "-")
wcharge <- symb__charge_1
# 32, two-character charge symbol, obsolete:
symb__charge_2 <- c("++", "--")
scharge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
echarge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 33, multicharacter charge symbol:
symb__class_2 <- map(c(":"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_3 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_4 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(symb__isotope_2, \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class <- c(symb__class_2, symb__class_3, symb__class_4)
sclass <- map(symb__class, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rclass <- map(symb__class, \(x) str_sub(x, 2)) |> unlist() |> unique()
```

Now it is possible to assess the intersection between the character sets describing different entities in SMILES (*a little deviation from the point of this part of the study and latter will be transferred to the appropriate position*):

-   Anything

-   Atoms

-   Bonds

-   Bond modifiers

-   Cis/Trans things

-   Atom properties written in square brackets

![**Figure 3.** Intersection between the classes of characters describing the different entities in SMILES.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_charClasses.png)
