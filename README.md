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

## General SMILES parsing strategy

The main idea is as follows:

1.  Computer program accepts the SMILES string, i.e., sequence of characters constituting symbols having chemical meaning.

2.  Computer program process this string from left to right one character at time.

    **What is meant by "computer program process"?**

    | - Computer program has default state.
    | - Every time computer program encounters new (next) character, state of the computer program changes accordingly (taking into account program's current state and what character it encounters).
    | - At each step computer program takes some action to produce an output.

3.  Computer program produces an output.

    **What is "output"?**

    | Data structure appropriate for the further computer processing and filled with the chemical data encoded by the input SMILES string.

    To get an insight into what kinds of state switching will be needed and possible for the program, it will be useful to check, which pairs of characters are possible in SMILES.

    | Topic of parsers is well developed, please, see the [Grune, D., & Jacobs, C. J. (2008). Introduction to parsing. In *Parsing techniques: A practical guide* (pp. 61-102). New York, NY: Springer New York.] and documentation and theory associated with the widely accepted parser generator software for the reference.

    To construct such a parser understanding of SMILES is needed, this understanding is about knowing symbols and characters, which could be present in the SMILES string; and rules of their arrangement. Thus, at the first stage all symbols and characters allowed in SMILES will be enumrated and classified.

## Enumeration and classification of the symbols and characters allowed in SMILES

| *This part may require some further adjustments and corrections, but should be OK in general.*

The following types (upper level of classification) of SMILES symbols (sequence of characters having specific meaning) could be enumerated:

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

-   construct symbol classes using meaningful combinations of facets mentioned above

-   construct corresponding classes of characters providing some convenience for parsing

-   assess their intersections and frequency in the available data, which will be useful while selecting particular parsing approach

And then, select particular parsing approach and set of rules within it and set of technologies for implementation to hopefully finally come up with the pretty normal SMILES parser.

### Atom symbols type and corresponding character classes

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two facets by their length:

-   symbols consisting of the single character

-   symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two facets by their grammatical requirements:

-   symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on the nature of their bonding:

-   symbols of the aromatic atoms

-   symbols of the aliphatic atoms

Thus, the following classes of atom symbols allowed in SMILES could be enumerated and labeled:

1.  Single character atom symbols of organic (from so called *organic* subset) aromatic atoms lacking the additional grammatical requirements and features (**atom_oar**):

| b, c, n, o, s, p

Corresponding characters could be designated as distinct character class, **w_atom_oar**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

2.  Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal**):

| B, C, N, O, S, P, F, I

Corresponding characters could be designated as distinct character class, **w_atom_oal**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

3.  Single character atom symbols of aromatic atoms enclosed within brackets (**atom_bar**):

| b, c, n, o, s, p

Corresponding characters could be designated as **w_atom_bar**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic. As it can be seen, this class contains the same symbols as w_atom_oar, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

4.  Single character atom symbols of bracket aliphatic atoms (**atom_bal**):

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **w_atom_bal**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic. As it can be seen, this class contains the same symbols as w_atom_oal, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

5.  Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal_2**):

| Cl, Br

Corresponding character classes could be designated as **s_atom_oal** & **e_atom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic. Further division of the characters describing symbol into two classes could be useful if the resulting parser will operate one character at time.

6.  Two character atom symbols of in-bracket aromatic atoms (**atom_bar_2**):

| se, as, te

Corresponding characters could be designated as **s_atom_bar** & **e_atom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

7.  Two character atom symbols of in-bracket aliphatic atoms (**atom_bal_2**):

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **s_atom_bal** & **e_atom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

### Symbol of anything and corresponding character class

**\*** is an allowed symbol in SMILES, it corresponds to any (atom) symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms :

8.  Single character symbol of any atom or basically **anything**:

| \*

Corresponding character could be designated as **w_any**, since it could mean basically any character or symbol.

### Square bracket symbols and corresponding character classes

#### What are they?

9.  Square bracket symbols **[, ]** is the SMILES way to mark the start of the atom record including its various properties and is the way to designate the end of the atom record.

| [, ]

Corresponding character classes could be designated as **s_bracket** & **e_bracket**.

### Bond symbols and corresponding character classes

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, all of them are single character and do not have other peculiar aspects, five of them correspond to the conventional type of chemical bond:**

10. Single character bond symbol corresponding to the single bond (**single_bond**):

| -

This single character symbol could be and typically is omitted, since by default all the atoms, which symbols are written side by side in SMILES string, are presumed to be connected by this type of bond. Corresponding character class will be designated as **w_single_bond**.

11. Single character bond symbol corresponding to the double bond (**double_bond**):

| =

Corresponding character class will be designated as **w_double_bond**.

12. Single character bond symbol corresponding to the triple bond (**triple_bond**):

| \#

Corresponding character class will be designated as **w_triple_bond**.

13. Single character symbol corresponding to the quadruple bond (**quadruple_bond**):

| \$

Corresponding character class will be designated as **w_quadruple_bond**.

14. Single character symbol corresponding to the aromatic bond (**aromatic_bond**):

| :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications. Corresponding character class will be designated as **w_aromatic_bond**.

15. Single character symbol corresponding to the absence of the bond between the two specific atoms (**no_bond**):

| .

### Bond modifying (multiplying) symbols and corresponding characters

#### What are they?

Bond modifying (multiplying) symbols, i.e. **modifiers**, are used in SMILES to extend the number of atoms, for which connections to the current atom could be written using linear notation (SMILES).

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

##### Thus, the following classes of bond multiplying symbols could be found in SMILES:

16. Single character bond multiplying symbols initiators of branching with implicit bond (**bm_ibi**):

| (

Corresponding character class could be designated as **w_bm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

17. Single character bond multiplying symbols initiators of rings with implicit bond (**bm_iri**):

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_bm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

18. Single character bond multiplying symbols terminators of branching with implicit bond (**bm_tbi**):

| )

Corresponding characters could be designated as **w_bm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

19. Single character bond multiplying symbols terminators of rings with implicit bond (**bm_tri**):

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **w_bm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

20. Two-character bond multiplying symbols initiators of branching with explicit bond (**bm_ibe**):

| ([-=#\$:.]

Corresponding character classes could be designated as **s_bm_ibe** & **e_bm_ibe** where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, second suffix **e** - for explicit.

21. Two-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_2**):

| [-=#\$:.][0-9]

Corresponding characters could be designated as **s_bm_ire** & **e_bm_ire**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit.

22. Three-character bond multiplying symbols initiators of rings with implicit bond (**bm_iri_3**):

| %[0-9][1-9], %[1-9][0-9]

Corresponding character classes could be designated as **s_bm_iri** & **r_bm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, next suffix **i** - for implicit.

23. Four-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_4**):

| [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding character classes could be designated as **s_bm_ire_4, n_bm_ire** & **r_bm_ire**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit; **4** in s_bm_ire_4 is used to differentiate this class from the s_bm_ire corresponding to the bm_ire_2, n_bm_ire and r_bm_ire are already unique.

24. Two-character bond multiplying symbols terminators of branching with explicit bond (**bm_tbe_2**):

| )[-=#\$:.]

Corresponding character classes could be designated as **s_bm_tbe** & **e_bm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

25. Two-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_2**):

|  [-=#\$:.][0-9]

Corresponding character classes could be designated as **s_bm_tre** & **e_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit.

26. Four-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_4**):

| [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding characters could be designated as **s_bm_tre, n_bm_tre** & **r_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit. **4** in s_bm_tre_4 is used to differentiate this class from the s_bm_tre corresponding to the bm_tre_2, n_bm_tre and r_bm_tre are already unique.

27. Three-character bond multiplying symbols terminators of rings with implicit bond (**bm_tri_3**):

| %[0-9][1-9], %[1-9][0-9]

Corresponding characters could be designated as **s_bm_tri** & **r_bm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, next suffix **i** - for implicit.

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two facets of cis/trans symbols are allowed in SMILES:

28. Cis/trans single character symbols on the left side of the rotary non-permissive bond (**l_ct**):

| /, \\

Corresponding single character characters could be designated as **l_ct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

29. Cis/trans symbols on the right side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **r_ct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\\F)=C/F are equivalent, since

| The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.
| *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****
| **\*** <http://opensmiles.org/opensmiles.html>

### All the symbols and corresponding character classes inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating specific isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

30. Single character isotope symbols (**isotope**):

| 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_isotope**, where prefix **w** stands for the whole symbol.

31. Multicharacter (from 2 to 3 characters) isotope symbols (**isotope_m**):

| [0-9][1-9], [1-9][0-9], [0-9][0-9][1-9], [0-9][1-9][0-9], [1-9][0-9][0-9]

Corresponding characters could be designated as **s_isotope & r_isotope**, where prefix **s** stands for the start and prefix **r** stands for the rest of the symbol.

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

32. Single character chirality symbol (**chiral**):

| \@

Corresponding character class could be designated as **w_chiral**, where prefix **w** stands for the whole symbol.

33. Two-character chirality symbols (**chiral_2**):

| [\@][\@]

Corresponding character classes could be designated as **s_chiral & e_chiral**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

34. Multicharacter (four or five character) chirality symbols (**chiral_m**):

| [\@]TH[1-2], [\@]AL[1-2], [\@]SP[1-3], [\@]TB[1-20], [\@]OH[1-30]

Corresponding character classes could be designated as **s_chiral_m, m_chiral & r_chiral**, where prefix **s** stands for the start, *prefix* **m** stands for the medium (two characters), prefix **r** stands for the rest of the symbol and *suffix* **m** stands for the multi, where it is needed.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 facets by their length:

35. Single character hydrogen symbol (**hydro**):

| H

Corresponding character classes could be designated as **w_hydro**, where prefix **w** stands for the whole symbol.

36. Two-character hydrogen symbols (**hydro_2**):

| H[2-9]

Corresponding character classes could be designated as **s_hydro & ehydro**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom (**charge**).

Charge symbols allowed in SMILES could be divided into 2 facets by their length:

37. Single character charge symbols (**charge**):

| [+-]

Corresponding character classes could be designated as **w_charge**, where prefix **w** stands for the whole symbol.

38. Two-character charge obsolete symbols (**charge_2**):

| [+][+], [-][-]

Corresponding characters could be designated as **s_charge & e_charge**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

39. Multicharacter (two or three characters) charge symbols (**charge_m**):

| [+-][1-9], [+-]1[0-5]

Corresponding characters could be designated as **s_charge_m & r_charge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi where it is needed.

##### Class symbols

Class symbols designate user-defined class of the atom.

Class symbols allowed in SMILES may have variable length, but there is no point to divide them into facets based on this aspect, so there is only:

40. Multicharacter (from 2 to 4 characters) class symbols (**class**):

| :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]

Corresponding character classes could be designated as **s_class & r_class**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

Information on symbols is summarized in **symbols.tsv**, all symbols are provided in **symbols_all.tsv**.

## Q1: are described symbols unique, i.e. is it possible to identify each SMILES symbols based only on characters constituting it?

No, as it can be seen from **symbols.tsv**. For example,

-   terminators and initiators of branching could be and often are the same by design

-   symbols of in-bracket atoms and bracket-free atoms could be the same (organic subset)

Also, several classes of symbols could be described or are described partially by the patterns [0-9] and [1-9], which makes parsing without consideration of the environment questionable.

In the previous version the attempt was taken to divide the whole symbols into the smaller subsets of characters, the result is as follows: this strategy does not pay off, character classes probably could be useful to construct the SMILES strings, but they provide no clear benefits for parsing:

> [!NOTE]
>
> It seems to be easier (and fast enough) to read the characters one by one until the longest possible sequence describing symbol is gathered (5 characters, I guess) and decide on the actual symbol afterwards, considering matches in **current sub-string** and **state** deduced from the previous symbols and length of the remaining SMILES string.
> 

## General SMILES parsing strategy revised

1.  Computer program initializes with the SMILES string having default **state** and empty **accumulator** of characters and empty **result**.

2.  Every time computer program encounters new (next) character it accumulates this character.

3.  Every time accumulator reaches its limits (longest symbol in SMILES or end of the string) its content is being evaluated, state changes accordingly, accumulator's content is being trimmed from left to right to delete all the characters evaluated as the whole symbol at this step.

4.  Every time state changes computer program takes some action to build up an output.

5.  When the end of the string is reached, computer program outputs the result.

> [!NOTE]
>
> The following text will be rewritten accordingly.
>



## General SMILES parsing strategy

The main idea is as follows:

1.  Computer program accepts the SMILES string, i.e., an ordered sequence of characters having chemical meaning as described previously.

2.  Computer program process this string from left to right one character at time.

    **What is meant by the "computer program process"?**

    | - Computer program has default state.
    | - Every time computer program encounters new (next) character, state of the computer program changes accordingly (taking into account program's current state and what character it encounters).
    | - At each step computer program takes some action to produce an output.

3.  Computer program produces an output.

    **What is the "output"?**

    | Data structure appropriate for the further computer processing and filled with the chemical data encoded by the input SMILES string.

    To get an insight into what kinds of state switching will be needed and possible for the program, it will be useful to check, which pairs of characters are possible in SMILES.

    | Topic of parsers is well developed, please, see the [Grune, D., & Jacobs, C. J. (2008). Introduction to parsing. In *Parsing techniques: A practical guide* (pp. 61-102). New York, NY: Springer New York.] for example.

## Pairs of SMILES characters, which are theoretically possible

**71** distinct characters were previously enumerated in SMILES, while parsing computer program can encounter only these characters:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -, #, \$, %, (, ), \*, ., /, :, \@, [, \\, ], +, =, a, A, b, B, c, C, d, D, e, E, f, F, g, G, h, H, i, I, k, K, l, L, m, M, n, N, o, O, p, P, r, R, s, S, t, T, u, U, v, V, W, X, y, Y, Z

Thus, it its quite possible to generate all the possible pairs of these characters to evaluate pairs' validity:

```         
library(tidyverse)
# Input chars, classes and types
data <- read_tsv("C:/.../chars_&_symbs.tsv") |>
                select(chars, charClass, symbClass, symbType) |>
                separate_longer_delim(chars, delim = ", ")
# Get unique characters
unique_chars <- data |> pull(chars) |> unique()
# Generate all the pairs possible in theory
pairs <- expand.grid(left_char = unique_chars, right_char = unique_chars)
```

**5041** distinct pairs are theoretically possible on the level of characters.

However, as it can be seen from Figures 1 and 2, many characters allowed in SMILES can belong to the several classes and types.

Thus, number of pairs is effectively larger considering classes and types:

```         
# Get all the theoretically possible pairs considering classes and types
pairs_labeled <- pairs |>
  inner_join(data, by = c("left_char" = "chars"), relationship = "many-to-many") |>
  rename(left_charClass = charClass, left_symbClass = symbClass, left_symbType = symbType) |>
  inner_join(data, by = c("right_char" = "chars"), relationship = "many-to-many") |>
  rename(right_charClass = charClass, right_symbClass = symbClass, right_symbType = symbType)
```

**99225** pairs of characters are theoretically possible in SMILES considering classes and types.

This number is huge.

# The following is the preliminary version, updates may be needed.

## Pairs of SMILES characters, which are allowed

Given the large number of the theoretically possible character to character transitions in SMILES, it should be reasonable to assess their viability, considering hierarchical nature of the labels.

## Pairs of symbol types, which are allowed in SMILES

```         
# Get the distinct theoretically possible pairs of symbol types
pairs_symbType <- pairs_labeled |> select(left_symbType, right_symbType) |> distinct()
```

**49** pairs of symbol types are possible in SMILES in theory, there is a need to check their viability using available community-driven SMILES specification and SMILES Parser Demo (<https://doc.gdb.tools/smilesDrawer/sd/example/index_light.html>) for testing:

| left_symbType       | right_symbType      | Allowed |
|---------------------|---------------------|---------|
| atom\_\_symbol      | atom\_\_symbol      | yes     |
| features\_\_symbol  | atom\_\_symbol      | yes     |
| anything\_\_symbol  | atom\_\_symbol      | yes     |
| bracket\_\_symbol   | atom\_\_symbol      | yes     |
| bond\_\_symbol      | atom\_\_symbol      | yes     |
| modifier\_\_symbol  | atom\_\_symbol      | yes     |
| cis_trans\_\_symbol | atom\_\_symbol      | yes     |
| atom\_\_symbol      | features\_\_symbol  | yes     |
| features\_\_symbol  | features\_\_symbol  | yes     |
| anything\_\_symbol  | features\_\_symbol  | yes     |
| bracket\_\_symbol   | features\_\_symbol  | yes     |
| bond\_\_symbol      | features\_\_symbol  | no      |
| modifier\_\_symbol  | features\_\_symbol  | no      |
| cis_trans\_\_symbol | features\_\_symbol  | no      |
| atom\_\_symbol      | anything\_\_symbol  | yes     |
| features\_\_symbol  | anything\_\_symbol  | yes     |
| anything\_\_symbol  | anything\_\_symbol  | yes     |
| bracket\_\_symbol   | anything\_\_symbol  | yes     |
| bond\_\_symbol      | anything\_\_symbol  | yes     |
| modifier\_\_symbol  | anything\_\_symbol  | yes     |
| cis_trans\_\_symbol | anything\_\_symbol  | yes     |
| atom\_\_symbol      | bracket\_\_symbol   | yes     |
| features\_\_symbol  | bracket\_\_symbol   | yes     |
| anything\_\_symbol  | bracket\_\_symbol   | yes     |
| bracket\_\_symbol   | bracket\_\_symbol   | no      |
| bond\_\_symbol      | bracket\_\_symbol   | yes     |
| modifier\_\_symbol  | bracket\_\_symbol   | yes     |
| cis_trans\_\_symbol | bracket\_\_symbol   | yes     |
| atom\_\_symbol      | bond\_\_symbol      | yes     |
| atom\_\_symbol      | modifier\_\_symbol  | yes     |
| features\_\_symbol  | bond\_\_symbol      | no      |
| features\_\_symbol  | modifier\_\_symbol  | no      |
| anything\_\_symbol  | bond\_\_symbol      | yes     |
| anything\_\_symbol  | modifier\_\_symbol  | yes     |
| bracket\_\_symbol   | bond\_\_symbol      | yes     |
| bracket\_\_symbol   | modifier\_\_symbol  | yes     |
| bond\_\_symbol      | bond\_\_symbol      | no      |
| bond\_\_symbol      | modifier\_\_symbol  | no      |
| modifier\_\_symbol  | bond\_\_symbol      | no      |
| modifier\_\_symbol  | modifier\_\_symbol  | yes     |
| cis_trans\_\_symbol | bond\_\_symbol      | no      |
| cis_trans\_\_symbol | modifier\_\_symbol  | yes     |
| atom\_\_symbol      | cis_trans\_\_symbol | yes     |
| features\_\_symbol  | cis_trans\_\_symbol | no      |
| anything\_\_symbol  | cis_trans\_\_symbol | yes     |
| bracket\_\_symbol   | cis_trans\_\_symbol | yes     |
| bond\_\_symbol      | cis_trans\_\_symbol | no      |
| modifier\_\_symbol  | cis_trans\_\_symbol | yes     |
| cis_trans\_\_symbol | cis_trans\_\_symbol | no      |

**Table 2.** List of symbol type to symbol type transitions, which are theoretically possible in SMILES.

From Table 2 it follows that among **49** symbol type to symbol type transitions theoretically possible in SMILES, there are **13** transitions, which are not allowed by the rules of the SMILES language as described in OpenSMILES documentation:

1.  bond\_\_symbol to features\_\_symbol

    | these symbols should always be separated by the bracket_symbol

2.  modifier\_\_symbol to features\_\_symbol

    | these symbols should always be separated by the bracket_symbol

3.  cis_trans\_\_symbol to features\_\_symbol

    | these symbols should always be separated by the bracket_symbol

4.  bracket\_\_symbol to bracket\_\_symbol

    | these symbols should always be separated by the atom_symbol

5.  features\_\_symbol to bond\_\_symbol

    | these symbols should always be separated by the bracket_symbol

6.  features\_\_symbol to modifier\_\_symbol

    | these symbols should always be separated by the bracket_symbol

7.  bond\_\_symbol to bond\_\_symbol

    | these symbols should always be separated by the atom_symbol

8.  bond\_\_symbol to modifier\_\_symbol

    | these symbols should always be separated by the bracket_symbol

9.  modifier\_\_symbol to bond\_\_symbol

    | these symbols should always be separated by the bracket_symbol

10. cis_trans\_\_symbol to bond\_\_symbol

    | these symbols should always be separated by the atom_symbol

11. features\_\_symbol to cis_trans\_\_symbol

    | these symbols should always be separated by the bracket_symbol

12. bond\_\_symbol to cis_trans\_\_symbol

    | these symbols should always be separated by the bracket_symbol

13. cis_trans\_\_symbol to cis_trans\_\_symbol

    | these symbols should always be separated by the atom_symbol

Thus, corresponding pairs of characters could be safely excluded from the further work:

```         
pairs_symbTypes_not <- read_tsv("C:/.../theory_pairs_symbType.tsv") |> filter(allowed == "no")
pairs_labeled <- pairs_labeled |> anti_join(pairs_symbTypes_not)
```

At this stage, **72709** pairs of characters are theoretically possible in SMILES considering classes, types and viability of symbol types pairing.

## Pairs of symbol classes, which are allowed in SMILES

[Here is the preliminary version, updates may be needed]{style="color:#8A350C"}.

Using the remaining set of character pairs, it is possible to do the same thing (assess viability) for the symbol classes:

```         
# Get the distinct theoretically possible pairs of symbol types
pairs_symbClass <- pairs_labeled |> select(left_symbClass, right_symbClass) |> distinct()
# Assess the number of theoretically possible symbol pairs
nrow(pairs_symbClass)
```

**909** pairs of symbol classes are possible in SMILES at this stage, there is a possibility to check their viability using available community-driven SMILES specification and SMILES Parser Demo (<https://doc.gdb.tools/smilesDrawer/sd/example/index_light.html>) for testing, however, the number of records is large this time, thus, it will be beneficial to filter out some of the not allowed pairings using the follwoing rules on the level of whole symbols:

1.  Aromatic bond, which is already deprecated, should not be paired with the aliphatic atoms.
2.  Isotope symbols could only be paired with the bracket on the left and bracket atom on the right.
3.  Features symbols, besides isotope, being on the left side could only be paired with the other feature symbols in the following order (from left to right): chiral -\> hydro -\> charge -\> class (gaps are allowed) and bracket.
4.  Features symbols, besides isotope, could only be paired with the other feature symbols or bracket symbols on the right and with the other feature symbols (in reverse order from 3) or bracket atom symbols on the left.
5.  Bracket atoms cannot be paired with the symbols contained outside the brackets and they can only be paired with the isotope symbols if those symbols are on the left and they can only be paired with the other bracket atom symbols on condition that those symbols has the same length, which is greater than 1.
6.  Organic atom symbols can only be paired with the symbols found outside the brackets.
7.  Bond symbols and bond modifying symbols can not be paired with the symbols contained inside the brackets.
8.  Bond modifying initiators and terminators of branching should be separated by at least one atom.
9.  Two character bond modifying initiator and terminator symbols could not precede other bond modifying initiator symbols.

Also, characters could be paired with the characters from the same symbol classes given that they belong to the one class, which members consist of more than one character considering the right order.

After applying this set of simple rules (with the latest update), **462** pairs of symbol classes are left, the following check-ups will be conducted on the level of character classes.

> [!NOTE]

> Toy parser shown that the decision making process concerning the characters of the current symbol mainly depends on the previous symbol.

> Thus, it will be alright to analyze one possible pair of symbols at time describing the appropriate course of action.

> After that the appropriate actions for each state possible will be gathered to describe the whole process.

## Checking pairs of symbol classes

It seems to be reasonable to check the remaining pairs of symbol classes as follows:

1. Simulate the parsing procedure for each case described by the pair of symbol characters, considering that pair is already known and the state is appropriate.

> [!NOTE]
> **Here is an idea of much simpler parsing strategy**:
> Read the SMILES string by character
> Accumulate characters of the current symbol
> At each step of accumulation search for the accumulated string in the list of the symbols allowed in SMILES
> Terminate the search and move to accumulating next symbol, when the longest symbol possible is exceeded by the accumulated string.

> [!NOTE]
> **Considering the previous note, additional symbol classes should be forbidden: pairs of identical classes, which were previously allowed because of their multicharacter nature.

## Checking pairs of symbol classes

It seems to be reasonable to check the remaining pairs of symbol classes as follows:

1. Simulate the parsing procedure for each case described by the pair of symbol characters, considering that pair is already known and the state is appropriate.
2. Check the results.
3. Analyze the results.

> [!NOTE]
> **Here is an idea of much simpler parsing procedure**:  
> - Read the SMILES string character by character.  
> - Accumulate current symbol.  
> - At each step of accumulation search for the accumulated string in the list of the symbols allowed in SMILES.  
> - Terminate the search and move to accumulating next symbol with regards to the last match, when the longest symbol possible is exceeded by the accumulated string.

So, it is time to build toy parser V2, which will parse first two symbols in the triples of symbols allowed in SMILES using this new cool strategy.

**This is much simpler than the current approach and should work. Also, it is a fully verifiable approach, since for each symbol pair it is possible to append each symbol allowed in SMILES and just check for mistakes.**

## Current status and further objectives

So, at this point symbols and characters allowed in SMILES are enumerated and classified, their relations are somewhat described, general parsing strategy is proposed.

Thus, it is possible to

-   design the data structure to store the parsed results

-   create basic realization of the parser

After that it will be possible to adjust this realization using

-   stats already obtained

-   stats, which could be gathered using publicly available sources

## Data structure to store the results of SMILES parsing

Probably, any in-computer representation of something could be characterized by the following aspects:

-   Human readability

-   Time and space (computational and storing) efficiency

-   Interchangeability

First two aspects are often in conflict, for example: one can read string and can not read the raw binary data, but human readable string usually takes more space to be stored and more time for processing than corresponding piece of the binary data.

Interchangeability here is defined as the ability for the data structure to be submitted elsewhere without further processing or with minimal processing. Interchangeability rather depends on the existing conventions, established practices.

Given that the basic realization of the parser is going to be designed and developed, it is safe to say that human readability and interchangeability are of priority.

From this, JSON (JavaScript Object Notation, <https://www.json.org/json-en.html>) seems to be a good choice: text format, which is completely language independent.

Substantially, as it was mentioned earlier, SMILES string contains information on atoms and bonds between them. Atoms and bonds could be used to describe the chemical structure further. Thus, it will be useful to have distinct, but related, records for them.

Thus, JSON object containing records on atoms and bonds between them is the desired output of the SMILES parser.

### Description of atoms

As it became clear from the OpenSMILES documentation and RSF 23-73-01058 (<https://github.com/RSF-23-73-01058/test_get-generic-scaffold-from-smiles>) the following parameters of atoms seem to be important to generate / parse and store:

-   Atom ID

-   Symbol

-   is organic

-   is aromatic

-   is in brackets

-   isotope

-   chirality

-   number of hydrogens

-   charge

-   class

### Description of bonds

As it became clear from the OpenSMILES documentation and RSF 23-73-01058 (<https://github.com/RSF-23-73-01058/test_get-generic-scaffold-from-smiles>) the following parameters of bonds seem to be important to generate / parse and store:

-  Bond ID

-  Atom_1 ID

-  Atom_2 ID

-  Type

-  is ring initiator/terminator

-  is branch initiator/terminator

 ## Basic parser

Parser should be written in some programming language, as of early 2026 C++ (generally speaking) seems to be a good option:

-   well established language, which has standards, reference documentation and examples

-   performance is considered to be high

-   is +/- known to the people outside of its community (like me), since syntax is influential: many other programming languages have somewhat similar syntax

-   compiled code could be called from many different environments relatively easy

-   code could be compiled to a WebAssembly module and be used in browser

The last two points are the main reasons in this case.

As it was mentioned, the following strategy of parsing is proposed:

1.  Computer program accepts the SMILES string, i.e., an ordered sequence of characters having chemical meaning as described previously.

2.  Computer program process this string from left to right one character at time.

    **What is meant by "computer program process"?**

    | - Computer program has default state.
    | - Every time computer program encounters new (next) character, state of the computer program changes accordingly (taking into account program's current state and what character it encounters).
    | - At each step computer program takes some action to produce an output.

3.  Computer program produces an output.

So, for starters, some code should be written to process each and every character in the SMILES string from left to right. As it can be seen from <https://en.cppreference.com/w/cpp/language/while.html>

**while** is an option to do just that.

Thus, the code of the SMILES parser in general could look like this:

```         
// not to run
while (i < smiles_str.length()) {
  // process i-th character:
  //    check state
  //    find character among the lableled characters
  //    If needed:
  //      update output in accordance with the current state and character's label
  //      update state
  i++;
}
```

From this code block it could be seen that some variables should be initialized prior to the main loop:

-   variable (object) to store the results

-   variable (objects) to store the state

-   variable(s) (objects) to store characters' labels and rules for processing

### Rust

On the other hand, Rust also seems to be an appropriate option for the task:

-   modern language: less legacy things both in terms of the code and documentation, which makes the search for information and actual coding much easier for beginners; ecosystem seems to be pretty compact and well structured

-   performance is considered to be high

-   venerated memory safety (probably not a deal breaker in the context of this task, but people could consider the resulting tool more reliable because of that, yes, they definately will)

-   compiled code could be called from many different environments

-   code could be compiled to a WebAssembly module and used in browser and is, as I can see, language is being developed with this option in mind

-   build with the ideas of functional programming in mind (opinionated)

### C++ vs Rust for the task according to the author's opinion

On the first glance for me it will be easier to use Rust considering the task, since it is more function-oriented, compact in the terms of documentation and standards, has memory safety.

It is an opinionated choice: no problem with C++, it allows to make things differently and safely, but C++ way will be harder generally due to the fact that there are a lot of things on it, many of which are object-oriented and are way too much for the task of a simple SMILES parser.

#### What is important

It is assumed that the results of this work will be used to feed the tools from R, JavaScript and Python ecosystems. How is it possible?

**C++ to R**, for example there is Rcpp package for that (<https://cran.r-project.org/web/packages/Rcpp/index.html>) it works fine.

**C++ to JavaScript (via Wasm)**, it surely works: <https://developer.mozilla.org/en-US/docs/WebAssembly/Guides/C_to_Wasm>

**C++ to Python**, definitely there are ways to do that, for example: <https://docs.python.org/3/extending/building.html#building>

**Rust to R**, there is extendR (<https://cran.r-project.org/web/packages/rextendr/index.html>)

**Rust to JavaScript (via Wasm)**, it surely works: <https://developer.mozilla.org/en-US/docs/WebAssembly/Guides/Rust_to_Wasm>

**Rust to Python,** there are ways to do that, for example: <https://github.com/pyo3/pyo3>

So, no surprise: both languages considered could be used for the task, however, using Rust is probably a more straightforward way. Time to read and check (<https://doc.rust-lang.org/book/>)

#### Decision

From what I have read so far it is quite clear that Rust is more than enough for me to proceed with the task: it has `structs` for the desired output, `mutables` in general to deal with the state, `references` to deal with the distinct characters from SMILES string in a safe and efficient way, `enums` or `HashSet` to check multiple conditions and compact ecosystem.

**Thus, Rust.**

> It should be noted that this decision is not so spontaneous as it may seems from the text above: over the years I have read news and discussions on the topic and still read, when I have time  
> <https://habr.com/>  
> <https://news.ycombinator.com/>  
> various books and other sources (I did not mention publications in the scientific peer-reviewed journals specifically, since reading them is like cheating, I mean, it may require some additional, not directly related to the task, skills, but probably all the concepts, which are in active use these days, could be traced back to such publications)  
> What is even more important I had a chance to talk with the people, colleagues, having experience and achievements in the field and did some work helping to understand what is what.

#### Proof-of-concept, toy SMILES parser

Still, it is hard to deduce all the possible obstacles from selective reading of the fascinating things only. To mitigate existing limitations it is possible to craft a toy parser with Rust, make Wasm module and check it alive before the deeper dive into analysis of SMILES strings.

Lets assume the SMILES parser is needed to parse strings describing structures consisting only of aliphatic carbon (C) atoms without additional properties where only single and double bonds are allowed and the whole structure is linear (no rings), like this:

| CCC=CCC

Many alkanes / alkenes (important organic molecules) are like this, thus, such a parser is not a joke.

The code in Rust will be as follows (I used online Rust playground to write and test it, https://play.rust-lang.org/?version=stable&mode=debug&edition=2024):

> [!NOTE]
> This code is not intended to be used in production

``` {.Rust .Rust}
fn main() {
use std::collections::HashSet; // data structures to store the conditions
use serde::Serialize;          // to get JSON from Struct
// define empty structs to store the results in  the desired format
    #[derive(Serialize)]
    #[derive(Clone)]
    struct Atom {
        atom_id: u32, // because, who knows?
        symbol: String,
        is_organic: bool,
        is_aromatic: bool,
        isin_bracket: bool,
        isotope: u8,
        chirality: String,
        n_hydrogens: u8,
        charge: i8,
        class: String
    }
    #[derive(Serialize)]
    #[derive(Clone)]
    struct Bond {
        bond_id: u32,
        atom_one: u32,
        atom_two: u32,
        bond_type: char
    }
    #[derive(Serialize)]
    struct Structure {
        atoms: Vec<Atom>,
        bonds: Vec<Bond>
    }
    // define the functions to parse limited SMILES
    fn parse_smiles(smiles_string: String) -> Structure {
        // define hash sets to store the allowed characters
        let mut allowed_atom: HashSet<char> = HashSet::new();
        let mut allowed_bond: HashSet<char> = HashSet::new();
        allowed_atom.insert('C');
        allowed_bond.insert('=');
        // define mutable to store the state
        let mut current_state: String = "".to_string();
        let mut prev_state: String = "".to_string();
        let mut current_atom_id: u32 = 1;
        let mut current_bond_id: u32 = 1;
        // Traverse the SMILES string filling the output
        let mut structure = Structure {
            atoms: Vec::new(),
            bonds: Vec::new()
        };
        for smiles_char in smiles_string.chars() {
            // update state
            if allowed_atom.contains(&smiles_char) {
                current_state = "atom".to_string();
            }
            if allowed_bond.contains(&smiles_char) {
                current_state = "bond".to_string();
            }
            println!("---");
            println!("previous state: {}", &prev_state.to_string());
            println!("current state: {}", &current_state.to_string());
            // Update the output
            if prev_state == "" && current_state == "atom" {
                // prepare on the first atom
                let first_atom = Atom {
                    atom_id: current_atom_id,
                    symbol: "C".to_string(),
                    is_organic: true,
                    is_aromatic: false,
                    isin_bracket: false,
                    isotope: 0,
                    chirality: "".to_string(),
                    n_hydrogens: 0,
                    charge: 0,
                    class: "".to_string()
                };
                // push the results
                structure.atoms.push(first_atom.clone());
                // increase atom id
                current_atom_id = current_atom_id + 1;
            }
            if prev_state == "atom" && current_state == "atom" {
                println!("atom-atom");
                // Add the data on the previous bond and current atom
                let typical_atom = Atom {
                    atom_id: current_atom_id,
                    symbol: "C".to_string(),
                    is_organic: true,
                    is_aromatic: false,
                    isin_bracket: false,
                    isotope: 0,
                    chirality: "".to_string(),
                    n_hydrogens: 0,
                    charge: 0,
                    class: "".to_string()
                };
                let typical_bond = Bond {
                    bond_id: current_bond_id,
                    atom_one: current_atom_id - 1,
                    atom_two: current_atom_id,
                    bond_type: '-'
                };
                // Push the results to the structure
                structure.atoms.push(typical_atom.clone());
                structure.bonds.push(typical_bond.clone());
                current_bond_id = current_bond_id + 1;
                current_atom_id = current_atom_id + 1;
            }
            if prev_state == "bond" && current_state == "atom" {
                println!("bond-atom");
                // Add the data on the previous bond and current atom
                let right_atom = Atom {
                    atom_id: current_atom_id,
                    symbol: "C".to_string(),
                    is_organic: true,
                    is_aromatic: false,
                    isin_bracket: false,
                    isotope: 0,
                    chirality: "".to_string(),
                    n_hydrogens: 0,
                    charge: 0,
                    class: "".to_string()
                };
                let double_bond = Bond {
                    bond_id: current_bond_id,
                    atom_one: current_atom_id - 1,
                    atom_two: current_atom_id,
                    bond_type: '='
                };
                // Push the results to structure
                structure.atoms.push(right_atom.clone());
                structure.bonds.push(double_bond.clone());
                current_bond_id = current_bond_id + 1;
                current_atom_id = current_atom_id + 1;
            }
            // Update previous state
            prev_state = current_state.clone();
        }
        // Return the value
        structure
    }
    // Execute
    let structure_rslt = parse_smiles("CC=C".to_string());
    let structure_json = serde_json::to_string_pretty(&structure_rslt).unwrap();
    println!("{}", structure_json);
}
```

Basically, this is the whole thing. However, logic should be extended, which will require changes to handle complexity; additional check-up and optimisations are required.

Still, with the given input:

> CC=C

This toy parser produces the following output:

``` JSON
{
  "atoms": [
    {
      "atom_id": 1,
      "symbol": "C",
      "is_organic": true,
      "is_aromatic": false,
      "isin_bracket": false,
      "isotope": 0,
      "chirality": "",
      "n_hydrogens": 0,
      "charge": 0,
      "class": ""
    },
    {
      "atom_id": 2,
      "symbol": "C",
      "is_organic": true,
      "is_aromatic": false,
      "isin_bracket": false,
      "isotope": 0,
      "chirality": "",
      "n_hydrogens": 0,
      "charge": 0,
      "class": ""
    },
    {
      "atom_id": 3,
      "symbol": "C",
      "is_organic": true,
      "is_aromatic": false,
      "isin_bracket": false,
      "isotope": 0,
      "chirality": "",
      "n_hydrogens": 0,
      "charge": 0,
      "class": ""
    }
  ],
  "bonds": [
    {
      "bond_id": 1,
      "atom_one": 1,
      "atom_two": 2,
      "bond_type": "-"
    },
    {
      "bond_id": 2,
      "atom_one": 2,
      "atom_two": 3,
      "bond_type": "="
    }
  ]
}
```

Which is correct. Thus, this toy / Proof-of-Concept example demonstrates for me that Rust is aligned well with the task. Must admit, memory safety and object-oriented things should be taken seriously during further development.

##### Compiling toy parser module to Web Assembly

In short, compilation using very clear instructions from

<https://developer.mozilla.org/en-US/docs/WebAssembly/Guides/Rust_to_Wasm>

was successful

<https://github.com/pavelVPo/IBMC_LSFBD_test_smiles_parser_Rust>

Not for production, just a test: code is not mature.

##### Licensing

**Rust (1.94.1)**, <https://rust-lang.org/> : [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) OR [MIT license](http://opensource.org/licenses/MIT)

**wasm-pack**, <https://www.npmjs.com/package/wasm-pack> : [MIT license](http://opensource.org/licenses/MIT) OR [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)

**wasm-bindgen**, <https://github.com/wasm-bindgen/wasm-bindgen> : [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) OR [MIT license](http://opensource.org/licenses/MIT)

**serde**, <https://crates.io/crates/serde> : [MIT license](http://opensource.org/licenses/MIT) OR [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)

##### Lessons learned from the toy parser writing and compiling

1.  Rust does not have the native JSON support, thus, additional code is needed to serialize / deserialize the data. Whether to use serde crate (<https://crates.io/crates/serde>) for Rust or write some dedicated code is the question. Probably, the usage of the dedicated crate is the optimal way, since serialization / deserialization is the common task, serde with the large user base exists, thus, serde.

2.  Important part of the desired parser is the `state` variable, since it will store parameters essential for parsing. The example above shows that the previous symbol is important for the successful parsing. Thus, previous work on SMILES analysis will help a lot a would continued and finalized with this understanding.

#### Proof-of-concept, toy SMILES parser V2

So, here is an idea: use just allowed SMILES symbols (without further division into character classes, etc) to parse the whole SMILES string.

To check idea's viability the toy parser V2 using Rust will be developed to parse the first two symbols (pair of any allowed symbols) out of three symbol string (third symbol is each symbol possible).

This will show if it is possible to distinguish between the SMILES symbols this way on the example of all allowed pairs of symbols and all possible symbols following them.

Assuming the state is appropriate.

**Here is the code (not for production, without considering right methods of string processing, just to test the idea and sheer performance):**

``` Rust
use std::io::{self, Read, Write};
use std::fs::File;
use rayon::prelude::*;

fn main() {
    // symbols, ... to be short
    const SYMBS_AR: [&str;2107] = ["b", "c", "n", "o", ..., ":111"];
// generate triples", "parse pairs", "check and collect results if parsing fails
    let fails: Vec<String> = SYMBS_AR.par_iter_mut()
            .map(|i|{
                let mut fails_i: Vec<String> = Vec::new();
                // mutable array for triple
                let mut query_ar: [&str;3] = ["";3];
                // mutable query string
                let mut current_str = String::from("n");
                // mutable string to store the first symbol
                let mut symbone_str = String::from("");
                // mutable string to store the second symbol
                let mut symbtwo_str = String::from("");
                query_ar[0] = i;
                // first inner loop
                for k in &SYMBS_AR {
                    query_ar[1] = k;
                    // second inner loop
                    for j in &SYMBS_AR {
                        query_ar[2] = j;
                        current_str = "".to_string();
                        let first_query_str = query_ar.join("");
                        let first_query: Vec<String> = first_query_str.chars().map(|c| c.to_string()).collect();
                        for i_q in &first_query {
                        current_str = [current_str, i_q.to_string()].join("");
                            // append to curent string
                            if SYMBS_AR.contains(&&current_str[..]) {
                                symbone_str = current_str.clone();
                            }
                        }
                        // find the second symbol in query
                        let second_query_str = first_query_str.replacen(&symbone_str, "", 1);
                        let second_query: Vec<String> = second_query_str.chars().map(|c| c.to_string()).collect();
                        current_str = "".to_string();
                        for i_q in &second_query {
                            current_str = [current_str, i_q.to_string()].join("");
                            // append to curent string
                            if SYMBS_AR.contains(&&current_str[..]) {
                                symbtwo_str = current_str.clone();
                                let tempi = format!("{}-{}-{}-{}-{}", query_ar[0].clone(), &k.clone(), &j.clone(), &symbone_str.clone(), &symbtwo_str.clone()).to_string();
                                if query_ar[0].to_string() != symbone_str.to_string() || query_ar[1].to_string() != symbtwo_str.to_string() {
                                    fails_i.push(tempi.to_string().clone());
                                }
                            }
                        }
                    }
                } 
                let result = fails_i.join("_qqq_").to_string().clone();
                result
            })
            .collect();
    let mut file = File::create("failed__triples-pairs.tsv").expect("Unable to create file");                                                                                                          
    for i in &fails{                                                                                                                                                                  
        file.write_all((i.as_bytes())).expect("Unable to write data");
    }
    // Memory consumption is huge; efficeint way to access / process strings should be considered 
}
```

Using the code above about 9 billions of SMILES symbols' triples were generated, first pair of symbols was parsed; no parsing errors were identified.

The whole process took several hours, so, compiled Rust code is fast before even the basic optimization.

Idea of more simple parser seems to be valid.
