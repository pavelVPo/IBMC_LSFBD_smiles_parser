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

    **What is an "output"?**

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

### Atom symbol type and corresponding symbol and character classes

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

> b, c, n, o, s, p

Corresponding characters could be designated as distinct character class, **w_atom_oar**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

2.  Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal**):

> B, C, N, O, S, P, F, I

Corresponding characters could be designated as distinct character class, **w_atom_oal**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

3.  Single character atom symbols of aromatic atoms enclosed within brackets (**atom_bar**):

> b, c, n, o, s, p

Corresponding characters could be designated as **w_atom_bar**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic. As it can be seen, this class contains the same symbols as w_atom_oar, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

4.  Single character atom symbols of bracket aliphatic atoms (**atom_bal**):

> H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **w_atom_bal**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic. As it can be seen, this class contains the same symbols as w_atom_oal, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

5.  Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal_2**):

> Cl, Br

Corresponding character classes could be designated as **s_atom_oal** & **e_atom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic. Further division of the characters describing symbol into two classes could be useful if the resulting parser will operate one character at time.

6.  Two character atom symbols of in-bracket aromatic atoms (**atom_bar_2**):

> se, as, te

Corresponding characters could be designated as **s_atom_bar** & **e_atom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

7.  Two character atom symbols of in-bracket aliphatic atoms (**atom_bal_2**):

> He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **s_atom_bal** & **e_atom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

### Symbol of anything and corresponding symbol and character class

**\*** is an allowed symbol in SMILES, it corresponds to any (atom) symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms :

8.  Single character symbol of any atom or basically **anything**:

> \*

Corresponding character could be designated as **w_any**, since it could mean basically any character or symbol.

### Square bracket symbol type and corresponding symbol and character classes

#### What are they?

9, 10. Square bracket symbols **[, ]** is the SMILES way to mark the start of the atom record including its various properties and is the way to designate the end of the atom record.

> [, ]

Corresponding symbol and character classes could be designated as **s_bracket** & **e_bracket**.

### Bond symbols and corresponding character classes

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, all of them are single character and do not have other peculiar aspects, five of them correspond to the conventional type of chemical bond:**

11. Single character bond symbol corresponding to the single bond (**single_bond**):

> -   

This single character symbol could be and typically is omitted, since by default all the atoms, which symbols are written side by side in SMILES string, are presumed to be connected by this type of bond. Corresponding character class will be designated as **w_single_bond**.

12. Single character bond symbol corresponding to the double bond (**double_bond**):

> =

Corresponding character class will be designated as **w_double_bond**.

13. Single character bond symbol corresponding to the triple bond (**triple_bond**):

> \#

Corresponding character class will be designated as **w_triple_bond**.

14. Single character symbol corresponding to the quadruple bond (**quadruple_bond**):

> \$

Corresponding character class will be designated as **w_quadruple_bond**.

15. Single character symbol corresponding to the aromatic bond (**aromatic_bond**):

> :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications. Corresponding character class will be designated as **w_aromatic_bond**.

16. Single character symbol corresponding to the absence of the bond between the two specific atoms (**no_bond**):

> .

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

17. Single character bond multiplying symbols initiators of branching with implicit bond (**bm_ibi**):

> (

Corresponding character class could be designated as **w_bm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

18. Single character bond multiplying symbols initiators of rings with implicit bond (**bm_iri**):

> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_bm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

19. Single character bond multiplying symbols terminators of branching with implicit bond (**bm_tbi**):

> )

Corresponding characters could be designated as **w_bm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

20. Single character bond multiplying symbols terminators of rings with implicit bond (**bm_tri**):

> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **w_bm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

21. Two-character bond multiplying symbols initiators of branching with explicit bond (**bm_ibe**):

> ([-=#\$:.]

Corresponding character classes could be designated as **s_bm_ibe** & **e_bm_ibe** where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, second suffix **e** - for explicit.

22. Two-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_2**):

> [-=#\$:.][0-9]

Corresponding characters could be designated as **s_bm_ire** & **e_bm_ire**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit.

23. Three-character bond multiplying symbols initiators of rings with implicit bond (**bm_iri_3**):

> \%[0-9][1-9], %[1-9][0-9]

Corresponding character classes could be designated as **s_bm_iri** & **r_bm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, next suffix **i** - for implicit.

24. Four-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_4**):

> [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding character classes could be designated as **s_bm_ire_4, n_bm_ire** & **r_bm_ire**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit; **4** in s_bm_ire_4 is used to differentiate this class from the s_bm_ire corresponding to the bm_ire_2, n_bm_ire and r_bm_ire are already unique.

25. Two-character bond multiplying symbols terminators of branching with explicit bond (**bm_tbe_2**):

> )[-=#\$:.]

Corresponding character classes could be designated as **s_bm_tbe** & **e_bm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

26. Two-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_2**):

> [-=#\$:.][0-9]

Corresponding character classes could be designated as **s_bm_tre** & **e_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit.

27. Four-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_4**):

> [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding characters could be designated as **s_bm_tre, n_bm_tre** & **r_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit. **4** in s_bm_tre_4 is used to differentiate this class from the s_bm_tre corresponding to the bm_tre_2, n_bm_tre and r_bm_tre are already unique.

28. Three-character bond multiplying symbols terminators of rings with implicit bond (**bm_tri_3**):

> \%[0-9][1-9], %[1-9][0-9]

Corresponding characters could be designated as **s_bm_tri** & **r_bm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, next suffix **i** - for implicit.

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two facets of cis/trans symbols are allowed in SMILES:

29. Cis/trans single character symbols on the left side of the rotary non-permissive bond (**l_ct**):

> /, \\

Corresponding single character characters could be designated as **l_ct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

30. Cis/trans symbols on the right side of the rotary non-permissive bond:

> /, \\

Corresponding characters could be designated as **r_ct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\\F)=C/F are equivalent, since

> The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.

> *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****

> **\*** <http://opensmiles.org/opensmiles.html>

### All the symbols and corresponding character classes inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating specific isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

31. Single character isotope symbols (**isotope**):

> 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_isotope**, where prefix **w** stands for the whole symbol.

32. Multicharacter (from 2 to 3 characters) isotope symbols (**isotope_m**):

> [0-9][1-9], [1-9][0-9], [0-9][0-9][1-9], [0-9][1-9][0-9], [1-9][0-9][0-9]

Corresponding characters could be designated as **s_isotope & r_isotope**, where prefix **s** stands for the start and prefix **r** stands for the rest of the symbol.

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

33. Single character chirality symbol (**chiral**):

> \@

Corresponding character class could be designated as **w_chiral**, where prefix **w** stands for the whole symbol.

34. Two-character chirality symbols (**chiral_2**):

> [\@][\@]

Corresponding character classes could be designated as **s_chiral & e_chiral**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

35. Multicharacter (four or five character) chirality symbols (**chiral_m**):

> [\@]TH[1-2], [\@]AL[1-2], [\@]SP[1-3], [\@]TB[1-20], [\@]OH[1-30]

Corresponding character classes could be designated as **s_chiral_m, m_chiral & r_chiral**, where prefix **s** stands for the start, *prefix* **m** stands for the medium (two characters), prefix **r** stands for the rest of the symbol and *suffix* **m** stands for the multi, where it is needed.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 facets by their length:

36. Single character hydrogen symbol (**hydro**):

> H

Corresponding character classes could be designated as **w_hydro**, where prefix **w** stands for the whole symbol.

37. Two-character hydrogen symbols (**hydro_2**):

> H[2-9]

Corresponding character classes could be designated as **s_hydro & ehydro**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom (**charge**).

Charge symbols allowed in SMILES could be divided into 2 facets by their length:

38. Single character charge symbols (**charge**):

> [+-]

Corresponding character classes could be designated as **w_charge**, where prefix **w** stands for the whole symbol.

39. Two-character charge obsolete symbols (**charge_2**):

> [+][+], [-][-]

Corresponding characters could be designated as **s_charge & e_charge**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

40. Multicharacter (two or three characters) charge symbols (**charge_m**):

> [+-][1-9], [+-]1[0-5]

Corresponding characters could be designated as **s_charge_m & r_charge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi where it is needed.

##### Class symbols

Class symbols designate user-defined class of the atom.

Class symbols allowed in SMILES may have variable length, but there is no point to divide them into facets based on this aspect, so there is only:

41. Multicharacter (from 2 to 4 characters) atom class symbol type (**class**):

> :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]

Corresponding character classes could be designated as **s_class & r_class**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

Information on symbols is summarized in **symbols.tsv**, all symbols are provided in **symbols_all.tsv**.

## Q1: are described symbols unique, i.e. is it possible to identify each SMILES symbols based only on characters constituting it?

No, as it can be seen from **symbols.tsv**. For example,

-   terminators and initiators of branching could be and often are the same by design

-   symbols of in-bracket atoms and bracket-free atoms could be the same (organic subset)

Also, several classes of symbols could be described or are described partially by the patterns [0-9] and [1-9], which makes parsing without consideration of the environment questionable.

In the previous version the attempt was taken to divide the whole symbols into the smaller subsets of characters, the result is as follows: this strategy does not pay off, character classes probably could be useful to construct the SMILES strings, but they provide no clear benefits for parsing:

> [!NOTE]

> It seems to be easier to read the characters one by one until the longest possible sequence describing symbol is gathered (5 characters, I guess) and decide on the actual symbol afterwards, considering matches in **current sub-string** and **state** deduced from the previous symbols and length of the remaining SMILES string.

> The following text will be rewritten accordingly.

## General SMILES parsing strategy revised

1.  Computer program initializes with the SMILES string having default **state** and empty **accumulator** of characters and empty **result**.

2.  Every time computer program encounters new (next) character it accumulates this character.

3.  Every time accumulator reaches its limits (longest symbol in SMILES or end of the string) its content is being evaluated, state changes accordingly, accumulator's content is being trimmed from left to right to delete all the characters evaluated as the whole symbol at this step.

4.  Every time state changes computer program takes some action to build up an output.

5.  When the end of the string is reached, computer program outputs the result.

## Pairs of symbol types, which are not allowed in SMILES

> [!NOTE]

> Chemistry deals with the complex objects, which structures are under strict constrains, however, even in general nonsensical chemical structures may make sense under specific conditions, representing the intermediate state, or due to the technical reasons (data analysis).

> Thus, there is no objective to produce only the chemically valid structures as the result of SMILES parsing.

> However, some restrictions concerning the SMILES string validity should be enforced and to identify them, all the theoretically possible pairs of symbol types will be enumerated.

Among the 36 theoretically possible pairs of symbol types the following will be forbidden:

-   **bond - bond:** bonds connect atoms, not bonds

-   **bond - modifier:** modifiers include bonds, bonds connect atoms, not bonds

-   **bond - property:** properties are allowed only inside the square brackets, bonds are allowed only outside the square brackets

-   **modifier - property:** properties are allowed only inside the square brackets, modifiers are allowed only outside the square brackets

-   **property - bond:** properties are allowed only inside the square brackets, bonds are allowed only outside the square brackets

-   **property - modifier:** properties are allowed only inside the square brackets, modifiers are allowed only outside the square brackets

## Pairs of symbol classes, which are not allowed in SMILES

1.  Out-of-bracket and in-bracket atoms are not allowed to be paired directly, the following pairs of symbol classes are forbidden:

    > atom_oal -\> atom_bal, atom_oal -\> atom_bal_2, atom_oal -\> atom_bar_2, atom_oal_2 -\> atom_bal, atom_oal_2 -\> atom_bal_2, atom_oal_2 -\> atom_bar, atom_oal_2 -\> atom_bar_2, atom_oar -\> atom_bal, atom_oar -\> atom_bal_2, atom_oar -\> atom_bar, atom_oar -\> atom_bar_2, atom_bal -\> atom_oal, atom_bal -\> atom_oal_2, atom_bal -\> atom_oar, atom_bal_2 -\> atom_oal, atom_bal_2 -\> atom_oal_2, atom_bal_2 -\> atom_oar, atom_bar -\> atom_oal, atom_bar -\> atom_oal_2, atom_bar -\> atom_oar, atom_bar_2 -\> atom_oal, atom_bar_2 -\> atom_oal_2, atom_bar_2 -\> atom_oar

    It should be noted that application of this rule depends on the state, since atom_bar_2 and atom_oal_2, etc. are the same symbols in different contexts.

2.  In-bracket atoms are not allowed to be paired directly, the following pairs of symbol classes are forbidden:

    > atom_bal -\> atom_bal, atom_bal -\> atom_bal_2, atom_bal -\> atom_bar, atom_bal -\> atom_bar_2, atom_bal_2 -\> atom_bal, atom_bal_2 -\> atom_bal_2, atom_bal_2 -\> atom_bar, atom_bal_2 -\> atom_bar_2, atom_bar -\> atom_bal, atom_bar -\> atom_bal_2, atom_bar -\> atom_bar, atom_bar -\> atom_bar_2, atom_bar_2 -\> atom_bal, atom_bar_2 -\> atom_bal_2, atom_bar_2 -\> atom_bar, atom_bar_2 -\> atom_bar_2

    This rule depends on the context.

3.  In-bracket atoms are not allowed to be paired directly with the anything symbol (any atom), the following pairs of symbol classes are forbidden:

    > atom_bal -\> anything, atom_bal_2 -\> anything, atom_bar -\> anything, atom_bar_2 -\> anything, anything -\> atom_bal, anything -\> atom_bal_2, anything -\> atom_bar, anything -\> atom_bar_2

4.  In-bracket atoms are not allowed to be paired directly with bond symbols, the following pairs of symbol classes are forbidden:

    > atom_bal -\> aromatic_bond, atom_bal -\> double_bond, atom_bal -\> no_bond, atom_bal -\> quadruple_bond, atom_bal -\> single_bond, atom_bal -\> triple_bond, atom_bal_2 -\> aromatic_bond, atom_bal_2 -\> double_bond, atom_bal_2 -\> no_bond, atom_bal_2 -\> quadruple_bond, atom_bal_2 -\> single_bond, atom_bal_2 -\> triple_bond, atom_bar -\> aromatic_bond, atom_bar -\> double_bond, atom_bar -\> no_bond, atom_bar -\> quadruple_bond, atom_bar -\> single_bond, atom_bar -\> triple_bond, atom_bar_2 -\> aromatic_bond, atom_bar_2 -\> double_bond, atom_bar_2 -\> no_bond, atom_bar_2 -\> quadruple_bond, atom_bar_2 -\> single_bond, atom_bar_2 -\> triple_bond, aromatic_bond -\> atom_bal, aromatic_bond -\> atom_bal_2, aromatic_bond -\> atom_bar, aromatic_bond -\> atom_bar_2, double_bond -\> atom_bal, double_bond -\> atom_bal_2, double_bond -\> atom_bar, double_bond -\> atom_bar_2, no_bond -\> atom_bal, no_bond -\> atom_bal_2, no_bond -\> atom_bar, no_bond -\> atom_bar_2, quadruple_bond -\> atom_bal, quadruple_bond -\> atom_bal_2, quadruple_bond -\> atom_bar, quadruple_bond -\> atom_bar_2, single_bond -\> atom_bal, single_bond -\> atom_bal_2, single_bond -\> atom_bar, single_bond -\> atom_bar_2, triple_bond -\> atom_bal, triple_bond -\> atom_bal_2, triple_bond -\> atom_bar, triple_bond -\> atom_bar_2

    This rule depends on the context.

5.  In-bracket atoms are not allowed to be paired directly with the modifiers, the following pairs of symbol classes are forbidden:

    > atom_bal -\> bm_ibe, atom_bal -\> bm_ibi, atom_bal -\> bm_ire_2, atom_bal -\> bm_ire_4, atom_bal -\> bm_iri, atom_bal -\> bm_iri_3, atom_bal -\> bm_tbe_2, atom_bal -\> bm_tbi, atom_bal -\> bm_tre_2, atom_bal -\> bm_tre_4, atom_bal -\> bm_tri, atom_bal -\> bm_tri_3, atom_bal -\> l_ct, atom_bal -\> r_ct, atom_bal_2 -\> bm_ibe, atom_bal_2 -\> bm_ibi, atom_bal_2 -\> bm_ire_2, atom_bal_2 -\> bm_ire_4, atom_bal_2 -\> bm_iri, atom_bal_2 -\> bm_iri_3, atom_bal_2 -\> bm_tbe_2, atom_bal_2 -\> bm_tbi, atom_bal_2 -\> bm_tre_2, atom_bal_2 -\> bm_tre_4, atom_bal_2 -\> bm_tri, atom_bal_2 -\> bm_tri_3, atom_bal_2 -\> l_ct, atom_bal_2 -\> r_ct, atom_bar -\> bm_ibe, atom_bar -\> bm_ibi, atom_bar -\> bm_ire_2, atom_bar -\> bm_ire_4, atom_bar -\> bm_iri, atom_bar -\> bm_iri_3, atom_bar -\> bm_tbe_2, atom_bar -\> bm_tbi, atom_bar -\> bm_tre_2, atom_bar -\> bm_tre_4, atom_bar -\> bm_tri, atom_bar -\> bm_tri_3, atom_bar -\> l_ct, atom_bar -\> r_ct, atom_bar_2 -\> bm_ibe, atom_bar_2 -\> bm_ibi, atom_bar_2 -\> bm_ire_2, atom_bar_2 -\> bm_ire_4, atom_bar_2 -\> bm_iri, atom_bar_2 -\> bm_iri_3, atom_bar_2 -\> bm_tbe_2, atom_bar_2 -\> bm_tbi, atom_bar_2 -\> bm_tre_2, atom_bar_2 -\> bm_tre_4, atom_bar_2 -\> bm_tri, atom_bar_2 -\> bm_tri_3, atom_bar_2 -\> l_ct, atom_bar_2 -\> r_ct,bm_ibe -\> atom_bal, bm_ibe -\> atom_bal_2, bm_ibe -\> atom_bar, bm_ibe -\> atom_bar_2, bm_ibi -\> atom_bal, bm_ibi -\> atom_bal_2, bm_ibi -\> atom_bar, bm_ibi -\> atom_bar_2, bm_ire_2 -\> atom_bal, bm_ire_2 -\> atom_bal_2, bm_ire_2 -\> atom_bar, bm_ire_2 -\> atom_bar_2, bm_ire_4 -\> atom_bal, bm_ire_4 -\> atom_bal_2, bm_ire_4 -\> atom_bar, bm_ire_4 -\> atom_bar_2, bm_iri -\> atom_bal, bm_iri -\> atom_bal_2, bm_iri -\> atom_bar, bm_iri -\> atom_bar_2, bm_iri_3 -\> atom_bal, bm_iri_3 -\> atom_bal_2, bm_iri_3 -\> atom_bar, bm_iri_3 -\> atom_bar_2, bm_tbe_2 -\> atom_bal, bm_tbe_2 -\> atom_bal_2, bm_tbe_2 -\> atom_bar, bm_tbe_2 -\> atom_bar_2, bm_tbi -\> atom_bal, bm_tbi -\> atom_bal_2, bm_tbi -\> atom_bar, bm_tbi -\> atom_bar_2, bm_tre_2 -\> atom_bal, bm_tre_2 -\> atom_bal_2, bm_tre_2 -\> atom_bar, bm_tre_2 -\> atom_bar_2, bm_tre_4 -\> atom_bal, bm_tre_4 -\> atom_bal_2, bm_tre_4 -\> atom_bar, bm_tre_4 -\> atom_bar_2, bm_tri -\> atom_bal, bm_tri -\> atom_bal_2, bm_tri -\> atom_bar, bm_tri -\> atom_bar_2, bm_tri_3 -\> atom_bal, bm_tri_3 -\> atom_bal_2, bm_tri_3 -\> atom_bar, bm_tri_3 -\> atom_bar_2, l_ct -\> atom_bal, l_ct -\> atom_bal_2, l_ct -\> atom_bar, l_ct -\> atom_bar_2, r_ct -\> atom_bal, r_ct -\> atom_bal_2, r_ct -\> atom_bar, r_ct -\> atom_bar_2

    This rule depends on the context.

6.  In-bracket atoms are not allowed to be immediately followed by the isotope classes, the following pairs of symbol classes are forbidden:

    atom_bal -\> isotope, atom_bal -\> isotope_m, atom_bal_2 -\> isotope, atom_bal_2 -\> isotope_m, atom_bar -\> isotope, atom_bar -\> isotope_m, atom_bar_2 -\> isotope, atom_bar_2 -\> isotope_m

7.  In-bracket atoms could not be immediately preceded by the in-bracket properties, besides, isotopic number; the following pairs of symbol classes are forbidden:

    > charge -\> atom_bal, charge -\> atom_bal_2, charge -\> atom_bar, charge -\> atom_bar_2, charge_2 -\> atom_bal, charge_2 -\> atom_bal_2, charge_2 -\> atom_bar, charge_2 -\> atom_bar_2, charge_m -\> atom_bal, charge_m -\> atom_bal_2, charge_m -\> atom_bar, charge_m -\> atom_bar_2, chiral -\> atom_bal, chiral -\> atom_bal_2, chiral -\> atom_bar, chiral -\> atom_bar_2, chiral_2 -\> atom_bal, chiral_2 -\> atom_bal_2, chiral_2 -\> atom_bar, chiral_2 -\> atom_bar_2, chiral_m -\> atom_bal, chiral_m -\> atom_bal_2, chiral_m -\> atom_bar, chiral_m -\> atom_bar_2, class -\> atom_bal, class -\> atom_bal_2, class -\> atom_bar, class -\> atom_bar_2, hydro -\> atom_bal, hydro -\> atom_bal_2, hydro -\> atom_bar, hydro -\> atom_bar_2, hydro_2 -\> atom_bal, hydro_2 -\> atom_bal_2, hydro_2 -\> atom_bar, hydro_2 -\> atom_bar_2

8.  In-bracket atoms could not be immediately precede the start of the square brackets, the following pairs of symbol classes are forbidden:

    > atom_bal -\> s_bracket, atom_bal_2 -\> s_bracket, atom_bar -\> s_bracket, atom_bar_2 -\> s_bracket

9.  Out-of-bracket atoms are not allowed to be paired with the in-bracket properties, the following pairs of symbols are forbidden:

    > atom_oal -\> charge, atom_oal -\> charge_2, atom_oal -\> charge_m, atom_oal -\> chiral, atom_oal -\> chiral_2, atom_oal -\> chiral_m, atom_oal -\> class, atom_oal -\> hydro, atom_oal -\> hydro_2, atom_oal -\> isotope, atom_oal -\> isotope_m, atom_oal_2 -\> charge, atom_oal_2 -\> charge_2, atom_oal_2 -\> charge_m, atom_oal_2 -\> chiral, atom_oal_2 -\> chiral_2, atom_oal_2 -\> chiral_m, atom_oal_2 -\> class, atom_oal_2 -\> hydro, atom_oal_2 -\> hydro_2, atom_oal_2 -\> isotope, atom_oal_2 -\> isotope_m, atom_oar -\> charge, atom_oar -\> charge_2, atom_oar -\> charge_m, atom_oar -\> chiral, atom_oar -\> chiral_2, atom_oar -\> chiral_m, atom_oar -\> class, atom_oar -\> hydro, atom_oar -\> hydro_2, atom_oar -\> isotope, atom_oar -\> isotope_m, charge -\> atom_oal, charge -\> atom_oal_2, charge -\> atom_oar, charge_2 -\> atom_oal, charge_2 -\> atom_oal_2, charge_2 -\> atom_oar, charge_m -\> atom_oal, charge_m -\> atom_oal_2, charge_m -\> atom_oar, chiral -\> atom_oal, chiral -\> atom_oal_2, chiral -\> atom_oar, chiral_2 -\> atom_oal, chiral_2 -\> atom_oal_2, chiral_2 -\> atom_oar, chiral_m -\> atom_oal, chiral_m -\> atom_oal_2, chiral_m -\> atom_oar, class -\> atom_oal, class -\> atom_oal_2, class -\> atom_oar, hydro -\> atom_oal, hydro -\> atom_oal_2, hydro -\> atom_oar, hydro_2 -\> atom_oal, hydro_2 -\> atom_oal_2, hydro_2 -\> atom_oar, isotope -\> atom_oal, isotope -\> atom_oal_2, isotope -\> atom_oar, isotope_m -\> atom_oal, isotope_m -\> atom_oal_2, isotope_m -\> atom_oar

10. Out-of-bracket aliphatic atoms are not allowed to be paired directly with the symbols of aromatic bonds, the following pairs of symbol classes are forbidden:

    > aromatic_bond -\> atom_oal, aromatic_bond -\> atom_oal_2, atom_oal -\> aromatic_bond, atom_oal_2 -\> aromatic_bond

11. Bonds are not allowed to be immediately followed by the ending square brackets, the following pairs of symbol classes are forbidden:

    > aromatic_bond -\> e_bracket, double_bond -\> e_bracket, no_bond -\> e_bracket, quadruple_bond -\> e_bracket, single_bond -\> e_bracket, triple_bond -\> e_bracket

12. Branching initiators and terminators with explicit bonds could not be paired with the intiators / terminators of branching having explicit bonds, the following symbol classes pairs are forbidden:

    > bm_ibe -\> aromatic_bond, bm_ibe -\> double_bond, bm_ibe -\> no_bond, bm_ibe -\> quadruple_bond, bm_ibe -\> single_bond, bm_ibe -\> triple_bond, bm_tbe_2 -\> aromatic_bond, bm_tbe_2 -\> double_bond, bm_tbe_2 -\> no_bond, bm_tbe_2 -\> quadruple_bond, bm_tbe_2 -\> single_bond, bm_tbe_2 -\> triple_bond

13. Some pairs of the of the bond modifiers having explicit bonds should be avoided due to the intristinct uncertainties:

    > bm_ibe -\> bm_ibe, bm_ibe -\> bm_ire_2, bm_ibe -\> bm_ire_4, bm_ibe -\> bm_tbe_2, bm_ibe -\> bm_tbi, bm_ibe -\> bm_tre_2, bm_ibe -\> bm_tre_4, bm_ibi -\> bm_ire_2, bm_ibi -\> bm_ire_4, bm_ibi -\> bm_tbe_2, bm_ibi -\> bm_tbi, bm_ibi -\> bm_tre_2, bm_ibi -\> bm_tre_4, bm_tbe_2 -\> bm_ibe, bm_tbe_2 -\> bm_ibi, bm_tbe_2 -\> bm_ire_2, bm_tbe_2 -\> bm_ire_4, bm_tbe_2 -\> bm_tbe_2, bm_tbe_2 -\> bm_tbi, bm_tbe_2 -\> bm_tre_2, bm_tbe_2 -\> bm_tre_4

The latter forbidden pairs of classes were constructed from the leftovers, thus, they could be less consistent:

1.  Ending square bracket should not precede the in-bracket atoms, in-bracket properties, and other bracket ending symbol:

    > e_bracket -\> atom_bal, e_bracket -\> atom_bal_2, e_bracket -\> atom_bar, e_bracket -\> atom_bar_2, e_bracket -\> charge, e_bracket -\> charge_2, e_bracket -\> charge_m, e_bracket -\> chiral, e_bracket -\> chiral_2, e_bracket -\> chiral_m, e_bracket -\> class, e_bracket -\> e_bracket, e_bracket -\> hydro, e_bracket -\> hydro_2, e_bracket -\> isotope, e_bracket -\> isotope_m

2.  Ending square bracket should not be preceded by the in-bracket symbols, which were not mentioned earlier and starting bracket:

    > atom_oal -\> e_bracket, atom_oal_2 -\> e_bracket, atom_oar -\> e_bracket, bm_ibe -\> e_bracket, bm_ibi -\> e_bracket, bm_ire_2 -\> e_bracket, bm_ire_4 -\> e_bracket, bm_iri -\> e_bracket, bm_iri_3 -\> e_bracket, bm_tbe_2 -\> e_bracket, bm_tbi -\> e_bracket, bm_tre_2 -\> e_bracket, bm_tre_4 -\> e_bracket, bm_tri -\> e_bracket, bm_tri_3 -\> e_bracket, l_ct -\> e_bracket, r_ct -\> e_bracket, s_bracket -\> e_bracket

3.  Starting bracket should not be followed by the not in-bracket symbols not mentioned earlier including another starting bracket:

    > s_bracket -\> aromatic_bond, s_bracket -\> atom_oal, s_bracket -\> atom_oal_2, s_bracket -\> atom_oar, s_bracket -\> bm_ibe, s_bracket -\> bm_ibi, s_bracket -\> bm_ire_2, s_bracket -\> bm_ire_4, s_bracket -\> bm_iri, s_bracket -\> bm_iri_3, s_bracket -\> bm_tbe_2, s_bracket -\> bm_tbi, s_bracket -\> bm_tre_2, s_bracket -\> bm_tre_4, s_bracket -\> bm_tri, s_bracket -\> bm_tri_3, s_bracket -\> charge, s_bracket -\> charge_2, s_bracket -\> charge_m, s_bracket -\> chiral, s_bracket -\> chiral_2, s_bracket -\> chiral_m, s_bracket -\> class, s_bracket -\> double_bond, s_bracket -\> hydro, s_bracket -\> hydro_2, s_bracket -\> l_ct, s_bracket -\> no_bond, s_bracket -\> quadruple_bond, s_bracket -\> r_ct, s_bracket -\> s_bracket, s_bracket -\> single_bond, s_bracket -\> triple_bond

4.  Starting bracket should not be preceded by the in-bracket symbols, the following pairs of symbol classes are forbidden:

    > charge -\> s_bracket, charge_2 -\> s_bracket, charge_m -\> s_bracket, chiral -\> s_bracket, chiral_2 -\> s_bracket, chiral_m -\> s_bracket, class -\> s_bracket, hydro -\> s_bracket, hydro_2 -\> s_bracket, isotope -\> s_bracket, isotope_m -\> s_bracket

5.  In-bracket property symbols should not be followed by the other various inappropriate symbols:

    > charge -\> anything, charge -\> atom_oal, charge -\> atom_oal_2, charge -\> atom_oar, charge -\> charge, charge -\> charge_2, charge -\> charge_m, charge -\> chiral, charge -\> chiral_2, charge -\> chiral_m, charge -\> hydro, charge -\> hydro_2, charge -\> isotope, charge -\> isotope_m, charge_2 -\> anything, charge_2 -\> atom_oal, charge_2 -\> atom_oal_2, charge_2 -\> atom_oar, charge_2 -\> charge, charge_2 -\> charge_2, charge_2 -\> charge_m, charge_2 -\> chiral, charge_2 -\> chiral_2, charge_2 -\> chiral_m, charge_2 -\> hydro, charge_2 -\> hydro_2, charge_2 -\> isotope, charge_2 -\> isotope_m, charge_m -\> anything, charge_m -\> atom_oal, charge_m -\> atom_oal_2, charge_m -\> atom_oar, charge_m -\> charge, charge_m -\> charge_2, charge_m -\> charge_m, charge_m -\> chiral, charge_m -\> chiral_2, charge_m -\> chiral_m, charge_m -\> hydro, charge_m -\> hydro_2, charge_m -\> isotope, charge_m -\> isotope_m, chiral -\> anything, chiral -\> atom_oal, chiral -\> atom_oal_2, chiral -\> atom_oar, chiral -\> chiral, chiral -\> chiral_2, chiral -\> chiral_m, chiral -\> isotope, chiral -\> isotope_m, chiral_2 -\> anything, chiral_2 -\> atom_oal, chiral_2 -\> atom_oal_2, chiral_2 -\> atom_oar, chiral_2 -\> chiral, chiral_2 -\> chiral_2, chiral_2 -\> chiral_m, chiral_2 -\> isotope, chiral_2 -\> isotope_m, chiral_m -\> anything, chiral_m -\> atom_oal, chiral_m -\> atom_oal_2, chiral_m -\> atom_oar, chiral_m -\> chiral, chiral_m -\> chiral_2, chiral_m -\> chiral_m, chiral_m -\> isotope, chiral_m -\> isotope_m, class -\> anything, class -\> atom_oal, class -\> atom_oal_2, class -\> atom_oar, class -\> charge, class -\> charge_2, class -\> charge_m, class -\> chiral, class -\> chiral_2, class -\> chiral_m, class -\> class, class -\> hydro, class -\> hydro_2, class -\> isotope, class -\> isotope_m, hydro -\> anything, hydro -\> atom_oal, hydro -\> atom_oal_2, hydro -\> atom_oar, hydro -\> chiral, hydro -\> chiral_2, hydro -\> chiral_m, hydro -\> hydro, hydro -\> hydro_2, hydro -\> isotope, hydro -\> isotope_m, hydro_2 -\> anything, hydro_2 -\> atom_oal, hydro_2 -\> atom_oal_2, hydro_2 -\> atom_oar, hydro_2 -\> chiral, hydro_2 -\> chiral_2, hydro_2 -\> chiral_m, hydro_2 -\> hydro, hydro_2 -\> hydro_2, hydro_2 -\> isotope, hydro_2 -\> isotope_m, isotope -\> atom_oal, isotope -\> atom_oal_2, isotope -\> atom_oar, isotope -\> charge, isotope -\> charge_2, isotope -\> charge_m, isotope -\> chiral, isotope -\> chiral_2, isotope -\> chiral_m, isotope -\> class, isotope -\> hydro, isotope -\> hydro_2, isotope -\> isotope, isotope -\> isotope_m, isotope_m -\> atom_oal, isotope_m -\> atom_oal_2, isotope_m -\> atom_oar, isotope_m -\> charge, isotope_m -\> charge_2, isotope_m -\> charge_m, isotope_m -\> chiral, isotope_m -\> chiral_2, isotope_m -\> chiral_m, isotope_m -\> class, isotope_m -\> hydro, isotope_m -\> hydro_2, isotope_m -\> isotope, isotope_m -\> isotope_m

6.   Cis/trans single character symbols should not be paired with some of the following classes of bond-related symbols:

   > bm_tbe_2 -\> l_ct, bm_tbe_2 -\> r_ct, l_ct -\> bm_ibe, l_ct -\> bm_ibi, l_ct -\> bm_tbe_2, l_ct -\> bm_tbi, l_ct -\> l_ct, l_ct -\> no_bond, l_ct -\> r_ct, r_ct -\> bm_ibe, r_ct -\> bm_ibi, r_ct -\> bm_tbe_2, r_ct -\> bm_tbi, r_ct -\> l_ct, r_ct -\> no_bond, r_ct -\> r_ct
    
So, at this point all the pairs of classes, which somehow contradict the SMILES rules, are enumerated. The list of these pairs could be used during the parsing procedure to stop it in cases, where there is no possibility to obtain the technically correct unambiguous output.
    
List of pairs could be updated if needed.

**Thus, it is possible to proceed with the basic parser, i.e. computer program which reads SMILES string from left to right, identifies meaningful symbols in the correct order and produces rather machine the human readable representation of the chemical graph.**

## Basic parser

### Technology

**Rust** (<https://doc.rust-lang.org/book/>) seems to be an appropriate option for the task:

-   modern language: less legacy things both in terms of the code and documentation, which makes the search for information and actual coding much easier for beginners; ecosystem seems to be pretty compact and well structured

-   performance is considered to be high

-   venerated memory safety (probably not a deal breaker in the context of this task, but people could consider the resulting tool more reliable because of that, yes, they definately will)

-   compiled code could be called from many different environments

-   code could be compiled to a WebAssembly module and used in browser and is, as I can see, language is being developed with this option in mind

-   build with the ideas of functional programming in mind (opinionated)

-   Rust to R, there is extendR (<https://cran.r-project.org/web/packages/rextendr/index.html>)

-   Rust to JavaScript (via Wasm), it surely works: <https://developer.mozilla.org/en-US/docs/WebAssembly/Guides/Rust_to_Wasm>

-   Rust to Python, there are ways to do that, for example: <https://github.com/pyo3/pyo3>

#### Some licensing

-   Rust (1.94.1), <https://rust-lang.org/> : [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) OR [MIT license](http://opensource.org/licenses/MIT)

-   wasm-pack, <https://www.npmjs.com/package/wasm-pack> : [MIT license](http://opensource.org/licenses/MIT) OR [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)

-   wasm-bindgen, <https://github.com/wasm-bindgen/wasm-bindgen> : [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) OR [MIT license](http://opensource.org/licenses/MIT)

-   serde (to manage data interchange using JSON), <https://crates.io/crates/serde> : [MIT license](http://opensource.org/licenses/MIT) OR [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)

### Data structure to store the results

It seems to be reasonable for the parser to produce an object consisting of the atom array (atoms and their properties) and bonds array (bonds, their properties, kinks to the corresponding atoms), also the input SMILES string maybe useful.

### Algorithm, draft

**input example:**

```         
smiles = "CCCCCCC"
```

**desired_output:**

```         
chem_struct = [
atoms: [{id, symbol, is_in_ring, ring_ids, is_in_branch, branch_id, is_aromatic, is_in_bracket, has_hs, isotopic_number, hs, charge, chirality, class}, ..., {...}],
bonds: [{id, atom_one, atom_two}, ... {...}],
symbols: ["C", ..., "C"],
input: original_smiles_string
]
```

**elements of state:**

```         
branches = empty LIFO data structure rings = empty set n_all = length(smiles) n_remain = length(smiles) subs = "" symb = "" symb_length = 0
```

**procedure:**

```         
while n > 0 do:
  // Get the substring
  if (n == n_all):
    subs = smiles[0 ... 5]
  else:
    subs = smiles[symb_length + 1 ... min(symb_length+6, n_all)]
  // Get the symbol (dummy function this time)
  symb = find_longest_symb(subs) symb_length = length(symb)
  // update an output (dummy function this time)
  update_output(output, state elements, symb)
  // update state (dummy function this time)
  update_state(state elements, symb)
  // decrease n
  n = n - symb_length
```

**return:**

```         
return output
```

### Code for parsing

