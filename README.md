# IBMC_LSFBD_smiles_parser

Parsing SMILES is useful, since this procedure provides deeper understanding of the possibilities and limitations of SMILES format.

For the SMILES parser to be testable and verifiable it is important to follow the accepted practices of parsing given the specification of the input.

In this repo the OpenSMILES specification (http://opensmiles.org/opensmiles.html) will be explored to create the pretty normal bottom-up SMILES parser for the https://github.com/pavelVPo/IBMC_LSFBD_webApp_input

|1| classify the characters, which are OpenSMILES-valid.
|2| generate all the possible pairs of characters

Results:

|1|

characters could be classified as:

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

1. 		character 		-> bracket atom symbol										-> atom symbol 						-> atom_H 		 								OK
2. 		character 		-> bracket atom symbol start								-> start atom symbol 				-> atom_symbol_start 							OK
3. 		character 		-> bracket atom symbol end									-> end atom symbol 					-> atom_symbol_end 								OK
4. 		character 		-> organic atom symbol 										-> whole atom symbol 				-> atom 										OK
5. 		character 		-> aromatic organic atom symbol								-> whole atom symbol 				-> atom 										OK
6. 		character 		-> wildcard atom 											-> * 								-> * 											OK
7. 		character 		-> bracket atom_symbol_start								-> start atom 						-> [ 											OK
8. 		character 		-> bracket atom_symbol_end									-> end atom 						-> ] 											OK
9. 		character 		-> whole chirality											-> whole chirality					-> ▲ 											OK
10. 	character 		-> start chirality											-> start chirality					-> ▲. 											OK
11. 	character 		-> mid chirality											-> mid chirality					-> .▲. 											OK
12. 	character 		-> end_chirality											-> end chirality 					-> .▲											OK
13.		character 		-> whole hydrogen add in brackets		 					-> H 								-> h 											OK
14.		character 		-> start hydrogen add in brackets		 					-> H. 								-> h. 											OK
15.		character 		-> mid hydrogen add in brackets		 						-> .H. 								-> .h. 											OK
16.		character 		-> end hydrogen add in brackets			 					-> .H 								-> .h 											OK
17.		charcter 		-> whole charge 				 							-> ± 								-> ± 											OK
18.		charcter 		-> start charge 				 							-> ± 								-> ±. 											OL
19.		charcter 		-> mid charge 	 				 							-> ± 								-> .±. 											OK
20.		charcter 		-> end charge 	 				 							-> ± 								-> .± 											OK
21.		character 		-> start atom class 										-> :. 								-> :. 											OK
22.		character 		-> mid atom class 											-> .:. 								-> .:. 											OK
23.		character 		-> end atom class 	 										-> .: 								-> .: 											OK
24. 	character 		-> bond 													-> - 								-> - 											OK
25. 	character 		-> component delimeter 										-> . 								-> .											OK
26. 	character 		-> branch start 											-> (								-> ( 											OK
27. 	character 		-> branch end 												-> )								-> ) 											OK
28. 	character 		-> whole ring 												-> % 								-> % 											OK
29. 	character 		-> start ring 												-> %. 								-> %. 											OK
30. 	character 		-> mid ring 												-> .%. 								-> .%. 											OK
31. 	character 		-> end ring 												-> .% 								-> .% 											OK
32.		character 		-> whole isotope 											-> ^ 								-> ^ 											OK
33.		character 		-> start isotope 											-> ^. 								-> ^. 											OK
34.		character 		-> mid isotope 												-> .^. 								-> .^. 											OK
35.		character 		-> end isotope  											-> .^ 								-> .^ 											OK

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



|2| 

All the possible pairs of character' types

/////////////////////////////////////////////////////////////////////////

1. 			atom_H							]
2. 			atom_H							▲
3. 			atom_H							▲.
4. 			atom_H							h
5. 			atom_H							h.
6. 			atom_H							±
7. 			atom_H							±.
8. 			atom_H							:.

9. 			atom_symbol_start				atom_symbol_end

10. 		atom_symbol_end					]
11. 		atom_symbol_end					▲
12. 		atom_symbol_end					▲.
13. 		atom_symbol_end					h
14. 		atom_symbol_end					h.
15. 		atom_symbol_end					±
16. 		atom_symbol_end					±.
17. 		atom_symbol_end					:.

18. 		atom 							atom
19. 		atom 							*
20. 		atom 							[
21. 		atom							]
22. 		atom							▲
23. 		atom							▲.
24. 		atom							h
25. 		atom							h.
26. 		atom							±
27. 		atom							±.
28. 		atom							:.
29. 		atom							-
30. 		atom							.
31. 		atom							(
32. 		atom							)
33. 		atom							%
34. 		atom							%.

35. 		* 	 							atom
36. 		* 	 							*
37. 		* 	 							[
38. 		* 								]
39. 		* 								▲
40. 		* 								▲.
41. 		* 								h
42. 		* 								h.
43. 		* 								±
44. 		* 								±.
45. 		* 								:.
46. 		* 								-
47. 		* 								.
48. 		* 								(
49. 		* 								)
50. 		* 								%
51. 		* 								%.

52. 		[ 								atom_H
53. 		[ 								atom_symbol_start
54. 		[ 								atom
55. 		[ 								*
56. 		[ 								^
57. 		[ 								^.

58. 		] 								atom
59. 		] 								*
60. 		] 								[
61. 		] 								=
62. 		] 								.
63. 		] 								(
64. 		] 								)
65. 		] 								*
66. 		] 								%
67. 		] 								%.
68. 		] 								*

69. 		▲ 								]
70. 		▲ 								▲.
71. 		▲ 								h
72. 		▲ 								h.
73. 		▲ 								±
74. 		▲ 								±.
75. 		▲ 								:.
76. 		▲ 								±

77. 		▲. 								.▲.
78. 		▲. 								.▲

79. 		.▲. 							.▲
80. 		.▲. 							.▲

81. 		h 								±
82. 		h 								±.
83. 		h 								:.

84. 		h. 								.h.
85. 		h. 								.h

86. 		.h 								±
87. 		.h 								±.
88. 		.h 								:.
89. 		.h 								]

90. 		± 								:.
91. 		± 								]

92. 		±. 								.±.
93. 		±. 								.±

94. 		±. 								.±.
95. 		±. 								.±

96. 		.±. 							.±

97. 		.± 								:.
98. 		.± 								]

99. 		:. 								.:.
100. 		:. 								.:

101. 		.:. 							.:.
102. 		.:. 							.:

103. 		.: 								]

104. 		- 								atom
105. 		- 								*
106. 		- 								[
107. 		- 								)
108. 		- 								%
109. 		- 								%.

110. 		. 								atom
111. 		. 								*
112. 		. 								[

113. 		( 								atom
114. 		( 								*
115. 		( 								[
115. 		( 								-
116. 		( 								(

117. 		) 								atom
118. 		) 								*
119. 		) 								[
120. 		) 								-
121. 		) 								.
122. 		) 								(
123. 		) 								)
124. 		) 								%
125. 		) 								%.

126. 		% 								atom
127. 		% 								*
128. 		% 								[
129. 		% 								(
130. 		% 								)
131. 		% 								%
132. 		% 								%.

133. 		%. 								.%.

134. 		.%. 							.%

135. 		%. 								.%.

136. 		.%. 							.%

137. 		.% 								atom
138. 		.% 								*
139. 		.% 								[
140. 		.% 								(
141. 		.% 								)
142. 		.% 								%
143. 		.% 								%.

144. 		^ 								atom_H
145. 		^ 								atom_symbol_start
146. 		^ 								atom
147. 		^ 								*

148. 		^. 								.^.
149. 		^. 								.^

150. 		.^. 							.^

151. 		.^ 								atom_H
152. 		.^ 								atom_symbol_start
153. 		.^ 								atom
154. 		.^ 								*

/////////////////////////////////////////////////////////////////////////

