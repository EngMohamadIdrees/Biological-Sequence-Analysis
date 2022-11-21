# Biological-Sequence-Analysis
 Program that does the following:

 *Create a function that constructs the suffix array of any string using the naive slow method O(n^2 log n) time but in O(n) space.
  *Test it on the first 10,000  characters from the file: http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta Output the suffix array to a file.

* Create another function that constructs the suffix array of any string using the efficient predfix doubling method O(n (log n)^2) which calls the built-in C++ sort() function several times, and validate your implementation by making sure that the output is equivalent to the previous method when applied to the same 10,000 characters. Also, construct the suffix array of the whole genome in the previous file and save it to another text file. Report the construction time in seconds.

* Test the above two functions by calls from main() in the way described above.

* comments explaining the code before each code line, small test cases (at least 30 strings of length at least 10) that cover a lot of cases, and following the CodingStyle.pdf file.

