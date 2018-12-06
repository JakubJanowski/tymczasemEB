# tymczasemEB

A project realised for Elements of Bioinformatics subject. <br />
The program computes an optimal matching of 2 DNA sequences based either on edit distance or similiarity. <br />
The sequences and edit distance or similiarity criteria are given in .txt files in /data directory. <br />
 <br />
The edit distance and similiarity criteria are 5x5 matrices with row and column values corresponding to sequence letter change or indel ('_') as shown: <br />
 <br />
   A  C  G  T  _ <br />
A  .  .  .  .  . <br />
C  .  .  .  .  . <br />
G  .  .  .  .  . <br />
T  .  .  .  .  . <br />
_  .  .  .  .  . <br />
 <br />
It is also possible to translate given DNA sequence to RNA sequence or aminoacid sequence.
