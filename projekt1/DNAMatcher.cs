using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DNA {
    public class DNAMatcher {
        private readonly int n;     // sequence1 length
        private readonly int m;     // sequence2 length
        private readonly byte[] u;  // sequence1
        private readonly byte[] w;  // sequence2
        private readonly int[,] d;  // distanceMatrix
        private readonly int[,] s;  // similiarityMatrix

        private int[,] D;           // edit distances matrix
        private int[,] S;           // similiarities matrix

        private bool similiarityMatrixExists = false;
        private bool isRNA;

        private enum MatchingType { Minimal, Maximal }

        private enum AminoAcids {
            STOP,
            Met,
            Phe,
            Leu,
            Ser,
            Tyr,
            Cys,
            Trp,
            Pro,
            His,
            Gln,
            Arg,
            Ile,
            Thr,
            Asn,
            Lys,
            Val,
            Ala,
            Asp,
            Glu,
            Gly,
            None
        }

        private static Dictionary<string, AminoAcids> aminoAcidCodons = new Dictionary<string, AminoAcids>() {
                {"UUU", AminoAcids.Phe},
                {"UUC", AminoAcids.Phe},
                {"UUA", AminoAcids.Leu},
                {"UUG", AminoAcids.Leu},
                {"UCU", AminoAcids.Ser},
                {"UCC", AminoAcids.Ser},
                {"UCA", AminoAcids.Ser},
                {"UCG", AminoAcids.Ser},
                {"UAU", AminoAcids.Tyr},
                {"UAC", AminoAcids.Tyr},
                {"UAA", AminoAcids.STOP},
                {"UAG", AminoAcids.STOP},
                {"UGU", AminoAcids.Cys},
                {"UGC", AminoAcids.Cys},
                {"UGA", AminoAcids.STOP},
                {"UGG", AminoAcids.Trp},
                {"CUU", AminoAcids.Leu},
                {"CUC", AminoAcids.Leu},
                {"CUA", AminoAcids.Leu},
                {"CUG", AminoAcids.Leu},
                {"CCU", AminoAcids.Pro},
                {"CCC", AminoAcids.Pro},
                {"CCA", AminoAcids.Pro},
                {"CCG", AminoAcids.Pro},
                {"CAU", AminoAcids.His},
                {"CAC", AminoAcids.His},
                {"CAA", AminoAcids.Gln},
                {"CAG", AminoAcids.Gln},
                {"CGU", AminoAcids.Arg},
                {"CGC", AminoAcids.Arg},
                {"CGA", AminoAcids.Arg},
                {"CGG", AminoAcids.Arg},
                {"AUU", AminoAcids.Ile},
                {"AUC", AminoAcids.Ile},
                {"AUA", AminoAcids.Ile},
                {"AUG", AminoAcids.Met},
                {"ACU", AminoAcids.Thr},
                {"ACC", AminoAcids.Thr},
                {"ACA", AminoAcids.Thr},
                {"ACG", AminoAcids.Thr},
                {"AAU", AminoAcids.Asn},
                {"AAC", AminoAcids.Asn},
                {"AAA", AminoAcids.Lys},
                {"AAG", AminoAcids.Lys},
                {"AGU", AminoAcids.Ser},
                {"AGC", AminoAcids.Ser},
                {"AGA", AminoAcids.Arg},
                {"AGG", AminoAcids.Arg},
                {"GUU", AminoAcids.Val},
                {"GUC", AminoAcids.Val},
                {"GUA", AminoAcids.Val},
                {"GUG", AminoAcids.Val},
                {"GCU", AminoAcids.Ala},
                {"GCC", AminoAcids.Ala},
                {"GCA", AminoAcids.Ala},
                {"GCG", AminoAcids.Ala},
                {"GAU", AminoAcids.Asp},
                {"GAC", AminoAcids.Asp},
                {"GAA", AminoAcids.Glu},
                {"GAG", AminoAcids.Glu},
                {"GGU", AminoAcids.Gly},
                {"GGC", AminoAcids.Gly},
                {"GGA", AminoAcids.Gly},
                {"GGG", AminoAcids.Gly}
        };

        public DNAMatcher(string sequence1, string sequence2, int[,] distanceMatrix, int[,] similiarityMatrix, bool isRNA = false) {
            if (isRNA) {
                u = sequence1.SplitInParts(3).Select(s => RNAToByte(s)).ToArray();
                w = sequence2.SplitInParts(3).Select(s => RNAToByte(s)).ToArray();
            }
            else {
                u = sequence1.Select(c => StringToByte(c)).ToArray();
                w = sequence2.Select(c => StringToByte(c)).ToArray();
            }
            n = u.Length;
            m = w.Length;
            d = distanceMatrix;
            s = similiarityMatrix;
            this.isRNA = isRNA;
        }

        public int ComputeEditDistance(out string[] matching) {
            D = new int[n + 1, m + 1];
            D[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    D[0, j] += d[StringToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    D[i, 0] += d[u[k], StringToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    D[i, j] = Math.Min(
                        D[i - 1, j - 1] + d[u[i - 1], w[j - 1]], Math.Min(
                        D[i, j - 1] + d[StringToByte('_'), w[j - 1]],
                        D[i - 1, j] + d[u[i - 1], StringToByte('_')]));
                }
            }

            PrintMatrix(D);

            matching = GetMatching(D, MatchingType.Minimal);

            return D[n, m];
        }

        public int ComputeSimiliarity(out string[] matching) {
            //if (isRNA)
            //    return ComputeRNASimiliarity(out matching);

            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                        S[0, j] += s[StringToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                        S[i, 0] += s[u[k], StringToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                        S[i, j] = Math.Max(
                            S[i - 1, j - 1] + s[u[i - 1], w[j - 1]], Math.Max(
                                S[i, j - 1] + s[StringToByte('_'), w[j - 1]],
                                S[i - 1, j] + s[u[i - 1], StringToByte('_')]));

                }
            }

            // TODO for RNA
            similiarityMatrixExists = true;
            PrintMatrix(S);

            // TODO for RNA
            matching = GetMatching(S, MatchingType.Maximal);

            return S[n, m];
        }

        private int ComputeRNASimiliarity(out string[] matching) {
            matching = null;
            return 0;
        }

        public int ComputeLocalMatching(out string[] matching) {
            if (!similiarityMatrixExists)
                ComputeSimiliarity(out matching);

            int x = 0;
            int y = 0;
            int max = 0;

            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= m; j++) {
                    if (max < S[i, j]) {
                        max = S[i, j];
                        x = i;
                        y = j;
                    }
                }
            }

            matching = GetMatching(S, MatchingType.Maximal, true, x, y);

            return max;
        }



        private byte StringToByte(char c)
        {
            byte DNAtoByte(char c1)
            {
                switch (c1)
                {
                    case 'A':
                        return 0;
                    case 'C':
                        return 1;
                    case 'G':
                        return 2;
                    case 'T':
                        return 3;
                }

                return 4;
            }

            byte RNAtoByte(char c1)
            {
                switch (c1)
                {
                    case '_': //STOP in RNA
                        return 0;
                }

                return 4;
            }

            if (isRNA)
                return RNAtoByte(c);

            return DNAtoByte(c);
        }

        private static char ByteToDNA(byte b) {
            switch (b) {
                case 0:
                    return 'A';
                case 1:
                    return 'C';
                case 2:
                    return 'G';
                case 3:
                    return 'T';
            }
            return '_';
        }

        private byte RNAToByte(string codon) {
            return (byte)aminoAcidCodons[codon];
        }

        private void PrintMatrix(int[,] M) {
            int padLenght = M.Cast<int>().Select(v => v.ToString().Length).Max() + 1;

            Console.Write("".PadLeft(2 * padLenght - 1));
            for (int i = 0; i < n; i++)
                Console.Write(ByteToDNA(u[i]).ToString().PadLeft(padLenght));
            Console.WriteLine();
            Console.Write("".PadLeft(padLenght - 1));

            for (int j = 0; j <= m; j++) {
                if (j > 0)
                    Console.Write(ByteToDNA(w[j - 1]) + " ");
                for (int i = 0; i <= n; i++)
                    Console.Write(M[i, j].ToString().PadLeft(padLenght));
                Console.WriteLine();
            }
        }


        private string[] GetMatching(int[,] matrix, MatchingType matchingType, bool stopAtNonPositive = false, int? startX = null, int? startY = null) {
            int x = startX ?? n;
            int y = startY ?? m;
            int lastX = x;
            int lastY = y;

            StringBuilder sequence1 = new StringBuilder();
            StringBuilder sequence2 = new StringBuilder();

            while ((x != 0 && y != 0) || (stopAtNonPositive && matrix[x, y] > 0)) {
                GetNextField(matrix, matchingType, ref x, ref y);

                if (x == lastX) {
                    sequence1.Insert(0, '_');
                    sequence2.Insert(0, ByteToDNA(w[y]));    // z nierówności trójkąta będzie działać
                }
                else if (y == lastY) {
                    sequence1.Insert(0, ByteToDNA(u[x]));
                    sequence2.Insert(0, '_');
                }
                else {
                    sequence1.Insert(0, ByteToDNA(u[x]));
                    sequence2.Insert(0, ByteToDNA(w[y]));
                }

                lastX = x;
                lastY = y;
            }

            return new string[] { sequence1.ToString(), sequence2.ToString() };
        }

        private static void GetNextField(int[,] matrix, MatchingType matchingType, ref int x, ref int y) {
            if (x == 0) {
                y--;
            }
            else if (y == 0) {
                x--;
            }
            else {
                int val = matrix[x, y];
                int next;

                switch (matchingType) {
                    case MatchingType.Maximal:
                        next = Math.Max(matrix[x - 1, y - 1], Math.Max(matrix[x, y - 1], matrix[x - 1, y]));
                        break;
                    case MatchingType.Minimal:
                    default:
                        next = Math.Min(matrix[x - 1, y - 1], Math.Min(matrix[x, y - 1], matrix[x - 1, y]));
                        break;
                }

                if (matrix[x - 1, y - 1] == next) {
                    x--;
                    y--;
                }
                else if (matrix[x, y - 1] == next) {
                    y--;
                }
                else if (matrix[x - 1, y] == next) {
                    x--;
                }
            }
        }
    }
}
