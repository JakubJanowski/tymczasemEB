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
        private readonly string sequence1;
        private readonly string sequence2;

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
                u = sequence1.Select(c => CharToByte(c)).ToArray();
                w = sequence2.Select(c => CharToByte(c)).ToArray();
            }
            n = u.Length;
            m = w.Length;
            d = distanceMatrix;
            s = similiarityMatrix;
            this.isRNA = isRNA;
            this.sequence1 = sequence1;
            this.sequence2 = sequence2;
        }

        public int ComputeEditDistance(out string[] matching) {
            D = new int[n + 1, m + 1];
            D[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    D[0, j] += d[CharToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    D[i, 0] += d[u[k], CharToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    D[i, j] = Math.Min(
                        D[i - 1, j - 1] + d[u[i - 1], w[j - 1]], Math.Min(
                        D[i, j - 1] + d[CharToByte('_'), w[j - 1]],
                        D[i - 1, j] + d[u[i - 1], CharToByte('_')]));
                }
            }

            PrintMatrix(D);

            matching = GetMatching(D, MatchingType.Minimal);

            return D[n, m];
        }

        public int ComputeSimiliarity(out string[] matching) {
            if (!isRNA)
                return ComputeSequenceSimiliarity(out matching, u, w);

            //todo stop jest spacją

            bool found = false;
            int bestResult = int.MinValue;
            Sequence uBestSubsequence = null;
            Sequence wBestSubsequence = null;
            int[,] bestS = null;
            Sequence[] uSequences = u.FindAminoAcidSequences((byte)AminoAcids.Met, (byte)AminoAcids.STOP).ToArray();
            Sequence[] wSequences = w.FindAminoAcidSequences((byte)AminoAcids.Met, (byte)AminoAcids.STOP).ToArray();
            foreach (Sequence uSequence in uSequences) {
                foreach (Sequence wSequence in wSequences) {
                    int result = ComputeSequenceSimiliarityRNA(uSequence.data, wSequence.data);
                    if(result > bestResult) {
                        found = true;
                        bestResult = result;
                        uBestSubsequence = uSequence;
                        wBestSubsequence = wSequence;
                        bestS = S.Clone() as int[,];
                    }
                }
            }

            if(!found) {
                matching = null;
                return int.MinValue;
            }

            matching = GetMatchingRNA(bestS, MatchingType.Maximal, uBestSubsequence, wBestSubsequence);
            return bestResult;
        }

        private int ComputeSequenceSimiliarity(out string[] matching, byte[] u, byte[] w) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    S[0, j] += s[CharToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    S[i, 0] += s[u[k], CharToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    S[i, j] = Math.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]], Math.Max(
                            S[i, j - 1] + s[CharToByte('_'), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], CharToByte('_')]));

                }
            }

            similiarityMatrixExists = true;
            //TODO print or not for RNA?
            PrintMatrix(S);
            
            matching = GetMatching(S, MatchingType.Maximal);

            return S[n, m];
        }

        private int ComputeSequenceSimiliarityRNA(byte[] u, byte[] w) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    S[0, j] += s[CharToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    S[i, 0] += s[u[k], CharToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    S[i, j] = Math.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]], Math.Max(
                            S[i, j - 1] + s[CharToByte('_'), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], CharToByte('_')]));

                }
            }

            similiarityMatrixExists = true;

            return S[n, m];
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

            if (isRNA)
                matching = null;// GetMatchingRNA(S, MatchingType.Maximal, true, x, y);
            else
                matching = GetMatching(S, MatchingType.Maximal, true, x, y);

            return max;
        }



        private byte CharToByte(char c) {
            byte DNAtoByte(char c1) {
                switch (c1) {
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

            byte RNAtoByte(char c1) {
                if (c1 == '_')
                    return (byte)AminoAcids.STOP;

                throw new Exception();
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

            StringBuilder outSequence1 = new StringBuilder();
            StringBuilder outSequence2 = new StringBuilder();

            while ((x != 0 || y != 0) || (stopAtNonPositive && matrix[x, y] > 0)) {
                GetNextField(matrix, matchingType, ref x, ref y);

                if (x == lastX) {
                    outSequence1.Insert(0, '_');
                    outSequence2.Insert(0, ByteToDNA(w[y]));    // z nierówności trójkąta będzie działać
                }
                else if (y == lastY) {
                    outSequence1.Insert(0, ByteToDNA(u[x]));
                    outSequence2.Insert(0, '_');
                }
                else {
                    outSequence1.Insert(0, ByteToDNA(u[x]));
                    outSequence2.Insert(0, ByteToDNA(w[y]));
                }

                lastX = x;
                lastY = y;
            }

            return new string[] { outSequence1.ToString(), outSequence2.ToString() };
        }

        private string[] GetMatchingRNA(int[,] matrix, MatchingType matchingType, Sequence uSequence, Sequence wSequence) {
            int x = uSequence.data.Length;
            int y = wSequence.data.Length;
            int lastX = x;
            int lastY = y;

            int stopSeq1 = (x + uSequence.index) * 3; // 3 chars for amino
            int stopSeq2 = (y + wSequence.index) * 3;


            StringBuilder outSequence1 = new StringBuilder();
            StringBuilder outSequence2 = new StringBuilder();

            while ((x != 0 || y != 0) || matrix[x, y] > 0) {
                GetNextField(matrix, matchingType, ref x, ref y);

                if (x == lastX) {
                    stopSeq2 -= 3;
                    outSequence1.Insert(0, "___");
                    outSequence2.Insert(0, sequence2.Substring(stopSeq2, 3));
                }
                else if (y == lastY) {
                    stopSeq1 -= 3;
                    outSequence1.Insert(0, sequence1.Substring(stopSeq1, 3));
                    outSequence2.Insert(0, "___");
                }
                else {
                    stopSeq1 -= 3;
                    stopSeq2 -= 3;
                    outSequence1.Insert(0, sequence1.Substring(stopSeq1, 3));
                    outSequence2.Insert(0, sequence2.Substring(stopSeq2, 3));
                }

                lastX = x;
                lastY = y;
            }

            return new string[] { outSequence1.ToString(), outSequence2.ToString() };
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
