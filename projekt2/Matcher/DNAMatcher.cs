using Models;
using System;
using System.Linq;
using System.Text;
using Utilities;

namespace Matcher
{
    public class DNAMatcher
    {
        private readonly byte[] uGlobal;  // sequence1
        private readonly byte[] wGlobal;  // sequence2
        private readonly int[,] s;  // similiarityMatrix

        private int[,] S;           // similiarities matrix

        public bool Verbose { get; set; } = false;

        public DNAMatcher(string sequence1, string sequence2, int[,] similiarityMatrix)
        {
            uGlobal  = sequence1.Select(c => DNAToByte(c)).ToArray();
            wGlobal = sequence2.Select(c => DNAToByte(c)).ToArray();
            s = similiarityMatrix;
        }

        protected int PenaltyFunction(int n)
        {
            return -n * n;
        }

        public int ComputeSimiliarityWithPenalty(out string[] matching)
        {
            int n = uGlobal.Length;
            int m = wGlobal.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    S[0, j] += s[DNAToByte('_'), wGlobal[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    S[i, 0] += s[uGlobal[k], DNAToByte('_')];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(
                        S[i - 1, j - 1] + s[uGlobal[i - 1], wGlobal[j - 1]],
                            S[i, j - 1] + s[DNAToByte('_'), wGlobal[j - 1]],
                            S[i - 1, j] + s[uGlobal[i - 1], DNAToByte('_')]);

            if (Verbose)
                PrintMatrix(uGlobal, wGlobal, S);

            matching = GetMatching(uGlobal, wGlobal, S, s);

            return S[n, m];
        }

        public int Hirschberg(out string[] matching)
        {
            return Hirschberg(uGlobal, wGlobal, out matching);
        }

        public int Hirschberg(byte[] u, byte[] w, out string[] matching)
        {
            byte[] Z = new byte[u.Length + w.Length];
            byte[] W = new byte[u.Length + w.Length];
            if (u.Length == 0)
            {
                StringBuilder matchingBuilder1 = new StringBuilder();
                StringBuilder matchingBuilder2 = new StringBuilder();
                for (int i = 0; i < w.Length; i++)
                {
                    matchingBuilder1.Append(DNAToByte('_'));
                    matchingBuilder2.Append(w[i]);
                }
                matching = new string[] { matchingBuilder1.ToString(), matchingBuilder2.ToString() };
                return 0;
            }
            else if (w.Length == 0)
            {
                StringBuilder matchingBuilder1 = new StringBuilder();
                StringBuilder matchingBuilder2 = new StringBuilder();
                for (int i = 0; i < u.Length; i++)
                {
                    matchingBuilder1.Append(u[i]);
                    matchingBuilder2.Append(DNAToByte('_'));
                }
                matching = new string[] { matchingBuilder1.ToString(), matchingBuilder2.ToString() };
                return 0;
            }
            else if (u.Length == 1 || w.Length == 1)
            {
                return NeedlemanWunsch(u, w, out matching);
            }
            else
            {
                int xlen = u.Length;
                int xmid = u.Length / 2;
                int ylen = w.Length;

                int[] ScoreL = NWScore(u.SubArray(0, xmid), w);
                int[] ScoreR = NWScore(u.SubArray(xmid, xlen - xmid).Reverse().ToArray(), w.Reverse().ToArray());
                for (int i = 0; i < ScoreL.Length; i++)
                {
                    ScoreL[i] += ScoreR[ScoreL.Length - i - 1];
                }
                int ymid = Array.IndexOf(ScoreL, ScoreL.Max());

                string[] matching1, matching2;

                int similiarity = Hirschberg(u.SubArray(0, xmid), w.SubArray(0, ymid), out matching1) + Hirschberg(u.SubArray(xmid, xlen - xmid), w.SubArray(ymid, ylen - ymid), out matching2);

                matching = new string[] { matching1[0] + matching2[0], matching1[1] + matching2[1] };
                return similiarity;
            }
        }

        private int NeedlemanWunsch(byte[] u, byte[] w, out string[] matching)
        {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    S[0, j] += s[DNAToByte('_'), w[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    S[i, 0] += s[u[k], DNAToByte('_')];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]],
                            S[i, j - 1] + s[DNAToByte('_'), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], DNAToByte('_')]);

            if (Verbose)
                ;// PrintMatrix(S);

            matching = GetMatching(u, w, S, s);

            return S[n, m];
        }

        private int[] NWScore(byte[] u, byte[] w)
        {
            int[,] Score = new int[2, w.Length + 1];
            Score[0, 0] = 0;

            for (int i = 0; i < w.Length; i++)
            {
                Score[0, i + 1] = Score[0, i] + s[DNAToByte('_'), w[i]];
            }

            for (int i = 0; i < u.Length; i++)
            {
                Score[1, 0] = Score[0, 0] + s[u[i], DNAToByte('_')];
                for (int j = 0; j < w.Length; j++)
                {
                    int scoreSub = Score[0, j] + s[u[i], w[j]];
                    int scoreDel = Score[0, j + 1] + s[u[i], DNAToByte('_')];
                    int scoreIns = Score[1, j] + s[DNAToByte('_'), w[j]];
                    Score[1, j + 1] = MathUtils.Max(scoreSub, scoreDel, scoreIns);
                }

                for (int j = 0; j <= w.Length; j++)
                {
                    Score[0, j] = Score[1, j];
                }
            }

            int[] lastLine = new int[w.Length + 1];
            for (int i = 0; i < w.Length; i++)
            {
                lastLine[i] = Score[1, i];
            }
            return lastLine;
        }

        private string[] GetMatching(byte[] u, byte[] w, int[,] matrix, int[,] costMatrix, bool stopAtNonPositive = false, Field startField = null)
        {
            Field field = startField ?? new Field(u.Length, w.Length);
            Field lastField = field.Clone();

            StringBuilder outSequence1 = new StringBuilder();
            StringBuilder outSequence2 = new StringBuilder();

            while ((field.x != 0 || field.y != 0) && (!stopAtNonPositive || matrix[field.x, field.y] > 0))
            {
                field = GetNextField(u, w, matrix, costMatrix, field);

                if (field.x == lastField.x)
                {
                    outSequence1.Insert(0, '_');
                    outSequence2.Insert(0, ByteToDNA(w[field.y]));    // z nierówności trójkąta będzie działać
                }
                else if (field.y == lastField.y)
                {
                    outSequence1.Insert(0, ByteToDNA(u[field.x]));
                    outSequence2.Insert(0, '_');
                }
                else
                {
                    outSequence1.Insert(0, ByteToDNA(u[field.x]));
                    outSequence2.Insert(0, ByteToDNA(w[field.y]));
                }

                lastField.Copy(field);
            }

            return new string[] { outSequence1.ToString(), outSequence2.ToString() };
        }

        private Field GetNextField(byte[] u, byte[] w, int[,] matrix, int[,] costMatrix, Field startField)
        {
            Field nextField = null;
            if (startField.x == 0)
            {
                nextField = startField;
                nextField.y--;
            }
            else if (startField.y == 0)
            {
                nextField = startField;
                nextField.x--;
            }
            else
            {
                int value = matrix[startField.x, startField.y];
                int nextValue;
                int bestNextValue = int.MinValue;

                Field field = new Field(startField.x - 1, startField.y - 1);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[u[field.x], w[field.y]] + nextValue == value && nextValue > bestNextValue)
                {   // check if the step is legal
                    bestNextValue = nextValue;
                    nextField = field.Clone();
                }
                field = new Field(startField.x, startField.y - 1);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[DNAToByte('_'), w[field.y]] + nextValue == value && nextValue > bestNextValue)
                {
                    bestNextValue = nextValue;
                    nextField = field.Clone();
                }
                field = new Field(startField.x - 1, startField.y);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[u[field.x], DNAToByte('_')] + nextValue == value && nextValue > bestNextValue)
                {
                    bestNextValue = nextValue;
                    nextField = field;
                }
            }

            return nextField;
        }

        private static byte DNAToByte(char c)
        {
            switch (c)
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

        private static char ByteToDNA(byte b)
        {
            switch (b)
            {
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

        private void PrintMatrix(byte[] u, byte[] w, int[,] M)
        {
            int padLenght = M.Cast<int>().Select(v => v.ToString().Length).Max() + 1;

            Console.Write("  " + "".PadLeft(padLenght));
            for (int i = 0; i < u.Length; i++)
                Console.Write(ByteToDNA(u[i]).ToString().PadLeft(padLenght));
            Console.WriteLine();
            Console.Write("  ");

            for (int j = 0; j <= w.Length; j++)
            {
                if (j > 0)
                    Console.Write(ByteToDNA(w[j - 1]) + " ");
                for (int i = 0; i <= u.Length; i++)
                    Console.Write(M[i, j].ToString().PadLeft(padLenght));
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
