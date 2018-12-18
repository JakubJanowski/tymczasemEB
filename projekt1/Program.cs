using CommandLine;
using CommandLine.Text;
using Matcher;
using System;
using System.IO;
using System.Linq;

namespace projekt1 {
    class Program {
        private static bool _printHelp = false;

        static void Main(string[] args) {
            ParserResult<Options> p = Parser.Default.ParseArguments<Options>(args)
              .WithParsed(opts => Run(opts));

            if (_printHelp)
                Console.WriteLine(HelpText.AutoBuild(p));

            Console.ReadKey();

        }

        private static void Run(Options options) {

            if (!ValidateOptions(options)) {
                _printHelp = true;
                return;
            }

            string sequence1 = File.ReadAllText(options.Sequence1Path);
            string sequence2 = File.ReadAllText(options.Sequence2Path);

            int[,] distanceMatrix = null;
            int[,] similiarityMatrix = null;
            IMatcher matcher;

            if (options.IsRNA) {
                if (!string.IsNullOrWhiteSpace(options.DistanceMatrixPath))
                    distanceMatrix = GetMatrix(options.DistanceMatrixPath, 21);
                if (!string.IsNullOrWhiteSpace(options.SimiliarityMatrixPath))
                    similiarityMatrix = GetMatrix(options.SimiliarityMatrixPath, 21);
                matcher = new RNAMatcher(sequence1, sequence2, distanceMatrix, similiarityMatrix);
            }
            else {
                if (!string.IsNullOrWhiteSpace(options.DistanceMatrixPath))
                    distanceMatrix = GetMatrix(options.DistanceMatrixPath, 5);
                if (!string.IsNullOrWhiteSpace(options.SimiliarityMatrixPath))
                    similiarityMatrix = GetMatrix(options.SimiliarityMatrixPath, 5);
                matcher = new DNAMatcher(sequence1, sequence2, distanceMatrix, similiarityMatrix);
            }

            matcher.Verbose = options.Verbose;

            string[] matching;

            if (options.Verbose) {
                Console.WriteLine("Input sequences:");
                Console.WriteLine(sequence1);
                Console.WriteLine(sequence2);
                Console.WriteLine();
            }

            if (distanceMatrix != null) {
                int editDistance = matcher.ComputeEditDistance(out matching);
                Console.WriteLine($"Optimal edit distance: {editDistance}\n");
                PrintMatching(matching);
            }

            if (similiarityMatrix != null) {
                int similiarity = matcher.ComputeSimiliarity(out matching);
                Console.WriteLine($"Optimal global similiarity: {similiarity}\n");
                PrintMatching(matching);

                int localSimiliarity = matcher.ComputeLocalSimiliarity(out matching);
                Console.WriteLine($"Optimal local similiarity: {localSimiliarity}\n");
                PrintMatching(matching);
            }
        }

        private static bool ValidateOptions(Options options) {
            if (string.IsNullOrWhiteSpace(options.DistanceMatrixPath) && string.IsNullOrWhiteSpace(options.SimiliarityMatrixPath)) {
                Console.WriteLine("At least one of distance matrix or similiarity matrix intput files should be specified.");
                return false;
            }
            return true;
        }

        private static int[,] GetMatrix(string filepath, int matrixSize) {
            string[] matrixLines;
            int[,] matrix;

            try {
                 matrixLines = File.ReadAllLines(filepath);
            } catch (Exception) {
                Console.WriteLine($"[Error] Could not read file {filepath}\n");
                return null;
            }

            try {
                matrix = new int[matrixSize, matrixSize];
                for (int y = 0; y < matrixSize; y++) {
                    int[] row = matrixLines[y].Split(new char[] { ' ' }).Select(s => int.Parse(s)).ToArray();
                    for (int x = 0; x < matrixSize; x++)
                        matrix[y, x] = row[x];
                }
            } catch(Exception) {
                Console.WriteLine($"[Error] File {filepath} should contain {matrixSize}x{matrixSize} matrix. Each row should be in different line and each cell in row should be separated by space.\n");
                return null;
            }

            return matrix;

        }

        private static void PrintMatching(string[] matching) {
            if (matching != null) {
                Console.WriteLine("Matching:");
                Console.WriteLine(matching[0]);
                Console.WriteLine(matching[1]);
            }
            Console.WriteLine();
            Console.WriteLine();
        }
    }
}
