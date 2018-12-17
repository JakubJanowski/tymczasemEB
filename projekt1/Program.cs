using CommandLine;
using CommandLine.Text;
using DNA;
using System;
using System.IO;
using System.Linq;
using System.Runtime.Remoting.Lifetime;

namespace projekt1 {
    class Program {
        private static bool _printHelp = false;

        static void Main(string[] args) {
            ParserResult<Options> p = Parser.Default.ParseArguments<Options>(args)
              .WithParsed(opts => Run(opts));
            
            if(_printHelp)
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

            if (!string.IsNullOrWhiteSpace(options.DistanceMatrixPath))
                distanceMatrix = GetMatrix(options.DistanceMatrixPath, 5);

            if (!string.IsNullOrWhiteSpace(options.SimiliarityMatrixPath))
                similiarityMatrix = options.IsRNA
                    ? GetMatrix(options.SimiliarityMatrixPath, 21)
                    : GetMatrix(options.SimiliarityMatrixPath, 5);

            DNAMatcher matcher = new DNAMatcher(sequence1, sequence2, distanceMatrix, similiarityMatrix, options.IsRNA);
            string[] matching;

            if (distanceMatrix != null && !options.IsRNA) {
                int editDistance = matcher.ComputeEditDistance(out matching);
                PrintMatching(matching);
                Console.WriteLine($"Optimal edit distance: {editDistance}\n");
            }

            if (similiarityMatrix != null) {
                int similiarity = matcher.ComputeSimiliarity(out matching);
                PrintMatching(matching);
                Console.WriteLine($"Optimal global similiarity: {similiarity}\n");

                if (!options.IsRNA) {
                    int localSimiliarity = matcher.ComputeLocalMatching(out matching);
                    PrintMatching(matching);
                    Console.WriteLine($"Optimal local similiarity: {localSimiliarity}\n");
                }
            }

            Console.ReadKey();
        }

        private static bool ValidateOptions(Options options) {
            if (string.IsNullOrWhiteSpace(options.DistanceMatrixPath) && string.IsNullOrWhiteSpace(options.SimiliarityMatrixPath)) {
                Console.WriteLine("At least one of distance matrix or similiarity matrix intput files should be specified.");
                return false;
            }
            return true;
        }

        private static int[,] GetMatrix(string filepath, int matrixSize) {
            string[] matrixLines = File.ReadAllLines(filepath);
            int[,] matrix = new int[matrixSize, matrixSize];

            for (int y = 0; y < matrixSize; y++) {
                int[] row = matrixLines[y].Split(new char[] { ' ' }).Select(s => int.Parse(s)).ToArray();
                for (int x = 0; x < matrixSize; x++)
                    matrix[y, x] = row[x];
            }
            return matrix;
        }

        private static void PrintMatching(string[] matching) {
            Console.WriteLine();
            if (matching != null) {
                Console.WriteLine(matching[0]);
                Console.WriteLine(matching[1]);
            }
            Console.WriteLine();
        }
    }
}
