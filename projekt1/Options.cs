using CommandLine;

namespace projekt1 {
    class Options {
        [Option('a', "sequenceA", Required = true, HelpText = "First sequence input file to be processed.")]
        public string Sequence1Path { get; set; }

        [Option('b', "sequenceB", Required = true, HelpText = "Second sequence input file to be processed.")]
        public string Sequence2Path { get; set; }

        [Option('d', "distanceMatrix", Required = false, HelpText = "Distance matrix input file. Matrix size is 5x5 (or 21x21 if RNA option was specified).")]
        public string DistanceMatrixPath { get; set; }

        [Option('s', "similiarityMatrix", Required = false, HelpText = "Similiarity matrix input file. Matrix size is 5x5 (or 21x21 if RNA option was specified).")]
        public string SimiliarityMatrixPath { get; set; }

        [Option('r', "rna", Required = false, HelpText = "Set this flag if input files contain RNA sequence. If not set, the sequences are treated as DNA sequence.")]
        public bool IsRNA { get; set; } = false;
    }
}
