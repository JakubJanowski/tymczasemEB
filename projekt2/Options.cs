using CommandLine;

namespace projekt1 {
    class Options {
        [Option('a', "sequenceA", Required = true, HelpText = "First sequence input file to be processed.")]
        public string Sequence1Path { get; set; }

        [Option('b', "sequenceB", Required = true, HelpText = "Second sequence input file to be processed.")]
        public string Sequence2Path { get; set; }

        [Option('s', "similiarityMatrix", Required = false, HelpText = "Similiarity matrix input file. Matrix size is 5x5 (or 21x21 if RNA option was specified).")]
        public string SimiliarityMatrixPath { get; set; }

        [Option('p', "penalty", Required = false, HelpText = "Set this flag to use penalty function.")]
        public bool PenaltyFunctionEnabled { get; set; } = false;

        [Option('v', "verbose", Required = false, HelpText = "Set this flag if edit distance and optimal similiarity matrices should be printed.")]
        public bool Verbose { get; set; } = false;
    }
}
