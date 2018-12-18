namespace Matcher {
    interface IMatcher {
        bool Verbose { get; set; }

        int ComputeEditDistance(out string[] matching);
        int ComputeSimiliarity(out string[] matching);
        int ComputeLocalSimiliarity(out string[] matching);
    }
}
