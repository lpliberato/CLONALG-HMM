using System.Collections.Generic;

namespace AIS.ClonalgPR.Models
{
    public class Result
    {
        public int MaximumIterations { get; set; }
        public double Average { get; set; }
        public double StandardDeviation { get; set; }
        public double Variance { get; set; }
        public double GreaterAffinity { get; set; }
        public double PercentHighAffinity { get; set; }
        public double PercentLowAffinity { get; set; }
        public double Time { get; set; }
    }
}
