using System.Collections.Generic;

namespace AIS.ClonalgPR.Models
{
    public class State
    {
        public StateEnum Name { get; set; }
        public Dictionary<char, double> EmissionProbabilities { get; set; }
        public Dictionary<StateEnum, double> TransitionProbabilities { get; set; }

        public State()
        {
            EmissionProbabilities = new Dictionary<char, double>();
            TransitionProbabilities = new Dictionary<StateEnum, double>();
        }
    }
}
