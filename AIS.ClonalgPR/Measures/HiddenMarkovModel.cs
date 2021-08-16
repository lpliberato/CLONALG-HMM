using AIS.ClonalgPR.Models;
using System;
using System.Collections.Generic;
using System.Linq;

namespace AIS.ClonalgPR.Measures
{
    public class HiddenMarkovModel : IDistance
    {
        private List<char[]> Observations { get; set; }
        private List<State> States { get; set; }
        private TypeBioSequence TypeBioSequence;

        public HiddenMarkovModel(List<char[]> observations, TypeBioSequence typeBioSequence)
        {
            Observations = observations;
            TypeBioSequence = typeBioSequence;
            States = new List<State>();
        }

        public void Train()
        {
            CreateStates();
            CreateProbabilities();
            CreateTransitions();
        }

        private int GetObservationsAmount()
        {
            return Observations.Count;
        }

        private int GetObservationSize()
        {
            return Observations[0].Length;
        }

        private void AddMatchState()
        {
            States.Add(new State()
            {
                Name = StateEnum.Match,
                EmissionProbabilities = InitEmissionProbabilities(),
                TransitionProbabilities = InitTransitionProbabilities()
            });
        }

        private void AddInsertState()
        {
            States.Add(new State()
            {
                Name = StateEnum.Insert,
                EmissionProbabilities = InitEmissionProbabilities(),
                TransitionProbabilities = InitTransitionProbabilities()
            });
        }

        private void AddDeleteState()
        {
            States.Add(new State()
            {
                Name = StateEnum.Delete,
                EmissionProbabilities = InitEmissionProbabilities(),
                TransitionProbabilities = InitTransitionProbabilities()
            });
        }

        private bool IsLastIndexObservation(int index)
        {
            var observationSize = Observations[0].Length;
            return index >= observationSize;
        }

        private void CreateStates()
        {
            var observationSize = GetObservationSize();

            for (int index = 0; index < observationSize; index++)
            {
                var nextIndex = index + 1;
                var observationsAmount = GetObservationsAmount();
                var symbolsAmount = GetAmountOfMatches(index);
                var nextSymbolsAmount = GetAmountOfMatches(nextIndex);

                if (symbolsAmount == observationsAmount)
                    AddMatchState();
                else if (symbolsAmount > 0 && (nextSymbolsAmount == observationsAmount || nextSymbolsAmount == 0))
                    AddInsertState();
            }
        }

        private void UpdateInsertStateCounterForProbabilities(int index, ref int insertStateCounter)
        {
            var observationsAmount = GetObservationsAmount();
            var symbolsAmount = GetAmountOfMatches(index);

            if (symbolsAmount != observationsAmount)
            {
                var currentIdex = index + 1;
                if (IsLastIndexObservation(currentIdex)) return;

                insertStateCounter++;
                var nextSymbolsAmount = GetAmountOfMatches(currentIdex);
                if (nextSymbolsAmount == observationsAmount)
                    insertStateCounter--;
            }
        }

        private char GetSymbol(int lineIndex, int columnIndex)
        {
            char[] sequence = Observations[lineIndex];
            return sequence[columnIndex];
        }

        private State GetCurrentState(int index)
        {
            return States[index];
        }

        private bool IsGap(char symbol)
        {
            return Constants.Gaps.Contains(symbol);
        }

        private int IncreaseObservationAmount(char symbol, Dictionary<char, int> symbols)
        {
            return symbols[symbol] + 1;
        }

        private double CalculateProbabilityBySymbol(int insertStateCounterColumn, char symbol, Dictionary<char, int> symbols)
        {
            var observationsAmount = GetObservationsAmount();
            return (double)symbols[symbol] / observationsAmount / (insertStateCounterColumn == 0 ? 1d : Convert.ToDouble(insertStateCounterColumn));
        }

        private void UpdateTheAmountOfSymbols(char symbol, ref Dictionary<char, int> symbols)
        {
            symbols[symbol] = IncreaseObservationAmount(symbol, symbols);
        }

        private void UpdateSymbolProbability(int insertStateCounterColumn, char symbol, Dictionary<char, int> symbols, ref Dictionary<char, double> probabilities)
        {
            probabilities[symbol] = CalculateProbabilityBySymbol(insertStateCounterColumn, symbol, symbols);
        }

        private void UpdateProbabilityAndAmountOfSymbols(int lineIndex, int columnIndex, int insertStateCounterColumn, ref Dictionary<char, int> symbols, ref Dictionary<char, double> probabilities)
        {
            var symbol = GetSymbol(lineIndex, columnIndex);
            if (!IsGap(symbol))
            {
                UpdateTheAmountOfSymbols(symbol, ref symbols);
                UpdateSymbolProbability(insertStateCounterColumn, symbol, symbols, ref probabilities);
            }
        }

        private void UpdateAllProbabilityAndAmountOfSymbols(int columnIndex, int insertStateCounterColumn, ref Dictionary<char, int> symbols, ref Dictionary<char, double> probabilities)
        {
            var observationsAmount = GetObservationsAmount();
            for (int lineIndex = 0; lineIndex < observationsAmount; lineIndex++)
                UpdateProbabilityAndAmountOfSymbols(lineIndex, columnIndex, insertStateCounterColumn, ref symbols, ref probabilities);
        }

        private void UpdateInsertStateCounterColumn(StateEnum state, ref int insertStateCounterColumn)
        {
            if (state == StateEnum.Insert) insertStateCounterColumn++;
            else insertStateCounterColumn = 0;
        }

        private void CreateProbabilities()
        {
            var observationSize = GetObservationSize();
            var insertStateCounter = 0;
            var insertStateCounterColumn = 0;

            for (int index = 0; index < observationSize; index++)
            {
                var currentIndex = index - insertStateCounter;
                var symbols = InitSymbols();
                var probabilities = InitEmissionProbabilities();
                State currentState = GetCurrentState(currentIndex);
                UpdateInsertStateCounterColumn(currentState.Name, ref insertStateCounterColumn);
                UpdateAllProbabilityAndAmountOfSymbols(index, insertStateCounterColumn, ref symbols, ref probabilities);
                UpdateEmissionProbabilities(currentState.EmissionProbabilities, probabilities);
                UpdateInsertStateCounterForProbabilities(index, ref insertStateCounter);
            }
        }

        private void UpdateEmissionProbabilities(Dictionary<char, double> emissionProbabilities, Dictionary<char, double> probabilities)
        {
            foreach (var item in probabilities)
                emissionProbabilities[item.Key] += item.Value;
        }

        private Dictionary<char, int> InitSymbols()
        {
            switch (TypeBioSequence)
            {
                case TypeBioSequence.DNA:
                    return new Dictionary<char, int>()
                    {
                        { 'A', 0 },
                        { 'C', 0 },
                        { 'G', 0 },
                        { 'T', 0 }
                    };
                case TypeBioSequence.RNA:
                    return new Dictionary<char, int>()
                    {
                        { 'A', 0 },
                        { 'C', 0 },
                        { 'G', 0 },
                        { 'U', 0 }
                    };
                case TypeBioSequence.PROTEIN:
                    return new Dictionary<char, int>()
                    {
                        { 'A', 0 },
                        { 'C', 0 },
                        { 'D', 0 },
                        { 'E', 0 },
                        { 'F', 0 },
                        { 'G', 0 },
                        { 'H', 0 },
                        { 'I', 0 },
                        { 'K', 0 },
                        { 'L', 0 },
                        { 'M', 0 },
                        { 'N', 0 },
                        { 'P', 0 },
                        { 'Q', 0 },
                        { 'R', 0 },
                        { 'S', 0 },
                        { 'T', 0 },
                        { 'V', 0 },
                        { 'W', 0 },
                        { 'Y', 0 }
                    };
                default:
                    return new Dictionary<char, int>()
                    {
                        { 'A', 0 },
                        { 'C', 0 },
                        { 'G', 0 },
                        { 'T', 0 }
                    };
            }
        }

        private Dictionary<char, double> InitEmissionProbabilities()
        {
            switch (TypeBioSequence)
            {
                case TypeBioSequence.DNA:
                    return new Dictionary<char, double>()
                    {
                        { 'A', 0d },
                        { 'C', 0d },
                        { 'G', 0d },
                        { 'T', 0d }
                    };
                case TypeBioSequence.RNA:
                    return new Dictionary<char, double>()
                    {
                        { 'A', 0d },
                        { 'C', 0d },
                        { 'G', 0d },
                        { 'U', 0d }
                    };
                case TypeBioSequence.PROTEIN:
                    return new Dictionary<char, double>()
                    {
                        { 'A', 0d },
                        { 'C', 0d },
                        { 'D', 0d },
                        { 'E', 0d },
                        { 'F', 0d },
                        { 'G', 0d },
                        { 'H', 0d },
                        { 'I', 0d },
                        { 'K', 0d },
                        { 'L', 0d },
                        { 'M', 0d },
                        { 'N', 0d },
                        { 'P', 0d },
                        { 'Q', 0d },
                        { 'R', 0d },
                        { 'S', 0d },
                        { 'T', 0d },
                        { 'V', 0d },
                        { 'W', 0d },
                        { 'Y', 0d }
                    };
                default:
                    return new Dictionary<char, double>()
                    {
                        { 'A', 0d },
                        { 'C', 0d },
                        { 'G', 0d },
                        { 'T', 0d }
                    };
            }
        }

        private Dictionary<StateEnum, double> InitTransitionProbabilities()
        {
            return new Dictionary<StateEnum, double>()
            {
                { StateEnum.Match, 0d },
                { StateEnum.Insert, 0d },
                { StateEnum.Delete, 0d }
            };
        }

        private void UpdateInsertStateCounterForTransitions(ref int insertStateCounter, int index, int observationsAmount)
        {
            var symbolsAmount = GetAmountOfMatches(index);
            var nextSymbolsAmount = GetAmountOfMatches(index + 1);

            if (symbolsAmount != observationsAmount && nextSymbolsAmount != observationsAmount)
                insertStateCounter++;
        }

        private bool IsLastIndexStates(int index)
        {
            return index == States.Count;
        }

        private int GetCurrentIndex(int index, int insertStateCounter)
        {
            return index - insertStateCounter;
        }

        private int GetIndexToTransition(int currentIndex, int previousIndex)
        {
            return previousIndex >= Observations[0].Length - 2 ? currentIndex : previousIndex;
        }

        private void CreateTransitions()
        {
            var insertStateCounter = 0;
            var observationSize = GetObservationSize();

            for (int previousIndex = 0; previousIndex < observationSize; previousIndex++)
            {
                var currentIndex = GetCurrentIndex(previousIndex, insertStateCounter);
                if (IsLastIndexStates(currentIndex)) return;

                State currentState = GetCurrentState(currentIndex);
                CreateTransitionsStates(previousIndex, currentIndex, currentState);

                var indexToTransition = GetIndexToTransition(currentIndex, previousIndex);
                UpdateInsertStateCounterForTransitions(ref insertStateCounter, indexToTransition, Observations.Count);
            }
        }

        private int GetAmountOfMatches(int index)
        {
            if (IsLastIndexObservation(index)) return 0;

            var count = 0;
            for (int i = 0; i < Observations.Count; i++)
            {
                var symbol = GetSymbol(i, index);

                if (!IsGap(symbol))
                    count++;
            }

            return count;
        }

        private int GetAmountOfInserts(int index)
        {
            var indexes = new List<int>();
            var observationSize = GetObservationSize();
            var observationsAmount = GetObservationsAmount();
            var matchCount = 0;
            var gapCount = 0;

            for (int i = index; i < observationSize; i++)
            {
                for (int j = 0; j < observationsAmount; j++)
                {
                    var symbol = GetSymbol(j, i);

                    if (!IsGap(symbol))
                    {
                        matchCount++;
                        if (!indexes.Contains(j))
                            indexes.Add(j);
                    }
                    else
                        gapCount++;
                }
                if (matchCount == observationsAmount) break;
            }

            return indexes.Count;
        }

        private State GetNextState(int currentIndex)
        {
            return currentIndex < States.Count - 1 ? States[currentIndex + 1] : null;
        }

        private double GetPercentageOfMachState(int symbolsAmount)
        {
            return Convert.ToDouble(symbolsAmount) / Convert.ToDouble(Observations.Count);
        }

        private double GetPercentageOfInsertState(double percentageOfMachState)
        {
            return 1d - percentageOfMachState;
        }

        private double CalculatePercentageOfMachState(int index, int currentIndex, State currentState, State nextState)
        {
            double percentageOfMachState = 0d;
            var symbolsAmount = 0;

            if (currentState.Name == StateEnum.Match && (nextState == null || nextState.Name == StateEnum.Match))
            {
                symbolsAmount = GetAmountOfMatches(index);
                percentageOfMachState = GetPercentageOfMachState(symbolsAmount);
            }
            else if (currentState.Name == StateEnum.Match && nextState.Name == StateEnum.Insert)
            {
                symbolsAmount = GetAmountOfInserts(currentIndex + 1);
                var percentageOfInsertState = GetPercentageOfMachState(symbolsAmount);
                percentageOfMachState = GetPercentageOfInsertState(percentageOfInsertState);
            }
            else if (currentState.Name == StateEnum.Insert && (nextState == null || nextState.Name == StateEnum.Match))
            {
                symbolsAmount = GetAmountOfInserts(currentIndex);
                percentageOfMachState = GetPercentageOfMachState(symbolsAmount);
            }

            return percentageOfMachState;
        }

        private double CalculatePercentageOfInsertState(double percentageOfMachState)
        {
            return GetPercentageOfInsertState(percentageOfMachState);
        }

        private void UpdateTransitions(State currentState, double percentageOfMachState, double percentageOfInsertState)
        {
            currentState.TransitionProbabilities = new Dictionary<StateEnum, double>()
                        {
                            { StateEnum.Match, percentageOfMachState },
                            { StateEnum.Insert, percentageOfInsertState },
                            { StateEnum.Delete, 0.0 }
                        };
        }

        private void CreateTransitionsStates(int previousIndex, int currentIndex, State currentState)
        {
            if (IsLastIndexStates(currentIndex)) return;

            var nextState = GetNextState(currentIndex);
            var percentageOfMachState = CalculatePercentageOfMachState(previousIndex, currentIndex, currentState, nextState);
            var percentageOfInsertState = CalculatePercentageOfInsertState(percentageOfMachState);

            UpdateTransitions(currentState, percentageOfMachState, percentageOfInsertState);
        }

        //public double CalculateLogOdds(char[] observations)
        //{
        //    var probability = CalculateTotalProbability(observations);
        //    return Math.Log(probability) - (observations.Length * Math.Log(1d / Constants.DNA.Length));
        //}

        private double ForwardViterbi(char[] observation, int index, int length)
        {
            var prob = 1d;
            var lastIndexState = index + length;
            var observationIndex = 0;

            for (int i = index; i < lastIndexState; i++)
            {
                var symbol = observation[observationIndex];
                var currentState = GetCurrentState(i);

                if (currentState.Name == StateEnum.Delete) continue;

                var emissionProbability = 0d;
                var transitionProbabilities = 0d;

                if (currentState.TransitionProbabilities != null)
                    transitionProbabilities = currentState.TransitionProbabilities.Max(m => m.Value);

                if (currentState.EmissionProbabilities != null)
                    emissionProbability = currentState.EmissionProbabilities[symbol];

                prob *= emissionProbability;
                prob *= transitionProbabilities;
                observationIndex++;

                if (prob == 0) break;
            }

            return prob;
        }

        public double CalculateCloneRate(double affinity, int length)
        {
            return affinity * length;
        }

        public double CalculateMutationRate(double affinity, int length)
        {
            return 1 - affinity;
        }

        public bool IsBetterAffinity(double affinityAB, double affinityM)
        {
            return affinityAB > affinityM;
        }

        public IEnumerable<Antibody> Order(List<Antibody> population)
        {
            return population.OrderByDescending(o => o.Affinity);
        }

        public int SequenceSize()
        {
            return States.Count();
        }

        public double Calculate(char[] sequence, int index, int length)
        {
            return ForwardViterbi(sequence, index, length);
        }
    }
}