using AIS.ClonalgPR.Measures;
using AIS.ClonalgPR.Models;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace AIS.ClonalgPR
{
    public class ClonalgPR
    {
        private List<Antibody> _memoryCells = new List<Antibody>();
        private List<string> _memoryCellsStr = new List<string>();
        private Result _results = new Result();
        private IDistance _distance = null;
        private List<Antigen> _antigens = new List<Antigen>();
        private TypeBioSequence _typeBioSequence;
        private int _antibodySize = 0;

        public Result Results
        {
            get { return _results; }
        }

        public List<string> MemoryCells
        {
            get { return _memoryCellsStr; }
        }

        public ClonalgPR(IDistance distance, List<Antigen> antigens, TypeBioSequence typeBioSequence, int antibodySize = Constants.MAX_SIZE_ANTIBODY)
        {
            _distance = distance;
            _antigens = antigens;
            _typeBioSequence = typeBioSequence;
            _antibodySize = antibodySize;
        }

        private void Affinity(List<Antibody> antibodies)
        {
            var antigenSize = _distance.SequenceSize();

            for (int i = 0; i < antibodies.Count(); i++)
            {
                var sequence = antibodies[i].Sequence;
                var length = sequence.Length;
                var index = Constants.Random.Next(0, antigenSize - length);
                antibodies[i].Affinity = _distance.Calculate(sequence, index, length);

            }
        }

        private List<Antibody> Clone(List<Antibody> antibodies)
        {
            var clones = new List<Antibody>();
            foreach (var antibody in antibodies)
            {
                var rate = _distance.CalculateCloneRate(antibody.Affinity, antibody.Length);
                var clonesAmount = Math.Round(rate * antibody.Length);
                for (int j = 0; j < clonesAmount; j++)
                    clones.Add(antibody.Clone());
            }

            return clones;
        }

        private List<Antibody> Initialize(int antibodyAmount = 0)
        {
            var antibodies = new List<Antibody>();

            if (antibodyAmount == 0)
                antibodyAmount = Constants.Random.Next(Constants.MINIMAL_AMOUNT_OF_ANTIBODIES, Constants.MAXIMUM_AMOUNT_OF_ANTIBODIES);

            //Constants.Random.Next(Constants.MIN_SIZE_ANTIBODY, Constants.MAX_SIZE_ANTIBODY);

            for (int i = 0; i < antibodyAmount; i++)
            {

                var sequence = new char[_antibodySize];

                for (int j = 0; j < _antibodySize; j++)
                    sequence[j] = GenerateSequences();

                antibodies.Add(new Antibody
                {
                    Sequence = sequence,
                    Length = sequence.Where(w => !Constants.Gaps.Contains(w)).Count()
                });
            }

            return antibodies;
        }

        private char GenerateSequences()
        {
            switch (_typeBioSequence)
            {
                case TypeBioSequence.DNA:
                    return Constants.DNA[Constants.Random.Next(0, Constants.DNA.Length)];
                case TypeBioSequence.RNA:
                    return Constants.RNA[Constants.Random.Next(0, Constants.RNA.Length)];
                case TypeBioSequence.PROTEIN:
                    return Constants.Aminoacids[Constants.Random.Next(0, Constants.Aminoacids.Length)];
                default:
                    return Constants.DNA[Constants.Random.Next(0, Constants.DNA.Length)];
            }
        }

        private void Insert(List<Antibody> antibodies)
        {
            if (antibodies == null || antibodies.Count() == 0) return;

            if (_memoryCells.Count() > 0)
            {
                for (int i = 0; i < antibodies.Count(); i++)
                {
                    for (int j = 0; j < _memoryCells.Count(); j++)
                    {
                        var memoryCell = _memoryCells[j];
                        var antibody = antibodies[i];

                        if (_distance.IsBetterAffinity(antibody.Affinity, memoryCell.Affinity))
                        {
                            memoryCell.Affinity = antibody.Affinity;
                            memoryCell.Length = antibody.Length;
                            memoryCell.Sequence = antibody.Sequence;
                            break;
                        }
                    }
                }
            }
            else
                antibodies.ForEach(antibody => _memoryCells.Add(antibody));
        }

        private List<Antibody> Mutation(List<Antibody> antibodies)
        {
            foreach (var antibody in antibodies)
            {
                var rate = _distance.CalculateMutationRate(antibody.Affinity, antibody.Length);
                var mutateAmount = (int)(antibody.Length * rate);
                var sequenceIndex = Constants.Random.Next(0, antibody.Length);
                var sequence = antibody.Sequence;

                for (int j = 0; j < mutateAmount; j++)
                {
                    var nucleotideOrAminoacid = GenerateSequences();
                    sequence[sequenceIndex] = nucleotideOrAminoacid;
                    antibody.Sequence = sequence;
                    sequenceIndex = Constants.Random.Next(0, antibody.Length);
                }
            }

            return antibodies;
        }

        private List<Antibody> Replace(List<Antibody> antibodies, int inferiorLimit)
        {
            if (antibodies == null || antibodies.Count() == 0) return antibodies;

            var antibodiesReplaced = antibodies.OrderBy(o => o.Affinity).Take(inferiorLimit).ToList();
            antibodiesReplaced.ForEach(antibody => antibodies.Remove(antibody));

            var amountOfAntibodiesReplaced = antibodiesReplaced.Count();
            if (amountOfAntibodiesReplaced > 0)
            {
                var _antibodies = Initialize(antibodyAmount: amountOfAntibodiesReplaced);
                _antibodies.ForEach(p => antibodies.Add(p));
            }
            return antibodies;
        }

        private List<Antibody> Select(List<Antibody> population, int numberHighAffinity)
        {
            return _distance.Order(population).Take(numberHighAffinity).ToList();
        }

        public void Execute(int maximumIterations, double percentHighAffinity, double percentLowAffinity, int index = 0)
        {
            var i = 1;
            var antibodies = Initialize();
            var numberHighAffinity = (int)Math.Round(percentHighAffinity * antibodies.Count());
            var numberLowAffinity = (int)Math.Round(percentLowAffinity * antibodies.Count());

            while (i < maximumIterations)
            {
                Affinity(antibodies);
                var selectedPopulation = Select(antibodies, numberHighAffinity);
                var clonedPopulation = Clone(selectedPopulation);
                var mutatedPopulation = Mutation(clonedPopulation);
                Affinity(mutatedPopulation);
                selectedPopulation = Select(mutatedPopulation, numberHighAffinity);
                Insert(selectedPopulation);
                antibodies = Replace(selectedPopulation, numberLowAffinity);
                i++;
            }
            SetMemoryCells();
            SaveMemoryCells(index);
        }

        private void SetMemoryCells()
        {
            _memoryCells
                .Where(memoryCell => memoryCell.Sequence.Length > 0)
                .ToList()
                .ForEach(memoryCell => _memoryCellsStr.Add(new string(memoryCell.Sequence)));
        }

        private void SaveMemoryCells(int index = 0)
        {
            if (_memoryCellsStr.Count() == 0) return;

            var filePath = Path.Combine(Helpers.GetPath(), $"memoryCells{(index == 0 ? _antibodySize : index)}.json");
            using (StreamWriter file = File.CreateText(filePath))
            {
                JsonSerializer serializer = new JsonSerializer();
                serializer.Serialize(file, _memoryCellsStr);
            }
        }
    }
}
