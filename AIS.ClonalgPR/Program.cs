using AIS.ClonalgPR.Measures;
using AIS.ClonalgPR.Models;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace AIS.ClonalgPR
{
    class Program
    {
        private const int START_REGION_SPIKE_PROTEIN = 21300;
        private const int END_REGION_SPIKE_PROTEIN = 25400;

        private static List<Antigen> antigens = new List<Antigen>();
        private static TypeBioSequence TypeBioSequence
        {
            get
            {
                return IsDNA() ?
                    TypeBioSequence.DNA :
                     IsRNA() ?
                    TypeBioSequence.RNA :
                    IsProtein() ?
                    TypeBioSequence.PROTEIN :
                    TypeBioSequence.UNKNOWN;
            }
        }
        private static bool IsSigle
        {
            get
            {
                antigens = GetAntigens();
                return antigens != null && antigens.Count > 0;
            }
        }
        private static List<string> PatternsProof
        {
            get 
            {
                return new List<string>() { "TAAA", "ATG", "TTT", "TCT", "GAT", "TGG", "ACT", "ATTTTG", "ATA", "CGCT", "CAA" };
            }
        }

        static void Main()
        {
            if (IsSigle)
                ExecuteSingle(antigens);
            else
            {
                var text = GetTextFile();
                var antigensMultiple = GetAntigensByFastaFileAligned(text);

                ExecuteMultiple(antigensMultiple);
            }
        }

        static void ExecuteSingle(List<Antigen> antigens)
        {
            var sequences = GetSequencesByAntigens(antigens);
            var markov = new HiddenMarkovModel(sequences, TypeBioSequence);

            markov.Train();

            for (int antibodySize = Constants.MIN_SIZE_ANTIBODY; antibodySize <= Constants.MAX_SIZE_ANTIBODY; antibodySize++)
            {
                var clonalgPR = new ClonalgPR(distance: markov, antigens: antigens, typeBioSequence: TypeBioSequence, antibodySize: antibodySize);
                clonalgPR.Execute(maximumIterations: 1000, percentHighAffinity: 0.6, percentLowAffinity: 0.4);
            }

            ReadAllFiles(antigens.Count());
        }
        static void ExecuteMultiple(List<List<Antigen>> antigens)
        {
            for (int index = 0; index < antigens.Count; index++)
            {
                var _antigens = antigens[index];
                var sequences = _antigens.Select(s => s.Sequence).ToList();
                var markov = new HiddenMarkovModel(sequences, TypeBioSequence);

                markov.Train();

                var clonalgPR = new ClonalgPR(distance: markov, antigens: _antigens, typeBioSequence: TypeBioSequence);
                clonalgPR.Execute(maximumIterations: 10000, percentHighAffinity: 0.9, percentLowAffinity: 0.1, index);
            }

            ReadAllFiles(antigens.Count());
        }

        private static bool IsDNA()
        {
            var sequence = antigens.Select(antigen => antigen.Sequence).FirstOrDefault();
            var regex = new Regex(@"^[AaCcTtGg\W]+$", RegexOptions.Multiline);
            var matches = regex.Matches(new string(sequence));
            return matches.Count() > 0;
        }

        private static bool IsRNA()
        {
            var sequence = antigens.Select(antigen => antigen.Sequence).FirstOrDefault();
            var regex = new Regex(@"^[AaCcUuGg\W]+$", RegexOptions.Multiline);
            var matches = regex.Matches(new string(sequence));
            return matches.Count() > 0;
        }

        private static bool IsProtein()
        {
            var sequence = antigens.Select(antigen => antigen.Sequence).FirstOrDefault();
            var regex = new Regex(@"^[AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy\W]+$", RegexOptions.Multiline);
            var matches = regex.Matches(new string(sequence));
            return matches.Count() > 0;
        }

        private static List<char[]> GetSequencesByAntigens(List<Antigen> antigens)
        {
            return antigens.Select(s => s.Sequence).ToList();
        }

        private static string[] GetCharacteristics(string antigen)
        {
            var splitOptions = new char[3] { '|', '\n', '\r' };
            return (!string.IsNullOrEmpty(antigen)) ? string.Join("", antigen).Split(splitOptions, StringSplitOptions.RemoveEmptyEntries) : new string[] { };
        }

        private static string GetNameAntigen(string[] characteristics)
        {
            return (characteristics.Length > 0) ? characteristics[0] : "";
        }

        private static string GetHostAntigen(string[] characteristics)
        {
            return (characteristics.Length > 1) ? characteristics[1] : "";
        }

        private static string GetSequenceAntigen(string[] characteristics)
        {
            var length = characteristics.Length - 3;
            var targetArray = new string[length];
            Array.Copy(characteristics, 3, targetArray, 0, length);
            return string.Join("", targetArray);
        }

        private static string GetTextFile()
        {
            var path = Helpers.GetPathFile();
            return File.ReadAllText(path);
        }

        private static List<Antigen> GetAntigens()
        {
            var text = GetTextFile();

            if (IsMsfFile(text))
                return GetAntigensByMsfFile(text);

            if (IsFastaFile(text))
                return GetAntigensByFastaFile(text);

            return new List<Antigen>();
        }

        private static List<List<Antigen>> GetAntigensByFastaFileAligned(string text)
        {
            var antigensMultiple = new List<List<Antigen>>();
            var antigens = new List<Antigen>();
            var sequences = text.Split("tr", StringSplitOptions.RemoveEmptyEntries);

            foreach (var item in sequences)
            {
                var antigenName = item.Split("|", StringSplitOptions.RemoveEmptyEntries);
                var antigenHost = antigenName[1].Split(" ", StringSplitOptions.RemoveEmptyEntries);
                var partOfAligment = antigenHost[1].Trim();

                if (!antigens.Any(antigen => antigen.Name == antigenName[0]))
                {
                    antigens.Add(
                        new Antigen()
                        {
                            Name = antigenName[0],
                            Host = antigenHost[0],
                            Sequence = partOfAligment.Trim().ToCharArray(),
                            Length = partOfAligment.Length
                        });
                }
                else
                {
                    antigensMultiple.Add(antigens);
                    antigens = new List<Antigen>();
                    antigens.Add(
                        new Antigen()
                        {
                            Name = antigenName[0],
                            Host = antigenHost[0],
                            Sequence = partOfAligment.Trim().ToCharArray(),
                            Length = partOfAligment.Length
                        });
                }
            }
            antigensMultiple.Add(antigens);

            return antigensMultiple;
        }

        private static bool IsFastaFile(string text)
        {
            var regex = new Regex(@"^>[A-Za-z1-9\/\|\s].+$", RegexOptions.Multiline);
            var matches = regex.Matches(text);
            return matches.Count > 0;
        }

        private static bool IsMsfFile(string text)
        {
            var regex = new Regex(@"[^Name: ]\w.+(?=\soo)", RegexOptions.Multiline);
            var matches = regex.Matches(text);
            return matches.Count > 0;
        }

        private static List<Antigen> GetAntigensByMsfFile(string text)
        {
            var antigens = new List<Antigen>();
            var regex = new Regex(@"[^Name: ]\w.+(?=\soo)", RegexOptions.Multiline);
            var matches = regex.Matches(text);

            foreach (var item in matches)
            {
                var sequenceName = item.ToString();
                var expression = @"(?<=" + sequenceName + @")\s{2,}[a-zA-Z].+$";
                regex = new Regex(expression, RegexOptions.Multiline);
                var partsOfSequence = regex.Matches(text);
                var sequence = string.Join("", partsOfSequence).Replace(" ", "");
                antigens.Add(new Antigen()
                {
                    Name = sequenceName,
                    Sequence = sequence.Trim().ToCharArray(),
                    Length = sequence.Length
                });
            }

            return antigens;
        }

        private static List<Antigen> GetAntigensByFastaFile(string text)
        {
            var antigens = new List<Antigen>();
            var sequences = text.Split(">", StringSplitOptions.RemoveEmptyEntries);

            foreach (var sequence in sequences)
            {
                var characteristics = GetCharacteristics(sequence);
                var _sequence = GetSequenceAntigen(characteristics);
                var sequenceSource = _sequence.Trim().ToUpper().ToCharArray();
                var charactersAmount = END_REGION_SPIKE_PROTEIN - START_REGION_SPIKE_PROTEIN;
                char[] sequenceDestination = new char[charactersAmount];
                Array.Copy(sequenceSource, START_REGION_SPIKE_PROTEIN, sequenceDestination, 0, charactersAmount);

                antigens.Add(new Antigen
                {
                    Name = GetNameAntigen(characteristics),
                    Host = GetHostAntigen(characteristics),
                    Sequence = sequenceDestination,
                    Length = _sequence.Length
                });
            }

            return antigens;
        }

        private static void ReadResults(int antigensCount)
        {
            var results = new List<Result>();

            if (IsSigle)
            {
                for (int index = Constants.MIN_SIZE_ANTIBODY; index <= Constants.MAX_SIZE_ANTIBODY; index++)
                {
                    try
                    {
                        var fileName = $"results{index}.json";
                        var filePath = Path.Combine(Helpers.GetPath(), fileName);
                        JsonSerializer serializer = new JsonSerializer();

                        using (StreamReader streamReader = new StreamReader(filePath))
                        {
                            JsonTextReader reader = new JsonTextReader(streamReader);
                            var result = serializer.Deserialize<Result>(reader);
                            results.Add(result);
                        }
                    }
                    catch (Exception)
                    {
                        continue;
                    }
                }
            }
            else
            {
                for (int index = 0; index < antigensCount; index++)
                {
                    try
                    {
                        var fileName = $"results{index}.json";
                        var filePath = Path.Combine(Helpers.GetPath(), fileName);
                        JsonSerializer serializer = new JsonSerializer();

                        using (StreamReader streamReader = new StreamReader(filePath))
                        {
                            JsonTextReader reader = new JsonTextReader(streamReader);
                            var result = serializer.Deserialize<Result>(reader);
                            results.Add(result);
                        }
                    }
                    catch (Exception)
                    {
                        continue;
                    }
                }
            }
        }

        private static void EvaluatePatterns(List<string> patterns)
        {
            var intersection = patterns.Intersect(PatternsProof);
            if (intersection.Count() == PatternsProof.Count())
                Console.WriteLine("Foi possível encontrar todos os padrões!");
            else if (intersection.Count() > 0)
                Console.WriteLine($"Foi possível encontrar {intersection.Count()} de {PatternsProof.Count()} padrões!");

            Console.WriteLine("");
            intersection.ToList().ForEach(pattern => Console.WriteLine(pattern));
        }

        private static List<string> ReadSingleFile()
        {
            var patterns = new List<string>();
            JsonSerializer serializer = new JsonSerializer();

            for (int index = Constants.MIN_SIZE_ANTIBODY; index <= Constants.MAX_SIZE_ANTIBODY; index++)
            {
                try
                {
                    var fileName = $"memoryCells{index}.json";
                    var filePath = Path.Combine(Helpers.GetPath(), fileName);

                    using (StreamReader streamReader = new StreamReader(filePath))
                    {
                        JsonTextReader reader = new JsonTextReader(streamReader);
                        var fragments = serializer.Deserialize<List<string>>(reader);
                        fragments.ForEach(fragment => patterns.Add(fragment));
                    }
                }
                catch (Exception)
                {
                    continue;
                }
            }

            return patterns;
        }

        private static List<string> ReadMultipleFiles(int antigensCount)
        {
            var patterns = new List<string>();
            JsonSerializer serializer = new JsonSerializer();

            for (int index = 0; index < antigensCount; index++)
            {
                try
                {
                    var fileName = $"memoryCells{index}.json";
                    var filePath = Path.Combine(Helpers.GetPath(), fileName);

                    using (StreamReader streamReader = new StreamReader(filePath))
                    {
                        JsonTextReader reader = new JsonTextReader(streamReader);
                        var fragments = serializer.Deserialize<List<string>>(reader);

                        for (int j = 0; j < fragments.Count(); j++)
                        {
                            if (patterns.Count() <= j)
                                patterns.Add(fragments[j]);
                            else
                                patterns[j] += fragments[j];
                        }
                    }
                }
                catch (Exception)
                {
                    continue;
                }
            }

            return patterns;
        }

        private static void ReadAllFiles(int antigensCount)
        {
            var patterns = new List<string>();

            if (IsSigle)
                patterns = ReadSingleFile();
            else
                patterns = ReadMultipleFiles(antigensCount);

            EvaluatePatterns(patterns);
        }
    }
}
