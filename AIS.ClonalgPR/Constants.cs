using System;

namespace AIS.ClonalgPR
{
    public static class Constants
    {
        private static char[] _aminoacids = null;
        private static char[] _dna = null;
        private static char[] _rna = null;
        private static Random _random = null;
        private static char[] _gaps = null;

        public const int MINIMAL_AMOUNT_OF_ANTIBODIES = 1000;
        public const int MAXIMUM_AMOUNT_OF_ANTIBODIES = 10000;
        public const int MIN_SIZE_ANTIBODY = 3;
        public const int MAX_SIZE_ANTIBODY = 9;

        public static char[] Aminoacids
        {
            get
            {
                if (_aminoacids == null)
                    _aminoacids = new char[20] { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };

                return _aminoacids;
            }
        }

        public static char[] DNA
        {
            get
            {
                if (_dna == null)
                    _dna = new char[4] { 'A', 'C', 'T', 'G' };

                return _dna;
            }
        }

        public static char[] RNA
        {
            get
            {
                if (_rna == null)
                    _rna = new char[4] { 'A', 'C', 'U', 'G' };

                return _rna;
            }
        }

        public static Random Random
        {
            get
            {
                if (_random == null)
                    _random = new Random();

                return _random;
            }
        }

        public static char[] Gaps
        {
            get
            {
                if (_gaps == null)
                    _gaps = new char[4] { '.', '-', 'x', 'X' };

                return _gaps;
            }
        }
    }
}
