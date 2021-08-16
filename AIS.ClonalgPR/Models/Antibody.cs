namespace AIS.ClonalgPR.Models
{
    public class Antibody
    {
        public double Affinity { get; set; }
        public string Antigen { get; set; }
        public int Length { get; set; }
        public char[] Sequence { get; set; }

        public Antibody Clone()
        {
            return new Antibody
            {
                Affinity = Affinity,
                Antigen = Antigen,
                Length = Length,
                Sequence = Sequence
            };
        }
    }
}
