using System;
using System.IO;

namespace AIS.ClonalgPR
{
    public class Helpers
    {
        public static string GetPath()
        {
            return Path.Combine(Environment.CurrentDirectory, "Data");
        }

        public static string GetPathFile()
        {
            return Path.Combine(Environment.CurrentDirectory, "Data", "Proteins.fa");
        }
    }
}
