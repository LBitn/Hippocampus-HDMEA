using System;
using BW;
namespace BW
{
  
  class Program
  {
      public int* ReadFrame(){
	
      }//ends ReadFrame
      
        static void Main(string[] args)
        {
            ChCoord ch = new ChCoord(53, 03);
            Console.WriteLine("The channel has coordinates: row: {0}; column: {1}", ch.Row, ch.Col);
            BrwRdr brwRdr = new BrwRdr();
            brwRdr.Open("G:\\290414raw.brw");
            //double[] datoscrudos;
            double[] datoscrudos;
            datoscrudos = brwRdr.GetRawData(ch,0,7000);
            foreach (double x in datoscrudos)
            {
                short z;
                double y = x;
                y *= 100;
                z = (short) y;
                Console.WriteLine(z);
            }
            //double nSec = brwRdr.RecDuration / 1000;
            long nFrames = brwRdr.RecNFrames;
            Console.WriteLine(nFrames);
            
            /*Console.WriteLine(string.Join(" ", datoscrudos));*/
            
            Console.ReadKey();
        }
    }
}
