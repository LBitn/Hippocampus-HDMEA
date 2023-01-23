using System;
using System.Text;
using System.IO;
using BW;

namespace ProbandoFunciones{
  class Programa{
    
    static public short[,] ReadFrame(int argumento){
      /*funcion que a partir de un argumento devuelve el frame
	numero argumento de un archivo brw
       */
      short[,] result;
      result= new short[,]{{1,2},{3,4}};
      return result;
      
    }

    static void Main(){
      short[,] z=ReadFrame(2);
      //string[,] puchas=Convert.ToString(z);
      Console.WriteLine("Tu abuelita en vinagre");
      
      foreach (short j in z){
      Console.WriteLine(j);
      }
      //System.IO.File.WriteAllLines("Test.dat",puchas);
      using (StreamWriter outfile = new StreamWriter("datos.dat"))
	{
	  for (int x = 0; x < 2; x++)
	    {
	      string content = "";
	      for (int y = 0; y < 2; y++)
		{
		  content += z[x, y].ToString("0") + "\t";
		}
	      outfile.WriteLine(content);
	    }
	}
 
    }
    
  }
}

