using System;
using System.Text;
using System.IO;
//using System.Random;

namespace ProbandoFunciones
{
  class Programa{


    static void Main(){
      short[,] prueba=new short[64,64];
      
      using (BinaryReader reader = new BinaryReader(File.Open("datos.bin", FileMode.Open)))
            {
	      for (int j=0; j<64; j++){
		 for (int k=0; k<64; k++){
		   prueba[j,k]=reader.ReadInt16();
		   Console.WriteLine(prueba[j,k]); 
		 }
	      }
            }


    }

  }
  
}