using System;
using System.Text;
using System.IO;
using BW;

namespace ProbandoFunciones{
  class Programa{
    
    static public short[,] ReadFrame(int argumento, BrwRdr archivo ){
      /*funcion que a partir de un argumento devuelve el frame
	numero argumento de un archivo brw que es el segundo argumento
       */
      //canales a medir (TODOS)
      chs=archivo.GetRecChsUnion();
      int aux1, aux2;
      short[,] result;
      //Un cuadrito mide a lo mas 64x64 electrodos
      result= new short[64,64];
      short auxz;
      double[] auxdato;
      double y;

      foreach (ChCoord  k in chs){
	//este son las cordeadas del Canal k
	aux1=k.Row;
	aux2=k.Col;	
	auxdato=archivo.GetRawData(k,argumento,1);//un solo dato de ese canal	
	//Hacemos la conversion horrible de hace rato
	y=auxdato[0];
	y *= 100;
	z = (short) y;
	//metemos el valor como nos gusta (C numera desde 0)
        result[aux1-1,aux2-1]=z;      

      }
      return result;
      
    }

    static void Main(){
      
      //abrir un archivo BrainWave
       BrwRdr chanfles = new BrwRdr();
       chanfles.Open("G:\\290414raw.brw");
           
       //probar la funcion: deberia de regresar una matriz de shorts
       short[,] z=ReadFrame(chanfles, 2);
       //string[,] puchas=Convert.ToString(z);
       Console.WriteLine("Tu abuelita en vinagre");
       
       /* foreach (short j in z){
	 Console.WriteLine(j);
	 }*/
       //System.IO.File.WriteAllLines("Test.dat",puchas);
       using (StreamWriter outfile = new StreamWriter("datos.dat"))
	 {
	   for (int x = 0; x < 64; x++)
	    {
	      string content = "";
	      for (int y = 0; y < 64; y++)
		{
		  content += z[x, y].ToString("0") + "\t";
		}
	      outfile.WriteLine(content);
	    }
	}
 
    }
    
  }
}

