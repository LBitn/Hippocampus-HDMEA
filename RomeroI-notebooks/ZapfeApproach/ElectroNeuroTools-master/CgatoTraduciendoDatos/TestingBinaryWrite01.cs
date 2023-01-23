using System;
using System.Text;
using System.IO;
//using System.Random;

namespace ProbandoFunciones
{
    class Programa
    { //El programa ejemplo que muestra como usar ReadFrame
        static public short[,] MatrizShortAleatorea(int argumento=2)
        {
            /* 
	       Funcion que crea una matriz aleatorea de argumento*argumento elementos
	       tipo short para checar BinaryWrite
	     */
	  short[,] result=new short[argumento,argumento];
	  Random aux= new Random();
	  short m;
	  
	  for (int j=0; j<argumento; j++){
	    for (int k=0; k< argumento; k++){
	       //creacion del numero arbitrario
	      m=(short)aux.Next(-100,100);		 
	       result[j, k] = m;
	       
            }
	  }
	  
	  return result;
	  
        }

        static void Main()
        { //Rutina "main" de prueba a ver si funciona


	  int argporom=64;
	  short[,] z = MatrizShortAleatorea(argporom);
	  Console.WriteLine("Tu abuelita en vinagre"); //mensaje de aviso, todo va bien hasta ahora.
       
            
            using (StreamWriter outfile = new StreamWriter("datos.dat"))
	      {//escribir en texto crudo
                for (int x = 0; x < argporom; x++) 
                {
                    string content = "";
                    for (int y = 0; y < argporom; y++) 
                    {
                        content += z[x, y].ToString("0") + "\t";
                    }
                    outfile.WriteLine(content);
                }
            }

	    using (BinaryWriter writer = new BinaryWriter(File.Open("TestBlabla.bin", FileMode.Create)))
	      {//escribir en binario
		writer.Write(1.250F);
		writer.Write(@"c:\Temp");
		writer.Write(10);
		writer.Write(true);
	      }

	    using (BinaryWriter chura = new BinaryWriter(File.Open("datos.bin", FileMode.Create)))
	      {//escribir en binario

		  for (int x = 0; x < argporom; x++) 
                {
                    
                    for (int y = 0; y < argporom; y++) 
                    {
		      chura.Write(z[x,y]);
                    }
		    
                    
                }
		
	      }
	    

        }

    }
}