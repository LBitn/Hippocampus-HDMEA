using System;
using System.Text;
using System.IO;
using BW;

namespace ProbandoFunciones
{
    class Programa
    { //El programa ejemplo que muestra como usar ReadFrame
        static public short[,] ReadFrame(BrwRdr archivo, int argumento=1)
        {
            /*funcion que a partir de un argumento devuelve el cuadro
          numero "argumento" de un "archivo" tipo brw.  */
            //canales a medir (TODOS). 
            /* Tal vez podamos acelerar un poco el proceso descartando los canales 
             que esten saturados o que marquen puro 0 con un condicional. Las funciones de BrwRdr 
             tienen algo para verificar eso, al parecer. */
            ChCoord[] chs = archivo.GetRecChsUnion();
            int aux1, aux2; //las cordenadas de cada canal
            short[,] result; //la matriz imagen, escrita como una coleccion de short integer
            //Un cuadrito mide a lo mas 64x64 electrodos
            result = new short[64, 64]; 
            short auxz; //variable auxiliar para el casting
            double[] auxdato; //variable auxiliar para el casting
            double y;//variable auxiliar para el casting

            foreach (var k in chs)
            {
                //este son las cordeadas del Canal k, funciones intrinsecas de la clase ChCoord
                aux1 = k.Row; 
                aux2 = k.Col;
                auxdato = archivo.GetRawData(k, argumento, 1);//un solo dato de ese canal	
                //Hacemos la conversion horrible de hace rato
                y = auxdato[0];
                y *= 100;
                auxz = (short)y;
                //metemos el valor como nos gusta (C numera desde 0)
                result[aux1 - 1, aux2 - 1] = auxz;

            }
            return result;

        }

        static void Main()
        { //Rutina "main" de prueba a ver si funciona

            //abrir un archivo BrainWave
            BrwRdr chanfles = new BrwRdr();
            chanfles.Open("G:\\290414raw.brw");

            //probar la funcion: deberia de regresar una matriz de shorts
            short[,] z = ReadFrame(chanfles, 5000); //Esta parte se tarda bastante al parecer
            //string[,] puchas=Convert.ToString(z);
            Console.WriteLine("Tu abuelita en vinagre"); //mensaje de aviso, todo va bien hasta ahora.

            /* foreach (short j in z){
          Console.WriteLine(j);
          }*/
            //System.IO.File.WriteAllLines("Test.dat",puchas);
            using (StreamWriter outfile = new StreamWriter("datos.dat"))
            {
                for (int x = 0; x < 64; x++) //es probable que sea mejor convertir esto a alguna funcion del archivo.
                {
                    string content = "";
                    for (int y = 0; y < 64; y++) //igual que el comentario anterior.
                    {
                        content += z[x, y].ToString("0") + "\t";
                    }
                    outfile.WriteLine(content);
                }
            }

        }

    }
}