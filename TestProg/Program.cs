using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Mathcad;
using ODIS.AMM;


namespace CalculateFunction

{
    class Program
    {

        static void Main(string[] args)
        {
            IMathcadApplication mc = new Application();
            IMathcadWorksheets mwk = mc.Worksheets;
            IMathcadWorksheet ws = mwk.Open(@"C:\temp_1\TestProg.mcdx");

            //1.1.параметров входящего потока

            int N = 1; //Задается параметр высокой интенсивности входящего потока

            Console.WriteLine("N = " + N);

            //Задаются параметры λ и κ входящего потока

            double Alpha = 0.25;
            //Console.WriteLine("Alpha = " + (double)Alpha);

            double Beta = 0.25;
            //Console.WriteLine("Beta = " + (double)Beta);

            ws.SetValue("Alpha", 0.25);
            ws.SetValue("Beta", 0.25);

            ws.Recalculate();

            double Lambda = ((double)Beta / (double)Alpha);
            Console.WriteLine("Lambda = " + (double)Lambda);

            double Sigma2 = Alpha / Math.Pow(Beta, 2);
            Console.WriteLine("Sigma2 = " + (double)Sigma2);

            double a1 = (1 / (double)Lambda);
            Console.WriteLine("a1 = " + a1);

            double Kappa = ((Math.Pow(Lambda, 3)) * (double)Sigma2 - (Math.Pow(a1, 2)));
            Console.WriteLine("Kappa = " + Kappa);


            //1.2.параметров маршрутизации

            int K = 4; //Число узлов сети

            ws.SetValue("K", 4);
            ws.Recalculate();

            Console.WriteLine("K = " + K);

            ComplexMatrix v = new ComplexMatrix(1, 4);
            v[1, 1] = 0.2;
            v[1, 2] = 0.2;
            v[1, 3] = 0.1;
            v[1, 4] = 0.5;

            ws.SetValue("v", v);
            ws.Recalculate();

            //Console.WriteLine(String.Format("{0:0}", "v = [" + v[1, 1].Real + ", " + v[1, 2].Real + ", " + v[1, 3].Real + ", " + v[1, 4].Real) + "]");
            //System.Console.WriteLine();

            // ComplexMatrix M = new ComplexMatrix(4, 4);
            ComplexMatrix M = new ComplexMatrix(4, 4);
            M[1, 1] = 0.1; M[1, 2] = 0.4; M[1, 3] = 0.4; M[1, 4] = 0;
            M[2, 1] = 0; M[2, 2] = 0.4; M[2, 3] = 0.1; M[2, 4] = 0.2;
            M[3, 1] = 0.1; M[3, 2] = 0; M[3, 3] = 0.2; M[3, 4] = 0.1;
            M[4, 1] = 0.3; M[4, 2] = 0.2; M[4, 3] = 0.5; M[4, 4] = 0;

            for (int i = 1; i <= 4; i++)
            {
                for (int j1 = 1; j1 <= 4; j1++)
                {
                    if (i == 1 && j1 == 1)
                        System.Console.Write("M = \t[" +
                            String.Format("{0:0.0}", M[i, j1].Real) + " ");
                    else
                        if (i == 4 && j1 == 4)
                        System.Console.Write(String.Format("{0:0.0}", M[i, j1].Real) + "]");
                    else
                        System.Console.Write(String.Format("{0:0.0}", M[i, j1].Real) + " ");
                }
                ws.SetValue("M", M);
                ws.Recalculate();

                //System.Console.WriteLine();
                System.Console.Write("\t");

            }
            //System.Console.WriteLine();

            //1.3.параметров обслуживания

            float[] AlphaS = new float[4];  //компоненты векторов αS 
            AlphaS[0] = 0.25f;
            AlphaS[1] = 0.25f;
            AlphaS[2] = 1.5f;
            AlphaS[3] = 4f;

            ws.SetValue("AlphaS", AlphaS);
            ws.Recalculate();
            
            Console.WriteLine("AlphaS = [" + AlphaS[0] + ", " + AlphaS[1] + ", " + AlphaS[2] + ", " + AlphaS[3] + "]");
            System.Console.WriteLine();

            float[] BetaS = new float[4];   //компоненты векторов βS
            BetaS[0] = 0.5f;
            BetaS[1] = 0.25f;
            BetaS[2] = 3f;
            BetaS[3] = 2f;

            ws.SetValue("BetaS", BetaS);
            ws.Recalculate();

            Console.WriteLine("BetaS = [" + BetaS[0] + ", " + BetaS[1] + ", " + BetaS[2] + ", " + BetaS[3] + "]");
            System.Console.WriteLine();

            var result = ws.GetValue("All SetValue");

            mc.Visible = true;

            // System.Console.WriteLine("Covar = ");
           
            Console.ReadKey();
        }
    }
}