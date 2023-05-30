using System;
using System.Collections.Generic;
using System.Text;

namespace ODIS.AMM
{
    public class AMMUtility
    {
        public const double NearZeroValue = 1E-15;

        /// <summary>
        /// Compute natural logarithm of Gamma(x) using the asymptotic Sterling's expansion.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static double logGamma(double x) // From Abramowitz & Stegun, Handbook of Mathematical Functions, 1964 [6.1.41]
        // The first 20 terms give the result with 50 digits. If x <= 0, rise exception.
        {
            int LGM_LIM = 7; // Нижняя граница аргумента, при котором требуется повышение точности в алгоритме вычисления log Г(x)

            int i;
            /* Negative argument: Error!         */
            if (x <= 0) return 0;

            if (x == 1 || x == 2) return 0;

            double z = 0;
            for (z = 0; x < LGM_LIM; x += 1)    /* Increase argument if necessary.   */
                z += Math.Log(x);

            double den = x;
            double x2 = x * x;                         /* Compute the asymptotic expansion  */
            double presum = (x - 0.5) * Math.Log(x) - x + 0.9189385332046727417803297364;
            decimal[] c = 
            {
              /* Asymtotic expansion coefficients             */
              1.0m / 12.0m, -1.0m / 360.0m, 1.0m / 1260.0m, -1.0m / 1680.0m, 1.0m / 1188.0m,
              -691.0m / 360360.0m, 1.0m / 156.0m, -3617.0m / 122400.0m, 43867.0m / 244188.0m,
              -174611.0m / 125400.0m, 77683.0m / 5796.0m, -236364091.0m / 1506960.0m,
              657931.0m / 300.0m, -3392780147.0m / 93960.0m, 1723168255201.0m / 2492028.0m,
              -7709321041217.0m / 505920.0m, 151628697551.0m / 396.0m,
              -26315271553053477373.0m / 2418179400.0m, 154210205991661.0m / 444.0m,
              - 261082718496449122051.0m / 21106800.0m
            };
            double sum = 0;
            for (i = 0; i < 20; i++)
            {
                sum = presum + (double)c[i] / den;
                if (sum == presum) break;
                den = den * x2;
                presum = sum;
            }
            return sum - z;                     /* Fit the increased argument if any  */

        }/*logGamma*/

        public static double GammaSeries(double a, double x)
        // См. Abramowitz & Stegun, Handbook of Mathematical Functions, 1964 [6.5.29]: М.Абрамовиц, И.Стиган. Справочник по специальным функциям (М: Мир, 1979)
        // Разложение в ряд Тейлора используется для вычисления неполной гамма функции P(a,x)
        {
            double sum = 0;
            double prev_sum = 0;
            double term = 0;
            double aa = a;
            long i = 0;
            term = sum = 1 / a;
            do
            {
                aa += 1;
                term *= x / aa;
                prev_sum = sum;
                sum += term;
                i++;
            } while (Math.Abs(prev_sum) != Math.Abs(sum));
            return sum * Math.Exp(-x + a * Math.Log(x) - logGamma(a));
        }

        public static double GammaFraction(double a, double x)
        // См. Abramowitz & Stegun, Handbook of Mathematical Functions, 1964 [6.5.31]: М.Абрамовиц, И.Стиган. Справочник по специальным функциям (М: Мир, 1979)
        // Разложение в цепную дробь Лежандра используется для вычисления неполной гамма функции P(a,x)
        //  P(a,x)=exp(-x +x*ln(a))*CF/logGamma(a),            
        //                                                                      
        //  где
        //        1    1-a   1   2-a   2
        //  CF = ---  ----- --- ----- --- ....
        //       x+    1+   x+   1+   x+
        //
        //  Используются подходящие дроби CF(n) = A(n) / B(n)
        //
        //  где
        //        A(n) = (s(n) * A(n-1) + r(n) * A(n-2)) * factor
        //        B(n) = (s(n) * B(n-1) + r(n) * B(n-2)) * factor
        //  причем                                                                 
        //        A(-1) = 1, B(-1) = 0, A(0) = s(0), B(0) = 1.                  
        //                                                                      
        //  Здесь
        //        s(0) = 0, s(1) = x, r(0) = 0, r(1) = 1,                       
        //                                                                      
        //  так что
        //        A(1) = one * factor, B(1) = r * factor                        
        //                                                                      
        //  и, следовательно,
        //        r(i) = k - a  if i = 2k,   k > 0                              
        //        r(i) = k      if i = 2k+1,                                    
        //        s(i) = 1      if i = 2k,                                      
        //        s(i) = x      if i = 2k+1                                     
        //                                                                      
        //  factor - шкалирующий множитель
        {
            double old_sum = 0; double factor = 1;
            double A0 = 0; double A1 = 1; double B0 = 1; double B1 = x;
            double sum = 1 / x; double z = 0; double ma = -a; double rfact;
            do
            {
                z += 1;
                ma += 1;
                /* two steps of recurrence replace A's & B's */
                A0 = (A1 + ma * A0) * factor;	/* i even */
                B0 = (B1 + ma * B0) * factor;
                rfact = z * factor;
                A1 = x * A0 + rfact * A1;	/* i odd, A0 already rescaled */
                B1 = x * B0 + rfact * B1;
                if (B1 != 0)
                {
                    factor = 1 / B1;
                    old_sum = sum;
                    sum = A1 * factor;
                }
            } while (Math.Abs(sum) != Math.Abs(old_sum));
            return Math.Exp(-x + a * Math.Log(x) - logGamma(a)) * sum;
        }/*GammaFraction*/

        static public double Gamma(double x)
        {
            // x>0
            return Math.Exp(logGamma(x));
        }

        static public double Beta(double a, double b)
        {
            //a,b>0
            return Gamma(a) * Gamma(b) / Gamma(a + b);
        }

        static double logBeta(double a, double b)
        {
            //a,b>0
            return logGamma(a) + logGamma(b) - logGamma(a + b);
        }

        public static double BetaFraction(double a, double b, double x)
        // См. Abramowitz & Stegun, Handbook of Mathematical Functions, 1964 [26.5.8]: М.Абрамовиц, И.Стиган. Справочник по специальным функциям (М: Мир, 1979)
        //
        //   Неполная бета-функция вычисляется с помощью разложения в цепную дробь
        //
        //         i_beta(a,b,x) = x^{a}*(1-x)^{b}*fraction / a * beta(a,b),
        //  где
        //                    1    d1   d2   d3   d4
        //        fraction = ---  ---- ---- ---- ---- ....
        //                    1+   1+   1+   1+   1+
        //
        //  Подходящие дроби: A(n) / B(n)
        //  где
        //        A(n) = (s(n) * A(n-1) + r(n) * A(n-2)) * factor
        //        B(n) = (s(n) * B(n-1) + r(n) * B(n-2)) * factor
        //  и
        //        A(-1) = 1, B(-1) = 0, A(0) = s(0), B(0) = 1.
        //
        //  Здесь s(0) = 0 и при n >= 1 s(n) = 1,
        //  а r(1) = 1 и при i >= 2
        //
        //        r(i) =  m(b-m)x / (a+i-1)(a+i)       когда i = 2m,
        //        r(i) = -(a+m)(a+b+m)x / (a+i-1)(a+i) когда i = 2m+1,
        //  factor - шкалирующий множитель, позволяющий избежать переполнения.
        //
        //  Итак, A(0) = 0 , B(0) = 1,
        //        r(1) = -(a+b)*x / (a+1)
        //        A(1) = A(0) + r(1)*A(-1) = r(1) = 1
        //        B(1) = B(0) + r(1)*B(-1) = 1
        {
            double old_bta = 0, factor = 1;
            double A0 = 0, A1 = 1, B0 = 1, B1 = 1;
            double bta = 1, am = a, ai = a;
            double iter = 0, r;

            do
            {
                // часть цикла, вычисляющая нечетные подходящие дроби
                // начинаем с i = 1, iter = 0
                ai += 1;              // = a+i
                r = -am * (am + b) * x / ((ai - 1) * ai);
                /* пересчет A и B в два шага           */
                A0 = (A1 + r * A0) * factor;  /* i НЕчетно */
                B0 = (B1 + r * B0) * factor;
                // часть цикла, вычисляющая нечетные подходящие дроби
                // начинаем с i = 2, iter = 1
                am += 1;
                iter += 1;
                ai += 1;
                r = iter * (b - iter) * x * factor / ((ai - 1) * ai);
                A1 = A0 + r * A1;     /* i четно, A0 и B0 уже шкалированы */
                B1 = B0 + r * B1;
                old_bta = bta;
                factor = 1 / B1;
                bta = A1 * factor;
            } while (Math.Abs(old_bta) != Math.Abs(bta));
            return bta * Math.Exp(a * Math.Log(x) + b * Math.Log(1 - x) - logBeta(a, b)) / a;
        }/*BetaFraction*/

        public delegate double MathFunction(double x);
        public static double integral(MathFunction function, double a, double b, double eps = 1E-6)
        {
            double result = 0;
            double dx = eps;
            double x = a;
            while (x < b)
            {
                result += dx * function(x + dx / 2);
                x += dx;
            }
            return result;
        }

        public static int doubleStringPositionAndDigitsMinDifference = 5;
        public static string DefaultDoubleStringFormat = "{0,8:f3} ";
        // исключаем инициализацию, поэтому изначально формат должен совпадать с "{0,"+doubleStringPositions+":f"+doubleStringDigits+"} "
        static int doubleStringPositions = 8;
        public static int DoubleStringPositions
        {
            get { return doubleStringPositions; }
            set
            {
                doubleStringPositions = value;
                if (doubleStringPositions < doubleStringDigits + doubleStringPositionAndDigitsMinDifference)
                    doubleStringPositions = doubleStringDigits + doubleStringPositionAndDigitsMinDifference;
                UpdateDefaultDoubleStringFormat();
            }
        }

        static int doubleStringDigits = 3;
        public static int DoubleStringDigits
        {
            get { return doubleStringDigits; }
            set
            {
                doubleStringDigits = value;
                if (doubleStringPositions < doubleStringDigits + doubleStringPositionAndDigitsMinDifference)
                    doubleStringDigits = doubleStringPositions - doubleStringPositionAndDigitsMinDifference;
                UpdateDefaultDoubleStringFormat();
            }
        }

        public static void ChangeDoubleStringFormat(int positions, int digits)
        {
            doubleStringPositions = positions;
            doubleStringDigits = digits;
            UpdateDefaultDoubleStringFormat();
        }

        private static void UpdateDefaultDoubleStringFormat()
        {
            DefaultDoubleStringFormat = "{0," + doubleStringPositions + ":f" + doubleStringDigits + "} ";
        }

        public static string DoubleToString(double x, string format = "default")
        {
            if (String.IsNullOrEmpty(format)) return x.ToString() + " ";
            else if (format == "default") format = DefaultDoubleStringFormat;
            return String.Format(format, x);
        }

        public static string DoubleToIntString(double x)
        {
            return DoubleToString(x, "{0," + doubleStringPositions + ":f0} ");
        }

        public static string IntToString(int i)
        {
            return DoubleToIntString(i);
        }

        public static double StringToDouble(string s)
        {
            s = CorrectDoubleSeparator(s);
            return double.Parse(s);
        }

        public static string CorrectDoubleSeparator(string text)
        {
            string separatorMustBe = System.Globalization.CultureInfo.CurrentCulture.NumberFormat.NumberDecimalSeparator;
            string separatorToChange = ".";
            if (separatorMustBe == separatorToChange) separatorToChange = ",";
            return text.Replace(separatorToChange, separatorMustBe);
        }

        public static string IntToSmartString(int i)
        {
            if (i < 1000) return i.ToString();
            else if (i < 1000000) return (i / 1000).ToString() + " тыс.";
            else if (i < 10000000)
            {
                int m = i / 1000000;
                return m.ToString() + "," + (i / 100000 - m * 10).ToString() + " млн.";
            }
            else return (i / 1000000).ToString() + " млн.";
        }
    }
}
