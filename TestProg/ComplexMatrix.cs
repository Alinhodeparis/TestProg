using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;

namespace ODIS.AMM
{
    class ComplexMatrixRow
    {
        //private List<double> Items = null;
        private Complex[] items = null;

        private int size = 1;
        public int Size { get { return size; } }

        public ComplexMatrixRow(int size)
        {
            SetSize(size);
            items = new Complex[size];
        }

        private void SetSize(int size)
        {
            if ((size > 0) && (size <= Matrix.MaxSize)) this.size = size;
            else throw new Exception("Invalid size of matrix row");
        }

        public Complex this[int index]
        {
            get
            {
                if ((index < 1) || (index > size)) throw new Exception("Invalid index in matrix row");
                else return items[index - 1];
            }
            set
            {
                if ((index < 1) || (index > size)) throw new Exception("Invalid index in matrix row");
                else items[index - 1] = value;
            }
        }


       /* internal MatrixRow Sort()
        {
            MatrixRow result = new MatrixRow(this.Size);
            var z = this.items.OrderBy(x => x);
            int i = 0;
            foreach (double x in z)
            {
                result.items[i] = x;
                i++;
            }
            return result;
        }*/
    }

    public class ComplexMatrix
    {
        public const int MaxSize = 1000;
        private ComplexMatrixRow[] rowsList = null;
        private int rows = 1;
        public int Rows { get { return rows; } }
        private int cols = 1;
        public int Cols { get { return cols; } }

        public ComplexMatrix(int rows, int cols)
        {
            SetSizes(rows, cols);
            rowsList = CreateRowsList(this.Rows, this.Cols);
        }

        /*public static Matrix LoadFromFile(string filepath)
        {
            if (File.Exists(filepath))
            {
                Matrix result = null;
                List<string> lines = new List<string>();
                string line="";
                StreamReader sr = new StreamReader(filepath);
                while ((line = sr.ReadLine()) != null) lines.Add(line);
                sr.Close();
                int rows = lines.Count;
                int cols = 0;
                for (int i = 0; i < rows; i++)
                {
                    string[] data = lines[i].Split(new char[]{'\t', ' '}, StringSplitOption.RemoveEmptyEntries);
                    if (i == 0)
                    {
                        cols = data.Length;
                        result = new Matrix(rows, cols);
                    }
                    for (int j = 0; j < cols; j++)
                        result[i + 1, j + 1] = StringToDouble(data[j]);
                }
                return result;
            }
            else throw new Exception("file not found:\r\n" + filepath);
        }*/

       /* private static double StringToDouble(string s)
        {
            s = CorrectDoubleSeparator(s);
            return double.Parse(s);
        }

        private static string CorrectDoubleSeparator(string text)
        {
            string separatorMustBe = System.Globalization.CultureInfo.CurrentCulture.NumberFormat.NumberDecimalSeparator;
            string separatorToChange = ".";
            if (separatorMustBe == separatorToChange) separatorToChange = ",";
            return text.Replace(separatorToChange, separatorMustBe);
        }
        */
       
        private ComplexMatrixRow[] CreateRowsList(int rows, int cols)
        {
            ComplexMatrixRow[] result = new ComplexMatrixRow[rows];
            for (int i = 0; i < rows; i++)
                result[i] = new ComplexMatrixRow(cols);
            return result;
        }

        private void SetSizes(int rows, int cols)
        {
            if ((rows < 1) || (cols < 1) || (rows > MaxSize) || (cols > MaxSize)) throw new Exception("Invalid matrix size");
            else
            {
                this.rows = rows;
                this.cols = cols;
            }
        }

        public Complex this[int index1, int index2]
        {
            get
            {
                if ((index1 < 1) || (index1 > Rows)) throw new Exception("Invalid index of matrix row");
                else return rowsList[index1 - 1][index2];
            }
            set
            {
                if ((index1 < 1) || (index1 > Rows)) throw new Exception("Invalid index of matrix row");
                else rowsList[index1 - 1][index2] = value;
            }
        }


        public void Resize(int rows, int cols)
        {
            ComplexMatrixRow[] newRowsList = CreateRowsList(rows, cols);
            int minRows = (rows > this.Rows) ? this.Rows : rows;
            int minCols = (cols > this.Cols) ? this.Cols : cols;
            for (int i = 1; i <= minRows; i++)
                for (int j = 1; j <= minCols; j++)
                    newRowsList[i - 1][j] = this[i, j];
            SetSizes(rows, cols);
            rowsList = newRowsList;
        }

        /// <summary>
        /// Обращение матрицы
        /// </summary>
        /// <returns>Возвращает обратную матрицу</returns>
        public ComplexMatrix Inversion() // здесь также автоматически вычисляется определитель - но это будем использовать позже: detA: Real
        {
            /* i,j,k,irev: Integer;
              m,c,s: Real;
              ind: Boolean;*/
            if (Rows != Cols) return null;
            int size = Rows;
            ComplexMatrix original = Copy();
            ComplexMatrix result = Get1Matrix(size);
            bool ind = true; // какой-то индикатор
            int i = 1;
            int irev = 0;
            int j, k; Complex c, detA;
            while (ind && (i < size))
            {
                k = i;
                double m = Complex.Abs(original[i, i]);
                for (j = i + 1; j <= size; j++)
                    if (m < Complex.Abs(original[j, i]))
                    {
                        m = Complex.Abs(original[j, i]);
                        k = j;
                    }
                if (m < 1.0E-30) ind = false;
                if (ind)
                {
                    if (i != k)
                    {
                        irev++;
                        for (j = i; j <= size; j++)
                        {
                            c = original[i, j];
                            original[i, j] = original[k, j];
                            original[k, j] = c;
                        }
                        for (j = 1; j <= size; j++)
                        {
                            c = result[i, j];
                            result[i, j] = result[k, j];
                            result[k, j] = c;
                        }
                    }
                    for (j = i + 1; j <= size; j++)
                    {
                        for (k = i + 1; k <= size; k++)
                            original[j, k] = original[j, k] - original[i, k] * original[j, i] / original[i, i];
                        for (k = 1; k <= size; k++)
                            result[j, k] = result[j, k] - result[i, k] * original[j, i] / original[i, i];
                    }
                }
                i++;
            }
            if (Complex.Abs(original[size, size]) < 1.0E-30) ind = false;
            if ((irev % 2) == 0) detA = 1;
            else detA = -1;
            if (ind)
                for (i = 1; i <= size; i++) detA = detA * original[i, i];
            else detA = 0;
            if (ind)
            {
                for (i = 1; i <= size; i++) result[size, i] = result[size, i] / original[size, size];
                for (i = size - 1; i >= 1; i--)
                    for (j = 1; j <= size; j++)
                    {
                        Complex s = 0;
                        for (k = i + 1; k <= size; k++) s = s - original[i, k] * result[k, j];
                        result[i, j] = (result[i, j] + s) / original[i, i];
                    }
            }
            //if (!ind) OutMessage("ОБРАЩЕНИЕ ВЫРОЖДЕННОЙ МАТРИЦЫ");
            return result;
        }
        
        /// <summary>
        /// Создает и возвращает копию матрицы
        /// </summary>
        /// <returns></returns>
        public ComplexMatrix Copy()
        {
            ComplexMatrix result = new ComplexMatrix(Rows, Cols);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    result[i, j] = this[i, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает копию действительной матрицы в комплексном типе
        /// </summary>
        /// <returns></returns>
        public static ComplexMatrix GetFromMatrix(Matrix A)
        {
            ComplexMatrix result = new ComplexMatrix(A.Rows, A.Cols);
            for (int i = 1; i <= A.Rows; i++)
                for (int j = 1; j <= A.Cols; j++)
                    result[i, j] = A[i, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает единичную матрицу размера size
        /// </summary>
        /// <param name="size"></param>
        /// <returns></returns>
        public static ComplexMatrix Get1Matrix(int size)
        {
            ComplexMatrix result = new ComplexMatrix(size, size);
            for (int i = 1; i <= size; i++)
                for (int j = 1; j <= size; j++)
                    if (i == j) result[i, j] = 1;
                    else result[i, j] = 0;
            return result;
        }

        /// <summary>
        /// Создает и возвращает нулевую матрицу
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        /// <returns></returns>
        public static ComplexMatrix Get0Matrix(int rows, int cols)
        {
            ComplexMatrix result = new ComplexMatrix(rows, cols);
            for (int i = 1; i <= rows; i++)
                for (int j = 1; j <= cols; j++)
                    result[i, j] = 0;
            return result;
        }

        /// <summary>
        /// Создает и возвращает диагональную матрицу на основе вектора-столбца или строки.
        /// Если исходная матрица имеет более 1 строк и столбцов, то в качестве вектора берется первая строка этой матрицы
        /// </summary>
        /// <param name="A">исходный вектор</param>
        /// <returns></returns>
        public static ComplexMatrix GetDiagonalMatrix(Matrix A)
        {
            int size = A.Cols;
            bool AisRowVector = true;
            if (size == 1)
            {
                size = A.Rows;
                AisRowVector = false;
            }
            ComplexMatrix result = ComplexMatrix.Get0Matrix(size, size);
            for (int i = 1; i <= size; i++)
                if (AisRowVector) result[i, i] = A[1, i];
                else result[i, i] = A[i, 1];
            return result;
        }

        /// <summary>
        /// Создает и возвращает транспонированную матрицу
        /// </summary>
        /// <returns></returns>
        public ComplexMatrix Transp()
        {
            ComplexMatrix result = new ComplexMatrix(Cols, Rows);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    result[j, i] = this[i, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает матрицу, равную сумме двух матриц
        /// </summary>
        /// <returns></returns>
        public ComplexMatrix Add(ComplexMatrix A)
        {
            ComplexMatrix result = new ComplexMatrix(Rows, Cols);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    result[i, j] = this[i, j] + A[i, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает матрицу, равную разнице двух матриц
        /// </summary>
        /// <returns></returns>
        public ComplexMatrix Sub(ComplexMatrix A)
        {
            ComplexMatrix result = new ComplexMatrix(Rows, Cols);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    result[i, j] = this[i, j] - A[i, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает произведение матрицы на другую матрицу
        /// </summary>
        /// <param name="M"></param>
        /// <returns></returns>
        public ComplexMatrix Multiply(ComplexMatrix M)
        {
            ComplexMatrix result = Get0Matrix(Rows, M.Cols);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= M.Cols; j++)
                    for (int k = 1; k <= Cols; k++)
                        result[i, j] += this[i, k] * M[k, j];
            return result;
        }

        /// <summary>
        /// Создает и возвращает матрицу, равную произведению матрицы на число
        /// </summary>
        /// <returns></returns>
        public ComplexMatrix MultiplyOnComplex(Complex x)
        {
            ComplexMatrix result = new ComplexMatrix(Rows, Cols);
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    result[i, j] = this[i, j] * x;
            return result;
        }

        /*public override string ToString()
        {
            return ToString(false, false, "");
        }

        const int headerWidth = 10;

        public string ToString(bool ShowRowNumbers=false, bool ShowColNumbers=false, string Header="")
        {
            string result = ""; // ! сделать через StringBuilder!
            if (ShowColNumbers) result = GetNumbersString(Header) + "\r\n";
            for (int i = 1; i <= Rows; i++)
            {
                string s = "";
                if (ShowRowNumbers) s = DoubleToString(i);
                result += s.PadRight(headerWidth);
                for (int j = 1; j <= Cols; j++) result += DoubleToString(this[i, j]);
                if (i < Rows) result += "\r\n";
            }
            return result;
        }

        public string GetNumbersString(string Header = "")
        {
            string s = Header;
            string result = s.PadLeft(headerWidth);
            for (int i = 1; i <= Cols; i++)
                result += DoubleToString(i);
            return result;
        }

        public static string DoubleToString(double x)
        {
            // было   return String.Format("{0,7} ", x); // было:{0,8:G} 
            return String.Format("{0,8} ", Math.Round(x,3));
        }*/

        /// <summary>
        /// Создает и возвращает в виде матрицы (вектор-строки) указанную строку матрицы
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public ComplexMatrix GetRow(int rowNumber)
        {
            ComplexMatrix result = new ComplexMatrix(1, Cols);
            for (int j = 1; j <= Cols; j++)
                result[1, j] = this[rowNumber, j];
            return result;
        }

        /*
        /// <summary>
        /// Создает и возвращает копию матрицы, в которой значения в строках отсортированы
        /// </summary>
        /// <returns></returns>
        public Matrix SortInRows()
        {
            Matrix result = new Matrix(Rows, Cols);
            for (int i = 1; i <= Rows; i++)
                result.rowsList[i] = this.rowsList[i].Sort();
            return result;
        }

        /// <summary>
        /// Возвращает сумму элементов строки rowNumber
        /// </summary>
        /// <param name="rowNumber"></param>
        /// <returns></returns>
        public double SumOfRowItems(int rowNumber)
        {
            double result = 0;
            for (int j = 1; j <= Cols; j++)
                result += this[rowNumber, j];
            return result;
        }
        */
        /// <summary>
        /// Заполняет матрицу значениями value
        /// </summary>
        /// <param name="value"></param>
        public void Fill(Complex value)
        {
            for (int i = 1; i <= Rows; i++)
                for (int j = 1; j <= Cols; j++)
                    this[i, j] = value;
        }

       /* public void FillRowAsDistribution(bool full = true)
        {
            double zeroProb = 0.3;
            double delta = 0.7;
            if (full) delta = 0;
            int digits = (int)Math.Log10(Cols) + 1;
            int decimals = (int)Math.Pow(10, digits);
            Random R = new Random();
            for (int i = 1; i <= Rows; i++)
            {
                double sum = 0;
                for (int j = 1; j <= Cols; j++)
                {
                    if (R.NextDouble() < zeroProb) this[i, j] = 0;
                    else this[i, j] = R.NextDouble();
                    sum += this[i, j];
                }
                int zhertva = R.Next(Cols) + 1;
                if (sum == 0) this[i, zhertva] = Math.Round(1 - R.NextDouble() * delta, digits);
                else
                {
                    sum = sum / (1 - R.NextDouble() * delta);
                    double sum2 = 0;
                    for (int j = 1; j <= Cols; j++)
                    {
                        if (j != zhertva)
                        {
                            this[i, j] = Math.Truncate(this[i, j] / sum * decimals) / decimals;
                            sum2 += this[i, j];
                        }
                    }
                    this[i, zhertva] = Math.Round((1 - sum2) * (1 - R.NextDouble() * delta), digits);
                    //if ((!full) && (R.NextDouble() < zeroProb)) this[i, zhertva] = 0;
                    if (this[i, 1] < 0) this[i, zhertva] = 0;
                }
            }
        }

        public void FillRandomFromValues(double[] values)
        {
            if ((values != null) && (values.Length > 0))
            {
                Random R = new Random();
                for (int i = 1; i <= Rows; i++)
                    for (int j = 1; j <= Cols; j++)
                        this[i, j] = values[R.Next(values.Length)];
            }
        } */
    }
}
