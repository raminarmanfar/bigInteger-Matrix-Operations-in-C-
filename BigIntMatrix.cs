using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace PIRLatticeGUI
{
    public class BigIntMatrix
    {
        static Random rnd = new Random(DateTime.Now.Millisecond);

        private BigInteger[,] Mat;
        private long rowsCount { get; set; }
        private long colsCount { get; set; }
        public BigIntMatrix(BigInteger[,] c1)
        {
            this.rowsCount = c1.GetLongLength(0);
            this.colsCount = c1.GetLongLength(1);
            Mat = new BigInteger[rowsCount, colsCount];
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                    Mat[i, j] = c1[i, j];
        }
        public BigIntMatrix(double[,] c1)
        {
            this.rowsCount = c1.GetLongLength(0);
            this.colsCount = c1.GetLongLength(1);
            Mat = new BigInteger[rowsCount, colsCount];
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                    Mat[i, j] = (BigInteger)c1[i, j];
        }
        public BigIntMatrix(long rowsCount, long colsCount)
        {
            this.rowsCount = rowsCount;
            this.colsCount = colsCount;
            Mat = new BigInteger[rowsCount, colsCount];
        }
        public BigIntMatrix(long rowsCount, long colsCount, BigInteger initialValue)
        {
            this.rowsCount = rowsCount;
            this.colsCount = colsCount;
            Mat = new BigInteger[rowsCount, colsCount];
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                    Mat[i, j] = initialValue;
        }

        //public enum MatrixType { SoftNoise, HardNoise, Random, RandomInvertible, RandomVectorMatrix,RandomInvertibleVector }
        public void GenerateSoftNoiseMatrix()
        {
            for (long i = 0; i < rowsCount; i++)
            {
                for (long j = 0; j < colsCount; j++)
                {
                    Mat[i, j] = rnd.Next(0, 2);
                    if (Mat[i, j] == 0) Mat[i, j] = -1;
                }
            }
        }
        public void GenerateHardNoiseMatrix(BigInteger q)
        {
            for (long i = 0; i < rowsCount; i++)
            {
                for (long j = 0; j < colsCount; j++)
                {
                    if (i == j) Mat[i, j] = q;
                    else
                    {
                        Mat[i, j] = rnd.Next(0, 2);
                        if (Mat[i, j] == 0) Mat[i, j] = -1;
                    }
                }
            }
        }

        public void GenerateRandomMatrix(BigInteger p, double l0)
        {
            long len = (long)Math.Ceiling(3.0 * l0 / 8.0) + 1;
            byte[] rndBytes = new byte[len];

            for (long i = 0; i < rowsCount; i++)
            {
                for (long j = 0; j < colsCount; j++)
                {
                    rnd.NextBytes(rndBytes);
                    BigInteger num = new BigInteger(rndBytes);
                    Mat[i, j] = MathClass.mod(num, p);
                }
            }
        }

        public void GenerateRandomVector(BigInteger p, double l0)
        {
            long len = (long)Math.Ceiling(3.0 * l0 / 8.0) + 1;
            byte[] rndBytes = new byte[len];

            for (long i = 0; i < rowsCount; i++)
            {
                rnd.NextBytes(rndBytes);
                BigInteger num = new BigInteger(rndBytes);
                Mat[i, i] = MathClass.mod(num, p);
            }
        }

        public BigIntMatrix GenerateRandomInvertibleMatrix(BigInteger p, double l0, bool isVectorMatrix)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);
            bool inversible = false;
            BigIntMatrix matInv = new BigIntMatrix(this.rowsCount, this.colsCount);
            do
            {
                if (isVectorMatrix) result.GenerateRandomVector(p, l0);
                else result.GenerateRandomMatrix(p, l0);
                try
                {
                    matInv = result.GaussJordanModInverse(p);
                    inversible = true;
                }
                catch(Exception ex)
                {
                    inversible = false;
                }
            } while (!inversible);

            Mat = result.Mat;
            return matInv;
        }

        public BigInteger GetElement(long Row, long Col)
        {
            if (Row > this.rowsCount || Col > this.colsCount)
                throw new System.ArgumentException("Given indices are out of matrix range", "original");
            else return this.Mat[Row, Col];
        }
        public void SetElement(long Row, long Col, BigInteger Value)
        {
            if (Row > this.rowsCount || Col > this.colsCount)
                throw new System.ArgumentException("Given indices are out of matrix range", "original");
            else this.Mat[rowsCount, Col] = Value;
        }
        public void SumElement(long Row, long Col, BigInteger Value)
        {
            if (Row > this.rowsCount || Col > this.colsCount)
                throw new System.ArgumentException("Given indices are out of matrix range", "original");
            else this.Mat[Row, Col] += Value;
        }

        public BigIntMatrix GetSubMatrix(long FromRow, long RowCount, long FromCol, long ColCount)
        {
            BigIntMatrix result;
            if (FromRow < 0 || FromRow > this.rowsCount || FromCol < 0 || FromCol > this.colsCount || FromRow + RowCount > this.rowsCount || FromCol + ColCount > this.colsCount)
            {
                throw new System.ArgumentException("Given indices are out of matrix range", "original");
            }
            else
            {
                result = new BigIntMatrix(RowCount, ColCount);
                for (long i = 0; i < RowCount; i++)
                {
                    for (long j = 0; j < ColCount; j++)
                    {
                        result.Mat[i, j] = this.Mat[i + FromRow, j + FromCol];
                    }
                }
            }
            return result;
        }

        public string Print(string Title)
        {
            string result = "\r\n>>>>>>>>>>>>>>> Matrix " + Title + "(" + rowsCount.ToString() + ", " + colsCount.ToString() + "):";
            for (long i = 0; i < rowsCount; i++)
            {
                result += "\r\n";
                for (long j = 0; j < colsCount; j++)
                {
                    result += Mat[i, j].ToString();
                    if (j + 1 < colsCount) result += "\t";
                }
            }
            return result + "\r\n";
        }
        public string Print()
        {
            string result = "\r\n>>>>>>>>>>>>>>> Matrix " + "(" + rowsCount.ToString() + ", " + colsCount.ToString() + "):";
            for (long i = 0; i < rowsCount; i++)
            {
                result += "\r\n";
                for (long j = 0; j < colsCount; j++)
                {
                    result += Mat[i, j].ToString();
                    if (j + 1 < colsCount) result += "\t";
                }
            }
            return result + "\r\n";
        }

        public BigIntMatrix Modulus(BigInteger p)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                {
                    result.Mat[i, j] = Mat[i, j] % p;
                    if (result.Mat[i, j] < 0) result.Mat[i, j] += p;
                }
            return result;
        }

        public static BigIntMatrix operator %(BigIntMatrix c1, BigInteger p)
        {
            BigIntMatrix result = new BigIntMatrix(c1.rowsCount, c1.colsCount);

            for (long i = 0; i < c1.rowsCount; i++)
                for (long j = 0; j < c1.colsCount; j++)
                {
                    result.Mat[i, j] = c1.Mat[i, j] % p;
                    if (result.Mat[i, j] < 0) result.Mat[i, j] += p;
                }
            return result;
        }

        public static BigIntMatrix operator +(BigIntMatrix c1, BigIntMatrix c2)
        {
            BigIntMatrix result = new BigIntMatrix(c1.rowsCount, c1.colsCount);
            if (c1.rowsCount != c2.rowsCount || c1.colsCount != c2.colsCount)
            {
                throw new System.ArgumentException("Matrices dimentions should be the same", "original");
            }
            else
            {
                for (long i = 0; i < c1.rowsCount; i++)
                    for (long j = 0; j < c1.colsCount; j++)
                        result.Mat[i, j] = c1.Mat[i, j] + c2.Mat[i, j];
            }
            return result;
        }

        public static BigIntMatrix operator -(BigIntMatrix c1, BigIntMatrix c2)
        {
            BigIntMatrix result = new BigIntMatrix(c1.rowsCount, c1.colsCount);
            if (c1.rowsCount != c2.rowsCount || c1.colsCount != c2.colsCount)
            {
                throw new System.ArgumentException("Matrices dimentions should be the same", "original");
            }
            else
            {
                for (long i = 0; i < c1.rowsCount; i++)
                    for (long j = 0; j < c1.colsCount; j++)
                        result.Mat[i, j] = c1.Mat[i, j] - c2.Mat[i, j];
            }
            return result;
        }

        public static BigIntMatrix operator *(BigIntMatrix c1, BigIntMatrix c2)
        {
            BigIntMatrix result = new BigIntMatrix(c1.rowsCount, c2.colsCount);
            if (c1.colsCount == c2.rowsCount)
            {
                for (int i = 0; i < result.rowsCount; i++)
                {
                    for (int j = 0; j < result.colsCount; j++)
                    {
                        result.Mat[i, j] = 0;
                        for (int k = 0; k < c1.colsCount; k++) // OR k<c2.rowsCount
                            result.Mat[i, j] += c1.Mat[i, k] * c2.Mat[k, j];
                    }
                }
            }
            else
            {
                throw new System.ArgumentException("First matrix row should be equal to the second matrix col", "original");
            }
            return result;
        }

        public static BigIntMatrix operator *(BigIntMatrix c1, BigInteger n)
        {
            BigIntMatrix result = new BigIntMatrix(c1.rowsCount, c1.colsCount);
            if (c1.colsCount == c1.rowsCount)
                for (int i = 0; i < result.rowsCount; i++)
                    for (int j = 0; j < result.rowsCount; j++)
                        result.Mat[i, j] = c1.Mat[i, j] * n;
            else
            {
                throw new System.ArgumentException("First matrix row should be equal to the second matrix col", "original");
            }
            return result;
        }

        public BigIntMatrix dotMul(BigIntMatrix c2)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);
            if (rowsCount != c2.rowsCount || colsCount != c2.colsCount)
            {
                throw new System.ArgumentException("Matrices dimentions should be the same", "original");
            }
            else
            {
                for (long i = 0; i < rowsCount; i++)
                    for (long j = 0; j < colsCount; j++)
                        result.Mat[i, j] = Mat[i, j] * c2.Mat[i, j];
            }
            return result;
        }

        public BigIntMatrix MatrixDuplicate()
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                    result.Mat[i, j] = Mat[i, j];
            return result;
        }

        public double[,] ToDouble()
        {
            double[,] result = new double[this.rowsCount, this.colsCount];
            for (long i = 0; i < rowsCount; i++)
                for (long j = 0; j < colsCount; j++)
                    result[i, j] = (double)this.Mat[i, j];
            return result;
        }

        public BigIntMatrix SwapCol(List<int> pList)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);
            for (int j = 0; j < this.rowsCount; j++)
            {
                for (int k = 0; k < this.colsCount; k++)
                {
                    long k2 = pList[k];
                    result.Mat[j, k2] = this.Mat[j, k];
                }
            }
            return result;
        }
        public BigIntMatrix GetMinor(long r, long c)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount - 1, colsCount - 1);
            long p = 0;
            for (long s = 0; s < rowsCount - 1; s++)
            {
                if (p == r) p++;
                long q = 0;
                for (long t = 0; t < colsCount - 1; t++)
                {
                    if (q == c) q++;
                    result.Mat[s, t] = Mat[p, q];
                    q++;
                }
                p++;
            }
            return result;
        }
        public BigIntMatrix RemoveRow(long r)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount - 1, colsCount);
            long p = 0;
            for (long s = 0; s < rowsCount - 1; s++)
            {
                if (p == r) p++;
                for (long t = 0; t < colsCount; t++)
                {
                    result.Mat[s, t] = Mat[p, t];
                }
                p++;
            }
            return result;
        }
        public BigIntMatrix RemoveCol(long c)
        {
            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount - 1);
            for (long s = 0; s < rowsCount; s++)
            {
                long q = 0;
                for (long t = 0; t < colsCount - 1; t++)
                {
                    if (q == c) q++;
                    result.Mat[s, t] = Mat[s, q];
                    q++;
                }
            }
            return result;
        }

        //---- Gauss Jordan Matrix Mod-Inverse ---
        public BigIntMatrix GaussJordanModInverse(BigInteger p)
        {
            BigIntMatrix workMat = new BigIntMatrix(2 * rowsCount, 2 * colsCount);
            //--- build work matrix ---
            long n = this.rowsCount;
            for (long i = 0; i < n; i++)
            {
                for (long j = 0; j < n; j++)
                {
                        workMat.Mat[i, j] = this.Mat[i, j];
                }
                workMat.Mat[i, i + n] = 1;
            }

            // --- partial pivoting ---
            for (long i = n - 1; i > 0; i--)
            {
                if (workMat.Mat[i - 1, 0] < workMat.Mat[i, 0])
                {
                    for (long j = 0; j < 2 * n; j++)
                    {
                        BigInteger d = workMat.Mat[i, j];
                        workMat.Mat[i, j] = workMat.Mat[i - 1, j];
                        workMat.Mat[i - 1, j] = d;
                    }
                }
            }

            // --- reducing to diagnoal matrix ---
            for (int i = 0; i < n; i++)
            {
                BigInteger mi = MathClass.EEA(workMat.Mat[i, i], p);
                for (int j = 0; j < 2 * n; j++)
                {
                    if (i != j)
                    {
                        BigInteger d = workMat.Mat[j, i] * mi;
                        for (int k = 0; k < 2 * n; k++)
                        {
                            workMat.Mat[j, k] = MathClass.mod(workMat.Mat[j, k] - workMat.Mat[i, k] * d, p);
                        }
                    }
                }
            }

            BigIntMatrix result = new BigIntMatrix(rowsCount, colsCount);

            //--- reducing to unit matrix ---
            for (int i = 0; i < n; i++)
            {
                BigInteger d = MathClass.EEA(workMat.Mat[i, i], p);
                for (int j = 0; j < 2 * n; j++)
                {
                    workMat.Mat[i, j] = workMat.Mat[i, j] * d;
                }

                for (int j = 0; j < n; j++)
                {
                    //--- check whether matrix is inversible or not ---
                    if (i != j && workMat.Mat[i, j] != 0)
                    {
                        throw new System.ArgumentException("ERROR: Matrix is not inversible...", "Matrix inversion error");
                    }

                    result.Mat[i, j] = MathClass.mod(workMat.Mat[i, j + n], p);
                }
            }

            return result;
        }

        public enum AppendMode { Beside, Under }
        public BigIntMatrix Append(BigIntMatrix c2, AppendMode appendMode)
        {
            BigIntMatrix result;
            if (appendMode == AppendMode.Beside)
            {
                result = new BigIntMatrix(this.rowsCount, this.colsCount + c2.colsCount);
                if (this.rowsCount == c2.rowsCount)
                {
                    for (int i = 0; i < this.rowsCount; i++)
                    {
                        for (int j = 0; j < this.colsCount; j++)
                            result.Mat[i, j] = this.Mat[i, j];
                        for (int j = 0; j < c2.colsCount; j++)
                            result.Mat[i, j + c2.colsCount] = c2.Mat[i, j];
                    }
                }
                else
                {
                    throw new System.ArgumentException("Matrices have no equal rows", "original");
                }
            }
            else
            {
                result = new BigIntMatrix(this.rowsCount + c2.rowsCount, this.colsCount);
                if (this.colsCount == c2.colsCount)
                {
                    for (int j = 0; j < this.colsCount; j++)
                    {
                        for (int i = 0; i < this.rowsCount; i++)
                            result.Mat[i, j] = this.Mat[i, j];
                        for (int i = 0; i < c2.rowsCount; i++)
                            result.Mat[i + this.rowsCount, j] = c2.Mat[i, j];
                    }
                }
                else
                {
                    throw new System.ArgumentException("Matrices have no equal columns", "original");
                }
            }
            return result;
        }
    }
}
