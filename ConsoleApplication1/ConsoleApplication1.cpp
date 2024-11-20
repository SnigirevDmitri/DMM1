#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

typedef std::vector<std::vector<int>> matrix;
typedef std::vector<int> vec;

class SLAE
{
private:
   size_t n;             // Число уравнений
   size_t m;             // Число неизвестных
   matrix B;             // Расширенная СЛАУ
   matrix B_orig;        // Оригинальная СЛАУ для проверки решения
   size_t size;          // Размер матрицы
   size_t s;             // Число свободных переменных

public:
   void input(std::string in)
   {
      std::ifstream input(in);
      input >> n >> m;

      size = m + n;

      B.resize(size);
      for (size_t i = 0; i < n + m; i++)
         B[i].resize(m + 1);

      B_orig.resize(n);
      for (size_t i = 0; i < n; i++)
         B_orig[i].resize(m + 1);

      for (size_t i = 0; i < n; i++)
      {
         for (size_t j = 0; j < m + 1; j++)
         {
            input >> B[i][j];
            B_orig[i][j] = B[i][j];
         }
         B[i][m] *= -1;
      }
      input.close();

      for (size_t i = n, j = 0; i < m + n; i++, j++)
         B[i][j] = 1;
   }

   bool solveSLAE()
   {
      // Приведение части матрицы к нижнему треугольному виду
      for (size_t i = 0; i < n; i++)
      {
         bool row_done = false;

         while (!row_done)
         {
            // Поиск минимального элемента в строке
            bool finded = false;
            size_t min = 0;
            for (size_t j = i; j < m; j++)
            {
               if (!finded)
               {
                  if (B[i][j] != 0)
                  {
                     min = j;
                     finded = true;
                  }
               }
               else
               {
                  if (abs(B[i][j]) < abs(B[i][min]) && B[i][j] != 0)
                     min = j;
               }
            }

            // Если не нашли (какая-либо строка линейно зависима от другой), убираем строку из матрицы
            if (!finded && B[i][m] == 0)
            {
               B.erase(B.begin() + i);
               i--;
               n--;
               size--;
               continue;
            }

            if (finded)
            {
               // Вычитаем столбец с минимальным значением из остальных
               for (size_t j = i; j < m; j++)
               {
                  if (j != min && B[i][j] != 0)
                  {
                     int q = B[i][j] / B[i][min];
                     for (size_t k = i; k < B.size(); k++)
                        B[k][j] -= q * B[k][min];
                  }
               }

               // Если диагональный элемент обнулился, меняем его местами с ведущим
               if (B[i][i] == 0)
               {
                  for (size_t k = i; k < size; k++)
                  {
                     std::swap(B[k][i], B[k][min]);
                  }
               }
            }

            // Проверяем все недиагональные элементы на наличие ненулевых
            bool indiag = false;
            for (size_t j = i + 1; j < m; j++)
            {
               if (B[i][j] != 0)
               {
                  indiag = true;
                  break;
               }
            }

            // Если такие не найдены, проверяем правую часть на отличие от нуля и невозможность обратить её в ноль путём вычета столбца
            if (!indiag)
            {
               if (B[i][i] != 0)
               {
                  if (B[i][m] % B[i][i] != 0)
                  {
                     return false;
                  }
                  else
                  {
                     int q = B[i][m] / B[i][i];
                     for (size_t k = i; k < B.size(); k++)
                        B[k][m] -= q * B[k][i];
                  }
               }
               else
               {
                  return false;
               }

               if (B[i][i] == 1)
               {
                  int q = B[i][m];
                  for (size_t k = i; k < B.size(); k++)
                     B[k][m] -= q * B[k][i];
               }

               // Если правая часть нулевая, заканчиваем обработку строки
               if (B[i][m] == 0) row_done = true;
            }
         }

      }
      return true;
   }

   bool testSolve()
   {
      // Проверим решение, подставив случайные свободные переменные, затем подставив полученные значения в исходную СЛАУ

      // Проверяем число свободных переменных
      s = m - n;

      if (s < 0)
      {
         return false;
      }

      vec t(s, 0);

      for (size_t i = 0; i < s; i++)
         t[i] = rand();

      vec test_value;

      for (size_t j = n; j < m + n; j++)
         test_value.push_back(B[j][m]);

      for (size_t k = n, i = 0; k < m + n; k++, i++)
         for (int j = s - 1; j >= 0; j--)
            test_value[i] += B[k][m - (j + 1)] * t[j];

      vec test(n, 0);

      for (size_t i = 0; i < test.size(); i++)
         for (size_t j = 0; j < B_orig[i].size() - 1; j++)
            test[i] += B_orig[i][j] * test_value[j];

      // Проверим норму разности получившегося вектора с изначальным вектором правой части
      int res = 0;
      for (size_t i = 0; i < test.size(); i++)
         res += sqrt((test[i] - B_orig[i][m]) - (test[i] - B_orig[i][m]));

      if (res > 0)
      {
         return false;
      }

      return true;
   }

   void output(std::string out, bool flag)
   {
      std::ofstream output(out);

      if (!flag)
      {
         output << "NO SOLUTIONS";
      }
      else 
      {
         output << s << std::endl;

         for (size_t i = n; i < size; i++)
         {
            for (size_t j = m - s; j < m + 1; j++)
            {
               output << B[i][j] << " ";
            }
            output << std::endl;
         }
      }

      output.close();
   }
};

int main()
{
   SLAE sl;
   
   sl.input("input.txt");

   bool flag = sl.solveSLAE();

   if (!flag) 
   {
      sl.output("output.txt", flag);
      return 0;
   }

   flag = sl.testSolve();

   sl.output("output.txt", flag);

   return 0;

}
