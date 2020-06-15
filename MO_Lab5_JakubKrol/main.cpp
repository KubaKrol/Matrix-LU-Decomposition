//Program działa tylko na macierzach kwadratowych.

#include<iostream>

#define MatrixSize 4

using namespace std;

void PrintMatrix(double matrixToPrint[MatrixSize][MatrixSize])
{
    int n = MatrixSize;

    for (int i = 0; i < n; i++)
	{
        for (int j = 0; j < n; j++)
		{
            cout << matrixToPrint[i][j] << " ";
        }
        cout << endl;
    }
}

void ClearMatrix(double matrixToClear[MatrixSize][MatrixSize])
{
	int n = MatrixSize;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) 
		{
			matrixToClear[i][j] = 0;
		}
	}
}

void SwapMatrixRow(double a[MatrixSize][MatrixSize], int row1, int row2)
{
	double tmp;
	int n = MatrixSize;

	for (int j = 0; j < n; j++)
	{
		tmp = a[row1][j];
		a[row1][j] = a[row2][j];
		a[row2][j] = tmp;
	}
}

//Procedura częściowego wyboru elementu podstawowego (partial pivoting -  zamiana kolejności wierszy)
//Macierz jest diagonalnie słaba
void PartialPivoting(double a[MatrixSize][MatrixSize], int diagonal)
{
	double maxValue = 0;
	double currentValue;
	bool found = false;
	int iMax = 0;
	int n = MatrixSize;

	for (int i = diagonal; i < n; i++)
	{
		currentValue = fabs(a[i][diagonal]);

		if (currentValue > maxValue)
		{
			maxValue = currentValue;
			iMax = i;
			found = true;
		}
	}

	if (found)
	{
		SwapMatrixRow(a, iMax, diagonal);
	}
	else
	{
		cout << "No value larger than 0 in a column" << endl;
	}
}

void SolveXAndY(double l[MatrixSize][MatrixSize], double u[MatrixSize][MatrixSize], double vector[MatrixSize])
{
	double x[MatrixSize], y[MatrixSize];
	double temp;
	int n = MatrixSize;

	for (int i = 0; i < n; i++) 
	{
		x[i] = 0;
		y[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		temp = 0.0;

		for (int j = 0; j < i; j++)
		{
			temp += l[i][j] * y[j];
		}

		y[i] = vector[i] - temp;
	}

	for (int i = n - 1; i >= 0; i--)
	{
		temp = 0.0;

		for (int j = i + 1; j < n; j++)
		{
			temp += u[i][j] * x[j];
		}

		x[i] = (1 / u[i][i]) * (y[i] - temp);
	}

	cout << "Rozwiazanie" << endl;

	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
}

//Dekompozycja macierzy A na macierz L oraz U.
void LUDecomposition(double a[MatrixSize][MatrixSize], double l[MatrixSize][MatrixSize], double u[MatrixSize][MatrixSize])
{
	int i, j, k;
	int matrixSize = MatrixSize;
	double ax;

	for (k = 0; k < matrixSize - 1; k++)
	{
		for (i = k + 1; i < matrixSize; i++)
		{
			//Procedura wyboru elementu podstawowego
			if (a[k][k] == 0.0) 
			{
				PartialPivoting(a, k);
			}

			ax = a[i][k];

			for (j = k; j < matrixSize; j++)
			{
				a[i][j] = a[i][j] - ((a[k][j] * ax) / a[k][k]);
			}

			l[i][k] = ax / a[k][k];
		}
	}

	//Zbuduj macierz U
	//Okazuje się, że macierz U to dokładnie ta sama macierz, jaką uzyskujemy przy eliminacji Gaussa.
	for (int m = 0; m < matrixSize; m++)
	{
		for (int n = 0; n < matrixSize; n++)
		{
			if (n < m)
			{
				//Wypełnij dolną część macierzy U wartościami równymi 0;
				u[m][n] = 0;
			}
			else
			{
				u[m][n] = a[m][n];
			}
		}
	}

	//Wypełnij przekątną macierzy L wartościami równymi 1;
	for (int i = 0; i < matrixSize; ++i)
	{
		l[i][i] = 1;
	}
}

int main()
{
	double L[MatrixSize][MatrixSize], U[MatrixSize][MatrixSize];

	//Wyczyść macierze L oraz U
	ClearMatrix(L);
	ClearMatrix(U);

	double A[MatrixSize][MatrixSize]
    {
		{1, -20, 30, -4},
        {2, -40, -6, 50},
		{9, -180, 11, -12},
		{-16, 15, -140, 13}
    };

	double Vector[MatrixSize]{ 35, 104, -366, -354 };

	int n = MatrixSize;

    PrintMatrix(A);

    LUDecomposition(A, L, U);
    cout << "L Decomposition is as follows..." << endl;

    PrintMatrix(L);
    cout << "U Decomposition is as follows..." << endl;
    PrintMatrix(U);

	cout << "Solution:" << endl;

	SolveXAndY(L, U, Vector);
    
    system("pause");

    return 0;
}