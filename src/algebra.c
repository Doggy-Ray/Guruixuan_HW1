#include "algebra.h"
#include <stdio.h>
#include <math.h>

#define EPSILON 1e-6

void print_matrix(Matrix m) 
{
    for (int i = 0; i < m.rows; i++) {
        for (int j = 0; j < m.cols; j++) {
            double value = m.data[i][j];
            if (fabs(value) < EPSILON)
                value = 0.0;  // 将接近零的值设为绝对零
            printf("%-4.2f\t", m.data[i][j]); 
        }
        printf("\n");
    }
}

Matrix create_matrix(int row, int col) 
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    for (int i = 0; i < row; i++) 
    {
        for (int j = 0; j < col; j++) 
        {
            m.data[i][j] = 0.0;
        }
    }
    return m;
}

Matrix add_matrix(Matrix a, Matrix b) 
{
    if (a.rows != b.rows || a.cols != b.cols) 
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < a.cols; j++) 
        {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b) 
{
    if (a.rows != b.rows || a.cols != b.cols) 
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < a.cols; j++) 
        {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b) 
{
    if (a.cols != b.rows) 
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < b.cols; j++) 
        {
            double sum = 0.0;
            for (int k = 0; k < a.cols; k++) 
            {
                sum += a.data[i][k] * b.data[k][j];
            }
            result.data[i][j] = sum;
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k) 
{
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < a.cols; j++) 
        {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a) 
{
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < a.cols; j++) 
        {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a)
{
    int n = a.rows;
    if (n == 1) 
    {
        return a.data[0][0];
    }
    double det = 0.0;
    for (int j = 0; j < n; j++) 
    {
        Matrix sub = create_matrix(n - 1, n - 1);
        for (int row = 1; row < n; row++) 
        {
            int col_sub = 0;
            for (int col = 0; col < n; col++) 
            {
                if (col == j) continue;
                sub.data[row - 1][col_sub++] = a.data[row][col];
            }
        }
        //奇偶换个符号
        double if_negative = ((j % 2) == 0 ? 1 : -1) * det_matrix(sub);
        det += a.data[0][j] * if_negative;
    }
    return det;
}

Matrix inv_matrix(Matrix a) 
{
    //不知道有没有必要判断是否是方阵
    if (a.rows != a.cols) 
    {
        printf("Error: need square matrix!!!.\n");
        return create_matrix(0, 0);
    }
    double det = det_matrix(a);
    int n = a.rows;
    Matrix adj = create_matrix(n, n);
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            Matrix sub = create_matrix(n - 1, n - 1);
            for (int row = 0; row < n; row++) 
            {
                if (row == i) continue;
                int sub_row = row - 1;
                if (row < i) sub_row++;
                for (int col = 0; col < n; col++) 
                {
                    if (col == j) continue;
                    int sub_col = col - 1;
                    if (col < j) sub_col++;
                    sub.data[sub_row][sub_col] = a.data[row][col];
                }
            }
            double if_negetive = ((i + j) % 2 == 0 ? 1 : -1) * det_matrix(sub);
            adj.data[j][i] = if_negetive; 
        }
    }
    return scale_matrix(adj, 1.0 / det);
}

int rank_matrix(Matrix a) 
{
    int rank = 0;
    int rows = a.rows;
    int cols = a.cols;
    double mat[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
    for (int i = 0; i < rows; i++) 
    {
        for (int j = 0; j < cols; j++) 
        {
            mat[i][j] = a.data[i][j];
        }
    }
    for (int col = 0; col < cols && rank < rows; col++) 
    {
        int pivot_row = -1;
        for (int i = rank; i < rows; i++) 
        {
            if (fabs(mat[i][col]) > EPSILON) 
            {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;
        if (pivot_row != rank) 
        {
            for (int j = col; j < cols; j++) 
            {
                double temp = mat[rank][j];
                mat[rank][j] = mat[pivot_row][j];
                mat[pivot_row][j] = temp;
            }
        }
        for (int i = rank + 1; i < rows; i++) 
        {
            double factor = mat[i][col] / mat[rank][col];
            for (int j = col; j < cols; j++) 
            {
                mat[i][j] -= factor * mat[rank][j];
            }
        }
        rank++;
    }
    return rank;
}

double trace_matrix(Matrix a) 
{
    if (a.rows != a.cols) 
    {
        printf("Error: need square matrix!!!.\n");
        return 0;
    }
    double trace = 0.0;
    for (int i = 0; i < a.rows; i++) 
    {
        trace += a.data[i][i];
    }
    return trace;
}