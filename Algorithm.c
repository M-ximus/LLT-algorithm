#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

const double PRES = 0.000000000000001;

enum errors{
    E_ERROR   = -1,
    E_BADARGS = -2
};

typedef struct sq_matrix{
    double** data;
    int      metrics;
} sq_matrix;

typedef struct vector{
    double* arr;
    int     dim;
} vector;

int get_matrix_from_file(sq_matrix* matr, const char* file_name);
int matrix_destr(sq_matrix* matr);
int init_vector(vector* vec, int dim);
int vector_destr(vector* vec);
int mul_sq_matrix_to_vector(sq_matrix* matr, vector* vec, vector* res);
int my_sqrt(double num, double* res);
int calc_L(sq_matrix* matr, sq_matrix* res);
int sq_matrix_trans(sq_matrix* matr, sq_matrix* res);
int print_sq_matrix(sq_matrix* matr);
int solve_eq(sq_matrix* L, vector* f, vector* res);
int equal(double first, double second);
int div_vectors(vector* res, vector* src1, vector* src2);
int euclid_norm(vector* vec, double* res);

// A*x = f
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        printf("[main] Bad num of arguments\n");
        exit(EXIT_FAILURE);
    }

    sq_matrix A;
    int ret = get_matrix_from_file(&A, argv[1]);
    if (ret < 0)
    {
        printf("[main] Getting matrix A from file error\n");
        exit(EXIT_FAILURE);
    }

    vector x;
    ret = init_vector(&x, A.metrics);
    if (ret < 0)
    {
        printf("[main] creating x vector error\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < x.dim; i++)
        x.arr[i] = (double) (i + 1);

    vector f; // heterogeneity
    ret = mul_sq_matrix_to_vector(&A, &x, &f);
    if (ret < 0)
    {
        printf("[main] get heterogeneity error\n");
        exit(EXIT_FAILURE);
    }

    sq_matrix L;
    ret = calc_L(&A, &L);
    if (ret < 0)
    {
        printf("[main] Find L error\n");
        exit(EXIT_FAILURE);
    }

    //print_sq_matrix(&L);

    vector result;
    ret = solve_eq(&L, &f, &result);
    if(ret < 0)
    {
        printf("[main] Error in solving equation\n");
        exit(EXIT_FAILURE);
    }

    vector div;
    ret = div_vectors(&div, &x, &result);
    if (ret < 0)
    {
        printf("[main] x - result error\n");
        exit(EXIT_FAILURE);
    }

    //for (int i = 0; i < div.dim; i++)
        //printf("%lg\n", div.arr[i]);

    double norm = 0.0;
    ret = euclid_norm(&div, &norm);
    if (ret < 0)
    {
        printf("[main] calc euclidus norm of div error\n");
        exit(EXIT_FAILURE);
    }

    printf("norm of div = %lg\n", norm);

    matrix_destr(&L);
    matrix_destr(&A);
    vector_destr(&x);
    vector_destr(&f);
    vector_destr(&result);
    vector_destr(&div);

    return 0;
}

int euclid_norm(vector* vec, double* res)
{
    if (vec == NULL || res == NULL)
    {
        printf("[euclid_norm] Bad args\n");
        return E_BADARGS;
    }

    double sum = 0.0;
    for (int i = 0; i < vec->dim; i++)
        sum += vec->arr[i] * vec->arr[i];

    int ret = my_sqrt(sum, res);
    if (ret < 0)
    {
        printf("[euclid_norm] Bad sqrt operation\n");
        return E_ERROR;
    }

    return 0;
}

int div_vectors(vector* res, vector* src1, vector* src2)
{
    if (res == NULL || src1 == NULL || src2 == NULL)
    {
        printf("[div_vectors] Bad args\n");
        return E_BADARGS;
    }

    if (src1->dim != src2->dim)
    {
        printf("[div_vectors] Dims don't equal\n");
        return E_ERROR;
    }

    int ret = init_vector(res, src1->dim);
    if (ret < 0)
    {
        printf("[div_vectors] Initializing res vector error\n");
        return E_ERROR;
    }

    for (int i = 0; i < src1->dim; i++)
        res->arr[i] = src1->arr[i] - src2->arr[i];

    return 0;
}

int solve_eq(sq_matrix* L, vector* f, vector* res)
{
    if (L == NULL || f == NULL || res == NULL)
    {
        printf("[solve_eq] Bad input args\n");
        return E_BADARGS;
    }

    if (L->metrics != f->dim)
    {
        printf("[solve_eq] Metrics of L %d != dim of f %d", L->metrics, f->dim);
        return E_ERROR;
    }

    int dim = L->metrics;

    vector temp_vector;
    int ret = init_vector(&temp_vector, dim);
    if (ret < 0)
    {
        printf("[solve_eq] init temp vector error\n");
        return E_ERROR;
    }

    for (int i = 0; i < dim; i++)
    {
        if (equal(L->data[i][i], 0.0))
        {
            printf("[solve_eq] %d diag elem == 0(with precision) in temple stage\n", i);
            vector_destr(&temp_vector);
            return E_ERROR;
        }

        temp_vector.arr[i] = f->arr[i];
        for (int col = 0; col < i; col++)
            temp_vector.arr[i] -= temp_vector.arr[col] * L->data[i][col];

        temp_vector.arr[i] /= L->data[i][i];
    }

    //for (int i = 0; i < dim; i++)
    //    printf("%lg\n", temp_vector.arr[i]);

    ret = init_vector(res, dim);
    if (ret < 0)
    {
        printf("[solve_eq] Init result vector error\n");
        return E_ERROR;
    }

    for (int i = dim - 1; i >= 0; i--)
    {
        if (equal(L->data[i][i], 0.0))
        {
            printf("[solve_eq] %d diag elem == 0(with precision) in last stage\n", i);
            vector_destr(&temp_vector);
            vector_destr(res);
            return E_ERROR;
        }

        res->arr[i] = temp_vector.arr[i];

        for (int row = i + 1; row < dim; row++)
            res->arr[i] -= temp_vector.arr[row] * L->data[row][i];

        res->arr[i] /= L->data[i][i];
    }

    //printf("\n");

    //for (int i = 0; i < dim; i++)
        //printf("%lg\n", res->arr[i]);

    vector_destr(&temp_vector);

    return 0;
}

int print_sq_matrix(sq_matrix* matr)
{
    if (matr == NULL)
    {
        printf("[printf_sq_matrix] Bad args\n");
        return E_BADARGS;
    }

    for (int i = 0; i < matr->metrics; i++)
    {
        for (int j = 0; j < matr->metrics; j++)
            printf("%lg ", matr->data[i][j]);
        printf("\n");
    }
}

int sq_matrix_trans(sq_matrix* matr, sq_matrix* res)
{
    if (matr == NULL || res == NULL)
    {
        printf("[sq_matrix_trans] Bad args");
        return E_BADARGS;
    }

    int dim = matr->metrics;
    res->metrics = matr->metrics;

    errno = 0; // for all
    double* arr = (double*) calloc(dim * dim, sizeof(double));
    if (arr == NULL)
    {
        perror("[calc_L] Alloc buff for matrix error\n");
        return E_ERROR;
    }

    res->data = (double**) calloc(dim, sizeof(double*));
    if (res->data == NULL)
    {
        perror("[calc_L] Alloc control pointers array for matrix\n");
        free(arr);
        return E_ERROR;
    }

    for (int row = 0; row < dim; row++)
        res->data[row] = arr + dim * row;

    for (int row = 0; row < dim; row++)
    {
        for (int col = 0; col < dim; col++)
            res->data[col][row] = matr->data[row][col];
    }

    return 0;
}

int equal(double first, double second)
{
    if (first > second)
    {
        if (first - second >= PRES)
            return 0;
    }
    else
    {
        if (second - first >= PRES)
            return 0;
    }

    return 1;
}

int calc_L(sq_matrix* matr, sq_matrix* res)
{
    if (matr == NULL || res == NULL)
    {
        printf("[calc_L] Bad args\n");
        return E_BADARGS;
    }

    int dim = matr->metrics;
    res->metrics = matr->metrics;

    errno = 0; // for all
    double* arr = (double*) calloc(dim * dim, sizeof(double));
    if (arr == NULL)
    {
        perror("[calc_L] Alloc buff for matrix error\n");
        return E_ERROR;
    }

    res->data = (double**) calloc(dim, sizeof(double*));
    if (res->data == NULL)
    {
        perror("[calc_L] Alloc control pointers array for matrix\n");
        free(arr);
        return E_ERROR;
    }

    for (int row = 0; row < dim; row++)
        res->data[row] = arr + dim * row;

    for (int k = 0; k < dim; k++)
    {
        double diag = matr->data[k][k];
        for (int i = 0; i < k; i++)
            diag -= res->data[k][i] * res->data[k][i];

        double diag_temp = 0.0;
        int ret = my_sqrt(diag, &diag_temp);
        if (ret < 0)
        {
            printf("[calc_L] Bad sqrt of diag element\n");
            matrix_destr(res);
            return E_ERROR;
        }

        if (equal(diag_temp, 0.0)) // with pres?
        {
            printf("[calc_L] diag element is zero");
            matrix_destr(res);
            return E_ERROR;
        }

        res->data[k][k] = diag_temp;

        for (int i = k + 1; i < dim; i++)
        {
            double temp = matr->data[i][k];
            for (int col = 0; col < k; col++)
                temp -= res->data[i][col] * res->data[k][col];

            res->data[i][k] = temp / diag_temp;
        }
    }

    return 0;
}

static int not_equal(double first, double second)
{
    if (first > second)
    {
        if (first - second >= PRES)
            return 1;
    }
    else
    {
        if (second - first >= PRES)
            return 1;
    }

    return 0;
}

int my_sqrt(double num, double* res)
{
    if (res == NULL)
    {
        printf("[my_sqrt] Bad args\n");
        return E_BADARGS;
    }

    if (num < 0.0 || !(num == num))
    {
        printf("[my_sqrt] Bad number %lg\n", num);
        return E_ERROR;
    }

    double max, min;
    if (num > 1.0)
    {
        max = num;
        min = 1.0;
    }
    else
    {
        max = 1.0;
        min = num;
    }

    double temp_point = 0.0;
    double err = 0.0;
    do {
        temp_point = (max + min) / 2;
        err = temp_point * temp_point - num;

        //printf("err = %lg\n", err);

        if (err < 0.0)
            min = temp_point;
        else
            max = temp_point;

    } while (not_equal(err, 0.0));

    *res = temp_point;

    return 0;
}

int vector_destr(vector* vec)
{
    if (vec == NULL)
    {
        printf("[vector_destr] Bad args\n");
        return E_BADARGS;
    }

    free(vec->arr);

    vec->arr = NULL;
    vec->dim = 0;

    return 0;
}

int mul_sq_matrix_to_vector(sq_matrix* matr, vector* vec, vector* res)
{
    if (matr == NULL || vec == NULL || res == NULL)
    {
        printf("[mul_sq_matrix_to_vector] Bad args\n");
        return E_BADARGS;
    }

    if (matr->metrics != vec->dim)
    {
        printf("[mul_sq_matrix_to_vector] Bad format of parametrs matrix metr %d != vector dim %d\n", matr->metrics, vec->dim);
        return E_ERROR;
    }
    int dim = vec->dim; // more comfortable

    int ret = init_vector(res, dim);
    if (ret < 0)
    {
        printf("[mul_sq_matrix_to_vector] Initializing vector error\n");
        return E_ERROR;
    }

    for (int row = 0; row < dim; row++)
    {
        res->arr[row] = 0.0;

        for(int i = 0; i < dim; i++)
            res->arr[row] += matr->data[row][i] * vec->arr[i];
    }

    return 0;
}

int init_vector(vector* vec, int dim)
{
    if (vec == NULL)
    {
        printf("[init_vector] Bad args\n");
        return E_BADARGS;
    }

    errno = 0;
    double* arr = (double*) calloc(dim, sizeof(double*));
    if (arr == NULL)
    {
        perror("[init_vector] Allocation arr error\n");
        return E_ERROR;
    }

    vec->dim = dim;
    vec->arr = arr;

    return 0;
}

int matrix_destr(sq_matrix* matr)
{
    if (matr == NULL)
    {
        printf("[matrix_destr] Bad input args\n");
        return E_ERROR;
    }

    free(matr->data[0]);
    free(matr->data);

    return 0;
}

int get_matrix_from_file(sq_matrix* matr, const char* file_name)
{
    if (file_name == NULL || matr == NULL)
    {
        printf("[get_matrix_from_file] Bad input args\n");
        return E_BADARGS;
    }

    errno = 0; // one for all

    FILE* in = fopen(file_name, "r");
    if (in == NULL)
    {
        perror("[get_matrix_from_file] Opening input file error\n");
        return E_ERROR;
    }

    int dim = 0;
    int ret = fscanf(in, "%d\n", &dim);
    if (ret <= 0)
    {
        perror("[get_matrix_from_file] Scanning n from matrix error\n");
        return E_ERROR;
    }
    if (dim < 0)
    {
        printf("[get_matrix_from_file] Bad dimention of matrix\n");
        return E_ERROR;
    }

    matr->metrics = dim;

    double* temp_data = (double*) calloc(dim * dim, sizeof(double));
    if (temp_data == NULL)
    {
        perror("[get_matrix_from_file] Allocaion memomry for matrix buffer error\n");
        return E_ERROR;
    }

    matr->data = (double**) calloc(dim, sizeof(double*));
    if (matr->data == NULL)
    {
        perror("[get_matrix_from_file] Allocation for controlling matrix memory error\n");
        return E_ERROR;
    }

    for (int row = 0; row < dim; row++)
    {
        matr->data[row] = temp_data + dim * row;
        for (int col = 0; col < dim; col++)
        {
            double temp = 0.0;
            ret = fscanf(in, "%lg", &(matr->data[row][col]));
            if (ret <= 0)
            {
                printf("[get_matrix_from_file] Scan %d row and %d column element\n", row, col);
                return E_ERROR;
            }
        }
    }

    fclose(in);

    return 0;
}
