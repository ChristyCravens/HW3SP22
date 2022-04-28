def solution(A):
    """This function finds the solution of a given augmented matrix and returns the values for x.

    :param A: the augmented matrix A, including values for b
    :return: values for x after solving Ax=b
    """
    # Define parameters for loops to work successfully, based on the number of values in A
    m = len(A)
    n = m + 1
    # Create a storage for x
    x = []
    # Matrix rows are not uniform, must resolve
    assert all([len(row) == m+1 for row in A[1:]])

    # For loop to find the pivot points
    for k in range(m):
        pivot = [abs(A[i][k]) for i in range(k, m)]
        i_max = pivot.index(max(pivot)) + k
        # Singular matrix
        assert A[i_max][k] != 0
        # Swapping out rows
        A[k], A[i_max] = A[i_max], A[k]
        for i in range(k+1, m):
            f=A[i][k]/A[k][k]
            for j in range(k+1, n):
                A[i][j] -= A[k][j] * f
            # Lower triangle matrix first has to be zeros
            A[i][k] = 0

    # Use for loops to solve Ax=b
    for i in range(m-1, -1, -1):
        x.insert(0, round(A[i][m], 4) / round(A[i][i], 4))
        for k in range(i-1, -1, -1):
            A[k][m] -= round(A[k][i], 4)*round(x[0], 4)

    # Round the values of x into X
    X = []
    for num in x:
        X.append(round(num, 4))
    return X