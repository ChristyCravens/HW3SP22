import copy, Gauss_Seidel as GS, solveMatrix as SM


def LUFactorization(A):
    """This function uses the Doolittle algorithm to perform LU factorization on a given matrix (A).
    The function then outputs two matrices, one for the upper and one for the lower.

    :param A: the argument for the original matrix given
    :return: the lower and upper triangles after performing LU Factorization
    """
    n = len(A)     # Number of loops for values in matrices
    lower = [[0 for x in range(n)] for y in range(n)]
    upper = [[0 for x in range(n)] for y in range(n)]

    # Break apart the matrix into upper and lower triangles
    for i in range(n):

        # Upper triangle for loop
        for k in range(i, n):

            # Sum of L(i, j) * U(j, k)
            sum = 0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])

            # Define what the upper triangle is
            upper[i][k] = A[i][k] - sum

        # Lower triangle for loop
        for k in range(i, n):
            if i == k:
                lower[i][i] = 1  # Setting the diagonal = 1
            else:
                # Sum of L(k, j) * U(j, i)
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])

                # Define what the lower triangle is
                lower[k][i] = int((A[k][i] - sum) /
                                  upper[i][i])
    Aaug = lower, upper
    return Aaug

# LUFactorization(A)


def Aug(Aaug):
    """Using the augmented matrix Aaug to separate it and store values for later use in the Doolittle function for ease.

    :param Aaug: augmented matrix A
    :return: A
    """
    # Copy the augmented matrix A as A, b, since it also includes values of b.
    A, b = copy.copy(Aaug)
    # Create a starting point for values of x
    x = 0
    # Separating A with for loop
    for c in A:
        n = b[x]
        c.append(n)
        x += 1
    return A

def Doolittle(Aaug):
    """ Takes the arguments from the LUFactorization function and uses Doolittle's method where [L][U][x]=[b] and
    [y]=[U][x], so [L][y]=[b]. This function will solve for [L][b]=[y] and [U][y]=[x]

    :param Aaug: the augmented matrix A made up of multiplying the matrices [L]*[U]
    :return: solution for [x], solution for [y]
    """
    # Bring over the values of L, U, A, b from previous functions
    L, U = LUFactorization(Aaug)
    A, b = GS.separateAugmented(Aaug)

    # Create argument one as L, b, then solve for the values of b
    arg1 = L, b
    print("Values of B: ", b)
    # Create the matrix (m1) of L,b and solve for values of y, print y
    m1 = Aug(arg1)
    Y = SM.solution(m1)
    print("Values of Y: ", Y)

    # Create argument two as U, y, then solve for the values of x and print
    arg2 = U, Y
    m2 = Aug(arg2)
    X = SM.solution(m2)
    print("Values of X: ", X)

    # Check the answer to verify that b equals the calculated value of b from the program above. If it is,
    # it will tell the user that the solution is valid. If it is not, it tells the user the solution is invalid.
    # If the solution is invalid, the user will want to verify their input values within their matrix.
    checkAns = GS.checkMatrixSoln(Aaug, X, augmented=True)
    if checkAns == b:
        print("Valid solution.")
    else:
        print("Invalid solution. Please verify input.")

# Doolittle(Aaug)


def main():
    """This function calls the Doolittle function to solve for the augmented matrix outlined within this function.

    :return: solutions from the Doolittle function
    """

    # Enter the values for the augmented matrix A to be solved, having the last values as b
    A = [[3, 9, 6, 4.6],
         [18, 48, 39, 27.2],
         [9, -27, 42, 9.0]]

    Doolittle(A)

main()