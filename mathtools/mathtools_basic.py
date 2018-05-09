def factorial(n):
    ''' Calculates factorial of an integer, n. '''
    tot = 1
    for i in range(n, 1, -1):
        tot *= i
    if tot == 1 and n != 1:
        raise ValueError("argument must be a natural number.")
    else:
        return tot


def tot_perm(n, k):
    ''' Calculates total permutations of k items chosen from n items. '''
    tot = factorial(n) // factorial(n - k)
    if n < 0 or k < 0:
        raise ValueError("arguments must be natural numbers.")
    else:
        return tot
