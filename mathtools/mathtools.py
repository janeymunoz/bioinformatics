def factorial(n):
    ''' Calculates factorial of an integer, n. '''
    if n < 1 or not isinstance(n, int):
        return "argument must be a natural number"
    else:
        total = 1
        for i in range(n, 1, -1):
            total *= i
    return total


def tot_perm(n, k):
    ''' Calculates total permutations of k items chosen from n items. '''
    if not isinstance(n, int) or not isinstance(k, int):
        return "arguments must be natural numbers"
    else:
        tot = factorial(n) // factorial(n - k)
    return tot
