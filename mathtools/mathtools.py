def sd(data_set, type):
    ''' Takes data set and returns its standard deviation.
    Either calculates the population standard deviation, or the
    sample standard deviation, depending on the type given in the
    argument. "pop" returns the population standard deviation,
    while "sample" returns the sample standard deviation.
    '''
    # Calculate the mean of the data set.
    mean_sum = 0
    num_points = len(data_set)
    for point in data_set:
        mean_sum += point
    mean = mean_sum / num_points
    # Find the square of each data point's distance to the mean.
    sd_sum = 0
    for point in data_set:
        sd_sum += ((point - mean) ** 2)
    # Divide sd_sum by num_points and take square-root.
    if type == "pop":
        sd = (sd_sum / num_points) ** 0.5
    elif type == "sample":
        sd = (sd_sum / (num_points - 1)) ** 0.5
    else:
        raise ValueError('argument of either "pop" or "sample" must be \
                entered after data set to indicate standard deviation type.')
    return sd


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
