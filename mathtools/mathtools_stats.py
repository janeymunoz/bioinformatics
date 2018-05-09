def sd(data, type):
    ''' Takes data set and returns its standard deviation.
    Either calculates the population standard deviation, or the
    sample standard deviation, depending on the type given in the
    argument. "pop" returns the population standard deviation,
    while "sample" returns the sample standard deviation.
    '''
    import statistics as stat
    # Calculate the mean of the data set.
    mean = stat.mean(data)
    # Find the square of each data point's distance to the mean.
    sd_sum = 0
    num_points = 0
    for point in data:
        sd_sum += ((point - mean) ** 2)
        num_points += 1
    # Divide sd_sum by num_points and take square-root.
    if type == "pop":
        sd = (sd_sum / num_points) ** 0.5
    elif type == "sample":
        sd = stat.stdev(data)
    else:
        raise ValueError('argument of either "pop" or "sample" must be \
                entered after data set to indicate standard deviation type.')
    return sd


def var(data):
    '''Calculates sample variance of a data set.'''
    import statistics as stat
    num_points = len(data)
    mean = stat.mean(data)
    s2 = 0
    for point in data:
        s2 += (point - mean) ** 2
    s2 = s2 / (num_points - 1)
    return s2


def var_covar(data):
    import statistics as stat
    var = stat.pvariance(data)
    return var
print(var_covar((17, 15, 23, 7, 9, 13)))


def sw_test(data):
    '''Performs the Shapiro-Wilk test.
    Tests the null hypothesis that a sample came from a normally distributed
    population. Recommended to be used on small populations (< 20). Data should
    be a list or tuple.
    '''
    import statistics as stat
    # Calculate the mean of the data set
    x_bar = stat.mean(data)
    # Sort data set small to large
    xi = sorted(data)
    # Get coefficients, ai.
    ai = []
