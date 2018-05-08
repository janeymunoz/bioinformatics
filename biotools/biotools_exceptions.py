class SeqError(Exception):
    ''' Base class for exceptions in this module. '''
    pass


class DNAError(Exception):
    ''' Exception raised for errors in the input sequence.
    Attributes:
        expression -- input expression in which the error occurred.
        message -- explanation of the error
    '''
    def __init__(self, expression, message):
        self.expression = expression
        self.message = 'the string must only have characters \
                       "A", "T", "G", "C"'


class RNAError(Exception):
    ''' Exception raised for errors in the input sequence.
    Attributes:
        expression -- input expression in which the error occurred.
        message -- explanation of the error
    '''
    def __init__(self, expression, message):
        self.expression = expression
        self.message = 'the string must only have characters \
                       "A", "U", "G", "C"'
