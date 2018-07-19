import time                                                
from contextlib import contextmanager

@contextmanager
def time_block(name):
    """
    This generator measure the time taken by a step, and prints the result
    in second to the console.
    
    Arguments:
        name (str): prefix to print
    """
    ts = time.time()
    yield
    te = time.time()
 
    print('%s %2.2f sec' %  (name, te-ts))

def timeit(method):
    """This decorator add a timing request for each call to the decorated function.
    
    Arguments:
        method (function): the method to decorate
    """

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts))
        return result

    return timed

def timeit_quiet(method):
    """This decorator add a timing request for each call to the decorated function.
    Same as cosmotool.timeit_ but is quieter by not printing the values of the arguments.
    
    Arguments:
        method (function): the method to decorate
    """

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('%r %2.2f sec' % (method.__name__, te-ts))
        return result

    return timed

