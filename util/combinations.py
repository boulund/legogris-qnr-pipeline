#!/usr/bin/python
from types import GeneratorType

def combinations(items):
    """
    Return a combination of keys/values for all values for each key in a dictionary.

    >>> combinations({'a': 1, 'bc': [2, 3], 'def': [4, 5, 6]})
    [{'a': 1, 'bc': 2, 'def': 4},
    {'a': 1, 'bc': 2, 'def': 5},
    {'a': 1, 'bc': 2, 'def': 6},
    {'a': 1, 'bc': 3, 'def': 4},
    {'a': 1, 'bc': 3, 'def': 5},
    {'a': 1, 'bc': 3, 'def': 6}]
    """
    return [dict((x[0], x[1]) for x in xs) for xs in _combinations(items.items())]
    # return [{x[0]: x[1] for x in xs} for xs in _combinations(items.items())]

def _combinations(items, current=[]):
    if not items:
        return current
    (key, vals) = items[0]
    if isinstance(vals, list):
        return flatten((_combinations(items[1:], current + [(key, val)]) for val in vals))
    return flatten(_combinations(items[1:], current + [(key, vals)]))

def flatten(gen):
    """
    Flatten a recursive set of generator functions to one generator function yielding all nested items.
    """
    for elem in gen:
        if isinstance(elem, GeneratorType):
            for i in flatten(elem):
                yield i
        else:
            yield elem
