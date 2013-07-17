#!/usr/bin/python
from types import GeneratorType

#Returns a combination of keys/values for all sublists
#e.g: [('a', 1), ('bc', [2, 3]), ('def', [4, 5, 6])] ->
#[[('a', 1), ('bc', 2), ('def', 4)],
#[('a', 1), ('bc', 2), ('def', 5)],
#[('a', 1), ('bc', 2), ('def', 6)],
#[('a', 1), ('bc', 3), ('def', 4)],
#[('a', 1), ('bc', 3), ('def', 5)],
#[('a', 1), ('bc', 3), ('def', 6)]]
def combinations(items, current=[]):
    if not items:
        return current
    (key, vals) = items[0]
    if isinstance(vals, list):
        return flatten((combinations(items[1:], current + [(key, val)]) for val in vals))
    return flatten(combinations(items[1:], current + [(key, vals)]))

#Flattens a recursive set of generator functions to one generator function with all enclosed elements in a line
def flatten(gen):
    for elem in gen:
        if isinstance(elem, GeneratorType):
            for i in flatten(elem):
                yield i
        else:
            yield elem
