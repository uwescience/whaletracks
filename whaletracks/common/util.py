"""Utilities for Whale Project."""

def complexifyString(stg, separator=";"):
  """
  Converts a string representation of complex numbers into
  a list of complex numbers.
  """
  splits = stg.split(separator)
  if len(splits[0]) > 0:
    return [complex(eval(s)[0], eval(s)[1]) for s in splits]
  else:
    return []
