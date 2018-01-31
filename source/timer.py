#!/usr/bin/env python

import time                                                
import Messages

def timeit(method):
  
  """
  Measuring the run time of a routine.
  """

  def timed(*args, **kw):
    
    ts = time.time()
    result = method(*args, **kw)
    te = time.time()

    msg = 'Runtime: %2.2f sec' % (te-ts)
    
    messages.log(method.__name__, msg, verbose=0)
    
    return result

  return timed
