import sys
import math
import getopt
import itertools
import MF as mf

'''
Permutation class
'''

class Permutation:
  @classmethod
  def inverse(cls, pi):
    inv = [0 for i in range(len(pi))]
    for i in range(len(pi)):
      inv[pi[i]] = i
    return inv


'''
Instance class
'''

class Instance:
  def __init__(self, _pi, _bpi, _biota):
    # Data validation
    if len(_bpi) != len(_biota):
      raise ValueError('The total of intergenic regions in source and target genomes is different.')
    if sum(_bpi) != sum(_biota):
      raise ValueError('The total of nucleotides in source and target genomes is different.')
    if min(_bpi) < 0 or min(_biota) < 0:
      raise ValueError('There is a negative intergenic region size.')
    if len(_pi) + 1 != len(_bpi):
      raise ValueError('There is an incompatible number of genes and intergenic regions.')
    if sorted(_pi) != list(range(1, len(_pi) + 1)):
      raise ValueError('The order of the genes is not an unsigned permutation.')

    self.size = len(pi)
    # Extended permutation
    self.pi = [0] + _pi + [len(_pi) + 1]
    # Intergenic regions
    # self.bpi[0] and self.biota[0] are dummy intergenic regions to follow the paper definition
    self.bpi = [0] + _bpi
    self.biota = [0] + _biota

  def print(self):
    print(self.pi)
    print(self.bpi)
    print(self.biota)

'''
Util class
'''

class Util:
  @classmethod
  def reversal(cls, I, params):
    i, j, x, y = params
    # Stop case
    if i > j:
      return
    # Intergenic parameters constraints
    if (x < 0) or (x > I.bpi[i]):
      raise Exception('Intergenic Reversal: Invalid Parameter x.')
    if (y < 0) or (y > I.bpi[(j + 1)]):
      raise Exception('Intergenic Reversal: Invalid Parameter y.')
    # Genes
    I.pi[i], I.pi[j] =  I.pi[j], I.pi[i]
    # Intergenic regions
    x_prime = (I.bpi[i] - x)
    y_prime = (I.bpi[j + 1] - y)
    I.bpi[j + 1] = x_prime + y_prime
    I.bpi[i] = x + y
    # Recursion
    Util.reversal(I, (i + 1, j - 1, 0, I.bpi[j]))

  @classmethod
  def transposition(cls, I, params):
    i, j, k, x, y, z = params
    # Intergenic parameters constraints
    if (x < 0) or (x > I.bpi[i]):
      raise Exception('Intergenic Transposition: Invalid Parameter x.')
    if (y < 0) or (y > I.bpi[j]):
      raise Exception('Intergenic Transposition: Invalid Parameter y.')
    if (z < 0) or (z > I.bpi[k]):
      raise Exception('Intergenic Transposition: Invalid Parameter z.')
    # x' and y'
    x_prime = I.bpi[i] - x
    y_prime = I.bpi[j] - y
    # Transposition simulated using three reversals
    Util.reversal(I, (i, k - 1, x, z))
    Util.reversal(I, (i, i + k - j - 1, x, y_prime))
    Util.reversal(I, (i + k - j, k - 1, z, x_prime))

  @classmethod
  def is_block(cls, I, a, b):
    return (abs(a - b)) == 1 and (abs(I.pi[a] - I.pi[b]) == 1) and (I.bpi[max(a, b)] == I.biota[max(I.pi[a], I.pi[b])])
  
  @classmethod
  def is_breakpoint(cls, I, i):
    return not Util.is_block(I, i, (i + 1))

  @classmethod
  def breakpoints(cls, I):
    breakpoints = []
    overcharged = []
    undercharged = []
    soft = []
    for i in range(len(I.pi) - 1):
      if Util.is_breakpoint(I, i):
        breakpoints.append(i)
        if abs(I.pi[i + 1] - I.pi[i]) != 1:
          soft.append(i)
        elif I.bpi[max((i + 1), i)] > I.biota[max(I.pi[i + 1], I.pi[i])]:
          overcharged.append(i)
        else:
          undercharged.append(i)
    return (breakpoints, soft, overcharged, undercharged)


'''
Algorithm class
'''

class Algorithm:

  def __init__(self, _I, _greedy):
    self.I = _I
    self.greedy = _greedy
    self.breakpoints = []
    self.overcharged = []
    self.undercharged = []
    self.soft = []
    self.sequence = []
    self.update()
    self.lower_bound = math.ceil(len(self.breakpoints) / 3.0)

  def update(self):
    self.breakpoints, self.soft, self.overcharged, self.undercharged = Util.breakpoints(self.I)

  # Lemma 3.4 - There are two or more overcharged breakpoints in I
  def case_i(self):
    if len(self.overcharged) < 2:
      return False
    # If there are more than two overcharged breakpoints, the third breakpoint will receive the excess
    if len(self.overcharged) > 2:
      i, j, k = self.overcharged[:3]
      a, b, c = self.I.bpi[i + 1], self.I.bpi[j + 1], self.I.bpi[k + 1]
      x, y = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
      z = (a + b + c) - (x + y)
      t1, t2 = mf.redistribution(a, b, c, x, y, z)
      params = ((i + 1), (j + 1), (k + 1)) + t1
      self.sequence.append(params)
      Util.transposition(self.I, params)
      params = ((i + 1), ((i + 1) + (k + 1) - (j + 1)), (k + 1)) + t2
      self.sequence.append(params)
      Util.transposition(self.I, params)
      return True
    # Otherwise, a soft or an undercharged breakpoint will receive the excess
    else:
      i, j, k = self.overcharged + [None]
      a = b = c = 0
      x = y = z = 0

      # The soft breakpoint is prioritized to receive the excess
      if len(self.soft):
        k = self.soft[0]
      # The undercharged breakpoint is the last option
      else:
        k = self.undercharged[0]

      if k < i:
        i, j, k = k, i, j
        a, b, c = self.I.bpi[i + 1], self.I.bpi[j + 1], self.I.bpi[k + 1]
        y, z = self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])], self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])]
        x = (a + b + c) - (y + z)
      elif i < k and k < j:
        i, j, k = i, k, j
        a, b, c = self.I.bpi[i + 1], self.I.bpi[j + 1], self.I.bpi[k + 1]
        x, z = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])]
        y = (a + b + c) - (x + z)
      else:
        a, b, c = self.I.bpi[i + 1], self.I.bpi[j + 1], self.I.bpi[k + 1]
        x, y = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
        z = (a + b + c) - (x + y)

      t1, t2 = mf.redistribution(a, b, c, x, y, z)
      params = ((i + 1), (j + 1), (k + 1)) + t1
      self.sequence.append(params)
      Util.transposition(self.I, params)
      params = ((i + 1), ((i + 1) + (k + 1) - (j + 1)), (k + 1)) + t2
      self.sequence.append(params)
      Util.transposition(self.I, params)
      return True

  # Lemma 3.5 - Case 1 (i, j) or (i + 1, j + 1)
  def soft_case_i(self, i, j):
    i, j = sorted([i, j])
    if abs(self.I.pi[i] - self.I.pi[j]) == 1:
      goal = self.I.biota[max(self.I.pi[i], self.I.pi[j])]
      if (self.I.bpi[(i + 1)] + self.I.bpi[(j + 1)]) >= goal:
        if self.I.bpi[(i + 1)] >= goal:
          x = goal
          y = 0
        else:
          x = self.I.bpi[(i + 1)]
          y = goal - x
        Util.reversal(self.I, ((i + 1), j, x, y))
        self.sequence.append(((i + 1), j, x, y))
        return True
    if abs(self.I.pi[(i + 1)] - self.I.pi[(j + 1)]) == 1:
      goal = self.I.biota[max(self.I.pi[(i + 1)], self.I.pi[(j + 1)])]
      if (self.I.bpi[(i + 1)] + self.I.bpi[(j + 1)]) >= goal:
        if self.I.bpi[(j + 1)] >= goal:
          x = self.I.bpi[(i + 1)]
          y = self.I.bpi[(j + 1)] - goal
        else:
          x = self.I.bpi[(i + 1)] - (goal - self.I.bpi[(j + 1)])
          y = 0
        Util.reversal(self.I, ((i + 1), j, x, y))
        self.sequence.append(((i + 1), j, x, y))
        return True
    return False

  # Lemma 3.5 - Case 2 (i + 1, j)
  def soft_case_ii(self, i, j):
    i, j = sorted([i, j])
    if Util.is_block(I, (i + 1), j):
      return False
    if abs(self.I.pi[(i + 1)] - self.I.pi[j]) == 1:
      goal = self.I.biota[max(self.I.pi[(i + 1)], self.I.pi[j])]
      if (self.I.bpi[(i + 1)] + self.I.bpi[(j + 1)]) >= goal:
        i, j, k = i, self.breakpoints[self.breakpoints.index(i) + 1], j
        x = y = z = 0
        if self.I.bpi[(i + 1)] >= goal:
          x = self.I.bpi[(i + 1)] - goal
          y = 0
          z = 0
        elif self.I.bpi[(k + 1)] >= goal:
          x = self.I.bpi[(i + 1)]
          y = 0
          z = goal
        else:
          x = 0
          y = 0
          z = goal - self.I.bpi[(i + 1)]
        self.sequence.append(((i + 1), (j + 1), (k + 1), x, y, z))
        Util.transposition(self.I, ((i + 1), (j + 1), (k + 1), x, y, z))
        return True
    return False

  # Lemma 3.5 - Case 3 (i, j + 1)
  def soft_case_iii(self, i, j):
    i, j = sorted([i, j])
    if abs(self.I.pi[i] - self.I.pi[(j + 1)]) == 1:
      goal = self.I.biota[max(self.I.pi[i], self.I.pi[(j + 1)])]
      if (self.I.bpi[(i + 1)] + self.I.bpi[(j + 1)]) >= goal:
        k = None
        x = y = z = 0
        if self.soft.index(i):
          i, j, k = self.soft[0], i, j
          if self.I.bpi[j + 1] >= goal:
            x = 0
            y = goal
            z = self.I.bpi[k + 1]
          elif self.I.bpi[k + 1] >= goal:
            x = 0
            y = 0
            z = self.I.bpi[k + 1] - goal
          else:
            x = 0
            y = goal - self.I.bpi[k + 1]
            z = 0
        else:
          i, j, k = i, j, self.soft[self.soft.index(j) + 1]
          if self.I.bpi[i + 1] >= goal:
            x = goal
            y = self.I.bpi[j + 1]
            z = 0
          elif self.I.bpi[j + 1] >= goal:
            x = 0
            y = self.I.bpi[j + 1] - goal
            z = 0
          else:
            x = goal - self.I.bpi[j + 1]
            y = 0
            z = 0
        self.sequence.append(((i + 1), (j + 1), (k + 1), x, y, z))
        Util.transposition(self.I, ((i + 1), (j + 1), (k + 1), x, y, z))
        return True
    return False

  # Lemma 3.5 - There is a pair of soft connected breakpoints in I
  def case_ii(self):
    # Inverse permutation indicating the position of each element i in pi
    inverse = Permutation.inverse(self.I.pi)

    # A Logical list indicating if exist a soft breakpoint in (i, (i + 1))
    is_soft = [False for _ in range(self.I.size + 1)]
    for i in self.soft:
      is_soft[i] = True

    for i in self.soft:
      # For each soft breakpoint, it is necessary to check the adjacencies of at most four elements: (pi[i] - 1), (pi[i] + 1), (pi[(i + 1)] - 1), and (pi[(i + 1)] + 1)
      elements = [e for e in [(self.I.pi[i] - 1), (self.I.pi[i] + 1), (self.I.pi[(i + 1)] - 1), (self.I.pi[(i + 1)] + 1)] if e >= 0 and e <= (self.I.size + 1)]
      for e in elements:
        # Check if (j, (j + 1)) is a soft breakpoint and it is connected with (i, (i + 1))
        j = inverse[e]
        if j <= self.I.size and is_soft[j] and (self.soft_case_i(i, j) or self.soft_case_ii(i, j) or self.soft_case_iii(i, j)):
          return True
        # Check if ((j - 1), j) is a soft breakpoint and it is connected with (i, (i + 1))
        j -= 1
        if j >= 0 and is_soft[j] and (self.soft_case_i(i, j) or self.soft_case_ii(i, j) or self.soft_case_iii(i, j)):
          return True
    return False

  # Lemma 3.8 - There is one overcharged breakpoint and at least one undercharged breakpoint in I
  def case_iii(self):
    if len(self.overcharged) != 1 or len(self.undercharged) == 0:
      return False
    # If the instance has only two breakpoints, a sequence of two reversals is applied
    if len(self.undercharged) + len(self.soft) == 1:
      i, j = sorted(self.overcharged + self.undercharged)
      if abs(i - j) != 1:
        Util.reversal(self.I, ((i + 1), j, 0, 0))
        self.sequence.append(((i + 1), j, 0, 0))
        return self.soft_case_i(i, j)
      else:
        if i in self.overcharged:
          goal = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])]
          Util.reversal(self.I, ((i + 1), j, goal, 0))
          self.sequence.append(((i + 1), j, goal, 0))
        else:
          goal = self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
          Util.reversal(self.I, ((i + 1), j, self.I.bpi[(i + 1)], (self.I.bpi[(j + 1)] - goal)))
          self.sequence.append(((i + 1), j, self.I.bpi[(i + 1)], (self.I.bpi[(j + 1)] - goal)))
        return True
    # Otherwise, a sequence of two transposition is applied
    else:
      i = j = k = 0
      a = b = c = 0
      x = y = z = 0
      # The soft breakpoint is prioritized as the third breakpoint
      if len(self.soft) > 0:
        i, j, k = sorted([self.overcharged[0], self.undercharged[0], self.soft[0]])
        a, b, c = self.I.bpi[(i + 1)], self.I.bpi[(j + 1)], self.I.bpi[(k + 1)]
        if i in self.overcharged:
          if j in self.undercharged:
            x, y = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
            z = (a + b + c) - (x + y)
          else:
            x, z = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])]
            y = (a + b + c) - (x + z)
        elif j in self.overcharged:
          if i in self.undercharged:
            y, x = self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])], self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])]
            z = (a + b + c) - (x + y)
          else:
            y, z = self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])], self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])]
            x = (a + b + c) - (y + z)
        else:
          if i in self.undercharged:
            z, x = self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])], self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])]
            y = (a + b + c) - (x + z)
          else:
            z, y = self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
            x = (a + b + c) - (y + z)
      # The undercharged breakpoint is the last option as the third breakpoint
      else:
        i, j, k = sorted([self.overcharged[0], self.undercharged[0], self.undercharged[1]])
        a, b, c = self.I.bpi[(i + 1)], self.I.bpi[(j + 1)], self.I.bpi[(k + 1)]
        if i in self.overcharged:
          x, y = self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
          z = (a + b + c) - (x + y)
        elif j in self.overcharged:
          y, x = self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])], self.I.biota[max(self.I.pi[i], self.I.pi[(i + 1)])]
          z = (a + b + c) - (x + y)
        else:
          z, y = self.I.biota[max(self.I.pi[k], self.I.pi[(k + 1)])], self.I.biota[max(self.I.pi[j], self.I.pi[(j + 1)])]
          x = (a + b + c) - (y + z)
      t1, t2 = mf.redistribution(a, b, c, x, y, z)
      params = ((i + 1), (j + 1), (k + 1)) + t1
      self.sequence.append(params)
      Util.transposition(self.I, params)
      params = ((i + 1), ((i + 1) + (k + 1) - (j + 1)), (k + 1)) + t2
      self.sequence.append(params)
      Util.transposition(self.I, params)
      return True

  # If this method applies a reversal, then the instance remains with no soft connected breakpoints.
  def classic_reversal(self, i, j):
    i, j = sorted([i, j])
    # Moving all the nucleotides from bpi[(j + 1) to bpi[(i + 1)]
    if (abs(self.I.pi[i] - self.I.pi[j]) == 1):
      Util.reversal(self.I, (i + 1, j, self.I.bpi[(i + 1)], self.I.bpi[(j + 1)]))
      self.sequence.append((i + 1, j, self.I.bpi[(i + 1)], self.I.bpi[(j + 1)]))
      return True
    # Moving all the nucleotides from bpi[(i + 1) to bpi[(j + 1)]
    if (abs(self.I.pi[(i + 1)] - self.I.pi[(j + 1)]) == 1):
      Util.reversal(self.I, (i + 1, j, 0, 0))
      self.sequence.append((i + 1, j, 0, 0))
      return True
    return False

  # If this method applies a transposition, then a super hard breakpoint is created.
  def classic_transposition(self, i, j):
    i, j = sorted([i, j])
    # ((i + 1), j) are adjacent in iota and not in pi
    if (abs((i + 1) - j) != 1) and (abs(self.I.pi[(i + 1)] - self.I.pi[j]) == 1):
      i, j, k = i, self.soft[self.soft.index(i) + 1], j
      self.sequence.append(((i + 1), (j + 1), (k + 1), 0, 0, 0))
      Util.transposition(self.I, ((i + 1), (j + 1), (k + 1), 0, 0, 0))
      return True
    # (i, (j + 1)) are adjacent in iota and not in pi
    if abs(self.I.pi[i] - self.I.pi[(j + 1)]) == 1:
      k = 0
      if self.soft.index(i):
        i, j, k = self.soft[self.soft.index(i) - 1], i, j
      else:
        i, j, k = i, j, self.soft[self.soft.index(j) + 1]
      self.sequence.append(((i + 1), (j + 1), (k + 1), 0, 0, 0))
      Util.transposition(self.I, ((i + 1), (j + 1), (k + 1), 0, 0, 0))
      return True
    return False

  def create_undercharged(self):
    # Inverse permutation indicating the position of each element i in pi
    inverse = Permutation.inverse(self.I.pi)

    # A Logical list indicating if exist a soft breakpoint in (i, (i + 1))
    is_soft = [False for _ in range(self.I.size + 1)]
    for i in self.soft:
      is_soft[i] = True

    # If exists a decreasing strip, then a reversal is applied.
    for i in self.soft:
      # For each soft breakpoint, it is necessary to check the adjacencies of at most four elements: (pi[i] - 1), (pi[i] + 1), (pi[(i + 1)] - 1), and (pi[(i + 1)] + 1)
      elements = [e for e in [(self.I.pi[i] - 1), (self.I.pi[i] + 1), (self.I.pi[(i + 1)] - 1), (self.I.pi[(i + 1)] + 1)] if e >= 0 and e <= (self.I.size + 1)]
      for e in elements:
        # Check if (j, (j + 1)) is a soft breakpoint
        j = inverse[e]
        if j <= self.I.size and is_soft[j] and self.classic_reversal(i, j):
          return True
        # Check if ((j - 1), j) is a soft breakpoint
        j -= 1
        if j >= 0 and is_soft[j] and self.classic_reversal(i, j):
          return True

    # Otherwise, a transposition is applied.
    for i in self.soft:
      # For each soft breakpoint, it is necessary to check the adjacencies of at most four elements: (pi[i] - 1), (pi[i] + 1), (pi[(i + 1)] - 1), and (pi[(i + 1)] + 1)
      elements = [e for e in [(self.I.pi[i] - 1), (self.I.pi[i] + 1), (self.I.pi[(i + 1)] - 1), (self.I.pi[(i + 1)] + 1)] if e >= 0 and e <= (self.I.size + 1)]
      for e in elements:
        # Check if (j, (j + 1)) is a soft breakpoint
        j = inverse[e]
        if j <= self.I.size and is_soft[j] and self.classic_transposition(i, j):
          return True
        # Check if ((j - 1), j) is a soft breakpoint
        j -= 1
        if j >= 0 and is_soft[j] and self.classic_transposition(i, j):
          return True
    return False

  # Lemma 3.10 - There is one overcharged breakpoint and no undercharged breakpoint in I
  def case_iv(self):
    if not self.create_undercharged():
      return False
    self.update()
    if self.case_ii():
      self.update()
      return self.case_i() or self.case_ii() or self.case_iii()
    else:
      return self.case_iii()


###################################
####  GREEDY STRATEGY METHODS  ####
###################################


  def transp_remove_three_bp(self, i, j, k): #removes left, middle, and right breakpoints
    if (1 <= i < j < k <= self.I.size+1) : # at this point we already know that pairs of elements in positions (k-1, i) and (j-1, k) from self.I.pi are consecutive
        if (self.I.pi[i-1] != self.I.pi[i]+1 and self.I.pi[i-1] != self.I.pi[i]-1 and
            self.I.pi[j-1] != self.I.pi[j]+1 and self.I.pi[j-1] != self.I.pi[j]-1 and
            self.I.pi[k-1] != self.I.pi[k]+1 and self.I.pi[k-1] != self.I.pi[k]-1): #positions where the reversal takes place must be breakpoints
          if (self.I.pi[j] + 1 == self.I.pi[i-1] or self.I.pi[j] - 1 == self.I.pi[i-1]) : #check if the elements in positions (i-1) and (j) from self.I.pi are consecutive
            lp, mp, rp = self.I.bpi[i], self.I.bpi[j], self.I.bpi[k] 
            li, mi, ri = self.I.biota[max(self.I.pi[(i-1)], self.I.pi[j])], self.I.biota[max(self.I.pi[(k-1)], self.I.pi[i])], self.I.biota[max(self.I.pi[(j-1)], self.I.pi[k])]
            if (lp + mp + rp == li + mi + ri) and (mp <= li + ri) and (mp + lp >= li) and (mp + rp >= ri) and (lp + rp >= mi): # necessary conditions to remove three breakpoints using a distribution between intergenic regions
              m_to_r = max(0,(ri-rp)) # minimum value for y
              max_m_to_r = mp - max(0,li-lp) # maximum value for y, if something from middle must also goes to the left
              while m_to_r <= max_m_to_r :
                m_to_l = mp - m_to_r #value of y'
                if (m_to_l <= li) and (m_to_r <= ri) and (m_to_l + lp >= li) and (m_to_r + rp >= ri):
                    l_to_l = li - m_to_l # value of x
                    l_to_m = lp - l_to_l # value of x'
                    r_to_r = ri - m_to_r # value of z'
                    r_to_m = rp - r_to_r # value of z
                    if (l_to_m + r_to_m == mi) :
                      self.sequence.append((i, j, k, l_to_l, m_to_r, r_to_m))
                      Util.transposition(self.I, (i, j, k, l_to_l, m_to_r, r_to_m))
                      return True
                m_to_r += 1
    return False

  def transp_remove_two_bp(self, i, j, k, tipo): #tipo 12 if left and middle breakpoints are removed, 13 if left and right breakpoints are removed or 23 if middle and right breakpoints are removed
    if (1 <= i < j < k <= self.I.size+1) : # check if indices are valid for a transposition
        if (self.I.pi[i-1] != self.I.pi[i]+1 and self.I.pi[i-1] != self.I.pi[i]-1 and
            self.I.pi[j-1] != self.I.pi[j]+1 and self.I.pi[j-1] != self.I.pi[j]-1 and
            self.I.pi[k-1] != self.I.pi[k]+1 and self.I.pi[k-1] != self.I.pi[k]-1): #positions where the transposition takes place must be breakpoints
          lp, mp, rp = self.I.bpi[i], self.I.bpi[j], self.I.bpi[k] 
          li, mi, ri = self.I.biota[max(self.I.pi[(i-1)], self.I.pi[j])], self.I.biota[max(self.I.pi[(k-1)], self.I.pi[i])], self.I.biota[max(self.I.pi[(j-1)], self.I.pi[k])]
          if (tipo == 12) and (lp + mp >= li) and (lp + rp >= mi) and (lp <= li + mi): # necessary conditions to remove left and middle breakpoints using a distribution between intergenic regions
            m_to_l = max(0,li-lp) # minimum value for y'
            while m_to_l <= mp :
              m_to_r = mp - m_to_l #value of y
              if (li >= m_to_l) and (lp >= li - m_to_l):
                  l_to_l = li - m_to_l # value of x
                  l_to_m = lp - l_to_l # value of x'
                  if (mi >= l_to_m) and (rp >= mi - l_to_m) :
                    r_to_m = mi - l_to_m # value of z
                    r_to_r = rp - r_to_m # value of z'
                    
                    self.sequence.append((i, j, k, l_to_l, m_to_r, r_to_m))
                    Util.transposition(self.I, (i, j, k, l_to_l, m_to_r, r_to_m))
                    return True
              m_to_l += 1
          elif (tipo == 13) and (mp + rp >= ri) and (lp + mp >= li) and (rp <= mi + ri) and (lp <= li + mi): # necessary conditions to remove left and right breakpoints using a distribution between intergenic regions
            m_to_r = max(0,ri-rp) # minimum value for y
            max_m_to_r = mp - max(0,li-lp) # maximum value for y, if something from middle must also goes to the left
            while m_to_r <= max_m_to_r :
              m_to_l = mp - m_to_r #value of y'
              if (ri >= m_to_r) and (li >= m_to_l):
                  r_to_r = ri - m_to_r # value of z'
                  r_to_m = rp - r_to_r # value of z
                  if (mi >= r_to_m) and (lp >= mi - r_to_m) :
                    l_to_m = mi - r_to_m # value of x'
                    l_to_l = lp - l_to_m # value of x 
                    
                    self.sequence.append((i, j, k, l_to_l, m_to_r, r_to_m))
                    Util.transposition(self.I, (i, j, k, l_to_l, m_to_r, r_to_m))
                    return True
              m_to_r += 1
          elif (tipo == 23) and (mp + rp >= ri) and (lp + rp >= mi) and (rp <= mi + ri): # necessary conditions to remove middle and right breakpoints using a distribution between intergenic regions
            m_to_r = max(0,ri-rp) # minimum value for y
            while m_to_r <= mp :
              m_to_l = mp - m_to_r #value of y'
              if (ri >= m_to_r) and (rp >= ri - m_to_r):
                  r_to_r = ri - m_to_r # value of z'
                  r_to_m = rp - r_to_r # value of z
                  if (mi >= r_to_m) and (lp >= mi - r_to_m) :
                    l_to_m = mi - r_to_m # value of x'
                    l_to_l = lp - l_to_m # value of x 
                    
                    self.sequence.append((i, j, k, l_to_l, m_to_r, r_to_m))
                    Util.transposition(self.I, (i, j, k, l_to_l, m_to_r, r_to_m))
                    return True
              m_to_r += 1
    return False

  def rev_remove_two_bp(self, i, j) :
    if (1 <= i < j <= self.I.size) : # at this point we already know that elements in positions (i-1) and (j) from self.I.pi are consecutive
      if (self.I.pi[i] != self.I.pi[i-1] + 1 and self.I.pi[j] != self.I.pi[j+1] + 1 and self.I.pi[i] != self.I.pi[i-1] - 1 and self.I.pi[j] != self.I.pi[j+1] - 1) : #positions where the reversal takes place must be soft breakpoints
        if (self.I.pi[j+1] == self.I.pi[i] + 1 or self.I.pi[j+1] == self.I.pi[i] - 1) : #check if elements in positions (i) and (j+1) from self.I.pi are consecutive
          lp, rp = self.I.bpi[i], self.I.bpi[j+1] # check if it is j+1 or j
          li, ri = self.I.biota[max(self.I.pi[(i-1)], self.I.pi[j])], self.I.biota[max(self.I.pi[(j+1)], self.I.pi[i])]
          if (lp + rp == li +ri):
            l_to_l = min(lp, li)           # value of x
            r_to_l = min(rp, li - l_to_l)  # value of y
            Util.reversal(self.I, (i, j, l_to_l, r_to_l))
            self.sequence.append((i, j, l_to_l, r_to_l))
            return True
    return False


  def transp_remove_three(self):
    inverse = Permutation.inverse(self.I.pi)
    for i in range(1,self.I.size) :
      k = inverse[self.I.pi[i]+1]+1
      if (k <= self.I.size) :
        j = inverse[self.I.pi[k] + 1]+1
        if self.transp_remove_three_bp(i, j, k) :
            return True
        j = inverse[self.I.pi[k] - 1]+1
        if self.transp_remove_three_bp(i, j, k) :
            return True
      k = inverse[self.I.pi[i]-1]+1
      if (k <= self.I.size) :
        j = inverse[self.I.pi[k] + 1]+1
        if self.transp_remove_three_bp(i, j, k) :
            return True
        j = inverse[self.I.pi[k] - 1]+1
        if self.transp_remove_three_bp(i, j, k) :
            return True
    return False


  def rev_transp_remove_two(self):
    inverse = Permutation.inverse(self.I.pi)
    for i in range(1,self.I.size) :
      j = inverse[self.I.pi[i-1] + 1] # leftmost will be removed
      if (self.rev_remove_two_bp(i, j)):
        return True
      k = inverse[self.I.pi[i] + 1] + 1 # middle will be removed
      if (self.transp_remove_two_bp(i, j, k, 12)):
        return True
      k = inverse[self.I.pi[j-1] + 1] # rightmost will be removed
      if (self.transp_remove_two_bp(i, j, k, 13)):
        return True
      k = inverse[self.I.pi[i] - 1] + 1 # middle will be removed
      if (self.transp_remove_two_bp(i, j, k, 12)):
        return True
      if (j >= 1 and self.I.pi[j-1] >= 1):
          k = inverse[self.I.pi[j-1] - 1] # rightmost will be removed
          if (self.transp_remove_two_bp(i, j, k, 13)):
            return True
      if (i >= 1 and self.I.pi[i-1] >= 1):
          j = inverse[self.I.pi[i-1] - 1] # leftmost will be removed
          if (self.rev_remove_two_bp(i, j)):
            return True
          k = inverse[self.I.pi[i] + 1] + 1 # middle will be removed
          if (self.transp_remove_two_bp(i, j, k, 12)):
            return True
          if (j >= 1 and self.I.pi[j-1] >= 1):
              k = inverse[self.I.pi[j-1] + 1] # rightmost will be removed
              if (self.transp_remove_two_bp(i, j, k, 13)):
                return True
          k = inverse[self.I.pi[i] - 1] + 1 # middle will be removed
          if (self.transp_remove_two_bp(i, j, k, 12)):
            return True
          if (j >= 1 and self.I.pi[j-1] >= 1):
              k = inverse[self.I.pi[j-1] - 1] # rightmost will be removed
              if (self.transp_remove_two_bp(i, j, k, 13)):
                return True
      k = inverse[self.I.pi[i] + 1] + 1 # middle will be removed
      if (k <= self.I.size):
          j = inverse[self.I.pi[k] + 1] + 1 # rightmost will be removed
          if (self.transp_remove_two_bp(i, j, k, 23)):
            return True 
          j = inverse[self.I.pi[k] - 1] + 1 # rightmost will be removed
          if (self.transp_remove_two_bp(i, j, k, 23)):
            return True

      k = inverse[self.I.pi[i] - 1] + 1 # middle will be removed
      if (k <= self.I.size):
          j = inverse[self.I.pi[k] + 1] + 1 # rightmost will be removed
          if (self.transp_remove_two_bp(i, j, k, 23)):
            return True
          j = inverse[self.I.pi[k] - 1] + 1 # rightmost will be removed
          if (self.transp_remove_two_bp(i, j, k, 23)):
            return True

    return False


######################
####  RUN METHOD  ####
######################


  def run(self):
    self.update()
    while len(self.breakpoints):
      if self.greedy and (self.transp_remove_three() or self.rev_transp_remove_two()): # Greedy strategy 
        self.update()
      elif self.case_i() or self.case_ii() or self.case_iii() or self.case_iv(): # Algorithm cases
        self.update()
      else:
        raise ValueError('Algorithm: None of the cases I, II, III, or IV was performed.')
    print(self.lower_bound, len(self.sequence), str(self.sequence).replace(' ', ''))


######################
####     MAIN     ####
######################

if __name__ == "__main__":
  try :
    opts, args = getopt.getopt(sys.argv[1:], "hs:",["file=", "greedy="])
  except getopt.GetoptError:
    sys.exit()
  greedy = False
  file = None
  for opt,arg in opts:
    if(opt == "--file"):
      file = arg
    if(opt == "--greedy"):
      greedy = bool(arg)
  f = open(file, 'r')
  line = f.readline()
  while line:
    pi, bpi, biota = [eval('[' + i + ']') for i in line.split()]
    I = Instance(pi, bpi, biota)
    A = Algorithm(I, greedy)
    A.run()
    line = f.readline()
  f.close()