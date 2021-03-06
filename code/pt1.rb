#!/usr/bin/env ruby1.9

require 'common'

class Vector
  def find_householder_reflection
    a = self.to_a
    a = a[0] if a[0].is_a?(Array)
    a[0] = a[0] + sign(a[0]) * self.r
    u = Vector[*a]
    norm_u_sqrd = u.r**2
    uut = u.covector.transpose * u.covector
    h = Matrix.identity(uut.row_size) - (uut * (2 / norm_u_sqrd))
    return h
  end
end

class Matrix

  #
  # Perform an LU decomposition.
  #
  def lu_decomposition
    return nil unless self.square?
    n = self.row_size()
    a = self
    l_n = []
    cvs = a.column_vectors.map { |v| v.to_a }
    for k in 0...cvs.length
      for j in 0...cvs.length
        l_new = Matrix.identity(n).to_a
        if l_new[j][k] == 1 || j < k
          next
        end
        l_new[j][k] = - (cvs[k][j] / cvs[k][k])
        l_n << l_new
        a = Matrix[*l_new] * Matrix[*cvs.transpose]
        cvs = a.column_vectors.map { |v| v.to_a }
      end
    end
    l_final = l_n.map { |m| Matrix[*m].inverse }.inject(&:*)
    u_final = a
    return l_final,u_final
  end

  #
  # Performs a Householder reflection for self to find its QR decomposition.
  #
  def householder
    return nil unless self.square?
    current_iteration = self
    init_dim = self.row_size
    h_list = []
    cv = current_iteration.column_vectors[0]
    h = (cv.find_householder_reflection - Matrix.identity(cv.size)).expand_to_dimensions(init_dim,init_dim) + Matrix.identity(init_dim)
    h_list << h
    current_iteration = h * current_iteration
    for i in 0...self.row_size
      cv = current_iteration.get_column_vector(i+1)
      break if cv.size < 2 || current_iteration.is_upper_triangular?
      h = (cv.find_householder_reflection - Matrix.identity(cv.size)).expand_to_dimensions(init_dim,init_dim) + Matrix.identity(init_dim)
      h_list << h
      current_iteration = h * current_iteration
    end
    q,r = h_list.inject(&:*), current_iteration
    return q,r
  end

  #
  # Performs a Givens rotation for self to find its QR decomposition.
  #
  def givens
    return nil unless self.square?
    n = self.row_size
    a = self
    g_n = []
    cvs = a.column_vectors.map { |v| v.to_a }
    for i in 0...cvs.length
      for j in 0...cvs.length
        next unless j > i
        g = Matrix.identity(n).to_a
        c = cvs[i][i] / Math.sqrt(cvs[i][i]**2 + cvs[i][j]**2)
        s = -cvs[i][j] / Math.sqrt(cvs[i][i]**2 + cvs[i][j]**2)
        g[i][i], g[j][j] = c, c
        g[j][i], g[i][j] = s, -s
        g = Matrix[*g]
        g_n << g
        a = g * a
        cvs = a.column_vectors.map { |v| v.to_a }
      end
    end
    q,r = g_n.map { |m| m.t }.inject(&:*), a
    return q,r
  end
  #
  # Expands a matrix to x,y
  #
  def expand_to_dimensions(x,y)
    curr_x, curr_y, a = self.row_size, self.column_size, self.to_a
    a.each_index do |row|
      for i in 0...(y - curr_y)
        a[row] = a[row].insert(0,0)
      end
    end
    for i in 0...(x - curr_x)
      a = a.insert(0,Array.new(y){0})
    end
    return Matrix.rows(a)
  end
  #
  # Finds a shorter cv
  #
  def get_column_vector(x)
    return Vector.elements(self.column(x)[x..-1])
  end
end


def solve_lu_hilbert(size)
  h = Matrix.hilbert(size)
  l,u = h.lu_decomposition
  b = Vector.elements(Array.new(size){1})
  x = u.inverse * l.inverse * b.covector.transpose
  err1 = ((l * u) - h).inf_norm
  err2 = ((h * x) - b).column(0).r
  return x, err1, err2
end

def solve_hilbert_householder(size)
  h = Matrix.hilbert(size)
  q,r = h.householder
  b = Vector.elements(Array.new(size){1})
  x = r.inverse * q.inverse * b.covector.transpose
  err1 = ((q * r) - h).inf_norm
  err2 = ((h * x) - b).column(0).r
  return x, err1, err2
end

def solve_hilbert_givens(size)
  h = Matrix.hilbert(size)
  q,r = h.givens
  b = Vector.elements(Array.new(size){1})
  x = r.inverse * q.inverse * b.covector.transpose
  err1 = ((q * r) - h).inf_norm
  err2 = ((h * x) - b).column(0).r
  return x, err1, err2
end

output = File.new("data1.txt", "w+")
towrite = ""

# LU
towrite << "\n=========================================\n=========== LU Decomposition ============\n=========================================\n"
for i in 2..20
  sol,err1,err2 = solve_lu_hilbert(i)
  towrite << "N = #{i}\n"
  towrite << "sol = #{sol.transpose.row_vectors[0]}\n"
  towrite << "err1 = #{err1}\n"
  towrite << "err2 = #{err2}\n"
  towrite << "\n"
end

# Householder
towrite << "\n=========================================\n========= Householder Method  ===========\n=========================================\n"
for i in 2..20
  sol,err1,err2 = solve_hilbert_householder(i)
  towrite << "N = #{i}\n"
  towrite << "sol = #{sol.transpose.row_vectors[0]}\n"
  towrite << "err1 = #{err1}\n"
  towrite << "err2 = #{err2}\n"
  towrite << "\n"
end

# Givens
towrite << "\n=========================================\n========== Givens Rotations  ============\n=========================================\n"
for i in 2..20
  sol,err1,err2 = solve_hilbert_givens(i)
  towrite << "N = #{i}\n"
  towrite << "sol = #{sol.transpose.row_vectors[0]}\n"
  towrite << "err1 = #{err1}\n"
  towrite << "err2 = #{err2}\n"
  towrite << "\n"
end

output.write(towrite)
