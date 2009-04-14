#!/usr/bin/env ruby1.9

require 'mathn'
require 'fileutils'

class Vector
  private
  def sign(x)
    return 1 if x > 0
    return -1 if x < 0
    return 0
  end
end

class Matrix
  #
  # Prints out the Matrix in a nice format
  #
  def pretty_print
    str = ""
    self.to_a.each do |row|
      row.each do |i|
        if i.to_i >= 0
          str << " "
        end
        if ("%.3f" % i).to_f == i.to_i
          str << "#{i.to_i}     "
        else
          str << "%.3f " % i
        end
      end
      str << "\n"
    end
    puts str
  end
  #
  # Finds the norm where p = infinity
  #
  def inf_norm
    self.to_a.map do |a|
      a.map do |ar|
        ar.abs
      end.inject(&:+)
    end.sort[0]
  end
  #
  # Determines whether or not a matrix is lower triangular, where all values above
  # the diagonal are 0.
  #
  def is_lower_triangular?
    triangular(self.column_vectors)
  end
  #
  # Determines whether or not a matrix is upper triangular, where all values below
  # the diagonal are 0.
  #
  def is_upper_triangular?
    triangular(self.row_vectors)
  end
  #
  # Creates an +n+ by +n+ Hilbert matrix.
  #
  def self.hilbert(n)
    m = Matrix.zero(n).to_a
    m = m.each_index.map{|row| m[row].each_index.map{|col| 1 / (row + col + 1)}}
    return Matrix.rows(m)
  end

  private
  def triangular(vecs)
    for i in 0...vecs.length
      vec = vecs[i].to_a
      unless i <= 1
        return false unless vec[0...i].all? { |n| n == 0 } and vec[i..-1].all? { |n| n != 0 }
      end
    end
    return true
  end
end
