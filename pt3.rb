#!/usr/bin/env ruby1.9

require 'common'

class Matrix
  #
  # Compute largest eigenvalue of a square matrix
  #
  def power_method
    a = self
    b = Array.new(a.column_size) { 1 }
    w = Array.new(a.column_size) { |i| i == 0 ? 1 : 0 }
    r, u, l = 0, 0, 0
    q = a * Vector[*b].covector.transpose
    q_old = q.to_a
    q_new = []

    for i in 0...q_old.size
      q_new.push(q_old[i][0])
    end

    u = u + Vector[*w].inner_product(Vector[*q_new])
    l = l + Vector[*w].inner_product(Vector[*b])
    r_prev = r
    r = (u/l).to_f
    b = q_new
    q = a * Vector[*b].covector.transpose
    q_old = q.to_a
    q_new = []

    for i in 0...q_old.size
      q_new.push(q_old[i][0])
    end

    while true
      u = u + Vector[*w].inner_product(Vector[*q_new])
      l = l + Vector[*w].inner_product(Vector[*b])
      r_prev = r
      r = (u/l).to_f
      b = q
      b = q_new
      q = a * Vector[*b].covector.transpose
      q_old = q.to_a
      q_new = []
      for i in 0...q_old.size
        q_new.push(q_old[i][0])
      end
      break if (r - r_prev).abs <= 0.00000001
    end
    return r
  end
end
