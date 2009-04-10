#!/usr/bin/env ruby1.9

require 'common'

class Matrix
  #Compute largest eigenvalue of a square matrix
  def power_method
    a = self
    #bArr is an array of ones, the original guess
    bArr = Array.new(a.column_size) { 1 }
    #wArr is [1,0,0,0,0...]
    wArr = Array.new(a.column_size) { |i| i == 0 ? 1 : 0 }
    #r,u,l values used to compute eigenvalue
    r, u, l = 0, 0, 0
    q = a * Vector[*bArr].covector.transpose
    qArrOld = q.to_a
    qArrNew = []
    #Converts aArrNew into an array of numbers instead of arrays
    for i in 0...qArrOld.size
      qArrNew.push(qArrOld[i][0])
    end

    u = u + Vector[*wArr].inner_product(Vector[*qArrNew])
    l = l + Vector[*wArr].inner_product(Vector[*bArr])
    r_prev = r
    r = (u/l).to_f
    b = q
    bArr = qArrNew
    q = a * Vector[*bArr].covector.transpose
    qArrOld = q.to_a
    qArrNew = []
    for i in 0...qArrOld.size
      qArrNew.push(qArrOld[i][0])
    end

    #Iteration for eigenvalue
    while true
      u = u + Vector[*wArr].inner_product(Vector[*qArrNew])
      l = l + Vector[*wArr].inner_product(Vector[*bArr])
      r_prev = r
      r = (u/l).to_f
      b = q
      bArr = qArrNew
      q = a * Vector[*bArr].covector.transpose
      qArrOld = q.to_a
      qArrNew = []
      for i in 0...qArrOld.size
        qArrNew.push(qArrOld[i][0])
      end
      break if (r - r_prev).abs <= 0.00000001
    end
    return r
  end
end
