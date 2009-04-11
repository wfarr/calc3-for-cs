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

output = File.new("data3.txt", "w+")
content = ""

leslie = Matrix[[0,  1.2, 1.1,0.9,0.1 ,0  ,0   ,0   ,0],
                [0.7,0,   0,  0,  0,   0  ,0   ,0   ,0],
                [0,  0.85,0,  0,  0,   0  ,0   ,0   ,0],
                [0,  0,   0.9,0,  0,   0  ,0   ,0   ,0],
                [0,  0,   0,  0.9,0,   0  ,0   ,0   ,0],
                [0,  0,   0,  0,  0.88,0  ,0   ,0   ,0],
                [0,  0,   0,  0,  0,   0.8,0   ,0   ,0],
                [0,  0,   0,  0,  0,   0,  0.77,0   ,0],
                [0,  0,   0,  0,  0,   0,  0,   0.40,0]]

# population distributions

x_0 = Matrix[[2.1], [1.9], [1.8], [2.1], [2.0], [1.7], [1.2], [0.9], [0.5]]
x_1 = leslie * x_0 #=> Matrix[[6.35], [1.47], [1.615], [1.62], [1.89], [1.76], [1.36], [0.924], [0.36]]
x_2 = leslie * x_1 #=> Matrix[[5.1875], [4.445], [1.2495], [1.4535], [1.458], [1.6632], [1.408], [1.0472], [0.3696]]
x_3 = leslie * x_2 #=> Matrix[[8.1624], [3.63125], [3.77825], [1.12455], [1.30815], [1.28304], [1.33056], [1.08416], [0.41888]]
x_4 = leslie * x_3 #=> Matrix[[9.656485], [5.71368], [3.0865625], [3.400425], [1.012095], [1.151172], [1.026432], [1.0245312], [0.433664]]
x_5 = leslie * x_4 #=> Matrix[[13.41322675], [6.7595395], [4.856628], [2.77790625], [3.0603825], [0.8906436], [0.9209376], [0.79035264], [0.40981248]]

# total populations

populations = [
               x_0.to_a.flatten.inject(&:+) * 10000,
               x_1.to_a.flatten.inject(&:+) * 10000,
               x_2.to_a.flatten.inject(&:+) * 10000,
               x_3.to_a.flatten.inject(&:+) * 10000,
               x_4.to_a.flatten.inject(&:+) * 10000,
               x_5.to_a.flatten.inject(&:+) * 10000
              ]

content << "2010: #{populations[1]}
2020: #{populations[2]}
2030: #{populations[3]}
2040: #{populations[4]}
2050: #{populations[5]}
"

# percentage of population change

content << "2000 -> 2010: %.2f\%
2010 -> 2020: %.2f\%
2020 -> 2030: %.2f\%
2030 -> 2040: %.2f\%
2040 -> 2050: %.2f\%" % [((populations[1] - populations[0]) / populations[0] * 100),
                         ((populations[2] - populations[1]) / populations[1] * 100),
                         ((populations[3] - populations[2]) / populations[2] * 100),
                         ((populations[4] - populations[3]) / populations[3] * 100),
                         ((populations[5] - populations[4]) / populations[4] * 100)]

# question 4

x_2new = x_2.to_a
x_2new[1][0] = x_2new[1][0] / 2
x_2new = Matrix[*x_2new]

x_3new = leslie * x_2new
x_4new = leslie * x_3new
x_5new = leslie * x_4new



# total populations

populations_new = [
                   x_2new.to_a.flatten.inject(&:+) * 10000,
                   x_3new.to_a.flatten.inject(&:+) * 10000,
                   x_4new.to_a.flatten.inject(&:+) * 10000,
                   x_5new.to_a.flatten.inject(&:+) * 10000
                  ]

content << "\n\nWith 2020's second set modified.
"

content <<  "2020: #{populations_new[0]}
2030: #{populations_new[1]}
2040: #{populations_new[2]}
2050: #{populations_new[3]}
"

# percentage of population change

content << "2020 -> 2030: %.2f\%
2030 -> 2040: %.2f\%
2040 -> 2050: %.2f\%" % [((populations_new[1] - populations_new[0]) / populations_new[0] * 100),
                         ((populations_new[2] - populations_new[1]) / populations_new[1] * 100),
                         ((populations_new[3] - populations_new[2]) / populations_new[2] * 100)]

output.write(content)
