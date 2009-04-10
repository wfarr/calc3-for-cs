require 'pt3'

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
x_1 = leslie * x_0
x_2 = leslie * x_1
x_3 = leslie * x_2
x_4 = leslie * x_3
x_5 = leslie * x_4

# total populations

populations = [
               x_0.to_a.flatten.inject(&:+) * 10000,
               x_1.to_a.flatten.inject(&:+) * 10000,
               x_2.to_a.flatten.inject(&:+) * 10000,
               x_3.to_a.flatten.inject(&:+) * 10000,
               x_4.to_a.flatten.inject(&:+) * 10000,
               x_5.to_a.flatten.inject(&:+) * 10000
              ]

puts "2010: #{populations[1]}"
puts "2020: #{populations[2]}"
puts "2030: #{populations[3]}"
puts "2040: #{populations[4]}"
puts "2050: #{populations[5]}"

# percentage of population change

puts "2000 -> 2010: %.2f\%" % ((populations[1] - populations[0]) / populations[0] * 100)
puts "2010 -> 2020: %.2f\%" % ((populations[2] - populations[1]) / populations[1] * 100)
puts "2020 -> 2030: %.2f\%" % ((populations[3] - populations[2]) / populations[2] * 100)
puts "2030 -> 2040: %.2f\%" % ((populations[4] - populations[3]) / populations[3] * 100)
puts "2040 -> 2050: %.2f\%" % ((populations[5] - populations[4]) / populations[4] * 100)

# puts leslie.power_method

# assuming x(n) = Ax(n-1), and A never changes,
#
# x_1 => Matrix[[6.35], [1.47], [1.615], [1.62], [1.89], [1.76], [1.36], [0.924], [0.36]]
# x_2 => Matrix[[5.1875], [4.445], [1.2495], [1.4535], [1.458], [1.6632], [1.408], [1.0472], [0.3696]]
# x_3 => Matrix[[8.1624], [3.63125], [3.77825], [1.12455], [1.30815], [1.28304], [1.33056], [1.08416], [0.41888]]
# x_4 => Matrix[[9.656485], [5.71368], [3.0865625], [3.400425], [1.012095], [1.151172], [1.026432], [1.0245312], [0.433664]]
# x_5 => Matrix[[13.41322675], [6.7595395], [4.856628], [2.77790625], [3.0603825], [0.8906436], [0.9209376], [0.79035264], [0.40981248]]

