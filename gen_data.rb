#!/usr/bin/env ruby
# coding: utf-8

def main()
  nof_plot = 100
  f = open("plot.csv", "w+")
  nof_plot.times do
    x = rand(100)
    y = rand(100)
    f.write("#{x}, #{y}\n")
  end
  f.close()
end

main()
