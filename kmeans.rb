#!/usr/bin/env ruby

require "rubygems"
require "narray"

def k_cluster(matrix, cls_num)
end

def get_euclidean_similarity(critics, person1, person2)
	movies = critics[person1].keys & critics[person2].keys
	return 0 if movies.empty?

	critics_person1 = NArray.to_na(movies.map { |movie| critics[person1][movie] })
	critics_person2 = NArray.to_na(movies.map { |movie| critics[person2][movie] })

	# 映画ごとの批評差を求める
	distances = (critics_person1 - critics_person2) ** 2
	# 批評差の総和を正規化して返す
	normalize(distances.sum)
end

# ピアソン相関定義
def get_pearson(array_x, array_y)
	x = NArray.to_na(array_x)
	y = NArray.to_na(array_y)

	n = x.size
	sum_x = x.sum
	sum_y = y.sum
	sum_x_square = (x ** 2).sum
	sum_y_square = (y ** 2).sum
	sum_xy = (x * y).sum

	numerator = sum_xy - (sum_x * sum_y.quo(n))
	denominator = sqrt(
		(sum_x_square - (sum_x ** 2).quo(n)) * (sum_y_square - (sum_y ** 2).quo(n)))

		denominator == 0 ? 0 : numerator.quo(denominator)
end
