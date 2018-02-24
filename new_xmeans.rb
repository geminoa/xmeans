#!/usr/bin/env ruby

require "pp"
require "rubygems"
require "narray"

class MyCluster
  attr_accessor :cid, :matrix

  def self.gen_random_2d_dfloat(row_num=100, max=100, offset=10)
    ret = NArray.dfloat(2, row_num.to_i).random(max - 2*offset) + offset
    return ret
  end

  def initialize(cls_id, mat=nil, narray_flg=nil)
    if mat.class == NArray
      @cid = cls_id.to_i
      @matrix = mat
    else  # 行列を与えない場合は、2次元空間で100個の点をランダムに配置
      @cid = cls_id.to_i
      @matrix = NArray.dfloat(2, 100).random(100 - 2*10) + 10 
    end
  end

  def xmeans()
    kmeans(5)
  end

  def kmeans(centroid_num=3)
    cond = 1e-30  # 終了条件(重心の更新差分がこの値より小さい場合に計算を終了)

    # @matrixの各点がどのクラスタに所属するかをランダムに割り振り
    labels = NArray.int(matrix.sizes[1]).random(centroid_num)

    end_flag = false
    prev_centroids = nil
    centroids = nil

    while( !end_flag )
      centroids = NArray.dfloat(@matrix.sizes[0], centroid_num)
      centroid_num.times {|cent_num|
        cent = NArray.dfloat(@matrix.sizes[0])
        labels.size.times {|cnt|
          if labels[cnt] == cent_num
            cent += @matrix[true, cnt]
          end
        }
        cent = cent / (labels.eq(cent_num)).sum  # 平均をとって重心を算出
        @matrix.sizes[0].times {|i|
          centroids[i, cent_num] = cent[i]
        }
      }

      if !prev_centroids
        prev_centroids = centroids.dup
        labels = update_labels(labels, prev_centroids)
        next	
      end

      # @matrixの各点について、どの重心に最も近いのかを再計算
      # labelsの更新
      labels = update_labels(labels, centroids)

      # 終了判定
      centroids.sizes[1].times {|i|
        if cond < Math.sqrt( ( (prev_centroids[true, i] - centroids[true, i])**2 ).sum )
          next	
        end
        end_flag = true  # 全てのprev_centroidsとcentroidsの差がcond以内なら終了
      }

      prev_centroids = centroids.dup
    end  # of while
    return centroids
  end

  def output_data(centroids, dimension=2)
    node_file = "node_data_#{dimension}d.tsv"
    cent_file = "centroid_data_#{dimension}d.tsv"

    open(node_file, "w+") {|f|
      @matrix.sizes[1].times do |i|
        po = @matrix[true, i]
        dimension.times {|dim|
          f.write(po[dim].to_s + "\t")
        }
        f.write("\n\n")
      end
    }

    open(cent_file, "w+") {|f|
      centroids.sizes[1].times do |i|
        po = centroids[true, i]
        dimension.times {|dim|
          f.write(po[dim].to_s + "\t")
        }
        f.write("\n\n")
      end
    }
  end	

  private
  def update_labels(labels, centroids)
    labels.size.times {|cnt|	
      min_dist = nil
      label = nil  # 更新されたlabelを代入

      point = @matrix[true, cnt]
      centroids.sizes[1].times {|i|
        dist = ( (point - centroids[true, i])**2 ).sum  # 距離を計算
        if min_dist
          if (min_dist > dist)
            min_dist = dist 
            label = i
          end
        else
          min_dist = dist
          label = i 
        end
      }
      labels[cnt] = label
    }
    return labels
  end

end

if __FILE__ == $0
  centroid_num = 5

  # 2次元データのプロット
  data2d = NArray.dfloat(2, 100).random(100)
  mc1 = MyCluster.new(1, data2d)
  centroids = mc1.kmeans(centroid_num)
  mc1.output_data(centroids, 2)

  # 3次元データのプロット
  data3d = NArray.dfloat(3, 100).random(100)
  mc2 = MyCluster.new(1, data3d)
  centroids = mc2.kmeans(centroid_num)
  mc2.output_data(centroids, 3)
end
