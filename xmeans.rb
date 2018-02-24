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
      @cid = cls_id
      @matrix = mat
    else  # 行列を与えない場合は、2次元空間で100個の点をランダムに配置
      @cid = cls_id
      @matrix = NArray.dfloat(2, 100).random(100 - 2*10) + 10 
    end
  end

  # クラスタの重心座標
  def centroid
    cent = []
    @matrix.dim.times do |i|
      cent[i] = @matrix[i, true].sum / @matrix.sizes[1]
    end
    return NArray.to_na(cent)
  end

  # クラスタの分散
  def variance
    cent = centroid()
    return Math.sqrt( ((@matrix - cent)**2).sum) / @matrix.sizes[1]
  end

  # 対数尤度
  def likelihood(cluster_num)
    node_num = @matrix.sizes[1]
    result = -((node_num / 2.0) * Math.log(2.0 * Math::PI)) -
      ( ((node_num * @matrix.dim) / 2.0) * Math.log(variance()) ) -
      ( ((node_num - cluster_num) / 2.0) ) +
      ( node_num * Math.log(node_num) )
    return result
  end


  def xmeans(div_limit=10)
    # はじめに2分割して、2つのparent clusterを生成
    labels, centroids = kmeans(2)
    parents = []  # 分割後のノードを格納する配列を用意

    2.times do |i|
      cls_size = labels.eq(i).sum
      parents[i] = NArray.dfloat(@matrix.sizes[0], cls_size)
    end

    # クラスタを配列に格納
    2.times do |i|
      cnt = 0
      labels.size.times {|j|
        if labels[j] == i
          parents[i][0, cnt] = @matrix[0, j]
          parents[i][1, cnt] = @matrix[1, j]
          cnt += 1
        end
      }
    end

    # さらに分割を試みる
    new_centroids = nil
    end_flag = false
    cluster_num = 2
    while (!end_flag)
      new_centroids = []
      next_parents = []
      new_labels = labels.dup

      parents_loop_cnt = 0  # 分割確定後のラベル更新で用いる
      bic_flag = false
      parents.each do |parent|
        p_cls = MyCluster.new(parents_loop_cnt - 2, parent)
        pBIC = likelihood(1) - 
          ((parent.dim + 1) / 2.0) * Math.log(parent.sizes[1])

        ch_labels, ch_centroids = p_cls.kmeans(2)
        children = []  # 分割後のノードを格納
        2.times do |i|
          cls_size = ch_labels.eq(i).sum
          children[i] = NArray.dfloat(p_cls.matrix.sizes[0], cls_size)
        end

        # クラスタを配列に格納
        2.times do |i|
          cnt = 0
          ch_labels.size.times {|j|
            if ch_labels[j] == i
              children[i][0, cnt] = @matrix[0, j]
              children[i][1, cnt] = @matrix[1, j]
              cnt += 1
            end
          }
        end

        c_likelihood = 0.0
        children.each do |child|
          c_cls = MyCluster.new(0, child)
          c_likelihood += likelihood(2) -
            ((1.0 + 2*child.dim + 2) / 2.0) * Math.log(parent.sizes[1])
        end

        cBIC = c_likelihood -
          ((1.0 + 2*parent.dim + 2) / 2.0) * Math.log(parent.sizes[1])
        puts "p: " + pBIC.to_s
        puts "c: " + cBIC.to_s

        # 分割すべきかを判定
        if (cBIC/1 > pBIC)
          #p ch_labels
          cluster_num += 1
          children.size.times do |i|
            next_parents << children[i]
          end

          # ch_centroidsはNArrayなので、一旦配列にする
          2.times do |i|
            tmp = ch_centroids[true, i]
            new_centroids << [tmp[0], tmp[1]]
          end

          # labelsの更新
          label_cnt = 0
          ch_label_cnt = 0
          labels.eq(parents_loop_cnt).each do |flg|
            if flg == 1 
              #p ch_label_cnt
              if ch_labels[ch_label_cnt] == 1
                new_labels[label_cnt] = cluster_num - 1
                ch_label_cnt += 1
              end
            end
            label_cnt += 1
          end
          #puts "\nlabels:\n"
          #labels.each{|l| print l}
          #puts "\nparent_num: " + (parents_loop_cnt + 1).to_s + " / " + cluster_num.to_s
          bic_flag = true 
        else  # BICにより分割しないと判定しなかった場合
          next_parents << parent
          new_centroids << [p_cls.centroid[0], p_cls.centroid[1]]
        end
        parents_loop_cnt += 1
      end  # of parents loop

      # parentを更新
      parents = next_parents
      labels = new_labels

      # 一度も分割が起こらなかったら終了
      if bic_flag == false
        end_flag = true
      elsif div_limit < cluster_num  # もしくはクラスタ数が制限を超えたら終了
        end_flag = true
      end
    end  # of while
    return labels, NArray.to_na(new_centroids)
  end

  def kmeans(centroid_num=3)
    cond = 1e-20  # 終了条件(重心の更新差分がこの値より小さい場合に計算を終了)

    # @matrixの各点がどのクラスタに所属するかをランダムに割り振り
    labels = NArray.int(matrix.sizes[1]).random(centroid_num)
    # クラスタに必ず1つ以上のノードが割り振られることを保証
    centroid_num.times {|i|
      labels[i] = i.to_f if labels[i] != nil
    }

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
        if (labels.eq(cent_num)).sum > 0.0
          cent = cent / (labels.eq(cent_num)).sum  # 平均をとって重心を算出
        end
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
    return [labels, centroids]
  end

  def output_data(centroids, data_dir, dimension=2)
    `mkdir #{data_dir}/#{dimension}d/`

    node_file = "#{data_dir}/#{dimension}d/node_data.tsv"
    cent_file = "#{data_dir}/#{dimension}d/centroid_data.tsv"

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
  centroid_num = 15

  # 2次元データのプロット
  ary = ""
  open("data/plot.csv"){|f| ary += f.read}
  ary = ary.split("\n")
  ary2 = []
  ary.each do |a|
    tmp = a.split(',')
    ary2 << [tmp[0].to_f, tmp[1].to_f]
  end

  data2d = NArray.to_na(ary2)
  #data2d = NArray.dfloat(2, 300).random(100)
  mc1 = MyCluster.new(1, data2d)
  mc1.xmeans
  labels, centroids = mc1.kmeans(centroid_num)
  mc1.output_data(centroids, "data", 2)

  # 3次元データのプロット
  data3d = NArray.dfloat(3, 200).random(100)
  mc2 = MyCluster.new(1, data3d)
  labels, centroids = mc2.kmeans(centroid_num)
  mc2.output_data(centroids, "data", 3)
end
