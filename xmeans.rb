#!/usr/bin/env ruby

require "pp"
require "rubygems"
require "narray"

class MyCluster
  attr_accessor :cid, :matrix

  def initialize(cls_id, mat=nil, narray_flg=nil)
    if mat.class == NArray
      @cid = cls_id
      @matrix = mat
    else  # generate 100 points in 2d space if no array given
      @cid = cls_id
      @matrix = NArray.dfloat(2, 100).random(100 - 2*10) + 10 
    end
  end

  # centriod of the cluster
  def centroid
    cent = []
    @matrix.dim.times do |i|
      cent[i] = @matrix[i, true].sum / @matrix.sizes[1]
    end
    return NArray.to_na(cent)
  end

  # variance of the cluster
  def variance
    cent = centroid()
    return Math.sqrt( ((@matrix - cent)**2).sum) / @matrix.sizes[1]
  end

  def likelihood(cluster_num)
    node_num = @matrix.sizes[1]
    result = -((node_num / 2.0) * Math.log(2.0 * Math::PI)) -
      ( ((node_num * @matrix.dim) / 2.0) * Math.log(variance()) ) -
      ( ((node_num - cluster_num) / 2.0) ) +
      ( node_num * Math.log(node_num) )
    return result
  end

  def xmeans(div_limit=10)
    # First, generate two parent clusters
    labels, centroids = kmeans(2)
    parents = []  # contains divided nodes

    2.times do |i|
      cls_size = labels.eq(i).sum
      parents[i] = NArray.dfloat(@matrix.sizes[0], cls_size)
    end

    # add divided clusters into parents array
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

    # try to divide more
    new_centroids = nil
    end_flag = false
    cluster_num = 2
    while (!end_flag)
      new_centroids = []
      next_parents = []
      new_labels = labels.dup

      parents_loop_cnt = 0  # it's used for updating labels
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

        # add clusters to children array
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

        # decide if devide or not
        if (cBIC/1 > pBIC)
          #p ch_labels
          cluster_num += 1
          children.size.times do |i|
            next_parents << children[i]
          end

          # convert ch_centroids from NArray to arrary
          2.times do |i|
            tmp = ch_centroids[true, i]
            new_centroids << [tmp[0], tmp[1]]
          end

          # update labels
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
        else  # case of not divide with BIC
          next_parents << parent
          new_centroids << [p_cls.centroid[0], p_cls.centroid[1]]
        end
        parents_loop_cnt += 1
      end  # of parents loop

      # update parents
      parents = next_parents
      labels = new_labels

      # terminate if division is occured
      if bic_flag == false
        end_flag = true
      # terminate if the number of clusters reaches to the limit
      elsif div_limit < cluster_num
        end_flag = true
      end
    end  # of while
    return labels, NArray.to_na(new_centroids)
  end

  def kmeans(centroid_num=3)
    cond = 1e-20  # condition for termination kmeans

    # assign which cluster belongs to for each of points of given matrix
    labels = NArray.int(matrix.sizes[1]).random(centroid_num)
    # confirm that cluster has one or more nodes at least
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
          cent = cent / (labels.eq(cent_num)).sum  # calculate centriod
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

      # re-calculation for each point of matrix
      labels = update_labels(labels, centroids)

      centroids.sizes[1].times {|i|
      # terminate if difference between centroids and all of prev_centroids
      # is less than cond
        if cond < Math.sqrt(
            ((prev_centroids[true, i] - centroids[true, i])**2).sum
        )
          next	
        end
        end_flag = true
      }

      prev_centroids = centroids.dup
    end  # of while
    return [labels, centroids]
  end

  def output_data(centroids, data_dir, dimension=2)
    `mkdir -p #{data_dir}/#{dimension}d/`

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
      label = nil

      point = @matrix[true, cnt]
      centroids.sizes[1].times {|i|
        dist = ( (point - centroids[true, i])**2 ).sum
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

def main(input_file=nil)
  centroid_num = 15

  # plot 2d data
  #ary = ""
  #open(input_file){|f| ary += f.read}
  #ary = ary.split("\n")
  #ary2 = []
  #ary.each do |a|
  #  tmp = a.split(',')
  #  ary2 << [tmp[0].to_f, tmp[1].to_f]
  #end

  #data2d = NArray.to_na(ary2)
  data2d = NArray.dfloat(2, 300).random(100)
  cluster_id = 1
  mc1 = MyCluster.new(cluster_id, data2d)
  labels, centroids = mc1.kmeans(centroid_num)
  mc1.output_data(centroids, "data/kmeans", 2)
  labels, centroids = mc1.xmeans
  mc1.output_data(centroids, "data/xmeans", 2)

  # plot 3d data
  data3d = NArray.dfloat(3, 200).random(100)
  mc2 = MyCluster.new(cluster_id, data3d)
  labels, centroids = mc2.kmeans(centroid_num)
  mc2.output_data(centroids, "data/kmeans", 3)
end

if __FILE__ == $0
  main(ARGV[0])
end
