type KmedoidsResult
    medoids::Vector{Int} #indexes of medoids
    assignments::Vector{Int} #cluster assignments
end

function mean_dist(dist, j)
    n = size(dist)[1]

    function denom(i)
        sum(l -> dist[i,l], 1:n)
    end

    sum(i -> dist[i,j]/denom(i), 1:n)
end

function initial_medoids(dist, k)
    n = size(dist)[1]

    scores = sort([1:n], by=(j -> mean_dist(dist, j)))

    scores[1:k]
end

function find_clusters(dist, medoids)
    n = size(dist)[1]
    k = size(medoids)[1]
    total_dist = 0.0

    clusters = [Int[] for i = 1:k]

    for i = 1:n
        (distance, index) = findmin(dist[medoids,i])
        total_dist += distance

        push!(clusters[index], i) 
    end

    (total_dist, clusters)
end

function new_medoids(dist, clusters)    
    medoids = zeros(size(clusters)[1])

    for (i, cluster) in enumerate(clusters)
        dist_within_cluster = [sum([dist[i,j] for j in cluster]) for i in cluster]
        best = findmin(dist_within_cluster)[2]
        medoids[i] = best
    end

    medoids
end

function cluster_membership(clusters, n)
    assignments = zeros(n)
    for (i, cluster) in clusters
        for object in cluster
            assignments[object] = i
        end
    end
    assignments
end

function kmedoids{R <: FloatingPoint}(dist::Matrix{R}, k::Int)
    medoids = initial_medoids(dist, k)
    (total_dist, clusters) = find_clusters(dist, medoids)
    old_dist = Inf

    while total_dist < old_dist 
        old_dist = total_dist
        medoids = new_medoids(dist, clusters)
        (total_dist, clusters) = find_clusters(dist, medoids)
    end
    
    KmedoidsResult(medoids, cluster_membership(clusters, size(dist)[1]))
end

