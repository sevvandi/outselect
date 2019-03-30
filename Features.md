Features
================

We have categorized the features in two main classes, with the first being generic features which measure various aspects of a dataset but are not particularly tailored towards outlier detection. Generic features can be further subdivided into 3 categories: simple, statistical and information theoretic based.

### Generic features

1.  **Simple features** These are related to the basic structure of a dataset. For our study these include the number of observations, number of attributes, ratio of observations to attributes, number of binary attributes, number of numerical attributes and the ratio of binary to numerical attributes.

2.  **Statistical features** These include statistical properties of skewness, kurtosis, mean to standard deviation ratio, and IQR to standard deviation ratio for all attributes of a dataset, i.e. if a dataset has *d* attributes, then there are *d* values for each statistical property. For skewness and kurtosis we include the mean, median, maximum and the 95% percentile of these *d* values as features. For IQR to standard deviation ratio we include the maximum and the 95% percentile. We also include the average mean to standard deviation ratio. As a correlation measure, we compute the correlation between all attributes and include the mean absolute correlation. We also perform Principal Component Analysis and include the standard deviation explained by the first Principal Component.

3.  **Information theoretic features** for measures of the amount of information present in a dataset, we compute the entropy of each attribute and include the mean in our feature set. Also we include the entropy of the whole dataset and the mutual information.

### Outlier features

In order to make our feature set richer we include density-based, residual-based and graph-based features. We compute these features for different subsets of the dataset, namely outliers, non-outliers and proxi-outliers. We define proxi-outliers as data points that are either far away from \`\`normal'' data or residing in low density regions. If there is a significant overlap between proxi-outliers and actual outliers, then we expect outlier detection methods to perform well on such datasets. The density, residual and graph-based features we consider are all ratios. It is either a ratio between proxi-outliers and outliers, or a ratio between outliers and non-outliers. An example is the ratio between average density of non-outliers and average density of outliers. We explain the outlier features below:

Features based on outliers and proxi-outliers fall into the category of landmarking features, which has been popular in meta-learning studies (Pfahringer, Bensusan, and Giraud-Carrier 2000, Peng et al. (2002),Smith-Miles (2009)). The concept of *landmarking* is to obtain a rough picture of the problem space with the aid of some simple tools. The problem space consists of all the datasets that are available with respect to the set of algorithms. Landmarking obtains performance measurements by running simple algorithms for each available dataset. The motivation for using landmarking features in meta-learning was that running simple algorithms might be faster and less expensive than computing certain complex features. For example, doing a simple regression to ascertain the degree of separability of a dataset may provide insight as to how complex algorithms perform classification on that dataset. As such, landmarking features can greatly contribute to finding strengths and weaknesses of algorithms.

#### Density based features

The density based features are computed either using the density based clustering algorithm DBSCAN (Ester et al. 1996) or kernel density estimation as follows:

1.  **DBSCAN features** We perform principal component analysis (PCA) and use DBSCAN for clustering in a lower-dimensional space. We focus on data-points that either belong to very small clusters, or do not belong to any cluster. Let us call these points dbscan-proxi-outliers, henceforth named dbscan-proxies. Once again, if dbscan-proxies are outliers then we expect density based outlier algorithms to perform well on such datasets. As features we include ..1. percentage of dbscan-proxies that are outliers, ..2. percentage of dbscan-proxies that are outliers which do not belong to any cluster and ..3. percentage of dbscan-proxies that are outliers which belong to very small clusters.

2.  **Kernel density estimate (KDE) related features** Here too, we perform PCA and compute kernel density estimates (KDE) on two dimensional PC spaces to reduce computational burden. We compute KDE as detailed by (Duong 2018) for the first 10 principal component (PC) pairs, and for each PC pair we find proxi-outliers. We compute the mean, median, standard deviation, IQR, minimum, maximum, 5th, and 95th percentiles of KDE values for outliers, non-outliers, proxi-outliers and non-proxi-outliers in each PC space. Next we take the computed summary statistics ratios of outliers to non-outliers, and proxi-outliers to non-proxi-outliers; for example the ratio between the mean KDE of non-outliers and the mean KDE of outliers. These ratios are computed for the first 10 two-dimensional PC spaces. As features, we include the average of each ratio for the set of PC spaces. In addition we also include the percentage of proxi-outliers that are outliers in our feature set. ..1.Local density features: We compute all features explained in 2) using a local density measure based on KDE instead of using KDE itself. The local density is computed by dividing the KDE value of a point by the mean KDE value of its *k*-nearest neighbours.

#### Residual based features

These features include summary statistics of residuals from linear models. First, we fit a series of linear models by randomly choosing the dependent variable and treating the rest of attributes as independent variables. For each model, data-points which have the top 3% of absolute residual values are deemed as proxi-outliers. Similar to KDE features, the mean, median, standard deviation, IQR, minimum, maximum, 5th, and 95th percentiles of residual values for outliers, non-outliers, proxi-outliers and non-proxi-outliers are computed. Next the respective ratios are computed for each linear model. Finally, the average of each ratio for all the models is included in the feature set. We also include the percentage of proxi-outliers that are outliers as a feature.

#### Graph based features

These features are based on graph-theoretic measures such as vertex degree, shortest path and connected components. First, we convert the dataset to a directed-graph based on each data-points' *k*-nearest neighbours using the software \*igraph (Csardi and Nepusz 2006). Next, we compute the degree of the vertices and label ones with the lowest degree as proxi-outliers. Then, similar to residual based features, we find the summary statistics of degree values for outliers, non-outliers, proxi-outliers and non-proxi-outliers and include ratios of outliers to non-outliers and proxi-outliers to non-proxi-outliers in our feature set. We also include the percentage of proxi-outliers which are actual outliers. Another set of graph features come from connected components. We compute the number of vertices in each connected component and, similar to degree calculations, compute summary statistics of these values for outliers, non-outliers, proxi-outliers and non-proxi-outliers and include the ratios as above. We also compute the shortest distances from outlier vertices to non-outlier vertices. Here the shortest distance from vertex *a* to vertex *b* is the minimum number of edges that connect *a* and *b*. We include some statistics about these distances: the percentage of outliers which have infinite shortest distance to all non-outlier vertices, i.e. outlier vertices which are not connected. For outliers that have finite shortest distances to some non-outlying vertex, we compute the percentage of outliers for which the shortest distance is 1.

References
----------

Csardi, Gabor, and Tamas Nepusz. 2006. “The Igraph Software Package for Complex Network Research.” *InterJournal, Complex Systems* 1695 (5): 1–9.

Duong, Tarn. 2018. *Ks: Kernel Smoothing*. <https://CRAN.R-project.org/package=ks>.

Ester, Martin, Hans-Peter Kriegel, Jörg Sander, Xiaowei Xu, and others. 1996. “A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise.” In *KDD*, 96:226–31. 34.

Peng, Yonghong, Peter A Flach, Carlos Soares, and Pavel Brazdil. 2002. “Improved Dataset Characterisation for Meta-Learning.” In *International Conference on Discovery Science*, 141–52. Springer.

Pfahringer, Bernhard, Hilan Bensusan, and Christophe G Giraud-Carrier. 2000. “Meta-Learning by Landmarking Various Learning Algorithms.” In *International Conference on Machine Learning (Icml)*, 743–50.

Smith-Miles, Kate A. 2009. “Cross-Disciplinary Perspectives on Meta-Learning for Algorithm Selection.” *ACM Computing Surveys (CSUR)* 41 (1). ACM: 6.