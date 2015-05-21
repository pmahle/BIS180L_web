---
layout: lab
title: Genetic Networks 1 - Clustering
hidden: true
tags:
     - networks
     - R
     - clustering
     - kmeans
     - hierarchical
     - heatmap
     - RNA-seq 
---

# Assignment template
Please download the [assignment template](http://jnmaloof.github.io/BIS180L_web/Assignments/Assignment_6_template.Rmd).  Place you answers to the exercises in this template.  When you are done push the .Rmd and a knitted .html file.  Start and issue indicating that the assignment is ready to grade.

# Clustering Introduction
As you learned last week, when we are dealing with genome scale data it is hard to come up with very specific summaries of the data unless you know exactly the question you are trying to ask computationally. Today we will be talking about three different ways to cluster data and get visual summaries of the expression of all genes that had a significant GxE interaction. Once we have these clusters, it allows us to ask further, more detailed questions such as what GO categories are enriched in each cluster, or are there specific metabolic pathways contained in the clusters? While clustering can be used in an exploratory way, the basics you will be learning today have been extended to very sophisticated statistical/machine learning methods used across many disciplines. In fact, there are many different methods used for clustering in R outlined in this **[CRAN VIEW](http://cran.r-project.org/web/views/Cluster.html)**.

The two clustering methods that we will be exploring are hierarchical clustering and k-means. These have important similarities and differences that we will discuss in detail throughout today. The basic idea with clustering is to find how similar the rows and/or columns in the dataset are based on the values contained within the data frame. You have already used a similar technique last week when you produced the MDS plot. This visualization of the the samples in the dataset showed that samples from similar genotype and treatment combinations were plotted next to one another based on their Biological Coefficient of Variation calculated across all of the counts of genes. 

# Hierarchical Clustering
An intuitive example is clustering the distances between know geographic locations, such as US Cities. I took this example from this [website](http://www.analytictech.com/networks/hiclus.htm).
 
1. Calculate a distance measure between all the row and column combinations in data set (think geographic distances between cities)
2. Start all the items in the data out as belonging to its own cluster (every city is its own cluster)
3. In the distance matrix, find the two closest points (find shortest distance between any two cities) 
4. Merge the two closest points into one cluster (merge BOS and NY in our example dataset)
5. Repeat steps 3 and 4 until all the items belong to a single large cluster

A special note: all the clusters at each merge take on the shortest distance between any one member of the cluster and the remaining clusters. For example, the distance between BOS and DC is 433 miles, but the distance between NY and DC is 233. Because BOS/NY are considered one cluster after our first round, their cluster distance to DC is 233. All three of these cities are then merged into one cluster DC/NY/BOS.



```bash
wget https://raw.githubusercontent.com/rjcmarkelz/BIS180L_web/feature-networks/data/us_cities.txt
```


```r
# make sure to change the path to where you downloaded this using wget
cities <- read.table("../data/us_cities.txt", sep = "\t", header = TRUE)

rownames(cities) <- cities$X # make the first column the rownames
cities <- cities[,-1] #remove first column
head(cities)
```

```
##       BOS   NY   DC  MIA  CHI  SEA   SF   LA  DEN
## BOS     0  206  429 1504  963 2976 3095 2979 1949
## NY    206    0  233 1308  802 2815 2934 2786 1771
## DC    429  233    0 1075  671 2684 2799 2631 1616
## MIA  1504 1308 1075    0 1329 3273 3053 2687 2037
## CHI   963  802  671 1329    0 2013 2142 2054  996
## SEA  2976 2815 2684 3273 2013    0  808 1131 1307
```

```r
plot(hclust(dist(cities))) 
```

![plot of chunk unnamed-chunk-2]({{ site.baseurl }}/figure/unnamed-chunk-2-1.png) 
**Exercise 1:**
Extending the example that I gave for BOS/NY/DC, what are the distances that define each split in the West Coast side of the hclust plot? 
*Hint 1: Start with the distances between SF and LA. Then look at the difference between that cluster up to SEA*
*Hint 2: Print cities, you only need to look at the upper right triangle of data matrix.*

What is the city pair and distance the joins the East Coast and West Coast cities? Fill in the values.
Hint: Think Midwest.


Now that we have that example out of the way, lets start using this technique on biological data. This week will be a review of all the cool data manipulation steps you have learned in past weeks. I would like to emphasize that printing dataframes (or parts of them) is a really fast way to get an idea about the data you are working with. Visual summaries like printed data, or plotting the data are often times the best way to make sure things are working the way they should be and you are more likely to catch errors. I have included visual summaries at all of the points where I would want to check on the data. 

If you remember last week, you found some genes that had significant GxE in the internode tissue. We are going to be taking a look at those same genes again this week. I followed the same tutorial that you did last week. That dataset only had 12 RNA-seq libraries in it. However, that subset of data was part of a much larger study that we are going to explore today. This dataset consists of RNA-seq samples collected from 2 genotypes of **Brassica rapa** (R500 and IMB211) that were grown in either dense (DP) or non-dense planting (NDP) treatments. Upendra also collected tissue from multiple tissue types including: Leaf, Petiole, Internode (you worked with this last week), and silique (the plant seed pod). There were also 3 biological replicates of each combination (Rep 1, 2, 3). If your head is spinning thinking about all of this, do not worry, data visualization will come to the rescue here in a second.

Remember last week when we were concerned with what distribution the RNA-seq count data was coming from so that we could have a good statistical model of it? Well, if we want to perform good clustering we also need to think about this because most of the simplest clustering assumes data to be from a normal distribution. I have precalculated this for you. First you will load my output from going through the tutorial from last week. Then you will pull in the new large table containing all the samples.

Lets start by loading the two datasets. Then we can subset the larger full dataset to include only the genes that we are interested in from our analysis last week.


```bash
wget https://raw.githubusercontent.com/rjcmarkelz/BIS180L_web/feature-networks/data/DEgenes_GxE.csv
wget https://raw.githubusercontent.com/rjcmarkelz/BIS180L_web/feature-networks/data/voom_transform_brassica.csv
```


```r
# make sure to change the path to where you downloaded this using wget
DE_genes <- read.table("../data/DEgenes_GxE.csv", sep = ",")
head(DE_genes) #check out the data
```

```
##               logFC   logCPM        LR       PValue          FDR
## Bra010821  6.657165 6.040918 153.15544 3.542436e-35 8.513891e-31
## Bra033034 -4.639893 6.611694  71.26078 3.129960e-17 3.761273e-13
## Bra035334 -4.142456 4.665331  68.34729 1.370938e-16 1.098304e-12
## Bra003598 -4.893648 3.900122  54.16121 1.846974e-13 1.109755e-09
## Bra016182  7.714360 6.330237  49.61407 1.871665e-12 8.996719e-09
## Bra013164  4.098373 9.329089  47.94483 4.383823e-12 1.756013e-08
```

```r
DE_gene_names <- row.names(DE_genes)

# make sure to change the path to where you downloaded this using wget
brass_voom_E <- read.table("../data/voom_transform_brassica.csv", sep = ",", header = TRUE)

GxE_counts <- as.data.frame(brass_voom_E[DE_gene_names,])

GxE_counts <- as.matrix(GxE_counts) # some of the downstream steps require a data matrix
```

Now that we have a dataframe containing our 255 significant GxE genes from the internode tissue, we can take a look at how these genes are acting across all tissues.


```r
hr <- hclust(dist(GxE_counts))
plot(hr)
```

![plot of chunk unnamed-chunk-5]({{ site.baseurl }}/figure/unnamed-chunk-5-1.png) 
What a mess! We are too overplotted in this direction. How about if we cluster by column instead? Notice we have to transpose the data using *t()*. Also, make sure you stretch out the window so you can view it! 


```r
hc <- hclust(dist(t(GxE_counts)))
plot(hc)
```

![plot of chunk unnamed-chunk-6]({{ site.baseurl }}/figure/unnamed-chunk-6-1.png) 
**Exercise 2:**
What is the general pattern in the h-clustering data? 
Using what you learned from the city example, what is the subcluster that looks very different than the rest of the samples? 
*Hint: It is a group of 3 libraries. You will have to plot this yourself and stretch it out. The rendering on the website compresses the output.*

We have obtained a visual summary using h-clustering. Now what? Well we can go a little further and start to define some important sub-clusters within our tree. We can do this using the following function. Once again make sure you plot it so you can stretch the axes. 


```r
?rect.hclust
hc <- hclust(dist(t(GxE_counts)))
plot(hc) #redraw the tree everytime before adding the rectangles
rect.hclust(hc, k = 4, border = "red")
```

![plot of chunk unnamed-chunk-7]({{ site.baseurl }}/figure/unnamed-chunk-7-1.png) 
**Exercise 3:**
With k = 4 as one of the arguments to the rect.hclust() function, what is the largest and smallest group contained within the rectangles? 
Play with the k-values between 3 and 7. Describe how the size of the clusters change when changing between k-values.

You may have noticed that your results and potential interpretation of the data could change very dramatically based on the how many subclusters you choose! This is one of the main drawbacks with this technique. Fortunately there are other packages that can help us determine what sub-clusters have good support. 


```r
install.packages("pvclust")
```

```
## Error in install.packages : Updating loaded packages
```

```r
library(pvclust)
?pvclust #check out the documentation

set.seed(12456) #important to run this
fit <- pvclust(GxE_counts, method.hclust = "ward.D", method.dist = "euclidean", nboot = 50)
```

```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.6)... Done.
## Bootstrap (r = 0.7)... Done.
## Bootstrap (r = 0.8)... Done.
## Bootstrap (r = 0.9)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.1)... Done.
## Bootstrap (r = 1.2)... Done.
## Bootstrap (r = 1.3)... Done.
## Bootstrap (r = 1.4)... Done.
```

```r
plot(fit) # dendogram with p-values
```

![plot of chunk unnamed-chunk-8]({{ site.baseurl }}/figure/unnamed-chunk-8-1.png) 
Normally we would do 1000+ bootstrap samples, to get our support for each of the branches in the tree, but we do not have lots of time today. The red values are the "Approximate Unbiased" (AU) values with numbers closer to 100 providing more support. Bootstrapping is a popular resampling technique that you can read about more *[here](http://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29)*.

**Exercise 4:**
After running the 50 bootstrap samples, leave the plot open. Then change nboot up to 500. In general what happens to AU comparing the two plots by flipping between them?


We will be discussing more methods for choosing the number of clusters in the k-means section. Until then, we will expand what we learned about h-clustering to do a more sophisticated visualization of the data. 

# Heatmaps
Heatmaps are another way to visualize h-clustering results. Instead of looking at either the rows (genes) or the columns (samples) like we did with the h-clustering examples we can view the entire data matrix at once. Well how do we do this? We could just print the matrix, but that would just be a bunch of numbers. Heatmaps take all the values within the data matrix and convert them to a color value. The human eye is really good at picking out patterns so lets convert that data matrix to a color value AND do some h-clustering. Although we can do this really easily with the heatmap() function that comes preloaded in R, there is a small library that provides some additional functionality to the heatmaps. We will start with the cities example because it is small and easy to see what is going on. 

*A general programming tip: always have little test data sets. That allows you to figure out what the functions are doing. If you understand how it works on a small scale then you will be better able to troubleshoot when scaling to large datasets.*


```r
install.packages("gplots") #not to be confused with ggplot2!
```

```
## Error in install.packages : Updating loaded packages
```

```r
library(gplots)
head(cities) # city example
```

```
##       BOS   NY   DC  MIA  CHI  SEA   SF   LA  DEN
## BOS     0  206  429 1504  963 2976 3095 2979 1949
## NY    206    0  233 1308  802 2815 2934 2786 1771
## DC    429  233    0 1075  671 2684 2799 2631 1616
## MIA  1504 1308 1075    0 1329 3273 3053 2687 2037
## CHI   963  802  671 1329    0 2013 2142 2054  996
## SEA  2976 2815 2684 3273 2013    0  808 1131 1307
```

```r
cities_mat <- as.matrix(cities)
cityclust <- hclust(dist(cities_mat))
?heatmap.2 #take a look at the arguments
heatmap.2(cities_mat, Rowv=as.dendrogram(cityclust), scale="row", density.info="none", trace="none")
```

![plot of chunk unnamed-chunk-9]({{ site.baseurl }}/figure/unnamed-chunk-9-1.png) 
**Exercise 5:**
We used the scale rows option. This is necessary so that every *row* in the dataset will be on the same scale when visualized in the heatmap. This is to prevent really large values somewhere in the dataset dominating the heatmap signal. Remember if you still have this dataset in memory you can take a look at a printed version to the terminal. Compare the distance matrix that you printed with the colors of the heat map. See the advantage of working with small test sets? Take a look at your plot of the cities heatmap and interpret what a dark red value and a light yellow value in the heatmap would mean in geographic distance. Provide an example of of each in your explanation.


##Now for some gene expression data. 


```r
hr <- hclust(dist(GxE_counts))
plot(hr)
```

![plot of chunk unnamed-chunk-10]({{ site.baseurl }}/figure/unnamed-chunk-10-1.png) 

```r
heatmap.2(GxE_counts, Rowv = as.dendrogram(hr), scale = "row", density.info="none", trace="none")
```

![plot of chunk unnamed-chunk-10]({{ site.baseurl }}/figure/unnamed-chunk-10-2.png) 
**Exercise 6:** The genes are definitely overplotted and we cannot tell one from another. However, what is the most obvious pattern that you can pick out from this data? Describe what you see. Make sure you plot this in your own session so you can stretch it out.
*Hint It will be a similar pattern as you noticed in the h-clustering example.*

**Exercise 7:** In the similar way that you interpreted the color values of the heatmap for the city example, come up with a biological interpretation of the yellow vs. red color values in the heatmap. What would this mean for the pattern that you described in exercise 6? Discuss.

# k-means clustering
K-means tries to fit "centering points" to your data by chopping your data up into however many mean centers you specify. If you pick 3, then you are doing 3-means centering to your data or trying to find the best 3 centers that describe all of your data. To build some intuition about how these things work there is a cool R package called "animation". You do not have to run the following code, but I will use it to demonstrate a few examples with the projector. These examples are based directly on the **[documentation](cran.r-project.org/web/packages/animation/animation.pdf)** for the functions so if you wanted to look at how some other commonly used methods work with a visual summary I encourage you to check out this package on your own time.


```r
# you do not have to run this code chunk
# library(animation) 

# kmeans.ani(x = cbind(X1 = runif(50), X2 = runif(50)), centers = 3,
# hints = c("Move centers!", "Find cluster?"), pch = 1:3, col = 1:3)

# kmeans.ani(x = cbind(X1 = runif(50), X2 = runif(50)), centers = 10,
# hints = c("Move centers!", "Find cluster?"), pch = 1:3, col = 1:3)

# kmeans.ani(x = cbind(X1 = runif(50), X2 = runif(50)), centers = 5,
# hints = c("Move centers!", "Find cluster?"), pch = 1:3, col = 1:3)
```

Now that you have a sense for how this k-means works, lets apply what we know to our data.  Lets start with 9 clusters, but please play with the number of clusters and see how it changes the visualization. We will need to run a quick Principle Component Analysis (similar to MDS), on the data in order to visualize how the clusters are changing with different k-values.


```r
library(ggplot2)
prcomp_counts <- prcomp(t(GxE_counts)) #gene wise
scores <- as.data.frame(prcomp_counts$rotation)[,c(1,2)]

set.seed(25) #make this repeatable as kmeans has random starting positions
fit <- kmeans(GxE_counts, 9)
clus <- as.data.frame(fit$cluster)
names(clus) <- paste("cluster")

plotting <- merge(clus, scores, by = "row.names")
plotting$cluster <- as.factor(plotting$cluster)

# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 4, stat = "identity") 
```

![plot of chunk unnamed-chunk-12]({{ site.baseurl }}/figure/unnamed-chunk-12-1.png) 
**Exercise 8:** Pretty Colors! Describe what you see visually with 2, 5, 9, and 15 clusters using either method. Why would it be a bad idea to have to few or to many clusters? Discuss with a specific example comparing few vs. many k-means. Justify your choice of too many and too few clusters by describing what you see in each case.

The final thing that we will do today is try to estimate, based on our data, what the ideal number of clusters is. For this we will use something called the Gap statistic. 



```r
install.packages("cluster")
```

```
## Error in install.packages : Updating loaded packages
```

```r
library(cluster)
set.seed(125)
gap <- clusGap(GxE_counts, FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
```

```
## Clustering k = 1,2,..., K.max (= 20): .. done
## Bootstrapping, b = 1,2,..., B (= 50)  [one "." per sample]:
## .................................................. 50
```

```r
plot(gap, main = "Gap Statistic")
```

![plot of chunk unnamed-chunk-13]({{ site.baseurl }}/figure/unnamed-chunk-13-1.png) 
This is also part of the cluster package that you loaded earlier. It will take a few minutes to calculate this statistic. In the mean time, check out some more information about it in the ?clusGap documentation. We could imagine that as we increase the number of k-means to estimate, we are always going to increase the fit of the data. The extreme examples of this would be if we had k = 255 for the total number of genes in the data set or k = 1. We would be able to fit the data perfectly in the k = 255 case, but what has it told us? It has not really told us anything. Just like you played with the number of k-means in Exercise 8, we can also do this computationally! We want optimize the fewest number of clusters that we have that can explain the variance in our dataset.

**Exercise 9:** Based on this Gap statistic plot at what number of k clusters (x-axis) do you start to see diminishing returns? To put this another way, at what value of k does k-1 and k+1 start to look the same for the first time? Or yet another way, when are you getting diminishing returns for adding more k-means? See if you can make the tradeoff of trying to capture a lot of variation in the data as the Gap statistic increases, but you do not want to add too many k-means because your returns diminish as you add more. Explain your answer using the plot as a guide to help you interpret the data.

Now we can take a look at the plot again and also print to the terminal what clusGap() calculated. 


```r
with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
```

```
## [1] 4
```
**Exercise 10:** What did clusGap() calculate? How does this compare to your answer from Exercise 9? Make a plot using the combined autoplot() and kmeans functions as you did before, but choose the number of k-means you chose and the number of k-means that are calculated from the Gap Statistic. Describe the differences in the plots.

Good Job Today! There was a lot of technical stuff to get through. We are going to build on all of this Thursday to construct co-expression networks and study their properties using a few of the techniques that you learned today. 











