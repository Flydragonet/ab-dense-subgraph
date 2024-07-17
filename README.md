# (alpha, beta)-dense subgraph computation

# Datasets

Because some datasets used in the paper are too large to be uploaded to GitHub, we have summarized the download links for the dataset in the table below.

Datasets used for performance studies.

| Dataset | Link |
| --- | --- |
| DBpedia | http://www.konect.cc/networks/dbpedia-genre/ |
| Digg | http://www.konect.cc/networks/digg-votes/ |
| Enron | http://www.konect.cc/networks/bag-enron/ |
| IMDB | http://www.konect.cc/networks/actor2/ |
| Livejournal | http://www.konect.cc/networks/livejournal-groupmemberships/ |
| Yahoo | http://www.konect.cc/networks/yahoo-song/ |
| Orkut | http://www.konect.cc/networks/orkut-groupmemberships/ |
| Twitter | http://www.konect.cc/networks/twitter/ |

# Preprocess

The dataset needs to be preprocessed into a specific format file to be input by the algorithm, with the format as follows: The first line consists of three integers representing |E|, |U|, |V|. Each subsequent line contains two numbers representing an edge (u, v). Note that |U| should be equal to the maximum value of u among all edges, and |V| should be equal to the maximum value of v among all edges. An example is given below, representing a (2,2)-clique:

```
4 2 2
1 1
1 2
2 1
2 2
```

# Compile

```
g++ code.cpp -o code -std=c++11 -O3
```

# Usage

```
./code <dataset_address>
```