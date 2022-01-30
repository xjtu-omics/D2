# D2: Compute DNA Density and Distance to periphery.
D2 is a fast and accurate tool for computing **D**NA **D**ensity and **D**istance to periphrery (DisTP). In our paper, we show that the DNA density and DisTP are highly correlated with nuclear activities. For detailed tutorial, please visit [wiki](https://github.com/xjtu-omics/D2/wiki).

paper.

![D2-logo](https://user-images.githubusercontent.com/37327473/151689962-b02ea629-3d7a-40bb-ab64-68ebb81a2594.png)

# Requirements
D2 was tested on Python v3.9.1, with the following basic requirements:

 * NumPy (tested on v1.20.2)
 * SciPy (tested on v1.6.2)
 * Pandas (tested on v1.2.4)
 * Seaborn (tested on v0.11.1)
 * Matplotlib (tested on v3.3.4)
 
# Typical workflow
**D2** is a user-friendly scripts, which takes [3dg files](https://github.com/xjtu-omics/D2/wiki/File-Format#3dg-file) and [an index file](https://github.com/xjtu-omics/D2/wiki/File-Format#index-file) as input, and outputs [bed-like file](https://github.com/xjtu-omics/D2/wiki/File-Format#den_dtp) storing the DNA density and DisTP. 

Besides, we provided here the scripts for computing the **enrichments of genetic markers**, including the construction of density-DisTP matrix and calculation of enrichments. In order to achieve it, users should apply the [genetic markers data](https://github.com/xjtu-omics/D2/wiki/File-Format#marker-file) **of same reference genome**.

Below is a typicle workflow using the [**test_data**](https://github.com/xjtu-omics/D2/tree/main/test_data).
## Compute DNA Density and DisTP
**D2.py D2** and **D2.py D2s** compute the DNA density and DisTP. The resulted DNA and DisTP are stored in bed-like format file.
  ```
  cd PATH/WHERE/D2/AT
  mkdir -p test_result/den_dtp
  
  python D3.py D3s test_data/dg_files test_data/hg19_diplo_20k.window.bed test_result/den_dtp
  ```
## Construct Density-DisTP Matrix
**D2.py sta** gives the density and DisTP ranges, and a scatter plot as below.

**D2.py map** puts the bins on density-DisTP matrix, and stores the probability of genomic bins appearing at matrix bins (states).

**D2.py ave** computes the average and standard deviation (sd) of density and DisTP.
  ```
  python D2.py sta test_result/den_dtp/den_dtp test_data/hg19_diplo_20k.window.bed test_result/test_map_sta
  python D2.py map test_result/den_dtp/den_dtp test_data/hg19_diplo_20k.window.bed test_result/test_map
  python D2.py ave test_map.txt test_ave.txt
  ```
 ![test_map_sta](https://user-images.githubusercontent.com/37327473/133371032-8a9061b8-c91f-4b9b-a143-a850fcafa32f.png)

## Marker Enrichment
**D2.py marks** indexes and concatenates the markers.
**D2.py enrich** plotted the enrichments of markers individually.
**D2.py hiera** ranks the matrix bins (states) by hierarchy cluster.
```
python D2.py marks test_data/marks/ test_data/hg19_diplo_20k.window.bed test_mark.mark.txt
mkdir mark_enrich
python D2.py enrich test_map.txt test_mark.mark.txt mark_enrich/gm12878
python D2.py hiera test_map.txt test_mark.mark.txt test_hiera
```
Markers Enrichments Results:
[gm12878_1-Active-Promoter_histplot.pdf](https://github.com/cyz0315/D2/files/7166776/gm12878_1-Active-Promoter_histplot.pdf)
[gm12878_13-Heterochrom_histplot.pdf](https://github.com/cyz0315/D2/files/7166777/gm12878_13-Heterochrom_histplot.pdf)

Hierarchy cluster Results:
[test_hiera_value_hierarchy.pdf](https://github.com/cyz0315/D2/files/7166779/test_hiera_value_hierarchy.pdf)
[test_hiera_hierarchy_hist.pdf](https://github.com/cyz0315/D2/files/7166780/test_hiera_hierarchy_hist.pdf)


# License
D2 is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University. For more information, please contact with Yizhuo che (cyz0315@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).
