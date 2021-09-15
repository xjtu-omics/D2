# D3: Compute DNA Density and Distance to periphery.
D3 is a fast and accurate tool for computing **D**NA **D**ensity and **D**istance to periphrery (DisTP). In our paper, we show that the DNA density and DisTP are highly correlated with nuclear activities. 

paper.

![D3-logo](https://user-images.githubusercontent.com/37327473/133361303-b0526de0-b8d8-429a-a61d-260cd74e2111.png)

# Requirements
D3 was tested on Python v3.9.1 (macOS and DeepinOS), with the following basic requirements:

 * NumPy (tested on v1.20.2)
 * SciPy (tested on v1.6.2)
 * Pandas (tested on v1.2.4)
 * Seaborn (tested on v0.11.1)
 * Matplotlib (tested on v3.3.4)
 
#  Tipycal Workflow
Below is a typical workflow using the test data.
## Compute DNA Density and DisTP
D3.py D3 and D3.py D3s compute the DNA density and DisTP. The resulted DNA and DisTP are stored in bed-like format file.
  ```
  cd PATH/WHERE/D3/AT
  mkdir -p test_result/den_dtp
  
  python D3.py D3s test_data/dg_files test_data/hg19_diplo_20k.window.bed test_result/den_dtp
  ```
## Construct Density-DisTP Matrix
D3.py sta gives the density and DisTP ranges, and a scatter plot as below.

D3.py map puts the bins on density-DisTP matrix, and stores the probability of genomic bins appearing at matrix bins (states).

D3.py ave computes the average and standard deviation (sd) of density and DisTP.
  ```
  python D3.py sta test_result/den_dtp/den_dtp test_data/hg19_diplo_20k.window.bed test_result/test_map_sta
  python D3.py map test_result/den_dtp/den_dtp test_data/hg19_diplo_20k.window.bed test_result/test_map
  python D3.py ave test_map.txt test_ave.txt
  ```
 ![test_map_sta](https://user-images.githubusercontent.com/37327473/133371032-8a9061b8-c91f-4b9b-a143-a850fcafa32f.png)

## Marker Enrichment
D3.py marks indexes and concatenates the markers.
D3.py enrich plotted the enrichments of markers individually.
D3.py hiera ranks the matrix bins (states) by hierarchy cluster.
```
python D3.py marks test_data/marks/ test_data/hg19_diplo_20k.window.bed test_mark.mark.txt
mkdir mark_enrich
python D3.py enrich test_map.txt test_mark.mark.txt mark_enrich/gm12878
python D3.py hiera test_map.txt test_mark.mark.txt test_hiera
```
Markers Enrichments:
[gm12878_1-Active-Promoter_histplot.pdf](https://github.com/cyz0315/D3/files/7166776/gm12878_1-Active-Promoter_histplot.pdf)
[gm12878_13-Heterochrom_histplot.pdf](https://github.com/cyz0315/D3/files/7166777/gm12878_13-Heterochrom_histplot.pdf)

Hierarchy cluster:
[test_hiera_value_hierarchy.pdf](https://github.com/cyz0315/D3/files/7166779/test_hiera_value_hierarchy.pdf)
[test_hiera_hierarchy_hist.pdf](https://github.com/cyz0315/D3/files/7166780/test_hiera_hierarchy_hist.pdf)


# License
D3 is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University. For more information, please contact with Yizhuo che (cyz0315@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).
