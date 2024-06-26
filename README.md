# GaiaCatalogQuery
A high-performance tool for quick searches and cross-matching in the Gaia DR3 catalog (or other catalogs).

This tool is based on HEALPix and the split algorithms described in my paper. It efficiently partitions large (~TB sizes) catalogs into files of approximately equal size according to the varying object densities in different regions. During searches, only the files corresponding to the relevant regions are read, and data structures such as KD-Tree are utilized to achieve very fast catalog queries. 

The tool has already been applied in many astronomical projects and has demonstrated a high degree of stability.

## How to Use

### Step 1
* Download the Gaia catalog (full or partial) from [this link](https://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/) or other mirrors, and uncompress the files to a folder.

### Step 2
* Clone this project:
  ```shell
  git clone https://github.com/yourusername/GaiaCatalogQuery.git
  ```
* Compile the project with the following commands:
  ```shell
  mkdir build
  cd build
  cmake ..
  make -j
  ```

### Step 3
Due to memory limits, large catalogs such as Gaia DR3 need to be first split into a temporary folder and then merged to obtain the final partitioned results.
* Split the Gaia catalog into a temporary folder, for example:
  ```shell
  ./GaiaCatalogQuery split /home/user/Gaia_dr3/gaia_source_uncompressed/ /home/user/gaia_temp/ 1000 20000000
  ```
    * 1000 is the single file size threshold. Setting it too high or too low is not advisable. (For HDD, 1000 is recommended. For SSD, you can set it higher, like 3000)
    * 20000000 corresponds to approximately 30 GB of memory usage. You can adjust this based on your memory size (larger values are better).
    * This process may take some time.
* Merge the temporary files into the final folder:
  ```shell
  ./GaiaCatalogQuery merge /home/user/gaia_temp/ /home/user/gaia_final/ 1000 20000000
  ```

### Step 4
After merging, you can query the Gaia catalog using the following command:
  ```shell
  ./GaiaCatalogQuery search ra dec radius(degree)
  ```
The result will be saved in `/dev/shm/results.csv`.

---

If you find this tool helpful, please refer to the following paper:

Zhang, Qi-qian, Dong-wei Fan, and Chen-zhou Cui. "The Efficient Indexing and Fusion Algorithms for Large-scale Catalogs Based on File." Progress in Astronomy 41 (2023): 429-447. [https://ui.adsabs.harvard.edu/abs/2023PrA....41..429Z/abstract](https://ui.adsabs.harvard.edu/abs/2023PrA....41..429Z/abstract)