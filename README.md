
# GaiaCatalogQuery
A high-performance tool for quick searches and cross-matching in the Gaia DR3 catalog (or other catalogs).

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
* Split the Gaia catalog into a temporary folder, for example:
  ```shell
  ./GaiaCatalogQuery split /home/user/Gaia_dr3/gaia_source_uncompressed/ /home/user/gaia_temp/ 1000 20000000
  ```
  This process may take some time.
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
Zhang, Qi-qian, Dong-wei Fan, and Chen-zhou Cui. "The Efficient Indexing and Fusion Algorithms for Large-scale Catalogs Based on File." Progress in Astronomy 41 (2023): 429-447.