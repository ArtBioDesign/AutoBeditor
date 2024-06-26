  
# AutoBeditor
## Project Introduction  
**AutoBeditor** is an editing sequence design tool for CRISPR/Cas mediated base editing technology.
![AutoPMD](https://github.com/editSeqDesign/AutoBeditor/blob/main/img/AutoBeditor1.png) 

## The main application scenarios of this software tool include:
### 1. sgRNA desgin
- **Description**：Support the design of high-throughput sgRNA for base editing.
### 2. Full process editing sequence design
- **Description**：Support the design of full process editing Sequence for base editing.


  
## Installation
### python packages 
We suggest using Python 3.8 for AutoBeditor.

```shell
conda env create -f environment.yml
```


## Usage & Example

### 1. sgRNA desgin
- **Step 1:** Upload the target information(TSV) file 、the genome(FAN) file、 the genome annotation(GFF) file to be edited.
- **Step 2:** provide the necessary configuration information.
    - Example configuration (json):
      ```json
       param={ 
        "genome_release" : "GCF001",
        "host" : "yb01",
        "be_names" : ["Target-AID"],
        "pams" : ["NG"]
        } 
      ```  
**Execute:**

```shell
python ./beditor/main.py 
```
**Output:**
- `SgRNA Design Catalog(YB01)` 

These files will be generated in the `XXX/AutoBeditor/data/output/` directory.  



### 2. Full process editing sequence design

**Input:**
- **Step 1:** Upload the plasmid template(gb) file、the target information(TSV) file 、the genome(FAN) file、 the genome annotation(GFF) file、 the enzyme(CSV) file to be edited.
- **Step 2:** provide the necessary configuration information.
    - Example configuration (json):
      ```json
       param={ 
        "genome_release" : "GCF001",
        "host" : "yb01",
        "be_names" : ["Target-AID"],
        "pams" : ["NG"]
        } 
      ```   
      

**Execute:**

```shell
python main.py 
```
**Output:**
- `SgRNA Design Catalog(YB01)` 
- `result.xlsx ` 


These files will be generated in the `XXX/AutoBeditor/data/output/` directory.  
