# getRunSummary
Oxford Nanopore Promethion run summary
getRunSummary.py is a python script that takes Experiment ID as input and outputs PromethION run summary for each positions and for each barcodes in a table format. The script utilises [minknow_api](minknow_api) to fetch information.  
The script will create a .tsv & xlsx file with "User given Experiment ID"_summary.tsv naming convention.

## Dependecies Installation
```sh
pip install minknow_api pandas colorama tqdm
```

## usage: 
```sh
python getRunSummary.py -ho localhost -p 9502 -e 6123123
```
#### Parameter explanation

```sh
General options:  
  -h, --help        show the help and exit  
  -ho, --host       localhost OR hostname of the device
  -p, --port        grpc port : 9501 OR 9502
  -e, --exid        experiment ID
```
