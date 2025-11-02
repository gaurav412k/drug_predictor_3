README

This script predicts the Primary and Seocndary drug based on bacterias and genes detected.
Make sure python is installed (latest version) in your system.

All the dependencies will be handled by the script, so you don't have to worry about that.

First prepare the input csv file, you can check the input_data.csv
it should be as following:-

type,input_feature,ct_mean
organism,Escherichia coli,22.5
genes,CTX-M,30.1
genes,tetM,28.9

once your input_data.csv is ready, now you can run the script
command:-
    python drug_prediction_script.py input_data.csv

This command we generate the result in excel format as well as it will be printed on the terminal

Currently the script can only handle single patient at a time. So modify the scipt for multiple sample accordingly.